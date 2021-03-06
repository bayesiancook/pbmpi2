
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include <cassert>
#include "Parallel.h"
#include <string.h>
#include "StringStreamUtils.h"

#include "MultiGenePhyloProcess.h"


void MultiGenePhyloProcess::New(int unfold)	{


	CreateMPI(0);
	AllocateAlignments(datafile);
	SetProfileDim();
	SetTree(treefile);

	Create();

	if (! GetMyid())	{
		GlobalBroadcastTree();
		MPI_Barrier(MPI_COMM_WORLD);

		if (topobf)	{
			if (topobf == 1)	{
				bffrac = 0;
			}
			if (topobf == 2)	{
				bffrac = -bfnfrac;
			}
			/*
			if (topobf == 3)	{
				cerr << "fix topo bf\n";
				cerr << "bffrac : " << bffrac << '\n';
			}
			*/
			// SetTopoBF();
		}

		if (sis == 1)	{
			sisfrac = 0;
			SetSIS();
		}
		if (sis == 2)	{
			SetSIS();
		}
		if (sis == 1)	{
			PriorSample();
		}
		else	{
			Sample();
		}

		GlobalUpdateParameters();
		GlobalSample();

		if (topobf)	{
			GlobalSetTopoBF();
		}
		GlobalUnfold();
	}
	else	{
		SlavePostNew();
		SlaveBroadcastTree();
		MPI_Barrier(MPI_COMM_WORLD);
	}

	if (BPP)	{
		BPP->RegisterWithTaxonSet(GetData()->GetTaxonSet());
	}
}

void MultiGenePhyloProcess::Open(istream& is, int unfold)	{

	CreateMPI(0);
	AllocateAlignments(datafile);
	SetProfileDim();

	tree = new Tree(GetData()->GetTaxonSet());
	if (GetMyid() == 0)	{
		istringstream s(treestring);
		tree->ReadFromStream(s);
		GlobalBroadcastTree();
	}
	else	{
		PhyloProcess::SlaveBroadcastTree();
	}

	Create();

	if (! GetMyid())	{

		GlobalBroadcastTree();
		MPI_Barrier(MPI_COMM_WORLD);

        // global method
		FromStream(is);

		if (sis)	{
			SetSIS();
		}

		GlobalUpdateParameters();

		if (topobf)	{
			GlobalSetTopoBF();
		}
		GlobalReadSiteRankFromStream(is);
		GlobalUnfold();
	}
	else	{
		SlavePostOpen();
		SlaveBroadcastTree();
		MPI_Barrier(MPI_COMM_WORLD);
        // starting from here, slave will enter wait loop
	}

	if (BPP)	{
		BPP->RegisterWithTaxonSet(GetData()->GetTaxonSet());
	}
}

void MultiGenePhyloProcess::SlavePostNew()	{
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->New(0);
		}
	}
}

void MultiGenePhyloProcess::SlavePostOpen()	{
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			string s = "dummy";
			istringstream is(s);
			process[gene]->Open(is,0);
			// process[gene]->SetTopoBF();
		}
	}
}

void MultiGenePhyloProcess::Create()	{

	PhyloProcess::Create();
	MultiGeneRateProcess::Create();
	MultiGeneBranchProcess::Create();
	MultiGeneMPIModule::Create();
}

void MultiGenePhyloProcess::Delete()	{

	MultiGeneMPIModule::Delete();
	MultiGeneBranchProcess::Delete();
	MultiGeneRateProcess::Delete();
	PhyloProcess::Delete();
}

void MultiGenePhyloProcess::AllocateAlignments(string datafile)	{

	ifstream is(datafile.c_str());
	is >> Ngene;
	ifstream* tis = 0;
	genename = new string[Ngene];
	genesize = new int[Ngene];
	genealloc = new int[Ngene];
	int* geneweight = new int[Ngene];
	// genedata = new SequenceAlignment*[Ngene];
	for (int gene=0; gene<Ngene; gene++)	{
		is >> genename[gene];
		SequenceAlignment* localdata = new FileSequenceAlignment(genename[gene]);
		int nstate = localdata->GetNstate();
		if (! gene)	{
			data = localdata;
		}
		else	{
			if (nstate != data->GetNstate())	{
				cerr << "error: all data files do not have the same alphabet\n";
				cerr << nstate << '\t' << data->GetNstate() << '\n';
				exit(1);
			}
		}

		genesize[gene] = localdata->GetNsite();
		geneweight[gene] = localdata->GetNsite() * localdata->GetNtaxa();
		if (gene)	{
			delete localdata;
		}
	}
	delete tis;
	tis = 0;

	// sort alignments by decreasing size
	int permut[Ngene];
	for (int gene=0; gene<Ngene; gene++)	{
		permut[gene] = gene;
	}
	for (int i=0; i<Ngene; i++)	{
		for (int j=Ngene-1; j>i; j--)	{
			if (geneweight[permut[i]] < geneweight[permut[j]])	{
			// if (genesize[permut[i]] < genesize[permut[j]])	{
				int tmp = permut[i];
				permut[i] = permut[j];
				permut[j] = tmp;
			}
		}
	}

	int totsize[nprocs];
	for (int i=0; i<nprocs; i++)	{
		totsize[i] = 0;
	}

	for (int i=0; i<Ngene; i++)	{
		int gene = permut[i];
		int size = geneweight[gene];
		// int size = genesize[gene];

		int min = 0;
		int jmin = 0;
		for (int j=1; j<nprocs; j++)	{
			if ((j==1) || (min > totsize[j]))	{
				min = totsize[j];
				jmin = j;
			}
		}
		genealloc[gene] = jmin;
		totsize[jmin] += size;
	}

	if (totsize[0])	{
		cerr << "error in alloc\n";
		exit(1);
	}
	int total = 0;
	for (int i=1; i<nprocs; i++)	{
		int tot = 0;
		for (int gene=0; gene<Ngene; gene++)	{
			if (genealloc[gene] == i)	{
				tot += geneweight[gene];
				// tot += genesize[gene];
				total++;
			}
		}
		if (tot != totsize[i])	{
			cerr << "error in allocation\n";
			cerr << tot << '\t' << totsize[i] << '\n';
			exit(1);
		}
	}
	if (total != Ngene)	{
		cerr << "error in total allocation\n";
		exit(1);
	}

	globalnsite = new int[nprocs];
	for (int i=0; i<nprocs; i++)	{
		globalnsite[i] = 0;
	}
	for (int gene=0; gene<Ngene; gene++)	{
		if ((genealloc[gene] < 0) || (genealloc[gene] >= nprocs))	{
			cerr << "alloc : " << genealloc[gene] << '\t' << gene << '\n';
			exit(1);
		}
		globalnsite[0] += genesize[gene];
		globalnsite[genealloc[gene]] += genesize[gene];
	}
	if (! myid)	{
		GlobalNsite = 0;
		cerr << '\n';
		cerr << "number of sites allocated to each slave:\n";
		for (int i=1; i<nprocs; i++)	{
			cerr << i << '\t' << globalnsite[i] << '\n';
			GlobalNsite += globalnsite[i];
		}
		cerr << '\n';
		cerr << "total: " << GlobalNsite << '\n';
	}
	
	// check total size
	if (! myid)	{
		int tot = 0;
		for (int i=1; i<nprocs; i++)	{
			tot += globalnsite[i];
		}
		if (tot != globalnsite[0])	{
			cerr << "error in total size during gene allocation\n";
			exit(1);
		}
		int tot2 = 0;
		for (int gene=0; gene<Ngene; gene++)	{
			tot2 += genesize[gene];
		}
		if (tot2 != tot)	{
			cerr << "error during alloc: total size does not match\n";
			exit(1);
		}
	}
	delete[] geneweight;
}

void MultiGenePhyloProcess::GlobalWriteSiteRankToStream(ostream& os)	{

	MESSAGE signal = WRITESITERANK;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	MPI_Status stat;
	for(int i=1; i<GetNprocs(); i++) {
		for (int gene=0; gene<Ngene; gene++)	{
			if (genealloc[gene] == i)	{
				int* rank = new int[genesize[gene]];
				MPI_Recv(rank,genesize[gene],MPI_INT,i,TAG1,MPI_COMM_WORLD,&stat);
				for (int k=0; k<genesize[gene]; k++)	{
					os << rank[k] << '\t';
				}
				os << '\n';
			}
		}
	}
}

/*
void MultiGenePhyloProcess::GlobalWriteSiteRankToStream(ostream& os)	{

	MESSAGE signal = WRITESITERANK;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	MPI_Status stat;
	for (int gene=0; gene<Ngene; gene++)	{
		int* rank = new int[genesize[gene]];
		MPI_Recv(rank,genesize[gene],MPI_INT,genealloc[gene],TAG1,MPI_COMM_WORLD,&stat);
		for (int k=0; k<genesize[gene]; k++)	{
			os << rank[k] << '\t';
		}
		os << '\n';
		delete[] rank;
	}
}
*/

void MultiGenePhyloProcess::SlaveWriteSiteRankToStream()	{

	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			MPI_Send(process[gene]->globalrank,genesize[gene],MPI_INT,0,TAG1,MPI_COMM_WORLD);
		}
	}
}

void MultiGenePhyloProcess::GlobalReadSiteRankFromStream(istream& is)	{

	cerr << "read site ranks\n";
	MESSAGE signal = READSITERANK;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	for(int i=1; i<GetNprocs(); i++) {
		for (int gene=0; gene<Ngene; gene++)	{
			if (genealloc[gene] == i)	{
				int* rank = new int[genesize[gene]];
				for (int k=0; k<genesize[gene]; k++)	{
					is >> rank[k];
				}
				MPI_Send(rank,genesize[gene],MPI_INT,i,TAG1,MPI_COMM_WORLD);
				delete[] rank;
			}
		}
	}
}

/*
void MultiGenePhyloProcess::GlobalReadSiteRankFromStream(istream& is)	{

	if (myid)	{
		for (int gene=0; gene<Ngene; gene++)	{
			int* rank = new int[genesize[gene]];
			for (int k=0; k<genesize[gene]; k++)	{
				is >> rank[k];
			}
			if (genealloc[gene] == myid)	{
				int* grank = process[gene]->globalrank;
				for (int k=0; k<genesize[gene]; k++)	{
					grank[k] = rank[k];
				}
			}
			delete[] rank;
		}
	}
}
*/

void MultiGenePhyloProcess::SlaveReadSiteRankFromStream()	{

	MPI_Status stat;
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			MPI_Recv(process[gene]->globalrank,genesize[gene],MPI_INT,0,TAG1,MPI_COMM_WORLD,&stat);
		}
	}
}


void MultiGenePhyloProcess::ToStream(ostream& os)	{

	MESSAGE signal = TOSTREAM;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	MultiGeneBranchProcess::ToStream(os);
	MultiGeneRateProcess::ToStream(os);

	string* geneparam = new string[Ngene];
	int* fullgeneparamsize = new int[Ngene];
	int* geneparamsize = new int[Ngene];
	int fullparamtotsize = 0;
	int paramtotsize = 0;

	MPI_Status stat;
	for(int i=1; i<GetNprocs(); i++) {
		MPI_Recv(&paramtotsize,1,MPI_INT,i,TAG1,MPI_COMM_WORLD,&stat);
		fullparamtotsize += paramtotsize;
		MPI_Recv(geneparamsize,Ngene,MPI_INT,i,TAG1,MPI_COMM_WORLD,&stat);
		char* c = new char[paramtotsize];
		MPI_Recv(c,paramtotsize,MPI_CHAR,i,TAG1,MPI_COMM_WORLD,&stat);
		int index = 0;
		for (int gene=0; gene<Ngene; gene++)	{
			if (genealloc[gene] == i)	{
				fullgeneparamsize[gene] = geneparamsize[gene];
				ostringstream os;
				for (int j=0; j<geneparamsize[gene]; j++)	{
					os << c[index++];
				}
				geneparam[gene] = os.str();
			}
		}
		if (index != paramtotsize)	{
			cerr << "error in MultiGenePhyloProcess::ToStream: non matching total size : " << index << '\t' << paramtotsize << '\n';
			exit(1);
		}
		delete[] c;
	}

	for (int gene=0; gene<Ngene; gene++)	{
        // os << "BEGIN\n";
		os << fullgeneparamsize[gene] << '\n';
		os << geneparam[gene];
		os << '\n';
        // os << "END\n";
	}

	delete[] fullgeneparamsize;
	delete[] geneparamsize;
	delete[] geneparam;
}

void MultiGenePhyloProcess::SlaveToStream()	{

	string* geneparam = new string[Ngene];
	int* geneparamsize = new int[Ngene];
	int paramtotsize = 0;

	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			ostringstream os;
			process[gene]->ToStream(os);
            // os.flush();
			geneparam[gene] = os.str();
			geneparamsize[gene] = geneparam[gene].size();
		}
		else	{
			geneparamsize[gene] = 0;
		}
		paramtotsize += geneparamsize[gene];
	}	
	char* c = new char[paramtotsize];
	int index = 0;
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			for (int j=0; j<geneparamsize[gene]; j++)	{
				c[index++] = geneparam[gene][j];
			}
		}
	}
	MPI_Send(&paramtotsize,1,MPI_INT,0,TAG1,MPI_COMM_WORLD);
	MPI_Send(geneparamsize,Ngene,MPI_INT,0,TAG1,MPI_COMM_WORLD);
	MPI_Send(c,paramtotsize,MPI_CHAR,0,TAG1,MPI_COMM_WORLD);

	delete[] c;
	delete[] geneparamsize;
	delete[] geneparam;
}

void MultiGenePhyloProcess::FromStream(istream& is)	{

	MESSAGE signal = FROMSTREAM;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	MultiGeneBranchProcess::FromStream(is);
	MultiGeneRateProcess::FromStream(is);

	string* geneparam = new string[Ngene];
	int* geneparamsize = new int[Ngene];
	int paramtotsize = 0;

	for (int gene=0; gene<Ngene; gene++)	{
        /*
        string tmp;
        is >> tmp;
        if (tmp != "BEGIN") {
            cerr << "begin error when reading multi gene from stream\n";
            exit(1);
        }
        */
		is >> geneparamsize[gene];
		paramtotsize += geneparamsize[gene];

		geneparam[gene].resize(geneparamsize[gene]);
		is.read(&geneparam[gene][0],geneparamsize[gene]);
        /*
        is >> tmp;
        if (tmp != "END")  {
            cerr << "end error when reading multi gene from stream\n";
            exit(1);
        }
        */
	}

	char* c = new char[paramtotsize];
	int index = 0;
	for (int gene=0; gene<Ngene; gene++)	{
		for (int j=0; j<geneparamsize[gene]; j++)	{
			c[index++] = geneparam[gene][j];
		}
	}

	MPI_Bcast(&paramtotsize,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(geneparamsize,Ngene,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(c,paramtotsize,MPI_CHAR,0,MPI_COMM_WORLD);

	delete[] c;
	delete[] geneparamsize;
	delete[] geneparam;
}

void MultiGenePhyloProcess::SlaveFromStream()	{

	string* geneparam = new string[Ngene];
	int* geneparamsize = new int[Ngene];
	int paramtotsize = 0;

	MPI_Bcast(&paramtotsize,1,MPI_INT,0,MPI_COMM_WORLD);

	MPI_Bcast(geneparamsize,Ngene,MPI_INT,0,MPI_COMM_WORLD);
	char* c = new char[paramtotsize];
	MPI_Bcast(c,paramtotsize,MPI_CHAR,0,MPI_COMM_WORLD);

	int index = 0;
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			
			ostringstream os;
			for (int j=0; j<geneparamsize[gene]; j++)	{
				os << c[index++];
			}
			string s = os.str();
			if (s.size() != geneparamsize[gene])	{
				cerr << "error in slave update params\n";
				cerr << geneparamsize[gene] << '\t' << s.size() << '\n';
				exit(1);
			}
			istringstream is(s);
			process[gene]->FromStream(is);
		}
		else	{
			index += geneparamsize[gene];
		}
	}

	delete[] c;
	delete[] geneparamsize;
	delete[] geneparam;

}

int MultiGenePhyloProcess::SpecialSlaveExecute(MESSAGE signal)	{

	switch(signal) {
	case TOSTREAM:
		SlaveToStream();
		return 1;
		break;
	case FROMSTREAM:
		SlaveFromStream();
		return 1;
		break;
	case SAMPLE:
		SlaveSample();
		return 1;
		break;
	case GENE_MOVE:
		SlaveGeneMove();
		return 1;
		break;
	case MEANALPHA:
		SlaveGetMeanAlpha();
		return 1;
		break;
	case COLLECTALPHA:
		SlaveCollectGeneAlphas();
		return 1;
		break;
	case MEANTOTLENGTH:
		SlaveGetMeanTotalLength();
		return 1;
		break;
	case COLLECTLENGTHS:
		SlaveCollectGeneBranchLengths();
		return 1;
		break;
	case FULLLIKELIHOOD:
		SlaveCollectFullLogLikelihood();
		return 1;
		break;
	case COLLECTLENGTHSUFFSTAT:
		SlaveCollectGeneLengthMappingSuffStat();
		return 1;
		break;
	case WRITESITERANK:
		SlaveWriteSiteRankToStream();
		return 1;
		break;
	case READSITERANK:
		SlaveReadSiteRankFromStream();
		return 1;
		break;
	case ACTIVATECPO:
		SlaveActivateCPO();
		return 1;
		break;
	case INACTIVATECPO:
		SlaveInactivateCPO();
		return 1;
		break;
	case COLLECTCPO:
		SlaveCollectPseudoMarginalLogLikelihood();
		return 1;
		break;
	case WRITESITELOGL:
		SlaveWriteSiteLogLToStream();
		return 1;
		break;

	default:
		return 0;
		// PhyloProcess::SlaveExecute(signal);
	}
}

void MultiGenePhyloProcess::GlobalWriteSiteLogLToStream(ostream& os)	{

	MESSAGE signal = WRITESITELOGL;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	MPI_Status stat;
    double* sitelogl = new double[globalnsite[0]];
    double* tmpsitelogl = new double[globalnsite[0]];
    for (int k=0; k<globalnsite[0]; k++)    {
        sitelogl[k] = 1.0;
    }
	for(int i=1; i<GetNprocs(); i++) {
        MPI_Recv(tmpsitelogl,globalnsite[i],MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
        int fromindex = 0;
        int index = 0;
		for (int gene=0; gene<Ngene; gene++)	{
			if (genealloc[gene] == i)	{
                for (int j=0; j<genesize[gene]; j++)    {
                    sitelogl[index+j] = tmpsitelogl[fromindex+j];
                }
                fromindex += genesize[gene];
			}
            index += genesize[gene];
		}
        if (fromindex != globalnsite[i])    {
            cerr << "error in MultiGenePhyloProcess::GlobalWriteSiteLogL: fromindex does not match\n";
        }
        if (index != globalnsite[0])    {
            cerr << "error in MultiGenePhyloProcess::GlobalWriteSiteLogL: index does not match\n";
        }
	}
    double total = 0;
    for (int k=0; k<globalnsite[0]; k++)    {
        if (sitelogl[k] == 1.0)  {
            cerr << "error in MultiGenePhyloProcess::GlobalWriteSiteLogL: out of frame\n";
            exit(1);
        }
        total += sitelogl[k];
    }
    os << total;
    for (int k=0; k<globalnsite[0]; k++)    {
        os << '\t' << sitelogl[k];
    }
    os << '\n';
    delete[] tmpsitelogl;
    delete[] sitelogl;
}

void MultiGenePhyloProcess::SlaveWriteSiteLogLToStream()	{
    double* sitelogl = new double[globalnsite[myid]];
    for (int k=0; k<globalnsite[myid]; k++) {
        sitelogl[k] = 1.0;
    }
    int index = 0;
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			genelnL[gene] = process[gene]->GetFullLogLikelihood(sitelogl+index);
            index += genesize[gene];
		}
	}
    if (index != globalnsite[myid]) {
        cerr << "error in MultiGenePhyloProcess::SlaveWriteSiteLogLToStream\n";
        exit(1);
    }
    MPI_Send(sitelogl,globalnsite[myid],MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
    delete[] sitelogl;
}

/*
void MultiGenePhyloProcess::SlaveResetTree()	{

	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->ResetTree();
		}
	}
}
*/

void MultiGenePhyloProcess::SlaveSetBFFrac()	{

	MPI_Bcast(&bffrac,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->bffrac = bffrac;
		}
	}
}

void MultiGenePhyloProcess::SlaveSetSISFrac()	{

	MPI_Bcast(&sisfrac,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->sisfrac = sisfrac;
		}
	}
}

void MultiGenePhyloProcess::GlobalSample()	{
	assert(myid == 0);
	MESSAGE signal = SAMPLE;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
}

void MultiGenePhyloProcess::SlaveSample()	{
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->Sample();
		}
	}
}

void MultiGenePhyloProcess::GlobalGeneMove()	{
	assert(myid == 0);
	MESSAGE signal = GENE_MOVE;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
}

void MultiGenePhyloProcess::SlaveGeneMove()	{
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->AugmentedMove();
		}
	}
}

void MultiGenePhyloProcess::SlaveBroadcastTree()	{

	int len;
	MPI_Bcast(&len,1,MPI_INT,0,MPI_COMM_WORLD);
	unsigned char* bvector = new unsigned char[len];
	MPI_Bcast(bvector,len,MPI_UNSIGNED_CHAR,0,MPI_COMM_WORLD);
	ostringstream os;
	for (int i=0; i<len; i++)	{
		os << bvector[i];
	}
	istringstream is(os.str());
	tree->ReadFromStream(is);
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			istringstream is(os.str());
			process[gene]->GetTree()->ReadFromStream(is);
			process[gene]->GetTree()->RegisterWith(GetData()->GetTaxonSet());
		}
	}
	delete[] bvector;
}

void MultiGenePhyloProcess::SlaveUnfold()	{
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->Unfold();
		}
	}
}

void MultiGenePhyloProcess::SlaveCollapse()	{
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->Collapse();
		}
	}
}

void MultiGenePhyloProcess::SlaveActivateSumOverRateAllocations()	{
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->ActivateSumOverRateAllocations();
			// process[gene]->sumratealloc = 1;
		}
	}
	// sumratealloc = 1;
}

void MultiGenePhyloProcess::SlaveInactivateSumOverRateAllocations()	{
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->InactivateSumOverRateAllocations();
			// process[gene]->sumratealloc = 0;
		}
	}
	// sumratealloc = 0;
}

void MultiGenePhyloProcess::SlaveActivateZip()	{
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->ActivateZip();
		}
	}
}

void MultiGenePhyloProcess::SlaveInactivateZip()	{
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->InactivateZip();
		}
	}
}

void MultiGenePhyloProcess::SlaveUpdateConditionalLikelihoods()	{
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->UpdateConditionalLikelihoods();
		}
	}
}

void MultiGenePhyloProcess::SlaveCollectLogLikelihood()	{
	double totlogl = 0;
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			genelnL[gene] = process[gene]->GetLogLikelihood();
			totlogl += genelnL[gene];
		}
	}
	MPI_Send(&totlogl,1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

void MultiGenePhyloProcess::SlaveCollectFullLogLikelihood()	{
	double totlogl = 0;
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			genelnL[gene] = process[gene]->GetFullLogLikelihood();
			totlogl += genelnL[gene];
		}
	}
	MPI_Send(&totlogl,1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

void MultiGenePhyloProcess::SlaveActivateCPO()  {
    int size;
	MPI_Bcast(&size,1,MPI_INT,0,MPI_COMM_WORLD);
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->ActivateCPO(size);
		}
	}
}

void MultiGenePhyloProcess::SlaveInactivateCPO()  {
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->InactivateCPO();
		}
	}
}

void MultiGenePhyloProcess::SlaveCollectPseudoMarginalLogLikelihood()	{
    double* tmpgenelogl = new double[Ngene];
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			tmpgenelogl[gene] = process[gene]->GetPseudoMarginalLogLikelihood();
		}
        else    {
            tmpgenelogl[gene] = 0;
        }
	}
	MPI_Send(tmpgenelogl,Ngene,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
    delete[] tmpgenelogl;
}


void MultiGenePhyloProcess::SlaveComputeNodeLikelihood(int fromindex,int auxindex) {
	double totlogl = 0;
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			genelnL[gene] = process[gene]->LocalComputeNodeLikelihood(fromindex,auxindex);
			totlogl += genelnL[gene];
		}
	}
	MPI_Send(&totlogl,1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

void MultiGenePhyloProcess::SlaveGetFullLogLikelihood()	{

	double totlogl = 0;
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			genelnL[gene] = process[gene]->GetFullLogLikelihood();
			totlogl += genelnL[gene];
		}
	}
	double sum[2];
	sum[0] = totlogl;
	// normally, should be the likelihood conditional on allocations
	sum[1] = totlogl;
	// sum[1] = logL;
	MPI_Send(&totlogl,2,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

void MultiGenePhyloProcess::SlaveComputeTopoBFLogLikelihoodRatio()	{

	double frac[2];
	MPI_Bcast(frac,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
	double delta = 0;
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			double tmp = process[gene]->ComputeTopoBFLogLikelihoodRatio(frac[0],frac[1]);
			delta += tmp;
		}
	}
	MPI_Send(&delta,1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

void MultiGenePhyloProcess::SlaveReset(int n,bool v) {
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->SlaveReset(n,v);
		}
	}
}

void MultiGenePhyloProcess::SlaveMultiply(int n,int m,bool v) {
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->SlaveMultiply(n,m,v);
		}
	}
}

void MultiGenePhyloProcess::SlaveMultiplyByStationaries(int n,bool v) {
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->SlaveMultiplyByStationaries(n,v);
		}
	}
}

void MultiGenePhyloProcess::SlaveInitialize(int n,int m,bool v) {
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->SlaveInitialize(n,m,v);
		}
	}
}

void MultiGenePhyloProcess::SlavePropagate(int n,int m,bool v,double t) {
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->SlavePropagate(n,m,v,t);
		}
	}
}

void MultiGenePhyloProcess::SlaveProposeMove(int n,double x) {
	PhyloProcess::SlaveProposeMove(n,x);
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->SlaveProposeMove(n,x);
		}
	}
}

void MultiGenePhyloProcess::SlaveRestoreBranch(int n) {
	PhyloProcess::SlaveRestoreBranch(n);
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->SlaveRestoreBranch(n);
		}
	}
}

void MultiGenePhyloProcess::SlaveRoot(int n) {
	PhyloProcess::SlaveRoot(n);
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->SlaveRoot(n);
		}
	}
}

void MultiGenePhyloProcess::SlaveBackupTree()	{
	GetTree()->Backup();
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->GetTree()->Backup();
		}
	}
}

void MultiGenePhyloProcess::SlaveRestoreTree()	{
	GetTree()->Restore();
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->GetTree()->Restore();
		}
	}
}

void MultiGenePhyloProcess::SlaveSwapTree()	{
	GetTree()->Swap();
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->GetTree()->Swap();
		}
	}
}

void MultiGenePhyloProcess::SlaveGibbsSPRScan(int idown, int iup)	{
	for(int i=0; i<GetNbranch(); i++) {
		loglarray[i] = 0.0;
	}
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->LocalGibbsSPRScan(idown,iup);
			for(int i=0; i<GetNbranch(); i++) {
				loglarray[i] += process[gene]->loglarray[i];
			}
		}
	}
	MPI_Send(loglarray,GetNbranch(),MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

void MultiGenePhyloProcess::SlavePropagateOverABranch(int l)	{

	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->SlavePropagateOverABranch(l);
		}
	}
}

void MultiGenePhyloProcess::LocalTryNNI(int l, int n, int* br, double* m, double* loglikelihood, int mimick)	{

	PhyloProcess::LocalTryNNI(l,n,br,m,loglikelihood,1);
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->LocalTryNNI(l,n,br,m,loglikelihood,0);
		}
	}
}

void MultiGenePhyloProcess::LocalFinalizeNNI(int n, int* br, int choice, int mimick)	{

	PhyloProcess::LocalFinalizeNNI(n,br,choice,1);
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->LocalFinalizeNNI(n,br,choice,0);
		}
	}
}

void MultiGenePhyloProcess::UpdateBranchLengthSuffStat() {

	for(int i=0; i<GetNbranch(); i++) {
		branchlengthsuffstatcount[i] = 0;
		branchlengthsuffstatbeta[i] = 0.0;
	}

	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->UpdateBranchLengthSuffStat();
			const double* count = process[gene]->GetBranchLengthSuffStatCount();
			const double* beta = process[gene]->GetBranchLengthSuffStatBeta();
			for(int i=0; i<GetNbranch(); i++) {
				branchlengthsuffstatcount[i] += count[i];
				branchlengthsuffstatbeta[i] += beta[i];
			}
		}
	}
}

void MultiGenePhyloProcess::SlaveSetMinMax()	{

	double minmax[2];
	MPI_Bcast(minmax,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
	SetMinMax(minmax[0],minmax[1]);
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->SetMinMax(minmax[0],minmax[1]);
		}
	}
}

void MultiGenePhyloProcess::GlobalReshuffleSites()	{

	MESSAGE signal = RESHUFFLE;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
}


void MultiGenePhyloProcess::SlaveReshuffleSites()	{

	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->NonMPIReshuffleSites();
		}
	}
}

void MultiGenePhyloProcess::ReadPB(int argc, char* argv[])	{

	string name = "";

	int burnin = -1;
	int every = 1;
	int until = -1;
	int ppred = 0;
	int rateprior = 0;
	int profileprior = 0;
	int rootprior = 0;
    int fulllogl = 0;

	// 1 : plain ppred (outputs simulated data)
	// 2 : diversity statistic
	// 3 : compositional statistic

	try	{

		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];
			if (s == "-div")	{
				ppred = 2;
			}
			else if (s == "-comp")	{
				ppred = 3;
			}
			else if (s == "-nsub")	{
				ppred = 4;
			}
			else if (s == "-ppred")	{
				ppred = 1;
			}
			else if (s == "-ppredrate")	{
				i++;
				string tmp = argv[i];
				if (tmp == "prior")	{
					rateprior = 1;
				}
				else if ((tmp == "posterior") || (tmp == "post"))	{
					rateprior = 0;
				}
				else	{
					cerr << "error after ppredrate: should be prior or posterior\n";
					throw(0);
				}
			}
			else if (s == "-ppredprofile")	{
				i++;
				string tmp = argv[i];
				if (tmp == "prior")	{
					profileprior = 1;
				}
				else if ((tmp == "posterior") || (tmp == "post"))	{
					profileprior = 0;
				}
				else	{
					cerr << "error after ppredprofile: should be prior or posterior\n";
					throw(0);
				}
			}
			else if (s == "-ppredroot")	{
				i++;
				string tmp = argv[i];
				if (tmp == "prior")	{
					rootprior = 1;
				}
				else if ((tmp == "posterior") || (tmp == "post"))	{
					rootprior = 0;
				}
				else	{
					cerr << "error after ppredroot: should be prior or posterior\n";
					throw(0);
				}
			}
			else if (s == "-logl")	{
				fulllogl = 1;
			}
			else if ( (s == "-x") || (s == "-extract") )	{
				i++;
				if (i == argc) throw(0);
				burnin = atoi(argv[i]);
				i++;
				if (i == argc) throw(0);
				every = atoi(argv[i]);
				i++;
				if (i == argc) throw(0);
				string tmp = argv[i];
				if (IsFloat(tmp))	{
					until = atoi(argv[i]);
				}
				else	{
					i--;
				}
			}
			else	{
				if (i != (argc -1))	{
					throw(0);
				}
				name = argv[i];
			}
			i++;
		}
	}
	catch(...)	{
		cerr << "error in command\n";
		cerr << '\n';
		MESSAGE signal = KILL;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Finalize();
		exit(1);
	}

	if (until == -1)	{
		until = GetSize();
	}
	if (burnin == -1)	{
		burnin = GetSize() / 5;
	}

	if (ppred)	{
		// PostPred(ppred,name,burnin,every,until,rateprior,profileprior,rootprior);
	}
	else if (fulllogl)	{
		ReadFullLogL(name,burnin,every,until);
	}
	else	{
		Read(name,burnin,every,until);
	}
}

void MultiGenePhyloProcess::ReadFullLogL(string name, int burnin, int every, int until)	{

	ifstream is((name + ".chain").c_str());
	if (!is)	{
		cerr << "error: no .chain file found\n";
		exit(1);
	}

    int samplesize = (int) ((until - burnin)/every);

	cerr << '\n';
	cerr << "burnin : " << burnin << "\n";
	cerr << "every  : " << every << '\n'; 
	cerr << "until  : " << until << '\n';
    cerr << "size   : " << samplesize << '\n';
	cerr << '\n';

	MESSAGE signal = ACTIVATECPO;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&samplesize,1,MPI_INT,0,MPI_COMM_WORLD);

	cerr << "burnin\n";
	for (int i=0; i<burnin; i++)	{
		cerr << ".";
		FromStream(is);
	}
    cerr << '\n';
    
    double mean = 0;
    double var = 0;

	cerr << "sample\n";
    for (int i=0; i<samplesize; i++)    {
		cerr << ".";
		FromStream(is);
		QuickUpdate();
		MPI_Status stat;
		MESSAGE signal = FULLLIKELIHOOD;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

		double total = 0;
        double tmp = 0;
		for(int i=1; i<GetNprocs(); ++i) {
			MPI_Recv(&tmp,1,MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
            total += tmp;
		}
        mean += total;
        var += total*total;
		
        for (int rep=1; rep<every; rep++)   {
			FromStream(is);
		}
	}
    cerr << '\n';
    mean /= samplesize;
    var /= samplesize;
    var -= mean*mean;

	signal = COLLECTCPO;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Status stat;
    double* genelogl = new double[Ngene];
    double* tmpgenelogl = new double[Ngene];
    for (int gene=0; gene<Ngene; gene++)	{
        genelogl[gene] = 0;
    }
	for(int proc=1; proc<GetNprocs(); proc++) {
		MPI_Recv(tmpgenelogl,Ngene,MPI_DOUBLE,proc,TAG1,MPI_COMM_WORLD,&stat);
		for (int gene=0; gene<Ngene; gene++)	{
			if (genealloc[gene] == proc)	{
				genelogl[gene] = tmpgenelogl[gene];
			}
		}
	}

    double totcposcore = 0;
    for (int gene=0; gene<Ngene; gene++)	{
        if (! genelogl[gene])   {
            cerr << "error: did not get cpo for gene : " << gene << '\n';
        }
        totcposcore += genelogl[gene];
    }

	ofstream os((name + ".logl").c_str());
	os << "posterior mean  ln L : " <<  mean << " +/- " << sqrt(var) << '\n';
	os << "pseudo marginal ln L : " <<  totcposcore << '\n';
	cerr << "posterior mean ln L : " <<  mean << " +/- " << sqrt(var) << '\n';
	cerr << "pseudo marginal ln L : " <<  totcposcore << '\n';
    cerr << '\n';

	ofstream gos((name + ".genepseudologl").c_str());
    for (int gene=0; gene<Ngene; gene++)	{
        gos << genelogl[gene] << '\n';
    }

    // clean up
    delete[] tmpgenelogl;
    delete[] genelogl;

	signal = INACTIVATECPO;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

}
