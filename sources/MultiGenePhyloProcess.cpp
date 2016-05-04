
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

#include "MultiGenePhyloProcess.h"

void MultiGenePhyloProcess::New()	{


	CreateMPI(0);
	// SetMPI(myid,nprocs);
	AllocateAlignments(datafile,treefile);
	SetProfileDim();
	SetTree(treefile);

	Create();

	if (! GetMyid())	{
		GlobalBroadcastTree();
		Sample();
		GlobalUpdateParameters();
		GlobalSample();
		GlobalUnfold();
	}
	else	{
		SlaveBroadcastTree();
	}

	if (BPP)	{
		BPP->RegisterWithTaxonSet(GetData()->GetTaxonSet());
	}
}


void MultiGenePhyloProcess::Create()	{

	if (! genelnL)	{
		PhyloProcess::Create();
		genelnL = new double[Ngene];
		tmpgenelnL = new double[Ngene];
		process = new PhyloProcess*[Ngene];
		for (int gene=0; gene<Ngene; gene++)	{
			process[gene] = 0;
		}
	}
}

void MultiGenePhyloProcess::Delete()	{

	if (genelnL)	{
		delete[] genelnL;
		delete[] tmpgenelnL;
		genelnL = 0;
		for (int gene=0; gene<Ngene; gene++)	{
			delete process[gene];
			process[gene] = 0;
		}
		delete[] process;
		PhyloProcess::Delete();
	}
}

void MultiGenePhyloProcess::AllocateAlignments(string datafile, string treefile)	{

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
		cerr << '\n';
		cerr << "number of sites allocated to each slave:\n";
		for (int i=1; i<nprocs; i++)	{
			cerr << i << '\t' << globalnsite[i] << '\n';
		}
		cerr << '\n';
		// cerr << "total: " << GetGlobalNsite() << '\n';
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

int MultiGenePhyloProcess::SpecialSlaveExecute(MESSAGE signal)	{

	switch(signal) {
	case SAMPLE:
		SlaveSample();
		return 1;
		break;
	case GENE_MOVE:
		SlaveGeneMove();
		return 1;
		break;

	default:
		return 0;
		// PhyloProcess::SlaveExecute(signal);
	}
}

// should be defined at the level of the non specialized PhyloProcess class
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
			process[gene]->Move();
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
			process[gene]->CloneTree();
			process[gene]->GetTree2()->RegisterWith(GetData()->GetTaxonSet());
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
			process[gene]->sumratealloc = 1;
		}
	}
	sumratealloc = 1;
}

void MultiGenePhyloProcess::SlaveInactivateSumOverRateAllocations()	{
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->InactivateSumOverRateAllocations(process[gene]->ratealloc);
			process[gene]->sumratealloc = 0;
		}
	}
	sumratealloc = 0;
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

	cerr << "in multigenephyloprocess::slavegetfullloglikelihood\n";
	exit(1);
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

void MultiGenePhyloProcess::SlaveRestore(int n) {
	PhyloProcess::SlaveRestore(n);
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->SlaveRestore(n);
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
	GetTree2()->Backup();
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->GetTree()->Backup();
			process[gene]->GetTree2()->Backup();
		}
	}
}

void MultiGenePhyloProcess::SlaveRestoreTree()	{
	GetTree()->Restore();
	GetTree2()->Restore();
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->GetTree()->Restore();
			process[gene]->GetTree2()->Restore();
		}
	}
}

void MultiGenePhyloProcess::SlaveSwapTree()	{
	GetTree()->Swap();
	GetTree2()->Swap();
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->GetTree()->Swap();
			process[gene]->GetTree2()->Swap();
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

void MultiGenePhyloProcess::SlaveSetMinMax()	{

	double minmax[2];
	MPI_Bcast(minmax,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->SetMinMax(minmax[0],minmax[1]);
		}
	}
}

void MultiGenePhyloProcess::SlaveDetach(int n,int m) {

	LocalDetach(n,m);
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->LocalDetach(n,m);
		}
	}
}

void MultiGenePhyloProcess::SlaveAttach(int n,int m,int p,int q) {

	LocalAttach(n,m,p,q);
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->LocalAttach(n,m,p,q);
		}
	}
}

void MultiGenePhyloProcess::SlaveDetach1(int n,int m) {

	LocalDetach1(n,m);
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->LocalDetach1(n,m);
		}
	}
}

void MultiGenePhyloProcess::SlaveAttach1(int n,int m,int p,int q) {

	LocalAttach1(n,m,p,q);
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->LocalAttach1(n,m,p,q);
		}
	}
}

void MultiGenePhyloProcess::SlaveDetach2(int n,int m) {

	LocalDetach2(n,m);
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->LocalDetach2(n,m);
		}
	}
}

void MultiGenePhyloProcess::SlaveAttach2(int n,int m,int p,int q) {

	LocalAttach2(n,m,p,q);
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->LocalAttach2(n,m,p,q);
		}
	}
}

void MultiGenePhyloProcess::SlaveSwapRoot()	{

	SwapRoot();
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->SwapRoot();
		}
	}
}


void MultiGenePhyloProcess::SlaveUpdateSiteRateSuffStat()	{

	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->UpdateSiteRateSuffStat();
		}
	}
}
