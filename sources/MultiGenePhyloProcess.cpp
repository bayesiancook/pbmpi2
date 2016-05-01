
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

void MultiGeneMPIModule::AllocateAlignments(string datafile, string treefile)	{

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
		SequenceAlignment* data = new FileSequenceAlignment(genename[gene],0,myid);
		int nstate = data->GetNstate();
		if (! gene)	{
			Nstate = nstate;
			statespace = data->GetStateSpace();
		}
		else	{
			if (Nstate != nstate)	{
				cerr << "error: all data files do not have the same alphabet\n";
				exit(1);
			}
		}

		genesize[gene] = data->GetNsite();
		geneweight[gene] = data->GetNsite() * data->GetNtaxa();
		delete data;
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
		for (int i=1; i<nprocs; i++)	{
			cerr << i << '\t' << globalnsite[i] << '\n';
		}
		cerr << '\n';
		cerr << "total: " << GetGlobalNsite() << '\n';
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

void MultiGenePhyloProcess::Create()	{

	if (! genelnL)	{
		genelnL = new double[Ngene];
		tmpgenelnL = new double[Ngene];
	}
}

void MultiGenePhyloProcess::Delete()	{

	if (genelnL)	{
		delete[] genelnL;
		delete[] tmpgenelnL;
		genelnL = 0;
	}
}

/*
double MultiGenePhyloProcess::GetLogLikelihood()	{
	GlobalCollectGeneLikelihoods();
	lnL = 0;
	for (int gene=0; gene<Ngene; gene++)	{
		lnL += genelnL[gene];
	}
	return lnL;
}

void MultiGenePhyloProcess::GlobalCollectGeneLikelihoods()	{
	// send signal
	assert(myid == 0);
	MESSAGE signal = LIKELIHOOD;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	MPI_Status stat;
	for (int j=1; j<nprocs; j++)	{
		MPI_Recv(tmpgenelnL,Ngene,MPI_DOUBLE,j,TAG1,MPI_COMM_WORLD,&stat);
		for (int gene=0; gene<Ngene; gene++)	{
			if (genealloc[gene] == j)	{
				genelnL[gene] = tmpgenelnL[gene];
			}
		}
	}
}

void MultiGenePhyloProcess::SlaveSendGeneLikelihoods()	{
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			genelnL[gene] = process[gene]->GetLogLikelihood();
		}
	}
	MPI_Send(genelnL,Ngene,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}
*/

void MultiGenePhyloProcess::SlaveLikelihood(int fromindex,int auxindex) {
	double totlogl = 0;
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			genelnL[gene] = process[gene]->LocalLikelihood(fromindex,auxindex);
			totlogl += genelnL[gene];
		}
	}
	MPI_Send(&totlogl,1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

void MultiGenePhyloProcess::GlobalGeneMove()	{
	assert(myid == 0);
	MESSAGE signal = GENE_MOVE;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
}

// should be defined at the level of the non specialized PhyloProcess class
void MultiGenePhyloProcess::GlobalSample()	{
	assert(myid == 0);
	MESSAGE signal = SAMPLE;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
}

void MultiGenePhyloProcess::SlaveGeneMove()	{
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->Move();
		}
	}
}

void MultiGenePhyloProcess::SlaveSample()	{
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->Sample();
		}
	}
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


void MultiGenePhyloProcess::SlaveRoot(int n) {
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->SlaveRoot(n);
		}
	}
}

void MultiGenePhyloProcess::SlavePropose(int n,double x) {
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->SlavePropose(n,x);
		}
	}
}

void MultiGenePhyloProcess::SlaveRestore(int n) {
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->SlaveRestore(n);
		}
	}
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

void MultiGenePhyloProcess::SlaveSMultiply(int n,bool v) {
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->SlaveSMultiply(n,v);
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

void MultiGenePhyloProcess::SlaveExecute(MESSAGE signal)	{

	switch(signal) {
	case LIKELIHOOD:
		SlaveSendGeneLikelihoods();
		break;
	case SAMPLE:
		SlaveSample();
		break;
	case UNFOLD:
		SlaveUnfold();
		break;
	case COLLAPSE:
		SlaveCollapse();
		break;
	case PARAMETER_DIFFUSION:
		SlaveUpdateParameters();
		break;
	case GENE_MOVE:
		SlaveGeneMove();
		break;

	default:
		PhyloProcess::SlaveExecute(signal);
	}
}
