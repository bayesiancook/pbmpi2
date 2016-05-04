
#include "MultiGeneRASCATGTRSBDPGammaPhyloProcess.h"
#include "Parallel.h"

void MultiGeneRASCATGTRSBDPGammaPhyloProcess::Create()	{
	
	/*
	ExpoConjugateGTRProfileProcess::Create();
	RASCATGTRSBDPSubstitutionProcess::Create();
	GammaBranchProcess::Create();
	*/
	RASCATGTRSBDPGammaPhyloProcess::Create();
	MultiGenePhyloProcess::Create();
	if (GetMyid())	{
		for (int gene=0; gene<Ngene; gene++)	{
			if (genealloc[gene] == myid)	{
				process[gene] = new RASCATGTRSBDPGammaPhyloProcess(Ncat,rrtype,kappaprior);
				process[gene]->SetParameters(genename[gene],treefile,iscodon,codetype,fixtopo,NSPR,NMHSPR,NTSPR,topolambda,topomu,toponstep,NNNI,nspec,ntspec,bpp,nbpp,ntbpp,bppnstep,bppname,bppcutoff,bppbeta,profilepriortype,dc,fixbl,sumovercomponents,proposemode,allocmode,sumratealloc,fasttopo,fasttopofracmin,fasttoponstep,fastcondrate);
				process[gene]->SetName(name);
				process[gene]->SetMPI(0,1);
				process[gene]->New();
				GetProcess(gene)->SetFixAlpha(GlobalAlpha());
				GetProcess(gene)->SetFixBL(GlobalBranchLengths());
			}
		}
	}
}
	
void MultiGeneRASCATGTRSBDPGammaPhyloProcess::Delete()	{
	if (GetMyid())	{
		for (int gene=0; gene<Ngene; gene++)	{
			if (genealloc[gene] == myid)	{
				delete process[gene];
				process[gene] = 0;
			}
		}
	}
	MultiGenePhyloProcess::Delete();
	RASCATGTRSBDPGammaPhyloProcess::Delete();
	/*
	GammaBranchProcess::Delete();
	RASCATGTRSBDPSubstitutionProcess::Delete();
	ExpoConjugateGTRProfileProcess::Delete();
	*/
}

double MultiGeneRASCATGTRSBDPGammaPhyloProcess::Move(double tuning)	{

	chronototal.Start();
	propchrono.Start();

	BranchLengthMove(tuning);
	BranchLengthMove(0.1 * tuning);

	if (! fixtopo)	{
		MoveTopo();
	}

	propchrono.Stop();

	for (int rep=0; rep<5; rep++)	{

		GlobalCollapse();

		GlobalGeneMove();

		GlobalUpdateParameters();

		GammaBranchProcess::Move(tuning,10);

		GlobalUpdateParameters();

		MultiGeneRateProcess::Move(0.3*tuning,50);
		MultiGeneRateProcess::Move(0.03*tuning,50);

		GlobalUpdateParameters();

		if (! fixrr){
			LengthRelRateMove(1,10);
			LengthRelRateMove(0.1,10);
			LengthRelRateMove(0.01,10);
		}

		GlobalUpdateParameters();

		GlobalUnfold();
	}

	chronototal.Stop();

	return 1;
}

double MultiGeneRASCATGTRSBDPGammaPhyloProcess::GlobalGetMeanNcomponent()	{

	MESSAGE signal = MEANNCOMPONENT;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	int ncomp = 0;
	MPI_Status stat;
	for(int i=1; i<GetNprocs(); i++) {
		int tmp;
		MPI_Recv(&tmp,1,MPI_INT,i,TAG1,MPI_COMM_WORLD,&stat);
		ncomp += tmp;
	}
	return ((double) ncomp) / Ngene;
}

double MultiGeneRASCATGTRSBDPGammaPhyloProcess::GlobalGetMeanStatEnt()	{

	MESSAGE signal = MEANSTATENT;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	double tot = 0;
	MPI_Status stat;
	for(int i=1; i<GetNprocs(); i++) {
		double tmp;
		MPI_Recv(&tmp,1,MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
		tot += tmp;
	}
	return tot / Ngene;
}

double MultiGeneRASCATGTRSBDPGammaPhyloProcess::GlobalGetMeanStatAlpha()	{

	MESSAGE signal = MEANSTATALPHA;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	double tot = 0;
	MPI_Status stat;
	for(int i=1; i<GetNprocs(); i++) {
		double tmp;
		MPI_Recv(&tmp,1,MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
		tot += tmp;
	}
	return tot / Ngene;
}

double MultiGeneRASCATGTRSBDPGammaPhyloProcess::GlobalGetMeanKappa()	{

	MESSAGE signal = MEANKAPPA;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	double tot = 0;
	MPI_Status stat;
	for(int i=1; i<GetNprocs(); i++) {
		double tmp;
		MPI_Recv(&tmp,1,MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
		tot += tmp;
	}
	return tot / Ngene;
}

void MultiGeneRASCATGTRSBDPGammaPhyloProcess::SlaveGetMeanNcomponent() {
	int ncomp = 0;
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			ncomp += GetProcess(gene)->GetNDisplayedComponent();
		}
	}
	MPI_Send(&ncomp,1,MPI_INT,0,TAG1,MPI_COMM_WORLD);
}

void MultiGeneRASCATGTRSBDPGammaPhyloProcess::SlaveGetMeanStatEnt() {
	double tot = 0;
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			tot += GetProcess(gene)->GetStatEnt();
		}
	}
	MPI_Send(&tot,1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

void MultiGeneRASCATGTRSBDPGammaPhyloProcess::SlaveGetMeanStatAlpha() {
	double tot = 0;
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			tot += GetProcess(gene)->GetMeanDirWeight();
		}
	}
	MPI_Send(&tot,1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

void MultiGeneRASCATGTRSBDPGammaPhyloProcess::SlaveGetMeanKappa() {
	double tot = 0;
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			tot += GetProcess(gene)->GetKappa();
		}
	}
	MPI_Send(&tot,1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

void MultiGeneExpoConjugateGTRProfileProcess::UpdateRRSuffStat() {

	for(int i=0; i<GetNrr(); i++) {
		rrsuffstatcount[i] = 0;
		rrsuffstatbeta[i] = 0.0;
	}

	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			GetGTRProfileProcess(gene)->UpdateRRSuffStat();
			const int* count = GetGTRProfileProcess(gene)->GetRRSuffStatCount();
			const double* beta = GetGTRProfileProcess(gene)->GetRRSuffStatBeta();
			for(int i=0; i<GetNrr(); i++) {
				rrsuffstatcount[i] += count[i];
				rrsuffstatbeta[i] += beta[i];
			}
		}
	}
}

void MultiGeneRASCATGTRSBDPGammaPhyloProcess::GlobalUpdateParameters() {

	if (GetNprocs() > 1)	{

		int nbranch = GetNbranch();
		int nrr = GetNrr();
		int nd = 3 + 3*nbranch + nrr;
		double dvector[nd]; 
		MESSAGE signal = PARAMETER_DIFFUSION;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

		// GlobalBroadcastTree();

		// First we assemble the vector of doubles for distribution
		int index = 0;
		dvector[index] = GetAlpha();
		index++;
		dvector[index] = meanalpha;
		index++;
		dvector[index] = varalpha;
		index++;
		
		for(int i=0; i<nbranch; ++i) {
			dvector[index] = blarray[i];
			index++;
		}
		for(int i=0; i<nbranch; ++i) {
			dvector[index] = branchmean[i];
			index++;
		}
		for(int i=0; i<nbranch; ++i) {
			dvector[index] = branchrelvar[i];
			index++;
		}
		
		for(int i=0; i<nrr ; ++i) {
			dvector[index] = rr[i];
			index++;
		}

		// Now send out the doubles and ints over the wire...
		MPI_Bcast(dvector,nd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	}
	else	{

	}
}

void MultiGeneRASCATGTRSBDPGammaPhyloProcess::SlaveUpdateParameters() {

	int nbranch = GetNbranch();
	int nrr = GetNrr();
	int nd = 3 + 3*nbranch + nrr;
	double dvector[nd]; 

	MPI_Bcast(dvector,nd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	int index = 0;
	SetAlpha(dvector[index]);
	index++;
	meanalpha = dvector[index];
	index++;
	varalpha = dvector[index];
	index++;
	for(int i=0; i<nbranch; ++i) {
		blarray[i] = dvector[index];
		index++;
	}
	for(int i=0; i<nbranch; ++i) {
		branchmean[i] = dvector[index];
		index++;
	}
	for(int i=0; i<nbranch; ++i) {
		branchrelvar[i] = dvector[index];
		index++;
	}

	for(int i=0; i<nrr; ++i) {
		rr[i] = dvector[index];
		index++;
	}

	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			if (GlobalAlpha())	{
				GetProcess(gene)->SetAlpha(GetAlpha());
			}
			else	{
				GetProcess(gene)->meanalpha = meanalpha;
				GetProcess(gene)->varalpha = varalpha;
			}
			if (GlobalBranchLengths())	{
				GetProcess(gene)->SetBranchLengths(GetBranchLengths());
			}
			else	{
				GetProcess(gene)->SetHyperParameters(branchmean,branchrelvar);
			}
			GetProcess(gene)->SetRR(GetRR());
			GetProcess(gene)->UpdateMatrices();
		}
	}
}


void MultiGeneRASCATGTRSBDPGammaPhyloProcess::SlaveExecute(MESSAGE signal)	{

	int ret = SpecialSlaveExecute(signal);
	if (! ret)	{
		switch(signal) {

			case MEANNCOMPONENT:
				SlaveGetMeanNcomponent();
				break;
			case MEANKAPPA:
				SlaveGetMeanKappa();
				break;
			case MEANSTATENT:
				SlaveGetMeanStatEnt();
				break;
			case MEANSTATALPHA:
				SlaveGetMeanStatAlpha();
				break;
			default:
			RASCATGTRSBDPGammaPhyloProcess::SlaveExecute(signal);
		}
	}
}

