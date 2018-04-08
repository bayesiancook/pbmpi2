
#include "MultiGeneRASCATGTRSBDPGammaPhyloProcess.h"
#include "Parallel.h"

void MultiGeneRASCATGTRSBDPGammaPhyloProcess::Create()	{
	
	RASCATGTRSBDPGammaPhyloProcess::Create();
	MultiGenePhyloProcess::Create();
	MultiGeneExpoConjugateGTRProfileProcess::Create();
	if (GetMyid())	{
		for (int gene=0; gene<Ngene; gene++)	{
			if (genealloc[gene] == myid)	{
				process[gene] = new RASCATGTRSBDPGammaPhyloProcess(Ncat,0,rrtype,kappaprior);
				process[gene]->SetParameters(genename[gene],treefile,partitionfile,iscodon,codetype,sis,sisfrac,sisnfrac,sisnrep,siscutoff,fixtopo,fixroot,roottax1,roottax2,topoburnin,topobf,bfburnin,bffrac,bfnfrac,bfnrep,blfactor,blfile,NSPR,NMHSPR,NTSPR,temperedbl,temperedgene,temperedrate,topolambda,topomu,toponstep,NNNI,nspec,ntspec,taxon1,taxon2,taxon3,taxon4,bpp,nbpp,ntbpp,bppnstep,bppname,bppcutoff,bppbeta,profilepriortype,dc,fixbl,sumovercomponents,proposemode,allocmode,fasttopo,fasttopofracmin,fasttoponstep,fastcondrate,dirpriortype,Nstatcomp,priorempmix,priormixtype,fixstatweight,fixstatalpha,fixstatcenter,reshuffle);
				process[gene]->SetName(name);
				process[gene]->SetMPI(0,1);
				GetProcess(gene)->SetFixAlpha(GlobalAlpha());
				GetProcess(gene)->SetFixBL(GlobalBranchLengths());
				if (! GlobalBranchLengths())	{
					GetProcess(gene)->hierarchicallengthprior = 1;
				}
				GetProcess(gene)->SetFixRR(1);
				if (nmodemax)	{
					GetProcess(gene)->SetNmodeMax(nmodemax);
				}
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
	MultiGeneExpoConjugateGTRProfileProcess::Delete();
	MultiGenePhyloProcess::Delete();
	RASCATGTRSBDPGammaPhyloProcess::Delete();
}

double MultiGeneRASCATGTRSBDPGammaPhyloProcess::Move(double tuning)	{

	chronototal.Start();
	propchrono.Start();

	if (GlobalBranchLengths())	{
		if ((topobf != 3) && ((topobf != 1) || (size < bfburnin)))	{
			BranchLengthMove(tuning);
			BranchLengthMove(0.1 * tuning);
		}
	}

	if (! fixtopo)	{
		if (! GlobalBranchLengths())	{
			cerr << "error in multigene multibl move topo\n";
			exit(1);
		}
		MoveTopo();
	}

	propchrono.Stop();

	for (int rep=0; rep<5; rep++)	{

		GlobalCollapse();

		MultiGeneBranchProcess::Move(tuning,10);
		MultiGeneBranchProcess::Move(0.1*tuning,10);

		GlobalUpdateParameters();

		GlobalGeneMove();

		GlobalUpdateParameters();

		MultiGeneRateProcess::Move(tuning,10);
		MultiGeneRateProcess::Move(0.3*tuning,10);
		MultiGeneRateProcess::Move(0.03*tuning,10);

		GlobalUpdateParameters();

		if (! fixrr){
            MoveRR();
            GlobalUpdateParameters();
			LengthRelRateMove(1,10);
			LengthRelRateMove(0.1,10);
			LengthRelRateMove(0.01,10);
            GlobalUpdateParameters();
		}

		GlobalUnfold();
	}

	chronototal.Stop();

	return 1;
}

void MultiGeneExpoConjugateGTRProfileProcess::UpdateRRSuffStat() {

	for(int i=0; i<GetNrr(); i++) {
		rrsuffstatcount[i] = 0;
		rrsuffstatbeta[i] = 0.0;
	}

	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			GetGTRProfileProcess(gene)->UpdateRRSuffStat();
			const double* count = GetGTRProfileProcess(gene)->GetRRSuffStatCount();
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
		int nd = 3;
		if (! GlobalBranchLengths())	{
			nd += 3*nbranch;
		}
		else	{
			nd += nbranch;
		}
		nd += nrr;
		if (kappaprior == 2)	{
			nd += 2;
		}
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
		if (! GlobalBranchLengths())	{
			for(int i=0; i<nbranch; ++i) {
				dvector[index] = branchmean[i];
				index++;
			}
			for(int i=0; i<nbranch; ++i) {
				dvector[index] = branchrelvar[i];
				index++;
			}
		}
		
		for(int i=0; i<nrr ; ++i) {
			dvector[index] = rr[i];
			index++;
		}

		if (kappaprior == 2)	{
			dvector[index] = kappamean;
			index++;
			dvector[index] = kapparelvar;
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
	int nd = 3;
	if (! GlobalBranchLengths())	{
		nd += 3*nbranch;
	}
	else	{
		nd += nbranch;
	}
	nd += nrr;
	if (kappaprior == 2)	{
		nd += 2;
	}
	double dvector[nd]; 

	MPI_Bcast(dvector,nd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	int index = 0;
	// SetAlpha(dvector[index]);
	SetRateParams(dvector[index],0);
	index++;
	meanalpha = dvector[index];
	index++;
	varalpha = dvector[index];
	index++;
	for(int i=0; i<nbranch; ++i) {
		blarray[i] = dvector[index];
		index++;
	}
	if (! GlobalBranchLengths())	{
		for(int i=0; i<nbranch; ++i) {
			branchmean[i] = dvector[index];
			index++;
		}
		for(int i=0; i<nbranch; ++i) {
			branchrelvar[i] = dvector[index];
			index++;
		}
	}

	for(int i=0; i<nrr; ++i) {
		rr[i] = dvector[index];
		index++;
	}

	if (kappaprior == 2)	{
		kappamean = dvector[index];
		index++;
		kapparelvar = dvector[index];
		index++;
	}

	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			if (kappaprior == 2)	{
				GetProcess(gene)->SetKappaHyperParams(kappamean,kapparelvar);
			}
			if (GlobalAlpha())	{
				GetProcess(gene)->SetRateParams(GetAlpha(),0);
				// GetProcess(gene)->SetAlpha(GetAlpha());
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
			case COLLECTKAPPAS:
				SlaveCollectKappas();
				break;
			default:
			RASCATGTRSBDPGammaPhyloProcess::SlaveExecute(signal);
		}
	}
}

