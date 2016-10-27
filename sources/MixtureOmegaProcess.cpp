
#include "MixtureOmegaProcess.h"
#include "Parallel.h"

double MixtureOmegaProcess::MoveOmega(double tuning)	{

	int nrep = 3;
	double ret = 0;
	if (GetNprocs() > 1)	{
		ret = MPIMoveOmega(tuning,nrep);
	}
	else	{
		ret = NonMPIMoveOmega(tuning,nrep);
	}
	return ret;
}

double MixtureOmegaProcess::GlobalOmegaIncrementalFiniteMove(int nrep=1)	{

	// send command and arguments
	MESSAGE signal = REALLOCOMEGA_MOVE;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&nrep,1,MPI_INT,0,MPI_COMM_WORLD);

	int NAccepted = 0;
	for (int rep=0; rep<nrep; rep++)	{
		// resample component weights based on current site allocations
		// and send them to slaves
		//UpdateOccupancyNumbers();
		ResampleOmegaWeights();
		MPI_Bcast(omegaweight,Nomega,MPI_DOUBLE,0,MPI_COMM_WORLD);
		// receive new site allocations from slave
		MPI_Status stat;
		int tmpalloc[GetNsite()];
		for(int i=1; i<GetNprocs(); ++i) {
			MPI_Recv(tmpalloc,GetNsite(),MPI_INT,i,TAG1,MPI_COMM_WORLD,&stat);
			for(int j=GetProcSiteMin(i); j<GetProcSiteMax(i); ++j) {
				if (ActiveSite(j))	{
					omegaalloc[j] = tmpalloc[j];
					if ((omegaalloc[j] < 0) || (omegaalloc[j] >= Nomega))	{
						cerr << "alloc overflow\n";
						exit(1);
					}
				}
			}
		}
	}

	// final cleanup
	//UpdateOccupancyNumbers(); 
	// UpdateModeProfileSuffStat();

	// CHECK that: might be useful depending on the exact submodel
	// ResampleWeights();
	// UpdateMatrices();

	return 1;
}

double MixtureOmegaProcess::SlaveOmegaIncrementalFiniteMove()	{
	int nrep;
	MPI_Bcast(&nrep,1,MPI_INT,0,MPI_COMM_WORLD);

	int NAccepted = 0;
	int Ntried = 0;

	double* bigarray = new double[Nomega * GetNsite()];
	double* bigcumul = new double[Nomega * GetNsite()];

	for (int rep=0; rep<nrep; rep++)	{

		// receive weights sent by master
		MPI_Bcast(omegaweight,Nomega,MPI_DOUBLE,0,MPI_COMM_WORLD);
		// do the incremental reallocation move on my site range
		for (int site=GetSiteMin(); site<GetSiteMax(); site++)	{
			if (ActiveSite(site))	{
				double* mLogSamplingArray = bigarray + site * Nomega;
				double* cumul = bigcumul + site * Nomega;

				int bk = omegaalloc[site];

				double max = 0;
				for (int mode = 0; mode < Nomega; mode++)	{
					mLogSamplingArray[mode] =  OmegaLogStatProb(site,mode);
					if ((!mode) || (max < mLogSamplingArray[mode]))	{
						max = mLogSamplingArray[mode];
					}
				}

				double total = 0;

				for (int mode = 0; mode < Nomega; mode++)	{
					double p = omegaweight[mode] * exp(mLogSamplingArray[mode] - max);
					total += p;
					cumul[mode] = total;
				}

				double q = total * rnd::GetRandom().Uniform();
				int mode = 0;
				while ( (mode<Nomega) && (q > cumul[mode])) mode++;
				if (mode == Nomega)	{
					cerr << "error in switch mode: gibbs overflow\n";
					exit(1);
				}

				int Accepted = (mode != bk);
				if (Accepted)	{
					NAccepted ++;
				}
				Ntried ++;
				omegaalloc[site] = mode;
			}
		}

		// send new allocations to master
		MPI_Send(omegaalloc,GetNsite(),MPI_INT,0,TAG1,MPI_COMM_WORLD);
	}
	
	delete[] bigarray;
	delete[] bigcumul;
	return ((double) NAccepted) / Ntried;

}

double MixtureOmegaProcess::MoveCompOmega(double tuning, int component)	{
	
	double bkomega = omega[component];

	double deltalogprob = -LogOmegaPrior(component) - OmegaSuffStatLogProb(component);

	double h = tuning * (rnd::GetRandom().Uniform() -0.5);
	double e = exp(h);
	omega[component] *= e;

	deltalogprob += h;
	deltalogprob += LogOmegaPrior(component) + OmegaSuffStatLogProb(component);

	int accepted = (rnd::GetRandom().Uniform() < exp(deltalogprob));
	if (! accepted)	{
		omega[component] = bkomega;
	}
	return accepted;	
}

double MixtureOmegaProcess::MoveOmegas(double tuning, int nrep)	{

	double nacc = 0;
	double ntot = 0;
	for (int rep=0; rep<nrep; rep++)	{
		for (int k=0; k<GetNomega(); k++)	{
			nacc += MoveCompOmega(tuning,k);
			ntot++;
		}
	}
	return nacc/ntot;
}

double MixtureOmegaProcess::MoveOmegaHyper(double tuning, int nrep)	{

	double nacc = 0;
	double ntot = 0;
	for (int rep=0; rep<nrep; rep++)	{
		nacc += MoveOmegaAlpha(tuning);
		ntot++;
		nacc += MoveOmegaBeta(tuning);
		ntot++;
	}
	return nacc/ntot;
}

double MixtureOmegaProcess::MoveOmegaAlpha(double tuning)	{

	double bkomegaalpha = omegaalpha;
	double deltalogprob = - LogOmegaPrior() - LogOmegaHyperPrior();
	double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
	double e = exp(m);
	omegaalpha *= e;
	deltalogprob += LogOmegaPrior() + LogOmegaHyperPrior();
	deltalogprob += m;
	int accepted = (log(rnd::GetRandom().Uniform()) < deltalogprob);

	if (!accepted)	{
		omegaalpha = bkomegaalpha;
	}
	return accepted;
}

double MixtureOmegaProcess::MoveOmegaBeta(double tuning)	{

	double bkomegabeta = omegabeta;
	double deltalogprob = - LogOmegaPrior() - LogOmegaHyperPrior();
	double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
	double e = exp(m);
	omegabeta *= e;
	deltalogprob += LogOmegaPrior() + LogOmegaHyperPrior();
	deltalogprob += m;
	int accepted = (log(rnd::GetRandom().Uniform()) < deltalogprob);

	if (!accepted)	{
		omegabeta = bkomegabeta;
	}
	return accepted;
}

double MixtureOmegaProcess::LogOmegaPrior()        {

	double total=0;
	if (omegaprior == 0)	{ // ratio of exponential random variables
		for (int k=0; k<GetNomega(); k++)	{
			total += -2 * log(1 + omega[k]);
		}
	}
	else if (omegaprior == 1)	{ // gamma
		for (int k=0; k<GetNomega(); k++)	{		
			total += omegaalpha*log(omegabeta) - rnd::GetRandom().logGamma(omegaalpha) + (omegaalpha-1)*log(omega[k]) - omegabeta*omega[k];
		}
	}
	else	{	
		for (int k=0; k<GetNomega(); k++)	{
			// jeffreys
			total += -log(omega[k]);
		}
	}
	return total;
}

void MixtureOmegaProcess::SampleOmegaAlloc()	{
	//for (int i=0; i<GetNsite(); i++)	{
	//	omegaalloc[i] = 0;
	//}
	for (int i=0; i<GetNsite(); i++)	{
		if (ActiveSite(i))	{
			int choose = rnd::GetRandom().FiniteDiscrete(GetNomega(),omegaweight);
			//AddSite(i,choose);
			omegaalloc[i] = choose;
		}
	}
}
	
void MixtureOmegaProcess::UpdateOmegaSuffStat()	{
	// UpdateSiteOmegaSuffStat();
	for (int l=0; l<GetNomega(); l++)	{
		compomegasuffstatbeta[l] = 0;
		compomegasuffstatcount[l] = 0;
	}
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i))	{
			compomegasuffstatbeta[GetOmegaSiteAlloc(i)] += siteomegasuffstatbeta[i];	
			compomegasuffstatcount[GetOmegaSiteAlloc(i)] += siteomegasuffstatcount[i];	
		}
	}
}	

void MixtureOmegaProcess::GlobalUpdateOmegaSuffStat()	{
	if (GetNprocs() > 1)	{
		// MPI2
		// should ask the slaves to call their UpdateRateSuffStat
		// and then gather the statistics;
		MPI_Status stat;
		MESSAGE signal = UPDATE_OMEGA;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

		for (int l=0; l<Nomega; l++)	{
			compomegasuffstatcount[l] = 0;
			compomegasuffstatbeta[l] = 0;
		}

		int ivector[Nomega];
		double dvector[Nomega];
		for(int i=1; i<GetNprocs(); i++) {
                	MPI_Recv(ivector,Nomega,MPI_INT,MPI_ANY_SOURCE,TAG1,MPI_COMM_WORLD,&stat);
			for (int l=0; l<Nomega; l++)	{
				compomegasuffstatcount[l] += ivector[l];
			}
		}
       		MPI_Barrier(MPI_COMM_WORLD);
		for(int i=1; i<GetNprocs(); i++) {
			MPI_Recv(dvector,Nomega,MPI_DOUBLE,MPI_ANY_SOURCE,TAG1,MPI_COMM_WORLD,&stat);
			for (int l=0; l<Nomega; l++)	{
				compomegasuffstatbeta[l] += dvector[l];
			}
       		}
	}
	else	{
		UpdateOmegaSuffStat();	
	}
}

void MixtureOmegaProcess::SlaveUpdateOmegaSuffStat()	{
	UpdateOmegaSuffStat();
	MPI_Send(compomegasuffstatcount,Nomega,MPI_INT,0,TAG1,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Send(compomegasuffstatbeta,Nomega,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

