
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#include "Parallel.h"
#include "OmegaProcess.h"
#include "GeneralPathSuffStatMatrixProfileProcess.h"
#include "Random.h"

double SingleOmegaProcess::LogOmegaPrior()        {

	if (omegaprior == 0)	{
		// ratio of exponential random variables
		return -2 * log(1 + *omega);
	}
	else	{	
		// jeffreys
		return -log(*omega);
	}
}

void SingleOmegaProcess::SampleOmega()        {

	*omega = 1.0;
}

void SingleOmegaProcess::UpdateOmegaSuffStat()	{
	UpdateSiteOmegaSuffStat();
	omegasuffstatbeta = 0;
	omegasuffstatcount = 0;
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i))	{
			omegasuffstatbeta += siteomegasuffstatbeta[i];	
			omegasuffstatcount += siteomegasuffstatcount[i];	
		}
	}
}	

double SingleOmegaProcess::MoveOmega(double tuning)        {

	int naccepted = 0;
	double bkomega = *omega;
	double deltalogprob = -LogOmegaPrior() - OmegaSuffStatLogProb();

	double h = tuning * (rnd::GetRandom().Uniform() -0.5);
	double e = exp(h);
	*omega *= e;

	UpdateOmega();
	deltalogprob += h;
	deltalogprob += LogOmegaPrior() + OmegaSuffStatLogProb();

	int accepted = (rnd::GetRandom().Uniform() < exp(deltalogprob));
	if (accepted)	{
		naccepted++;	
	}
	else	{
		*omega = bkomega;
		UpdateOmega();
	}
	return naccepted;	
}

void SingleOmegaProcess::GlobalUpdateOmegaSuffStat()	{
	if (GetNprocs() > 1)	{
		// MPI2
		// should ask the slaves to call their UpdateRateSuffStat
		// and then gather the statistics;
		MPI_Status stat;
		MESSAGE signal = UPDATE_OMEGA;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

		omegasuffstatcount = 0;
		omegasuffstatbeta = 0;	

		int ivalue;
		double dvalue;
		for(int i=1; i<GetNprocs(); i++) {
                	MPI_Recv(&ivalue,1,MPI_INT,MPI_ANY_SOURCE,TAG1,MPI_COMM_WORLD,&stat);
			omegasuffstatcount += ivalue;
		}
       		MPI_Barrier(MPI_COMM_WORLD);
		for(int i=1; i<GetNprocs(); i++) {
			MPI_Recv(&dvalue,1,MPI_DOUBLE,MPI_ANY_SOURCE,TAG1,MPI_COMM_WORLD,&stat);
			omegasuffstatbeta += dvalue;
       		}
	}
	else	{
		UpdateOmegaSuffStat();	
	}
}

void SingleOmegaProcess::SlaveUpdateOmegaSuffStat()	{
	UpdateOmegaSuffStat();
	MPI_Send(&omegasuffstatcount,1,MPI_INT,0,TAG1,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Send(&omegasuffstatbeta,1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}
			

double MultipleOmegaProcess::GlobalOmegaIncrementalFiniteMove(int nrep=1)	{

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

double MultipleOmegaProcess::SlaveOmegaIncrementalFiniteMove()	{
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

void MultipleOmegaProcess::ResampleOmegaWeights()	{
	int omegaoccupancy[Nomega];
	for (int l=0; l<Nomega; l++)	{
		omegaoccupancy[l]=0;
		for (int i=0; i<GetNsite(); i++)	{
			if (omegaalloc[i] == l)	{
				omegaoccupancy[l]++;
			}
		}
	}
	double total = 0;
	for (int k=0; k<Nomega; k++)	{
		omegaweight[k] = rnd::GetRandom().sGamma(omegaweightalpha + omegaoccupancy[k]);
		if (omegaweight[k] < 1e-10)	{
			omegaweight[k] = 1e-10;
		}
		total += omegaweight[k];
	}
	for (int k=0; k<Nomega; k++)	{
		omegaweight[k] /= total;
	}
}

void MultipleOmegaProcess::SampleOmega()	{
	for (int i=0; i<GetNomega(); i++)	{
		omega[i] = 1.0;
	}
}
	
void MultipleOmegaProcess::SampleOmegaWeights()	{
	double total = 0;
	for (int k=0; k<GetNomega(); k++)	{
		omegaweight[k] = rnd::GetRandom().sGamma(omegaweightalpha);
		total += omegaweight[k];
	}
	for (int k=0; k<GetNomega(); k++)	{
		omegaweight[k] /= total;
	}

}

double MultipleOmegaProcess::LogOmegaPrior()        {
	double total=0;
	if (omegaprior == 0)	{
		for (int i=0; i<GetNomega(); i++)	{
			// ratio of exponential random variables
			total += -2 * log(1 + omega[i]);
		}
	}
	else	{	
		for (int i=0; i<GetNomega(); i++)	{
			// jeffreys
			total += -log(omega[i]);
		}
	}
	return total;
}

double MultipleOmegaProcess::MoveOmega(double tuning)        {
	double total=0;
	total += MoveOmegaValues(tuning);
	total += GlobalOmegaIncrementalFiniteMove();
	return total;
}

double MultipleOmegaProcess::MoveOmegaValues(double tuning)        {

	int naccepted = 0;
	for (int l=0; l<GetNomega(); l++)	{
		int naccepted = 0;
		double bkomega = omega[l];
		double deltalogprob = -LogOmegaPrior() - OmegaSuffStatLogProb();

		double h = tuning * (rnd::GetRandom().Uniform() -0.5);
		double e = exp(h);
		omega[l] *= e;

		UpdateOmega();
		deltalogprob += h;
		deltalogprob += LogOmegaPrior() + OmegaSuffStatLogProb();

		int accepted = (rnd::GetRandom().Uniform() < exp(deltalogprob));
		if (accepted)	{
			naccepted++;	
		}
		else	{
			omega[l] = bkomega;
			UpdateOmega();
		}
	}
	return naccepted;	
}
	
void MultipleOmegaProcess::UpdateOmegaSuffStat()	{
	UpdateSiteOmegaSuffStat();
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

void MultipleOmegaProcess::GlobalUpdateOmegaSuffStat()	{
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

void MultipleOmegaProcess::SlaveUpdateOmegaSuffStat()	{
	UpdateOmegaSuffStat();
	MPI_Send(compomegasuffstatcount,Nomega,MPI_INT,0,TAG1,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Send(compomegasuffstatbeta,Nomega,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

