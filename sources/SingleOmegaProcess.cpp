
#include "SingleOmegaProcess.h"
#include "Parallel.h"

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
	// UpdateSiteOmegaSuffStat();
	omegasuffstatbeta = 0;
	omegasuffstatcount = 0;
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i))	{
			omegasuffstatbeta += siteomegasuffstatbeta[i];	
			omegasuffstatcount += siteomegasuffstatcount[i];	
		}
	}
}	

double SingleOmegaProcess::MoveOmega(double tuning)	{

	int nrep = 10;
	GlobalUpdateOmegaSuffStat();
	for (int rep=0; rep<nrep; rep++)	{
		SimpleMoveOmega(tuning);
		SimpleMoveOmega(0.3*tuning);
	}
}

double SingleOmegaProcess::SimpleMoveOmega(double tuning)        {

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

