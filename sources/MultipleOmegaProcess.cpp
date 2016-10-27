#include "MultipleOmegaProcess.h"
#include "Parallel.h"

double MultipleOmegaProcess::MPIMoveOmega(double tuning, int nrep)	{

	double total=0;

	for (int rep=0; rep<nrep; rep++)	{

		total += GlobalOmegaIncrementalFiniteMove(1);

		GlobalUpdateOmegaSuffStat();

		total += MoveOmegaValuesAndHyperParameters(tuning,nrep);

		ResampleOmegaWeights();

		// perhaps some broadcast here, or just call GlobalUpdateParameters
	}

	return total;
}

double MultipleOmegaProcess::NonMPIMoveOmega(double tuning, int nrep)	{

	double total=0;

	for (int rep=0; rep<nrep; rep++)	{

		// yet to be implemented
		// total += OmegaIncrementalFiniteMove();

		UpdateOmegaSuffStat();

		total += MoveOmegaValuesAndHyperParameters(tuning,nrep);

		ResampleOmegaWeights();
	}

	return total;
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
	/*
	for (int k=0; k<GetNomega(); k++)	{
		omega[k] = 1.0;
	}
	*/
	SampleOmegas();
	omegaweightalpha = 1.0;
	SampleOmegaWeights();
	SampleOmegaAlloc();
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
