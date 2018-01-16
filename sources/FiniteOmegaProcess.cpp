#include "FiniteOmegaProcess.h"
#include "BiologicalSequences.h"
#include "Parallel.h"


double FiniteOmegaProcess::MPIMoveOmega(double tuning, int nrep)	{

	double total=0;

	for (int rep=0; rep<nrep; rep++)	{

		total += GlobalOmegaIncrementalFiniteMove(1);

		GlobalUpdateOmegaSuffStat();

		if (!fixomega)	{
			total += MoveOmegaValuesAndHyperParameters(tuning,nrep);
		}

		ResampleOmegaWeights();

		// perhaps some broadcast here, or just call GlobalUpdateParameters
	}

	return total;
}

double FiniteOmegaProcess::NonMPIMoveOmega(double tuning, int nrep)	{

	double total=0;

	for (int rep=0; rep<nrep; rep++)	{

		// yet to be implemented
		// total += OmegaIncrementalFiniteMove();

		UpdateOmegaSuffStat();

		if (!fixomega)	{
			total += MoveOmegaValuesAndHyperParameters(tuning,nrep);
		}

		ResampleOmegaWeights();
	}

	return total;
}

void FiniteOmegaProcess::ResampleOmegaWeights()	{
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

void FiniteOmegaProcess::SampleOmega()	{
	/*
	for (int k=0; k<GetNomega(); k++)	{
		omega[k] = 1.0;
	}
	*/
	SampleOmegaHyper();
	omegaweightalpha = 1.0;
	if (empomegamix)	{
		ReadOmegaFix(omegamixtype);
		SetOmegaFix();
	}
	else	{
		SampleOmegas();
		SampleOmegaWeights();
	}
	SampleOmegaAlloc();
}
	
void FiniteOmegaProcess::SampleOmegaWeights()	{
	double total = 0;
	for (int k=0; k<GetNomega(); k++)	{
		omegaweight[k] = rnd::GetRandom().sGamma(omegaweightalpha);
		total += omegaweight[k];
	}
	for (int k=0; k<GetNomega(); k++)	{
		omegaweight[k] /= total;
	}

}


void FiniteOmegaProcess::ReadOmegaFix(string filename)	{
	if ((filename == "o10") || (filename == "O10"))	{
		Nomegafixcomp = o10N;
		omegafix = new double[Nomegafixcomp];
		empomegaweight = new double[Nomegafixcomp];
		for (int i=0; i<Nomegafixcomp; i++)	{
			omegafix[i] = o10OmegaFix[i];
			empomegaweight[i] = o10OmegaWeight[i];
		}
	}
	else	{
		cerr << "nothing else yet implemented\n";
		exit(1);
	}
}

void FiniteOmegaProcess::SetOmegaFix()	{

	if (! Nomegafixcomp)	{
		cerr << "error in set omega fix\n";
		exit(1);
	}
	Nomega = Nomegafixcomp;
	for (int k=0; k<Nomega; k++)	{
		omegaweight[k] = empomegaweight[k];
		omega[k] = omegafix[k];
	}
	fixnomegacomp = true;
}
