
#include "MVNProfileProcess.h"
#include "Random.h"

void MVNProfileProcess::Create()	{

	if (! kappa)	{
		kappa = new double[GetDim()];
		for (int j=0; j<GetDim(); j++)	{
			kappa[j] = 1;
		}
		mean = new double[GetDim()];
		df = GetDim();
		covmatrix = new InverseWishartMatrix(kappa,GetDim(),df);
	}
}

void MVNProfileProcess::Delete()	{

	if (kappa)	{
		delete[] kappa;
		kappa = 0;
		delete[] mean;
		mean = 0;
		delete covmatrix;
		covmatrix = 0;
	}
}

void MVNProfileProcess::SampleFrequencyStat(double* logprof)	{

	covmatrix->drawVal(logprof);
}

double MVNProfileProcess::LogFrequencyStatPrior(double* logprof)	{

	for (int k=0; k<GetDim(); k++)	{
		logprof[k] -= mean[k];
	}

	double ret = covmatrix->logValProb(logprof);

	for (int k=0; k<GetDim(); k++)	{
		logprof[k] += mean[k];
	}

	return ret;
}

void MVNProfileProcess::UpdateFrequencyStat(double* logprof, double* prof)	{

	double tot = 0;
	for (int k=0; k<GetDim(); k++)	{
		prof[k] = exp(logprof[k]);
		if (prof[k] < stateps)	{
			prof[k] = stateps;
		}
		tot += prof[k];
	}
	for (int k=0; k<GetDim(); k++)	{
		prof[k] /= tot;
	}
}

double MVNProfileProcess::MoveKappa(double tuning, int nrep)	{

	double naccepted = 0;
	for (int rep=0; rep<nrep; rep++)	{
		for (int k=0; k<GetDim(); k++)	{
			double deltalogprob = - LogHyperPrior();
			double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
			double e = exp(m);
			kappa[k] *= e;
			deltalogprob += LogHyperPrior();
			deltalogprob += m;
			int accepted = (log(rnd::GetRandom().Uniform()) < deltalogprob);
			if (accepted)	{
				naccepted++;
			}
			else	{
				kappa[k] /= e;
			}
		}
	}
	return naccepted / nrep / GetDim();
}

void MVNProfileProcess::SampleHyper()	{

	for (int k=0; k<GetDim(); k++)	{
		kappa[k] = rnd::GetRandom().sExpo();
	}
	covmatrix->SetDiagonal();
	// covmatrix->Sample();
}

double MVNProfileProcess::LogHyperPrior()	{

	double total = 0;
	for (int k=0; k<GetDim(); k++)	{
		total -= kappa[k];
	}
	total += covmatrix->GetLogProb();
	return total;
}


double MVNProfileProcess::LogProfileProposeMove(double* logprofile, double tuning, int n)	{

	if (!n)	{
		n = GetDim();
	}
	int* indices = new int[n];
	rnd::GetRandom().DrawFromUrn(indices,n,GetDim());
	for (int i=0; i<n; i++)	{
		double h = tuning * (rnd::GetRandom().Uniform() - 0.5);
		logprofile[indices[i]] += h;
	}
	delete[] indices;
	return 0;
}
