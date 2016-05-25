
#include "DirichletProfileProcess.h"
#include "Random.h"

void DirichletProfileProcess::Create()	{

	if (! dirweight)	{
		dirweight = new double[GetDim()];
	}
}

void DirichletProfileProcess::Delete()	{

	if (dirweight)	{
		delete[] dirweight;
		dirweight = 0;
	}
}

double DirichletProfileProcess::GetMeanDirWeight()	{
	double total = 0;
	for (int k=0; k<GetDim(); k++)	{
		total += dirweight[k];
	}
	return total;
}

double DirichletProfileProcess::GetCenterStatEnt()	{
	double totalweight = 0;
	for (int k=0; k<GetDim(); k++)	{
		totalweight += dirweight[k];
	}
	double total = 0;
	for (int k=0; k<GetDim(); k++)	{
		double w = dirweight[k] / totalweight;
		total -= w * log(w);
	}
	return total;
}

void DirichletProfileProcess::SampleFrequencyStat(double* prof)	{

	double statmin = stateps;

	double total = 0;
	int infreached = 0;
	for (int k=0; k<GetDim(); k++)	{
		prof[k] = rnd::GetRandom().sGamma(dirweight[k]);
		if (prof[k] < statmin)	{
			prof[k] = statmin;
			infreached = 1;
		}
		total += prof[k];
	}
	for (int k=0; k<GetDim(); k++)	{
		prof[k] /= total;
	}
	if (infreached)	{
		statinfcount++;
	}
	totstatcount++;
}

double DirichletProfileProcess::LogFrequencyStatPrior(double* prof)	{
	double total = 0;
	double totalweight = 0;
	for (int k=0; k<GetDim(); k++)	{
		total += (dirweight[k] - 1) * log(prof[k]) - rnd::GetRandom().logGamma(dirweight[k]);
		totalweight += dirweight[k];
	}
	total += rnd::GetRandom().logGamma(totalweight);
	return total;
}

double DirichletProfileProcess::MoveDirWeights(double tuning, int nrep)	{
	double naccepted = 0;
	for (int rep=0; rep<nrep; rep++)	{
		for (int k=0; k<GetDim(); k++)	{
			double deltalogprob = - LogHyperPrior() - LogStatPrior();
			double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
			double e = exp(m);
			dirweight[k] *= e;
			deltalogprob += LogHyperPrior() + LogStatPrior();
			deltalogprob += m;
			int accepted = (log(rnd::GetRandom().Uniform()) < deltalogprob);
			if (accepted)	{
				naccepted++;
			}
			else	{
				dirweight[k] /= e;
			}
		}
	}
	return naccepted / nrep / GetDim();
}

void DirichletProfileProcess::SampleHyper()	{

	for (int i=0; i<GetDim(); i++)	{
		dirweight[i] = 1;
	}
}

void DirichletProfileProcess::PriorSampleHyper()	{

	for (int i=0; i<GetDim(); i++)	{
		dirweight[i] = rnd::GetRandom().sExpo();
	}
}

double DirichletProfileProcess::LogHyperPrior()	{

	double total = 0;
	for (int k=0; k<GetDim(); k++)	{
		total -= dirweight[k];
	}
	return total;
}

