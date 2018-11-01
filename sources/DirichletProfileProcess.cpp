
#include "DirichletProfileProcess.h"
#include "Random.h"

void DirichletProfileProcess::Create()	{

    if (! dirweight)	{
        if (priorempmix)   {
            SetEmpiricalPriorStatMix();
        }
        else    {
            Nstatcomp = 1;
            dirweight = new double*[1];
            dirweight[0] = new double[GetDim()];
            statweight = new double[1];
            statweight[0] = 1.0;
		}
	}
}

void DirichletProfileProcess::Delete()	{

    if (dirweight)	{
        for (int k=0; k<Nstatcomp; k++) {
            delete[] dirweight[k];
        }
        delete[] dirweight;
        dirweight = 0;
        delete[] statweight;
        statweight = 0;
    }
}

double DirichletProfileProcess::GetMeanDirWeight()	{
    double total = 0;
    for (int cat=0; cat<Nstatcomp; cat++)   {
        for (int k=0; k<GetDim(); k++)	{
            total += dirweight[cat][k];
        }
	}
    total /= Nstatcomp;
	return total;
}

double DirichletProfileProcess::GetCenterStatEnt()	{
    double mean = 0;
    for (int cat=0; cat<Nstatcomp; cat++)   {
        double total = 0;
        for (int k=0; k<GetDim(); k++)  {
            total += dirweight[cat][k];
        }
        for (int k=0; k<GetDim(); k++)  {
            double tmp = dirweight[cat][k] / total;
            mean -= tmp * log(tmp);
        }
    }
    mean /= Nstatcomp;
    return mean;
}

void DirichletProfileProcess::SampleFrequencyStat(double* prof)	{

    int cat = 0;
	if (Nstatcomp > 1)	{
		cat = rnd::GetRandom().FiniteDiscrete(Nstatcomp,statweight);	
    }
    for (int k=0; k<GetDim(); k++)	{
        prof[k] = rnd::GetRandom().sGamma(dirweight[cat][k]);
    }

    // check for numerical errors
    double statmin = stateps;
    int infreached = 0;
    double total = 0;
	for (int k=0; k<GetDim(); k++)	{
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
	if (Nstatcomp == 1)	{
		double totweight = 0;
		for (int l=0; l<GetDim(); l++)	{
			total += (dirweight[0][l] - 1) * log(prof[l]) - rnd::GetRandom().logGamma(dirweight[0][l]);
			totweight += dirweight[0][l];
		}
		total += rnd::GetRandom().logGamma(totweight);
	}
	else	{
		double tmp[Nstatcomp];
		double max = 0;
		for (int k=0; k<Nstatcomp; k++)	{
			double tot = 0;
            double totweight = 0;
			for (int l=0; l<GetDim(); l++)	{
				tot += (dirweight[k][l] - 1) * log(prof[l]) - rnd::GetRandom().logGamma(dirweight[k][l]);
                totweight += dirweight[k][l];
			}
			tot += rnd::GetRandom().logGamma(totweight);
			tmp[k] = tot;
			if ((!k) || (max < tot))	{
				max = tot;
			}
		}
		double tot = 0;
		for (int k=0; k<Nstatcomp; k++)	{
			tot += statweight[k] * exp(tmp[k] - max);
		}
		total = log(tot) + max;
	}
	return total;
}

double DirichletProfileProcess::GetPostStatProb(double* prof, double* post)	{

    if (Nstatcomp == 1) {
        cerr << "error: in DirichletProfileProcess::GetPostStatProb with Nstatcomp == 1\n";
        exit(1);
    }

	double max = 0;
	for (int k=0; k<Nstatcomp; k++)	{
		double tot = 0;
        double totweight = 0;
		for (int l=0; l<GetDim(); l++)	{
			tot += (dirweight[k][l] - 1) * log(prof[l]) - rnd::GetRandom().logGamma(dirweight[k][l]);
            totweight += dirweight[k][l];
		}
		tot += rnd::GetRandom().logGamma(totweight);
		post[k] = tot;
		if ((!k) || (max < tot))	{
			max = tot;
		}
	}
	double tot = 0;
	for (int k=0; k<Nstatcomp; k++)	{
		post[k] = statweight[k] * exp(post[k] - max);
		tot += post[k];
	}
	for (int k=0; k<Nstatcomp; k++)	{
		post[k] /= tot;
	}
	return log(tot) + max;
}

double DirichletProfileProcess::MoveHyper(double tuning, int nrep)	{
    if (((!fixstatcenter) && fixstatalpha) || (fixstatcenter && (!fixstatalpha)))   {
        cerr << "error in DirichletProfileProcess::Move: currently, should move all hyperparameters or all fixed\n";
        exit(1);
    }
    if (! fixstatcenter)    {
        if (Nstatcomp > 1)  {
            cerr << "in DirichletProfileProcess:: dirweight should be fixed with Nstatcomp > 1\n";
            exit(1);
        }
        for (int cat=0; cat<Nstatcomp; cat++)   {
            MoveDirWeights(cat,tuning,nrep);
        }
    }
    if (! fixstatweight)	{
        if (Nstatcomp > 1) {
            MoveStatWeight();
        }
    }
}

double DirichletProfileProcess::MoveStatWeight()	{

	int* count = new int[Nstatcomp];
	GetPostCount(count);
	double total = 0;
	for (int k=0; k<Nstatcomp; k++)	{
		statweight[k] = rnd::GetRandom().sGamma(1.0 + count[k]);
		total += statweight[k];
	}
	for (int k=0; k<Nstatcomp; k++)	{
		statweight[k] /= total;
	}
	delete[] count;
	return 1.0;
}

double DirichletProfileProcess::MoveDirWeights(int cat, double tuning, int nrep)	{

	double naccepted = 0;
	for (int rep=0; rep<nrep; rep++)	{
		for (int k=0; k<GetDim(); k++)	{
			double deltalogprob = - LogHyperPrior() - LogStatPrior();
			double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
			double e = exp(m);
			dirweight[cat][k] *= e;
			deltalogprob += LogHyperPrior() + LogStatPrior();
			deltalogprob += m;
			int accepted = (log(rnd::GetRandom().Uniform()) < deltalogprob);
			if (accepted)	{
				naccepted++;
			}
			else	{
				dirweight[cat][k] /= e;
			}
		}
	}
	return naccepted / nrep / GetDim();
}

void DirichletProfileProcess::SampleHyper()	{

    if (Nstatcomp == 1) {
		for (int i=0; i<GetDim(); i++)	{
			dirweight[0][i] = 1;
		}
	}
	else	{
        // nothing to do
    }
}

void DirichletProfileProcess::PriorSampleHyper()	{

	if (Nstatcomp == 1)	{
		for (int i=0; i<GetDim(); i++)	{
			dirweight[0][i] = rnd::GetRandom().sExpo();
		}
	}
	else	{
        // nothing to do
    }
}

void DirichletProfileProcess::SetEmpiricalPriorStatMix()	{

	ifstream is(priormixtype.c_str());
	if (!is)	{
		cerr << "error in DirichletProfileProcess::SetEmpiricalPriorStatMix\n";
		exit(1);
	}

	is >> Nstatcomp;

    dirweight = new double*[Nstatcomp];
    for (int k=0; k<Nstatcomp; k++) {
        dirweight[k] = new double[GetDim()];
    }
    statweight = new double[Nstatcomp];

	for (int k=0; k<Nstatcomp; k++)	{
		is >> statweight[k];
        double statalpha;
		is >> statalpha;
		for (int l=0; l<GetDim(); l++)	{
			is >> dirweight[k][l];
            dirweight[k][l] *= statalpha;
		}
	}
}

double DirichletProfileProcess::LogHyperPrior()	{
	double total = 0;
    for (int cat=0; cat<Nstatcomp; cat++)   {
		for (int k=0; k<GetDim(); k++)	{
			total -= dirweight[cat][k];
		}
	}
}

