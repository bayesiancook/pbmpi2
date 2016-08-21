
#include "DirichletProfileProcess.h"
#include "Random.h"

void DirichletProfileProcess::Create()	{

	if (dirpriortype)	{
		if (! dirweight)	{
			dirweight = new double[GetDim()];
		}
	}
	else	{
		if (! statalpha)	{
			statalpha = new double[Nstatcomp];
			statweight = new double[Nstatcomp];
			statcenter = new double*[Nstatcomp];
			for (int k=0; k<Nstatcomp; k++)	{
				statcenter[k] = new double[GetDim()];
			}
		}
	}
}

void DirichletProfileProcess::Delete()	{

	if (dirpriortype)	{
		if (dirweight)	{
			delete[] dirweight;
			dirweight = 0;
		}
	}
	else	{
		if (statalpha)	{
			delete[] statalpha;
			statalpha = 0;
			delete[] statweight;
			for (int k=0; k<Nstatcomp; k++)	{
				delete[] statcenter[k];
			}
			delete[] statcenter;
		}
	}
}

double DirichletProfileProcess::GetMeanDirWeight()	{
	double total = 0;
	if (dirpriortype)	{
		for (int k=0; k<GetDim(); k++)	{
			total += dirweight[k];
		}
	}
	else	{
		for (int k=0; k<Nstatcomp; k++)	{
			total += statalpha[k];
		}
		total /= Nstatcomp;
	}
	return total;
}

double DirichletProfileProcess::GetCenterStatEnt()	{
	double total = 0;
	if (dirpriortype)	{
		double totalweight = 0;
		for (int k=0; k<GetDim(); k++)	{
			totalweight += dirweight[k];
		}
		for (int k=0; k<GetDim(); k++)	{
			double w = dirweight[k] / totalweight;
			total -= w * log(w);
		}
	}
	else	{
		for (int k=0; k<Nstatcomp; k++)	{
			for (int l=0; l<GetDim(); l++)	{
				total -= statcenter[k][l] * log(statcenter[k][l]);
			}
		}
		total /= Nstatcomp;
	}
	return total;
}

void DirichletProfileProcess::SampleFrequencyStat(double* prof)	{

	double statmin = stateps;
	double total = 0;
	int infreached = 0;


	if (dirpriortype)	{
		for (int k=0; k<GetDim(); k++)	{
			prof[k] = rnd::GetRandom().sGamma(dirweight[k]);
		}
	}
	else	{
		int k = rnd::GetRandom().FiniteDiscrete(Nstatcomp,statweight);	
		for (int l=0; l<GetDim(); l++)	{
			prof[l] = rnd::GetRandom().sGamma(statalpha[k] * statcenter[k][l]);
		}
	}

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
	if (dirpriortype)	{
		double totalweight = 0;
		for (int k=0; k<GetDim(); k++)	{
			total += (dirweight[k] - 1) * log(prof[k]) - rnd::GetRandom().logGamma(dirweight[k]);
			totalweight += dirweight[k];
		}
		total += rnd::GetRandom().logGamma(totalweight);
	}
	else	{
		double tmp[Nstatcomp];
		double max = 0;
		for (int k=0; k<Nstatcomp; k++)	{
			double tot = 0;
			for (int l=0; l<GetDim(); l++)	{
				tot += (statalpha[k] * statcenter[k][l] - 1) * log(prof[k]) - rnd::GetRandom().logGamma(statalpha[k]*statcenter[k][l]);
			}
			tot += rnd::GetRandom().logGamma(statalpha[k]);
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

	if (dirpriortype)	{
		cerr << "error: in DirichletProfileProcess::GetPostStatProb; under simple Dirichlet base distribution\n";
		exit(1);
	}

	double max = 0;
	for (int k=0; k<Nstatcomp; k++)	{
		double tot = 0;
		for (int l=0; l<GetDim(); l++)	{
			tot += (statalpha[k] * statcenter[k][l] - 1) * log(prof[k]) - rnd::GetRandom().logGamma(statalpha[k]*statcenter[k][l]);
		}
		tot += rnd::GetRandom().logGamma(statalpha[k]);
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

double DirichletProfileProcess::MoveStatCenter(double tuning, int n, int nrep)	{

	double naccepted = 0;
	double bk[GetDim()];
	for (int rep=0; rep<nrep; rep++)	{
		for (int k=0; k<Nstatcomp; k++)	{
			for (int l=0; l<GetDim(); l++)	{
				bk[l] = statcenter[k][l];
			}
			double deltalogprob = - LogHyperPrior() - LogStatPrior();
			double logh = ProfileProposeMove(statcenter[k],tuning,n,0,0,0);
			deltalogprob += LogHyperPrior() + LogStatPrior();
			deltalogprob += logh;
			int accepted = (log(rnd::GetRandom().Uniform()) < deltalogprob);
			if (accepted)	{
				naccepted++;
			}
			else	{
				for (int l=0; l<GetDim(); l++)	{
					statcenter[k][l] = bk[l];
				}
			}
		}
	}
	return naccepted / nrep / Nstatcomp;
}

double DirichletProfileProcess::MoveStatAlpha(double tuning, int nrep)	{

	double naccepted = 0;
	for (int rep=0; rep<nrep; rep++)	{
		for (int k=0; k<Nstatcomp; k++)	{
			double deltalogprob = - LogHyperPrior() - LogStatPrior();
			double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
			double e = exp(m);
			statalpha[k] *= e;
			deltalogprob += LogHyperPrior() + LogStatPrior();
			deltalogprob += m;
			int accepted = (log(rnd::GetRandom().Uniform()) < deltalogprob);
			if (accepted)	{
				naccepted++;
			}
			else	{
				statalpha[k] /= e;
			}
		}
	}
	return naccepted / nrep / Nstatcomp;
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

	if (dirpriortype)	{
		for (int i=0; i<GetDim(); i++)	{
			dirweight[i] = 1;
		}
	}
	else	{
		if (priorempmix)	{
			SetEmpiricalPriorStatMix();
		}
		else	{
			if (! fixstatweight)	{
				double tot = 0;
				for (int k=0; k<Nstatcomp; k++)	{
					statweight[k] = rnd::GetRandom().sExpo();
					tot += statweight[k];
				}
				for (int k=0; k<Nstatcomp; k++)	{
					statweight[k] /= tot;
				}
			}
			if (! fixstatalpha)	{
				for (int k=0; k<Nstatcomp; k++)	{
					statalpha[k] = GetDim() * rnd::GetRandom().sExpo();
				}
			}
			if (! fixstatcenter)	{
				for (int k=0; k<Nstatcomp; k++)	{
					double tot = 0;
					for (int l=0; l<GetDim(); k++)	{
						statcenter[k][l] = rnd::GetRandom().sExpo();
						tot += statcenter[k][l];
					}
					for (int l=0; l<GetDim(); l++)	{
						statcenter[k][l] /= tot;
					}
				}
			}
		}
	}
}

void DirichletProfileProcess::PriorSampleHyper()	{

	if (dirpriortype)	{
		for (int i=0; i<GetDim(); i++)	{
			dirweight[i] = rnd::GetRandom().sExpo();
		}
	}
	else	{
		if (priorempmix)	{
			SetEmpiricalPriorStatMix();
		}
		else	{
			if (! fixstatweight)	{
				double tot = 0;
				for (int k=0; k<Nstatcomp; k++)	{
					statweight[k] = rnd::GetRandom().sExpo();
					tot += statweight[k];
				}
				for (int k=0; k<Nstatcomp; k++)	{
					statweight[k] /= tot;
				}
			}
			if (! fixstatalpha)	{
				for (int k=0; k<Nstatcomp; k++)	{
					statalpha[k] = GetDim() * rnd::GetRandom().sExpo();
				}
			}
			if (! fixstatcenter)	{
				for (int k=0; k<Nstatcomp; k++)	{
					double tot = 0;
					for (int l=0; l<GetDim(); k++)	{
						statcenter[k][l] = rnd::GetRandom().sExpo();
						tot += statcenter[k][l];
					}
					for (int l=0; l<GetDim(); l++)	{
						statcenter[k][l] /= tot;
					}
				}
			}
		}
	}
}

void DirichletProfileProcess::SetEmpiricalPriorStatMix()	{

	ifstream is(priormixtype.c_str());
	if (!is)	{
		cerr << "error in DirichletProfileProcess::SetEmpiricalPriorStatMix\n";
		exit(1);
	}
	is >> Nstatcomp;
	cerr << "prior stat mixture: number of components: " << Nstatcomp << '\n';

	for (int k=0; k<Nstatcomp; k++)	{
		is >> statweight[k];
		is >> statalpha[k];
		for (int l=0; l<GetDim(); l++)	{
			is >> statcenter[k][l];
		}
	}
}

double DirichletProfileProcess::LogHyperPrior()	{

	double total = 0;
	if (dirpriortype)	{
		for (int k=0; k<GetDim(); k++)	{
			total -= dirweight[k];
		}
	}
	else	{
		if (! fixstatalpha)	{
			for (int k=0; k<Nstatcomp; k++)	{
				total -= statalpha[k]/GetDim();;
			}
		}
	}
	return total;
}

