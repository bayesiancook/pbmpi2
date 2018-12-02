
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#include "SubstitutionProcess.h"
#include "Random.h"

#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* Substitution Processes
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
//	* allocations / deallocations
//-------------------------------------------------------------------------

void SubstitutionProcess::Create()	{
	RateProcess::Create();
	ProfileProcess::Create();
}

void SubstitutionProcess::Delete() {
	ProfileProcess::Delete();
	RateProcess::Delete();
};

void SubstitutionProcess::CreateCondSiteLogL()	{
	if (condsitelogL)	{
		cerr << "error in SubstitutionProcess::CreateSiteLogL\n";
		exit(1);
	}
	sitelogL = new double[GetNsite()];
	condsitelogL = new double*[GetNsite()];
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		condsitelogL[i] = new double[GetNrate(i)];
	}
}

void SubstitutionProcess::DeleteCondSiteLogL()	{
	if (condsitelogL)	{
		for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
			delete[] condsitelogL[i];
		}
		delete[] condsitelogL;
		delete[] sitelogL;
		condsitelogL = 0;
		sitelogL = 0;
	}
}

double*** SubstitutionProcess::CreateConditionalLikelihoodVector()	{

	double*** condl = new double**[GetNsite()];
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		condl[i] = new double*[GetNrate(i)];
		for (int j=0; j<GetNrate(i); j++)	{
			condl[i][j] = new double[GetNstate(i) + 1];
			double* tmp = condl[i][j];
			for (int k=0; k<GetNstate(i); k++)	{
				tmp[k] = 1.0;
			}
			tmp[GetNstate(i)] = 0;
		}
	}
	return condl;
}

void SubstitutionProcess::DeleteConditionalLikelihoodVector(double*** condl)	{
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		for (int j=0; j<GetNrate(i); j++)	{
			delete[] condl[i][j];
		}
		delete[] condl[i];
	}
	delete[] condl;
	condl = 0;
}

void SubstitutionProcess::SamplePaths(BranchSitePath** patharray, int* stateup, int* statedown, double time)	{
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i))	{
			if (patharray[i])	{
				cerr << "error: patharray already exists\n";
				exit(1);
				delete patharray[i];
			}
			patharray[i] = SampleSitePath(i,stateup[i],statedown[i],time);
		}
	}
}

void SubstitutionProcess::SampleRootPaths(BranchSitePath** patharray, int* rootstate)	{
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i))	{
			if (patharray[i])	{
				cerr << "error: patharray already exists\n";
				exit(1);
				delete patharray[i];
			}
			patharray[i] = SampleRootSitePath(i,rootstate[i]);
		}
	}
}

//-------------------------------------------------------------------------
//	* elementary computations on conditional likelihood vectors 
//	(CPU level 1)
//-------------------------------------------------------------------------

// set the vector uniformly to 1 
void SubstitutionProcess::Reset(double*** t, bool condalloc, int phase)	{
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i))	{
			for (int j=0; j<GetNrate(i); j++)	{
				if ((! condalloc) || (ratealloc[i] == j))	{
					double* tmp = t[i][j];
					int nstate = GetNstate(i);
					for (int k=0; k<nstate; k++)	{
						(*tmp++) = phase;
					}
					*tmp = 0;
					tmp -= nstate;
				}
			}
		}
	}
}
	
// initialize the vector according to the data observed at a given leaf of the tree (contained in const int* state)
// steta[i] == -1 means 'missing data'. in that case, conditional likelihoods are all 1
void SubstitutionProcess::Initialize(double*** t, const int* state, bool condalloc)	{
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i))	{
			for (int j=0; j<GetNrate(i); j++)	{
				if ((! condalloc) || (ratealloc[i] == j))	{
					double* tmp = t[i][j];
					int nstate = GetNstate(i);
					tmp[nstate] = 0;
					if (state[i] == -1)	{
						for (int k=0; k<nstate; k++)	{
							(*tmp++) = 1.0;
						}
						tmp -= nstate;
					}
					else	{
						for (int k=0; k<nstate; k++)	{
							(*tmp++) = 0;
						}
						tmp -= nstate;
						tmp[state[i]] = 1.0;
					}
				}
			}
		}
	}
}

// multiply two conditional likelihood vectors, term by term
void SubstitutionProcess::Multiply(double*** from, double*** to, bool condalloc)	{
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i))	{
			for (int j=0; j<GetNrate(i); j++)	{
				if ((! condalloc) || (ratealloc[i] == j))	{
					double* tmpfrom = from[i][j];
					double* tmpto = to[i][j];
					int nstate = GetNstate(i);
					for (int k=0; k<nstate; k++)	{
						(*tmpto++) *= (*tmpfrom++);
					}
					*tmpto += *tmpfrom;
					tmpto -= nstate;
					tmpfrom -= nstate;
				}
			}
		}
	}
}

// Add two conditional likelihood vectors, term by term
void SubstitutionProcess::Add(double*** from, double*** to, bool condalloc)	{
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i))	{
			for (int j=0; j<GetNrate(i); j++)	{
				if ((! condalloc) || (ratealloc[i] == j))	{
					double* tmpfrom = from[i][j];
					double* tmpto = to[i][j];
					int nstate = GetNstate(i);
					for (int k=0; k<nstate; k++)	{
						(*tmpto++) += (*tmpfrom++);
					}
					*tmpto += *tmpfrom;
					tmpto -= nstate;
					tmpfrom -= nstate;
				}
			}
		}
	}
}

// multiply a conditional likelihood vector by the (possibly site-specific) stationary probabilities of the process
void SubstitutionProcess::MultiplyByStationaries(double*** to, bool condalloc)	{
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i))	{
			const double* stat = GetStationary(i);
			for (int j=0; j<GetNrate(i); j++)	{
				if ((! condalloc) || (ratealloc[i] == j))	{
					double* tmpto = to[i][j];
					int nstate = GetNstate(i);
					for (int k=0; k<nstate; k++)	{	
						(*tmpto++) *= (*stat++);
					}
					tmpto -= nstate;
					stat -= nstate;
				}
			}
		}
	}
}

// to avoid numerical errors: all entries for a given site and a given rate
// are divided by the largest among them
// and the residual is stored in the last entry of the vector
void SubstitutionProcess::Offset(double*** t, bool condalloc)	{
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i))	{
			for (int j=0; j<GetNrate(i); j++)	{
				if ((! condalloc) || (ratealloc[i] == j))	{
					double* tmp = t[i][j];
					double max = 0;
					for (int k=0; k<GetNstate(i); k++)	{
						if (tmp[k] <0)	{
							cerr << "error in pruning: negative prob : " << tmp[k] << "\n";
							exit(1);
							// tmp[k] = 0;
						}
						if (max < tmp[k])	{
							max = tmp[k];
						}
					}
					if (max == 0)	{
						max = 1e-12;
					}
					if (max < 0)	{
						cerr << "error in pruning (offset function): null likelihood\n";
						exit(1);
					}
					if (max > 0)	{
						for (int k=0; k<GetNstate(i); k++)	{
							tmp[k] /= max;
						}
					}
					tmp[GetNstate(i)] += log(max);
				}
			}
		}
	}
}

//-------------------------------------------------------------------------
//	* likelihood computation (last step, once the recursion has proceeded throughout the entire tree) 
//	(CPU level 2)
//-------------------------------------------------------------------------

double SubstitutionProcess::ComputeLikelihood(double*** aux, bool condalloc)	{

	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i))	{
			if (condalloc)	{
				int j = ratealloc[i];
				double* t = aux[i][j];
				double tot = 0;
				int nstate = GetNstate(i);
				for (int k=0; k<nstate; k++)	{
					tot += (*t++);
				}
				if (tot == 0)	{
					// dirty !
					tot = 1e-12;
				}
				sitelogL[i] = log(tot) + (*t);
                if (isnan(sitelogL[i])) {
                    cerr << "in SubstitutionProcess::ComputeLikelihood: nan\n";
                    cerr << tot << '\n';
                    exit(1);
                }
				t -= nstate;
			}
			else	{
				double max = 0;
				double* logl = condsitelogL[i];
				for (int j=0; j<GetNrate(i); j++)	{
					double* t = aux[i][j];
					double tot = 0;
					int nstate = GetNstate(i);
					for (int k=0; k<nstate; k++)	{
						tot += (*t++);
					}
					if (tot < 0)	{
						cerr << "error in SubstitutionProcess::ComputeLikelihood: negative prob\n";
						cerr << tot << '\n';
						exit(1);
					}
					if (isnan(tot))	{
						cerr << "error in SubstitutionProcess::ComputeLikelihood: tot is nan\n";
						for (int k=0; k<nstate; k++)	{
							cerr << GetStationary(i)[k] << '\t';
						}
						cerr << '\n';
						cerr << GetMyid() << '\t' << GetNprocs() << '\n';
						exit(1);
					}
					if (tot == 0)	{
						// dirty !
						tot = 1e-12;
					}
					logl[j] = log(tot) + (*t);
                    if (isinf(logl[j])) {
                        cerr << "in SubstitutionProcess::ComputeLikelihood: inf\n";
                        cerr << tot << '\t' << (*t) << '\t' << log(tot) << '\t' << log(tot) + (*t) << '\n';
                        exit(1);
                    }
					t -= nstate;
					if ((!j) || (max < logl[j]))	{
						max = logl[j];
					}
				}
				double total = 0;
				double meanrate = 0;
				for (int j=0; j<GetNrate(i); j++)	{
					double tmp = GetRateWeight(i,j) * exp(logl[j] - max);
					total += tmp;
					meanrate += tmp * GetRate(i,j);
				}
				sitelogL[i] = log(total) + max;
                if (isnan(sitelogL[i])) {
                    cerr << "in SubstitutionProcess::ComputeLikelihood: nan\n";
                    cerr << total << '\t' << max << '\n';
                    exit(1);
                }
				meanrate /= total;
				meansiterate[i] = meanrate;
			}
		}
	}

	logL = 0;
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i))	{
			logL += sitelogL[i];
		}
	}
	return logL;
}
	

//-------------------------------------------------------------------------
//	* sample the allocation of each site to one of the available rate categories
// 	for each site, the rate caegory is chosen with probability proportional
//	to the likelihood for that site, under that rate category
//	is aux == 0, this likelihood is assumed to be already computed and stored in int** condsitelogL (a member of SubstitutionProcess)
//	otherwise, aux is assumed to be a conditional likelihood vector multiplied by stationaries, and thus can be used to compute those likelihoods 
//	(CPU level 2)
//-------------------------------------------------------------------------

void SubstitutionProcess::DrawAllocations(double*** aux)	{

	if (condflag)	{
		cerr << "error in SubstitutionProcess::DrawAllocations: condflag is true\n";
		cerr << GetMyid() << '\n';
		exit(1);
	}
	if (aux)	{
		ComputeLikelihood(aux,false);
	}
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i))	{
			if (GetNrate(i) == 1)	{
				ratealloc[i] = 0;
			}
			else	{
				double max = 0;
				double* logl = condsitelogL[i];
				for (int j=0; j<GetNrate(i); j++)	{
					if ((!j) || (max < logl[j]))	{
						max = logl[j];
					}
				}
				double cumul = 0;
				double* p = new double[GetNrate(i)];
				for (int j=0; j<GetNrate(i); j++)	{
					cumul += GetRateWeight(i,j) * exp(logl[j] - max);
					p[j] = cumul;
				}
				double u = rnd::GetRandom().Uniform() * cumul;
				int j = 0;
				while ((j<GetNrate(i)) && (p[j] < u)) j++;
				if (j == GetNrate(i))	{
					cerr << "error in SubstitutionProcess::SampleAlloc\n";
					exit(1);
				}
				ratealloc[i] = j;
				delete[] p;
			}
		}
	}
}

//-------------------------------------------------------------------------
//	* sample states for each site, based on vectors of posterior probability stored in double*** t
//	this is assumed to be conditional on rate allocations 
//	(CPU level 2)
//-------------------------------------------------------------------------

void SubstitutionProcess::ChooseStates(double*** t, int* states)	{
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i))	{
			int j = ratealloc[i];
			double* tmp = t[i][j];
			double total = 0;
			for (int k=0; k<GetNstate(i); k++)	{
				total += tmp[k];
			}
			double u = rnd::GetRandom().Uniform() * total;
			double tot = tmp[0];
			int k = 0;
			while ((k<GetNstate(i)) && (tot < u))	{
				k++;
				if (k==GetNstate(i))	{
					cerr << "error in SubstitutionProcess::ChooseState\n";
					exit(1);
				}
				tot += tmp[k];
			}
			states[i] = k;
			for (int l=0; l<GetNstate(i); l++)	{
				tmp[l] = 0;
			}
			tmp[k] =  1;
		}
	}
}

//-------------------------------------------------------------------------
//	* elementary computations on conditional likelihood vectors 
//	(CPU level 1)
//-------------------------------------------------------------------------

// set the vector uniformly to 1 
void SubstitutionProcess::SiteReset(int i, double** t, bool condalloc)	{
	for (int j=0; j<GetNrate(i); j++)	{
		if ((! condalloc) || (ratealloc[i] == j))	{
			double* tmp = t[j];
			int nstate = GetNstate(i);
			for (int k=0; k<nstate; k++)	{
				(*tmp++) = 1.0;
			}
			*tmp = 0;
			tmp -= nstate;
		}
	}
}
	
// initialize the vector according to the data observed at a given leaf of the tree (contained in const int* state)
// steta[i] == -1 means 'missing data'. in that case, conditional likelihoods are all 1
void SubstitutionProcess::SiteInitialize(int i, double** t, const int state, bool condalloc)	{

	for (int j=0; j<GetNrate(i); j++)	{
		if ((! condalloc) || (ratealloc[i] == j))	{
			double* tmp = t[j];
			int nstate = GetNstate(i);
			/*
			if ((state < -1) || (state >= nstate))	{
				cerr << "error: state out of bound\n";
				cerr << state << '\n';
				exit(1);
			}
			*/
			tmp[nstate] = 0;
			if (state == -1)	{
				for (int k=0; k<nstate; k++)	{
					(*tmp++) = 1.0;
				}
				tmp -= nstate;
			}
			else	{
				for (int k=0; k<nstate; k++)	{
					(*tmp++) = 0;
				}
				tmp -= nstate;
				tmp[state] = 1.0;
			}
		}
	}
}

// multiply two conditional likelihood vectors, term by term
void SubstitutionProcess::SiteMultiply(int i, double** from, double** to, bool condalloc)	{
	for (int j=0; j<GetNrate(i); j++)	{
		if ((! condalloc) || (ratealloc[i] == j))	{
			double* tmpfrom = from[j];
			double* tmpto = to[j];
			int nstate = GetNstate(i);
			for (int k=0; k<nstate; k++)	{
				(*tmpto++) *= (*tmpfrom++);
			}
			*tmpto += *tmpfrom;
			tmpto -= nstate;
			tmpfrom -= nstate;
		}
	}
}

// multiply a conditional likelihood vector by the (possibly site-specific) stationary probabilities of the process
void SubstitutionProcess::SiteMultiplyByStationaries(int i, double** to, bool condalloc)	{
	const double* stat = GetStationary(i);
	for (int j=0; j<GetNrate(i); j++)	{
		if ((! condalloc) || (ratealloc[i] == j))	{
			double* tmpto = to[j];
			int nstate = GetNstate(i);
			for (int k=0; k<nstate; k++)	{	
				(*tmpto++) *= (*stat++);
			}
			tmpto -= nstate;
			stat -= nstate;
		}
	}
}

// to avoid numerical errors: all entries for a given site and a given rate
// are divided by the largest among them
// and the residual is stored in the last entry of the vector
void SubstitutionProcess::SiteOffset(int i, double** t, bool condalloc)	{
	for (int j=0; j<GetNrate(i); j++)	{
		if ((! condalloc) || (ratealloc[i] == j))	{
			double* tmp = t[j];
			double max = 0;
			for (int k=0; k<GetNstate(i); k++)	{
				if (tmp[k] <0)	{
					cerr << "error in pruning: negative prob : " << tmp[k] << "\n";
					exit(1);
					tmp[k] = 0;
				}
				if (max < tmp[k])	{
					max = tmp[k];
				}
			}
			/*
			if (max == 0)	{
				// max = 1e-12;
				cerr << "error in pruning: null likelihood\n";
				exit(1);
			}
			*/
			if (max > 0)	{
				for (int k=0; k<GetNstate(i); k++)	{
					tmp[k] /= max;
				}
			}
			tmp[GetNstate(i)] += log(max);
		}
	}
}

//-------------------------------------------------------------------------
//	* likelihood computation (last step, once the recursion has proceeded throughout the entire tree) 
//	(CPU level 2)
//-------------------------------------------------------------------------

double SubstitutionProcess::SiteComputeLikelihood(int i, double** aux, bool condalloc)	{

	if (condalloc)	{
		int j = ratealloc[i];
		double* t = aux[j];
		double tot = 0;
		int nstate = GetNstate(i);
		for (int k=0; k<nstate; k++)	{
			tot += (*t++);
		}
		if (tot == 0)	{
			// dirty !
			tot = 1e-12;
			/*
			cerr << "pruning : 0 \n";
			for (int k=0; k<GetNstate(i); k++)	{
				cerr << t[k] << '\n';
			}
			exit(1);
			*/
		}
		sitelogL[i] = log(tot) + (*t);
		t -= nstate;
	}
	else	{
		double max = 0;
		double* logl = condsitelogL[i];
		for (int j=0; j<GetNrate(i); j++)	{
			double* t = aux[j];
			double tot = 0;
			int nstate = GetNstate(i);
			for (int k=0; k<nstate; k++)	{
				tot += (*t++);
			}
			if (tot < 0)	{
				cerr << "negative total in site likelihood\n";
				exit(1);
				// dirty !
				tot = 1e-12;
				/*
				cerr << "pruning : 0 \n";
				for (int k=0; k<GetNstate(i); k++)	{
					cerr << t[k] << '\n';
				}
				exit(1);
				*/
			}
			logl[j] = log(tot) + (*t);
			t -= nstate;
			if ((!j) || (max < logl[j]))	{
				max = logl[j];
			}
		}
		double total = 0;
		double meanrate = 0;
		for (int j=0; j<GetNrate(i); j++)	{
			double tmp = 0;
			if (! isinf(logl[j]))	{
				tmp = GetRateWeight(i,j) * exp(logl[j] - max);
			}
			total += tmp;
			meanrate += tmp * GetRate(i,j);
		}
		sitelogL[i] = log(total) + max;
		if (isnan(sitelogL[i]))	{
			cerr << "nan site logl\n";
			cerr << max << '\t' << total << '\n';
			for (int j=0; j<GetNrate(i); j++)	{
				cerr << '\t' << logl[j] << '\t' << exp(logl[j] - max) << '\n';
			}
			exit(1);
		}
		meanrate /= total;
		meansiterate[i] = meanrate;
	}

	return sitelogL[i];
}
	

//-------------------------------------------------------------------------
//	* sample the allocation of each site to one of the available rate categories
// 	for each site, the rate caegory is chosen with probability proportional
//	to the likelihood for that site, under that rate category
//	is aux == 0, this likelihood is assumed to be already computed and stored in int** condsitelogL (a member of SubstitutionProcess)
//	otherwise, aux is assumed to be a conditional likelihood vector multiplied by stationaries, and thus can be used to compute those likelihoods 
//	(CPU level 2)
//-------------------------------------------------------------------------

void SubstitutionProcess::SiteDrawAllocations(int i, double** aux)	{

	if (aux)	{
		SiteComputeLikelihood(i,aux);
	}

	if (GetNrate(i) == 1)	{
		ratealloc[i] = 0;
	}
	else	{
		double max = 0;
		double* logl = condsitelogL[i];
		for (int j=0; j<GetNrate(i); j++)	{
			if ((!j) || (max < logl[j]))	{
				max = logl[j];
			}
		}
		double cumul = 0;
		double* p = new double[GetNrate(i)];
		for (int j=0; j<GetNrate(i); j++)	{
			cumul += GetRateWeight(i,j) * exp(logl[j] - max);
			p[j] = cumul;
		}
		double u = rnd::GetRandom().Uniform() * cumul;
		int j = 0;
		while ((j<GetNrate(i)) && (p[j] < u)) j++;
		if (j == GetNrate(i))	{
			cerr << "error in SubstitutionProcess::SampleAlloc\n";
			exit(1);
		}
		ratealloc[i] = j;
		delete[] p;
	}
}

//-------------------------------------------------------------------------
//	* sample states for each site, based on vectors of posterior probability stored in double*** t
//	this is assumed to be conditional on rate allocations 
//	(CPU level 2)
//-------------------------------------------------------------------------

int SubstitutionProcess::SiteChooseState(int i, double** t)	{

	int j = ratealloc[i];
	double* tmp = t[j];
	double total = 0;
	for (int k=0; k<GetNstate(i); k++)	{
		total += tmp[k];
	}
	double u = rnd::GetRandom().Uniform() * total;
	double tot = tmp[0];
	int k = 0;
	while ((k<GetNstate(i)) && (tot < u))	{
		k++;
		if (k==GetNstate(i))	{
			cerr << "error in SubstitutionProcess::ChooseState\n";
			exit(1);
		}
		tot += tmp[k];
	}
	for (int l=0; l<GetNstate(i); l++)	{
		tmp[l] = 0;
	}
	tmp[k] =  1;
	return k;
}

