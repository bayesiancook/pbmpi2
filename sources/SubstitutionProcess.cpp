
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

void SubstitutionProcess::Create(int site, int dim, int insitemin, int insitemax)	{
	sitemin = insitemin;
	sitemax = insitemax;
	//cout << sitemin << "  " << sitemax << endl;
	if (! ratealloc)	{
		RateProcess::Create(site);
		ProfileProcess::Create(site, dim);
		// ratealloc = new int[sitemax - sitemin];
		ratealloc = new int[GetNsite()];
	}
}

void SubstitutionProcess::Delete() {
	if (ratealloc)	{
		delete[] ratealloc;
		ratealloc = 0;
		ProfileProcess::Delete();
		RateProcess::Delete();
	}
};

void SubstitutionProcess::CreateCondSiteLogL()	{
	if (condsitelogL)	{
		cerr << "error in SubstitutionProcess::CreateSiteLogL\n";
		exit(1);
	}
	// sitelogL = new double[sitemax - sitemin];
	// condsitelogL = new double*[sitemax - sitemin];
	sitelogL = new double[GetNsite()];
	meansiterate = new double[GetNsite()];
	condsitelogL = new double*[GetNsite()];
	for (int i=sitemin; i<sitemax; i++)	{
	// for (int i=0; i<GetNsite(); i++)	{
		condsitelogL[i] = new double[GetNrate(i)];
	}
}

void SubstitutionProcess::DeleteCondSiteLogL()	{
	if (condsitelogL)	{
		for (int i=sitemin; i<sitemax; i++)	{
		// for (int i=0; i<GetNsite(); i++)	{
			delete[] condsitelogL[i];
		}
		delete[] condsitelogL;
		delete[] sitelogL;
		delete[] meansiterate;
		condsitelogL = 0;
		sitelogL = 0;
	}
}

double*** SubstitutionProcess::CreateConditionalLikelihoodVector()	{
	//cout << "VECTOR ALLOCATION: " << sitemax << "  " << sitemin << endl;
	// double*** condl = new double**[sitemax - sitemin];
	double*** condl = new double**[GetNsite()];
	for (int i=sitemin; i<sitemax; i++)	{
	// for (int i=0; i<GetNsite(); i++)	{
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
	//cout << "Test element " << condl[0][1][0] << endl;
	return condl;
}

void SubstitutionProcess::DeleteConditionalLikelihoodVector(double*** condl)	{
	for (int i=sitemin; i<sitemax; i++)	{
	// for (int i=0; i<GetNsite(); i++)	{
		for (int j=0; j<GetNrate(i); j++)	{
			delete[] condl[i][j];
		}
		delete[] condl[i];
	}
	delete[] condl;
	condl = 0;
}

//-------------------------------------------------------------------------
//	* elementary computations on conditional likelihood vectors 
//	(CPU level 1)
//-------------------------------------------------------------------------

// set the vector uniformly to 1 
void SubstitutionProcess::Reset(double*** t, bool condalloc)	{
	for (int i=sitemin; i<sitemax; i++)	{
	// for (int i=0; i<GetNsite(); i++)	{
		for (int j=0; j<GetNrate(i); j++)	{
			if ((! condalloc) || (ratealloc[i] == j))	{
				double* tmp = t[i][j];
				int nstate = GetNstate(i);
				for (int k=0; k<nstate; k++)	{
					(*tmp++) = 1.0;
					// tmp[k] = 1.0;
				}
				*tmp = 0;
				tmp -= nstate;
				// tmp[GetNstate(i)] = 0;
			}
		}
	}
}
	
// initialize the vector according to the data observed at a given leaf of the tree (contained in const int* state)
// steta[i] == -1 means 'missing data'. in that case, conditional likelihoods are all 1
void SubstitutionProcess::Initialize(double*** t, const int* state, bool condalloc)	{
	for (int i=sitemin; i<sitemax; i++)	{
	// for (int i=0; i<GetNsite(); i++)	{
		for (int j=0; j<GetNrate(i); j++)	{
			if ((! condalloc) || (ratealloc[i] == j))	{
				double* tmp = t[i][j];
				int nstate = GetNstate(i);
				tmp[nstate] = 0;
				if (state[i] == -1)	{
					for (int k=0; k<nstate; k++)	{
						(*tmp++) = 1.0;
						// tmp[k] = 1.0;
					}
					tmp -= nstate;
				}
				else	{
					for (int k=0; k<nstate; k++)	{
						(*tmp++) = 0;
						// tmp[k] = 0;
					}
					tmp -= nstate;
					tmp[state[i]] = 1.0;
				}
			}
		}
	}
}

// multiply two conditional likelihood vectors, term by term
void SubstitutionProcess::Multiply(double*** from, double*** to, bool condalloc)	{
	for (int i=sitemin; i<sitemax; i++)	{
	// for (int i=0; i<GetNsite(); i++)	{
		for (int j=0; j<GetNrate(i); j++)	{
			if ((! condalloc) || (ratealloc[i] == j))	{
				double* tmpfrom = from[i][j];
				double* tmpto = to[i][j];
				int nstate = GetNstate(i);
				for (int k=0; k<nstate; k++)	{
					(*tmpto++) *= (*tmpfrom++);
					// tmpto[k] *= tmpfrom[k];
				}
				*tmpto += *tmpfrom;
				tmpto -= nstate;
				tmpfrom -= nstate;
				// tmpto[GetNstate(i)] += tmpfrom[GetNstate(i)];
			}
		}
	}
}

// multiply a conditional likelihood vector by the (possibly site-specific) stationary probabilities of the process
void SubstitutionProcess::MultiplyByStationaries(double*** to, bool condalloc)	{
	for (int i=sitemin; i<sitemax; i++)	{
	// for (int i=0; i<GetNsite(); i++)	{
		const double* stat = GetStationary(i);
		for (int j=0; j<GetNrate(i); j++)	{
			if ((! condalloc) || (ratealloc[i] == j))	{
				double* tmpto = to[i][j];
				int nstate = GetNstate(i);
				for (int k=0; k<nstate; k++)	{	
					(*tmpto++) *= (*stat++);
					// tmpto[k] *= stat[k];
				}
				tmpto -= nstate;
				stat -= nstate;
			}
		}
	}
}

// to avoid numerical errors: all entries for a given site and a given rate
// are divided by the largest among them
// and the residual is stored in the last entry of the vector
void SubstitutionProcess::Offset(double*** t, bool condalloc)	{
	for (int i=sitemin; i<sitemax; i++)	{
	// for (int i=0; i<GetNsite(); i++)	{
		for (int j=0; j<GetNrate(i); j++)	{
			if ((! condalloc) || (ratealloc[i] == j))	{
				double* tmp = t[i][j];
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
				if (max == 0)	{
					max = 1e-12;
					/*
					cerr << "error in pruning: null likelihood\n";
					exit(1);
					*/
				}
				for (int k=0; k<GetNstate(i); k++)	{
					tmp[k] /= max;
				}
				tmp[GetNstate(i)] += log(max);
			}
		}
	}
}

//-------------------------------------------------------------------------
//	* likelihood computation (last step, once the recursion has proceeded throughout the entire tree) 
//	(CPU level 2)
//-------------------------------------------------------------------------

double SubstitutionProcess::ComputeLikelihood(double*** aux, bool condalloc)	{
	for (int i=sitemin; i<sitemax; i++)	{
	// for (int i=0; i<GetNsite(); i++)	{
		if (condalloc)	{
			int j = ratealloc[i];
			double* t = aux[i][j];
			double tot = 0;
			int nstate = GetNstate(i);
			for (int k=0; k<nstate; k++)	{
				tot += (*t++);
				// tot += t[k];
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
			// sitelogL[i] = log(tot) + t[GetNstate(i)];
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
					// tot += t[k];
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
				logl[j] = log(tot) + (*t);
				// logl[j] = log(tot) + t[GetNstate(i)];
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
			meanrate /= total;
			meansiterate[i] = meanrate;
		}
	}

	logL = 0;
	for (int i=sitemin; i<sitemax; i++)	{
	// for (int i=0; i<GetNsite(); i++)	{
		logL += sitelogL[i];
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

	if (aux)	{
		ComputeLikelihood(aux);
	}
	for (int i=sitemin; i<sitemax; i++)	{
	// for (int i=0; i<GetNsite(); i++)	{
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
		}
	}
}

//-------------------------------------------------------------------------
//	* sample states for each site, based on vectors of posterior probability stored in double*** t
//	this is assumed to be conditional on rate allocations 
//	(CPU level 2)
//-------------------------------------------------------------------------

void SubstitutionProcess::ChooseStates(double*** t, int* states)	{
	for (int i=sitemin; i<sitemax; i++)	{
	// for (int i=0; i<GetNsite(); i++)	{
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

