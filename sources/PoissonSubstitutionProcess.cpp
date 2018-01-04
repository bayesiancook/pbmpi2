
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "PoissonSubstitutionProcess.h"

#include "Parallel.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* Poisson Substitution Processes
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

// same methods as above, but specialized for the CAT poisson model
// probably less interesting to parallelize, at least in a first step


//-------------------------------------------------------------------------
//	* conditional likelihood propagation
//	(CPU level 3)
//-------------------------------------------------------------------------

void PoissonSubstitutionProcess::Propagate(double*** from, double*** to, double time, bool condalloc)	{
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i))	{
			const double* stat = GetStationary(i);
			for (int j=0; j<GetNrate(i); j++)	{
				if ((! condalloc) || (ratealloc[i] == j))	{
					double* tmpfrom = from[i][j];
					double* tmpto = to[i][j];
					double expo = exp(-GetRate(i,j) * time);
					double tot = 0;
					int nstate = GetNstate(i);
					for (int k=0; k<nstate; k++)	{
						tot += (*tmpfrom++) * (*stat++);
						// tot += tmpfrom[k] * stat[k];
					}
					tmpfrom -= nstate;
					stat -= nstate;
					tot *= (1-expo);
					for (int k=0; k<nstate; k++)	{	
						(*tmpto++) = expo * (*tmpfrom++) + tot;
						// tmpto[k] = expo * tmpfrom[k] + tot;
					}
					(*tmpto) = (*tmpfrom);
					tmpto -= nstate;
					tmpfrom -= nstate;
					// tmpto[GetNstate(i)] = tmpfrom[GetNstate(i)];
				}
			}
		}
	}
}

void PoissonSubstitutionProcess::SitePropagate(int i, double** from, double** to, double time, bool condalloc)	{

	const double* stat = GetStationary(i);
	for (int j=0; j<GetNrate(i); j++)	{
		if ((! condalloc) || (ratealloc[i] == j))	{
			double* tmpfrom = from[j];
			double* tmpto = to[j];
			double expo = exp(-GetRate(i,j) * time);
			double tot = 0;
			int nstate = GetNstate(i);
			for (int k=0; k<nstate; k++)	{
				tot += (*tmpfrom++) * (*stat++);
				// tot += tmpfrom[k] * stat[k];
			}
			tmpfrom -= nstate;
			stat -= nstate;
			tot *= (1-expo);
			for (int k=0; k<nstate; k++)	{	
				(*tmpto++) = expo * (*tmpfrom++) + tot;
				// tmpto[k] = expo * tmpfrom[k] + tot;
			}
			(*tmpto) = (*tmpfrom);
			tmpto -= nstate;
			tmpfrom -= nstate;
			// tmpto[GetNstate(i)] = tmpfrom[GetNstate(i)];
		}
	}
}


//-------------------------------------------------------------------------
//	* sample substitution mappings conditional on states at nodes 
//	(CPU level 3)
//-------------------------------------------------------------------------

// root version
BranchSitePath* PoissonSubstitutionProcess::SampleRootSitePath(int site, int state)	{

	return new BranchSitePath(0,state);
}

// general version
BranchSitePath* PoissonSubstitutionProcess::SampleSitePath(int site, int stateup, int statedown, double time) 	{

	const double* stat = GetStationary(site);
	double rate = GetRate(site);
	double l = rate * time;
	double pi = stat[statedown];

	int m = 0;
	int mmax = 1000;
	
	if (!l)	{
		if (stateup != statedown)	{
			cerr << "error in PoissonSubstitutionProcess::SampleSitePath: efflength == 0 but stateup != statedown\n";
			exit(1);
		}
		return new BranchSitePath(0,statedown);
	}

	if (stateup == statedown)	{
		double fact = pi * exp(-l);
		double total = exp(-l);
		double q = rnd::GetRandom().Uniform() * (exp(-l) * (1 - pi) + pi);
		while ((m<mmax) && (total < q))	{
			m++;
			fact *= l / m;
			total += fact;
		}
		if (m == mmax)	{
			suboverflowcount ++;
		}
	}
	else	{
		double fact = pi * exp(-l);
		double total = 0;
		double q = rnd::GetRandom().Uniform() * (1 - exp(-l)) * pi;
		while ((m<mmax) && (total < q))	{
			m++;
			fact *= l / m;
			total += fact;
		}
		if (m == mmax)	{
			suboverflowcount ++;
		}
	}
	return new BranchSitePath(m,statedown);
}

//-------------------------------------------------------------------------
//	* gather sufficient statistics 
//	(CPU level 3)
//-------------------------------------------------------------------------


void PoissonSubstitutionProcess::AddMeanSuffStat(double*** down, double*** up, double*** at, double branchlength, double* meanratecount, double& meanlengthcount, double** meanprofilecount, int* nonmissing)   {

    for (int i=GetSiteMin(); i<GetSiteMax(); i++)   {
        if (ActiveSite(i))  {
            if (nonmissing[i] == 2)    {
                
                if (branchlength)   {
                    cerr << "error: non root for missing 2\n";
                    exit(1);
                }

                int nstate = GetNstate(i);
                int j = ratealloc[i];
                if (j)  {
                    cerr << "error in PoissonSubstitutionProcess::AddMeanSuffStat: rate alloc is not 0\n";
                    exit(1);
                }

                double* a = at[i][j];
                double totz[nstate];
                double tot = 0;
                for (int k=0; k<nstate; k++)    {
                    totz[k] = a[k];
                    tot += totz[k];
                }
                for (int k=0; k<nstate; k++)    {
                    totz[k] /= tot;
                }
                AddZipToTrueMeanProfileSuffStat(i,totz,meanprofilecount[i]);

            }
            else if (nonmissing[i] == 1)  {

                if (! down) {
                    cerr << "error in PoissonSubstitutionProcess::AddMeanSuffStat: non root calculation called on root\n";
                    exit(1);
                }

                if (! branchlength) {
                    cerr << "error: root for missing 1\n";
                    exit(1);
                }

                int nstate = GetNstate(i);
                const double* stat = GetStationary(i);
                double rate = GetRate(i);
                double length = rate * branchlength;
                double e = exp(-length);

                int j = ratealloc[i];
                if (j)  {
                    cerr << "error in PoissonSubstitutionProcess::AddMeanSuffStat: rate alloc is not 0\n";
                    exit(1);
                }
                double* d = down[i][j];
                double* u = up[i][j];

                double totw = 0;
                double totn = 0;
                double totz[nstate];
                for (int k=0; k<nstate; k++)    {
                    totz[k] = 0;
                }

                for (int k=0; k<nstate; k++)    {
                    for (int l=0; l<nstate; l++)    {
                        if (k == l) {
                            double w = stat[k] * d[k] * u[l] * (e + (1-e)*stat[l]);
                            totw += w;
                            totn += w * length * stat[k] / (e + (1-e)*stat[k]);
                            totz[l] += w * (1-e) * stat[k] / (e + (1-e)*stat[k]);
                        }
                        else    {
                            double w = stat[k] * d[k] * u[l] * (1-e) * stat[l];
                            totw += w;
                            if (length < 1e-8)  {
                                totn += w;
                            }
                            else    {
                                totn += w * length / (1-e);
                            }
                            totz[l] += w;
                        }
                    }
                }
                totn /= totw;
                if (isnan(totn))    {
                    cerr << "totn is nan\n";
                    exit(1);
                }
                for (int k=0; k<nstate; k++)    {
                    totz[k] /= totw;
                }

                AddZipToTrueMeanProfileSuffStat(i,totz,meanprofilecount[i]);

                meanratecount[i] += totn;
                meanlengthcount += totn;
            }
        }
    }
}

void PoissonSubstitutionProcess::AddSiteRateSuffStat(double* siteratesuffstatcount, double* siteratesuffstatbeta, double branchlength, BranchSitePath** patharray, int* nonmissing)	{
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i) && (nonmissing[i] == 1))	{
			siteratesuffstatcount[i] += patharray[i]->GetNsub();
			siteratesuffstatbeta[i] += branchlength;
		}
	}
}


void PoissonSubstitutionProcess::AddBranchLengthSuffStat(double& count, double& beta, BranchSitePath** patharray, int* nonmissing)	{
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i) && (nonmissing[i] == 1))	{
			count += patharray[i]->GetNsub();
			beta += GetRate(i);
		}
	}
}


void PoissonSubstitutionProcess::AddSiteProfileSuffStat(double** siteprofilesuffstatcount, BranchSitePath** patharray, bool root)	{
	cerr << "error: in PoissonSubstitutionProcess::AddSiteProfileSuffStat: deprecated\n";
	exit(1);
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i))	{
			if (root || patharray[i]->GetNsub())	{
				siteprofilesuffstatcount[i][GetRandomStateFromZip(i,patharray[i]->GetFinalState())]++;
			}
		}
	}
}

void PoissonSubstitutionProcess::AddZipToTrueMeanProfileSuffStat(int site, const double* p, double* q)   {

    for (int k=0; k<GetOrbitSize(site); k++)    {
        q[GetStateFromZip(site,k)] += p[k];
    }
    if (GetZipSize(site) > GetOrbitSize(site))  {
        double* pi = GetProfile(site);
        int k = GetOrbitSize(site);
        double condp = p[k] / zipstat[site][k];
        for (int l=0; l<GetDim(); l++)  {
            if (! InOrbit(site,l))  {
                q[l] += condp * pi[l];
            }
        }
    }
}

void PoissonSubstitutionProcess::ChooseTrueStates(BranchSitePath** patharray, int* nodestateup, int* nodestatedown, bool root)	{
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i))	{
			int tmp = nodestateup[i];
			if (root || patharray[i]->GetNsub())	{
				tmp = GetRandomStateFromZip(i,patharray[i]->GetFinalState());
			}
			nodestatedown[i] = tmp;
		}
	}
}

//-------------------------------------------------------------------------
//	* recomputing the equilibrium frequency profiles of the recoded process
//
//-------------------------------------------------------------------------


void PoissonSubstitutionProcess::CreateZip()	{
	if (! zipstat)	{
		zipstat = new double*[GetNsite()];
		for (int i=0; i<GetNsite(); i++)	{
			zipstat[i] = new double[GetDim()];
		}
	}
}

void PoissonSubstitutionProcess::DeleteZip()	{
	for (int i=0; i<GetNsite(); i++)	{
		delete[] zipstat[i];
	}
	delete[] zipstat;
	zipstat = 0;
}

void PoissonSubstitutionProcess::UpdateZip()	{
	
	for (int i=0; i<GetNsite(); i++)	{
		if (ActiveSite(i))	{
			UpdateZip(i);
		}
	}
}

void PoissonSubstitutionProcess::UpdateZip(int i)	{
		double total = 0;
		double* pi = GetProfile(i);
		for (int k=0; k<GetOrbitSize(i); k++)	{
			int n = GetStateFromZip(i,k);
			zipstat[i][k] = 0;
			zipstat[i][k] = pi[GetStateFromZip(i,k)];
			total += zipstat[i][k];
		}
		if (GetZipSize(i) > GetOrbitSize(i))	{
			zipstat[i][GetOrbitSize(i)] = 1-total;
		}
}

int PoissonSubstitutionProcess::GetRandomStateFromZip(int site, int zipstate)	{
	int truestate = 0;
	if ((GetZipSize(site) != GetOrbitSize(site)) && (zipstate == GetOrbitSize(site)))	{
		double v = rnd::GetRandom().Uniform();
		double u = zipstat[site][GetOrbitSize(site)] * v;
		double total = 0;
		double* pi = GetProfile(site);

		int choose = -1;
		while ((choose < GetDim()) && (total < u))	{
			choose ++;
			if (choose == GetDim())	{
				cerr << "error in getstatefromzip\n";
				cerr << choose << '\t' << GetDim() << '\n';
				cerr << total << '\n';
				cerr << v << '\t' << zipstat[site][GetOrbitSize(site)] << '\t' << u << '\n';
				cerr << total - zipstat[site][GetOrbitSize(site)] << '\n';
				double newtotal = 0;
				cerr << '\n';
				for (int k=0; k<GetDim(); k++)	{
					cerr << pi[k] << '\t';
					cerr << InOrbit(site,k) << '\n';
					if (! InOrbit(site,k))	{
						newtotal += pi[k];
					}
				}
				cerr << '\n';
				for (int k=0; k<=GetOrbitSize(site); k++)	{
					cerr << zipstat[site][k] << '\n';
				}
				cerr << "new total : " << newtotal << '\t' << newtotal - total << '\t' << total - choose << '\n';
				exit(1);
			}
			if (!InOrbit(site,choose))	{
				total += pi[choose];
			}
		}
		truestate = choose;
	}
	else	{
		truestate = GetStateFromZip(site,zipstate);
	}
	return truestate;
}


void PoissonSubstitutionProcess::UnzipBranchSitePath(BranchSitePath** patharray, int* nodestateup, int* nodestatedown){
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i))	{
			int nsub = patharray[i]->GetNsub();
			patharray[i]->nsub=0;
			double* times = new double[nsub+1];
			for(int j = 0; j < nsub; j++){
				times[j] = rnd::GetRandom().Uniform();
				for(int k = 0; k < j; k++){
					if(times[k]>times[j]	){
						times[nsub] = times[k];
						times[k] = times[j];
						times[j] = times[nsub];
					}
				}
			}
			double mem = 0;
			for(int j = 0; j < nsub; j++){
				times[nsub] = times[j];
				times[j] = times[j] - mem;
				mem = times[nsub];
			}
			times[nsub]=1-mem;

			int previousstate = nodestateup[i];
			patharray[i]->Init()->SetState(previousstate);
			double* pi = GetProfile(i);
			for(int j = 0; j < nsub-1; j++){
				int newstate = rnd::GetRandom().DrawFromDiscreteDistribution(pi, GetDim());
				if(newstate != previousstate){
				      patharray[i]->Append(newstate, times[j]);
				      previousstate=newstate;
				}
				else{
				      times[j+1]+=times[j];
				}
			}
			if(previousstate == nodestatedown[i]){
				times[nsub] += times[nsub-1];
			}
			else{
				patharray[i]->Append(nodestatedown[i], times[nsub-1]);
			}
			patharray[i]->Last()->SetRelativeTime(times[nsub]);
			delete[] times;
		}
	}
}

