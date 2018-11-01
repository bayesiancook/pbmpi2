
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "PoissonMixtureProfileProcess.h"
#include "Random.h"
#include "Parallel.h"


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* PoissonMixtureProfileProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------


void PoissonMixtureProfileProcess::Create()	{
	if (! profilesuffstatcount)	{
		PoissonProfileProcess::Create();
		MixtureProfileProcess::Create();
		allocprofilesuffstatcount = new double[GetNmodeMax()*GetDim()];
		profilesuffstatcount  = new double*[GetNmodeMax()];
		alloctmpprofilesuffstatcount = new double[GetNmodeMax()*GetDim()];
		tmpprofilesuffstatcount  = new double*[GetNmodeMax()];
		for (int i=0; i<GetNmodeMax(); i++)	{
			profilesuffstatcount[i] = allocprofilesuffstatcount + i*GetDim();
			tmpprofilesuffstatcount[i] = alloctmpprofilesuffstatcount + i*GetDim();
		}
	}
}

void PoissonMixtureProfileProcess::Delete() {
	if (profilesuffstatcount)	{
		delete[] profilesuffstatcount;
		delete[] allocprofilesuffstatcount;
		delete[] tmpprofilesuffstatcount;
		delete[] alloctmpprofilesuffstatcount;
		profilesuffstatcount = 0;
		PoissonProfileProcess::Delete();
		MixtureProfileProcess::Delete();
	}
}

void PoissonMixtureProfileProcess::GlobalUpdateModeProfileSuffStat()	{

	if (GetNprocs() > 1)	{
		MESSAGE signal = UPDATE_MPROFILE;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Reduce(alloctmpprofilesuffstatcount,allocprofilesuffstatcount,Ncomponent*GetDim(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Bcast(allocprofilesuffstatcount,Ncomponent*GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	}

	else	{
		UpdateModeProfileSuffStat();
	}
}

void PoissonMixtureProfileProcess::SlaveUpdateModeProfileSuffStat()	{

	UpdateModeProfileSuffStat();
	for (int k=0; k<Ncomponent*GetDim(); k++)	{
		alloctmpprofilesuffstatcount[k] = allocprofilesuffstatcount[k];
	}
	MPI_Reduce(alloctmpprofilesuffstatcount,allocprofilesuffstatcount,Ncomponent*GetDim(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(allocprofilesuffstatcount,Ncomponent*GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);
}

void PoissonMixtureProfileProcess::UpdateModeProfileSuffStat()	{

	for (int i=0; i<GetNcomponent(); i++)	{
		for (int k=0; k<GetDim(); k++)	{
			profilesuffstatcount[i][k] = 0;
		}
	}
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i))	{
			const double* count = GetSiteProfileSuffStatCount(i);
			int cat = alloc[i];
			for (int k=0; k<GetDim(); k++)	{
				profilesuffstatcount[cat][k] += count[k];
			}
		}
	}
}

double PoissonMixtureProfileProcess::ProfileSuffStatLogProb(int cat)	{

    if (Nstatcomp > 1)  {
        cerr << "in PoissonMixtureProfileProcess::ProfileSuffStatLogProb\n";
        exit(1);
    }

	double total = 0;

    if (Nstatcomp == 1) {
        double priorweight = 0;
        double postweight = 0;
        for (int k=0; k<GetDim(); k++)	{
            total += rnd::GetRandom().logGamma(dirweight[0][k] + profilesuffstatcount[cat][k]) - rnd::GetRandom().logGamma(dirweight[0][k]);
            priorweight += dirweight[0][k];
            postweight += dirweight[0][k] + profilesuffstatcount[cat][k];
        }
        total += rnd::GetRandom().logGamma(priorweight) - rnd::GetRandom().logGamma(postweight);
    }
    else    {
    }
	return total;
}

double PoissonMixtureProfileProcess::DiffLogSampling(int cat, int site)	{

    if (Nstatcomp > 1)  {
        cerr << "error in PoissonMixtureProfileProcess::DiffLogsampling\n";
        exit(1);
    }

	const double* nsub = GetSiteProfileSuffStatCount(site);
	double* catnsub = profilesuffstatcount[cat];
	double totalsub = 0;
	double priorweight = 0;
	int grandtotal = 0;
	for (int k=0; k<GetDim(); k++)	{
		totalsub += nsub[k];
		priorweight += dirweight[0][k];
		grandtotal += catnsub[k];
	}
	
	double total = 0;
	for (int j=0; j< totalsub; j++)	{
		total -= log(priorweight + grandtotal + j);
	}
	for (int k=0; k<GetDim(); k++)	{
		for (int j=0; j< nsub[k]; j++)	{
			total += log(dirweight[0][k] + catnsub[k] + j);
		}
	}
	return total;
}

double PoissonMixtureProfileProcess::LogStatProb(int site, int cat)	{
	const double* nsub = GetSiteProfileSuffStatCount(site);
	double total = 0;
	for (int k=0; k<GetDim(); k++)	{
		total += nsub[k] * log(profile[cat][k]);
	}
	return total;
}

double PoissonMixtureProfileProcess::MoveProfile()	{
	for (int i=0; i<GetNcomponent(); i++)	{
		MoveProfile(i);
	}
	return 1;
}

int PoissonMixtureProfileProcess::SamplePriorAllocation(int cat)   {

    double logprob[Nstatcomp];
    int max = 0;
    for (int i=0; i<Nstatcomp; i++) {
        double priorweight = 0;
        double postweight = 0;
        double total = 0;
        for (int k=0; k<GetDim(); k++)	{
            total += rnd::GetRandom().logGamma(dirweight[i][k] + profilesuffstatcount[cat][k]) - rnd::GetRandom().logGamma(dirweight[i][k]);
            priorweight += dirweight[i][k];
            postweight += dirweight[i][k] + profilesuffstatcount[cat][k];
        }
        total += rnd::GetRandom().logGamma(priorweight) - rnd::GetRandom().logGamma(postweight);

        if ((!i) || (max < total))  {
            max = total;
        }
        logprob[i] = total;
    }

    double cumul[Nstatcomp];
    double total = 0;
    for (int i=0; i<Nstatcomp; i++) {
        double tmp = exp(logprob[i] - max);
        total += tmp;
        cumul[i] = total;
    }

    double u = total * rnd::GetRandom().Uniform();
    int i = 0;
    while ((i<Nstatcomp) && (u>cumul[i]))   {
        i++;
    }
    if (i == Nstatcomp) {
        cerr << "error in PoissonMixtureProfileProcess::SamplePriorCat: overflow\n";
        exit(1);
    }
    return i;
}

double PoissonMixtureProfileProcess::MoveProfile(int cat)	{
    int priorcat = 0;

    if (Nstatcomp > 1)  {
        priorcat = SamplePriorAllocation(cat);
    }
    
	double total = 0;
	for (int k=0; k<GetDim(); k++)	{
		profile[cat][k] = rnd::GetRandom().sGamma(dirweight[priorcat][k] + profilesuffstatcount[cat][k]);
		if (profile[cat][k] < stateps)	{
			profile[cat][k] = stateps;
		}
		total += profile[cat][k];
	}
	for (int k=0; k<GetDim(); k++)	{
		profile[cat][k] /= total;
	}
	return 1;
}

/*
double PoissonMixtureProfileProcess::LogStatIntPrior()	{
	double total = 0;
	for (int k=0; k<GetNcomponent(); k++)	{
		if (occupancy[k])	{
			total += LogStatIntPrior(k);
		}
	}
	return total;
}

double PoissonMixtureProfileProcess::LogStatIntPrior(int cat)	{

	double total = 0;
	double tot = 0;
	double totweight = 0;
	for (int k=0; k<GetDim(); k++)	{
		total += rnd::GetRandom().logGamma(dirweight[0][k] + profilesuffstatcount[cat][k]);
		total -= rnd::GetRandom().logGamma(dirweight[0][k]);
		totweight += dirweight[0][k];
		tot += profilesuffstatcount[cat][k];
	}
	total -= rnd::GetRandom().logGamma(totweight + tot);
	total += rnd::GetRandom().logGamma(totweight);
	return total;
}
	
double PoissonMixtureProfileProcess::MoveDirWeights(double tuning, int nrep)	{

	cerr << "in PoissonMixtureProfileProcess::MoveDirWeights\n";
	exit(1);
	UpdateOccupancyNumbers();
	UpdateModeProfileSuffStat();
	double naccepted = 0;
	for (int rep=0; rep<nrep; rep++)	{
		for (int k=0; k<GetDim(); k++)	{
			double deltalogprob = - LogHyperPrior() - LogStatIntPrior();
			double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
			double e = exp(m);
			dirweight[0][k] *= e;
			deltalogprob += LogHyperPrior() + LogStatIntPrior();
			deltalogprob += m;
			int accepted = (log(rnd::GetRandom().Uniform()) < deltalogprob);
			if (accepted)	{
				naccepted++;
			}
			else	{
				dirweight[0][k] /= e;
			}
		}
	}
	SampleStat();
	return naccepted / nrep / GetDim();
}
*/

void PoissonMixtureProfileProcess::RemoveSite(int site, int cat)	{
	if (cat != -1)	{
		alloc[site] = -1;
		occupancy[cat] --;
		if (activesuffstat)	{
			const double* nsub = GetSiteProfileSuffStatCount(site);
			double* catnsub = profilesuffstatcount[cat];
			for (int k=0; k<GetDim(); k++)	{
				catnsub[k] -= nsub[k];
			}
		}
	}
}

void PoissonMixtureProfileProcess::AddSite(int site, int cat)	{
	alloc[site] = cat;
	occupancy[cat] ++;
	UpdateZip(site);
	if (activesuffstat)	{
		const double* nsub = GetSiteProfileSuffStatCount(site);
		double* catnsub = profilesuffstatcount[cat];
		for (int k=0; k<GetDim(); k++)	{
			catnsub[k] += nsub[k];
		}
	}
}

void PoissonMixtureProfileProcess::SwapComponents(int cat1, int cat2)	{

	MixtureProfileProcess::SwapComponents(cat1,cat2);
	for (int k=0; k<GetDim(); k++)	{
		double tmp = profilesuffstatcount[cat1][k];
		profilesuffstatcount[cat1][k] = profilesuffstatcount[cat2][k];
		profilesuffstatcount[cat2][k] = tmp;
	}
}

