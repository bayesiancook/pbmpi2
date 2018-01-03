
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "PoissonSiteSpecificProfileProcess.h"
#include "Random.h"
#include "Parallel.h"


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* PoissonSiteSpecificProfileProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------


void PoissonSiteSpecificProfileProcess::Create()	{
	if (! logsum)	{
		PoissonProfileProcess::Create();
		SiteSpecificProfileProcess::Create();
        logsum = new double[GetDim()];
	}
}

void PoissonSiteSpecificProfileProcess::Delete() {
	if (logsum)	{
        delete[] logsum;
        logsum = 0;
		PoissonProfileProcess::Delete();
		SiteSpecificProfileProcess::Delete();
	}
}

double PoissonSiteSpecificProfileProcess::ProfileSuffStatLogProb(int site)	{
	double total = 0;
	double priorweight = 0;
	double postweight = 0;
    const int* count = GetSiteProfileSuffStatCount(site);
	for (int k=0; k<GetDim(); k++)	{
		total += rnd::GetRandom().logGamma(dirweight[k] + count[k]) - rnd::GetRandom().logGamma(dirweight[k]);
		priorweight += dirweight[k];
		postweight += dirweight[k] + count[k];
	}
	total += rnd::GetRandom().logGamma(priorweight) - rnd::GetRandom().logGamma(postweight);
	return total;
}

void PoissonSiteSpecificProfileProcess::GlobalResampleSiteProfiles()    {

    if (GetNprocs() > 1)    {
        MPI_Status stat;
        MESSAGE signal = MOVE_SPROFILE;
        MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
    }
    else    {
        ResampleSiteProfiles();
    }
}

void PoissonSiteSpecificProfileProcess::ResampleSiteProfiles()	{
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		ResampleSiteProfile(i);
	}
}

void PoissonSiteSpecificProfileProcess::ResampleSiteProfile(int site)	{
	double total = 0;
    const int* count = GetSiteProfileSuffStatCount(site);
	for (int k=0; k<GetDim(); k++)	{
		profile[site][k] = rnd::GetRandom().sGamma(dirweight[k] + count[k]);
		if (profile[site][k] < stateps)	{
			profile[site][k] = stateps;
		}
		total += profile[site][k];
	}
	for (int k=0; k<GetDim(); k++)	{
		profile[site][k] /= total;
	}
}

/*
double PoissonSiteSpecificProfileProcess::LogStatIntPrior()	{
	double total = 0;
	for (int k=0; k<GetNsite(); k++)	{
        total += LogStatIntPrior(k);
	}
	return total;
}

double PoissonSiteSpecificProfileProcess::LogStatIntPrior(int site)	{

	double total = 0;
	int tot = 0;
	double totweight = 0;
    const int* count = GetSiteProfileSuffStatCount(site);
	for (int k=0; k<GetDim(); k++)	{
		total += rnd::GetRandom().logGamma(dirweight[k] + count[site][k]);
		total -= rnd::GetRandom().logGamma(dirweight[k]);
		totweight += dirweight[k];
		tot += count[site][k];
	}
	total -= rnd::GetRandom().logGamma(totweight + tot);
	total += rnd::GetRandom().logGamma(totweight);
	return total;
}
*/

void PoissonSiteSpecificProfileProcess::GlobalUpdateProfileHyperSuffStat()  {

	if (GetNprocs() > 1)	{
        MPI_Status stat;
        MESSAGE signal = UPDATE_HYPERPROFILE;
        MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
        for (int k=0; k<GetDim(); k++)  {
            logsum[k] = 0;
        }
        double tmp[GetDim()];
        for(int i=1; i<GetNprocs(); i++) {
            MPI_Recv(tmp,GetDim(),MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
            for (int k=0; k<GetDim(); k++)  {
                logsum[k] += tmp[k];
            }
        }
    }
    else    {
        UpdateProfileHyperSuffStat();
    }
}

void PoissonSiteSpecificProfileProcess::SlaveUpdateProfileHyperSuffStat()   {

    UpdateProfileHyperSuffStat();
    double tmp[GetDim()];
    for (int k=0; k<GetDim(); k++)  {
        tmp[k] = logsum[k];
    }
	MPI_Send(tmp,GetDim(),MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

void PoissonSiteSpecificProfileProcess::UpdateProfileHyperSuffStat()   {

    for (int k=0; k<GetDim(); k++)  {
        logsum[k] = 0;
    }
    for (int i=GetSiteMin(); i<GetSiteMax(); i++)   {
        for (int k=0; k<GetDim(); k++)  {
            logsum[k] += log(profile[i][k]);
        }
    }
}

double PoissonSiteSpecificProfileProcess::ProfileHyperSuffStatLogProb() {

    double totlogprob = 0;
    double totweight = 0;
    for (int k=0; k<GetDim(); k++)  {
        totlogprob -= rnd::GetRandom().logGamma(dirweight[k]);
        totweight += dirweight[k];
    }
    totlogprob += rnd::GetRandom().logGamma(totweight);
    totlogprob *= GetNsite();
    for (int k=0; k<GetDim(); k++)  {
        totlogprob += (dirweight[k]-1)*logsum[k];
    }
    return totlogprob;
}

double PoissonSiteSpecificProfileProcess::MoveDirWeights(double tuning, int nrep)	{

	double naccepted = 0;
	for (int rep=0; rep<nrep; rep++)	{
		for (int k=0; k<GetDim(); k++)	{
			double deltalogprob = - LogHyperPrior() - ProfileHyperSuffStatLogProb();
			double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
			double e = exp(m);
			dirweight[k] *= e;
			deltalogprob += LogHyperPrior() + ProfileHyperSuffStatLogProb();
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

