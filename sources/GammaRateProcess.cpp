
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "GammaRateProcess.h"
#include "Random.h"

#include "Parallel.h"
#include <string.h>

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* GammaRateProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------


void GammaRateProcess::Create()	{
	if (! rate)	{
		RateProcess::Create();
		rate = new double[GetNsite()];
	}
}

void GammaRateProcess::Delete() 	{
	delete[] rate;
	rate = 0;
	RateProcess::Delete();
}

void GammaRateProcess::ToStream(ostream& os)	{
	os << alpha << '\n';
    for (int i=0; i<GetNsite(); i++)    {
        os << rate[i] << '\t';
    }
    os << '\n';
}

void GammaRateProcess::FromStream(istream& is)	{
	is >> alpha;
    for (int i=0; i<GetNsite(); i++)    {
        is >> rate[i];
    }
}

void GammaRateProcess::BackupRate() {
    bkalpha = alpha;
    for (int i=0; i<GetNsite(); i++)    {
        bkrate[i] = rate[i];
    }
}

void GammaRateProcess::RestoreRate() {
    alpha = bkalpha;
    for (int i=0; i<GetNsite(); i++)    {
        rate[i] = bkrate[i];
    }
}

void GammaRateProcess::SampleRate()	{
	if (! FixAlpha())	{
		alpha = 1.0;
	}
    SampleSiteRates();
}

void GammaRateProcess::SampleSiteRates()    {
    for (int i=GetSiteMin(); i<GetSiteMax(); i++)    {
        rate[i] = rnd::GetRandom().Gamma(alpha,alpha);
    }
}

void GammaRateProcess::ResampleSiteRates()    {
    for (int i=GetSiteMin(); i<GetSiteMax(); i++)    {
        rate[i] = rnd::GetRandom().Gamma(alpha + GetSiteRateSuffStatCount(i), alpha + GetSiteRateSuffStatBeta(i));
    }
}

void GammaRateProcess::GlobalResampleSiteRates()    {

    if (GetNprocs() > 1)    {
        MPI_Status stat;
        MESSAGE signal = MOVE_SRATE;
        MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
    }
    else    {
        ResampleSiteRates();
    }
}

void GammaRateProcess::PriorSampleRate()	{
	double tmpalpha = alpha;
	if (! FixAlpha())	{
		double a = meanalpha * meanalpha / varalpha;
		double b = meanalpha / varalpha;
		int count = 0;
		while ((count < 1000) && (tmpalpha < alphamin))	{
			tmpalpha = rnd::GetRandom().Gamma(a,b);
			count++;
		}
		if (count == 1000)	{
			cerr << "error in GammaRateProcess::PriorSampleRate\n";
			exit(1);
		}
	}
    SampleSiteRates();
}

void GammaRateProcess::GlobalUpdateRateHyperSuffStat()  {

	if (GetNprocs() > 1)	{
        MPI_Status stat;
        MESSAGE signal = UPDATE_HYPERRATE;
        MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
        ratesum = 0;
        logratesum = 0;
        double tmp[2];
        for(int i=1; i<GetNprocs(); i++) {
            MPI_Recv(tmp,2,MPI_DOUBLE,MPI_ANY_SOURCE,TAG1,MPI_COMM_WORLD,&stat);
            ratesum += tmp[0];
            logratesum += tmp[1];
        }
    }
    else    {
        UpdateRateHyperSuffStat();
    }
}

void GammaRateProcess::SlaveUpdateRateHyperSuffStat()   {

    UpdateRateHyperSuffStat();
    double tmp[2];
    tmp[0] = ratesum;
    tmp[1] = logratesum;
	MPI_Send(tmp,2,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

void GammaRateProcess::UpdateRateHyperSuffStat()   {
    logratesum = 0;
    ratesum = 0;
    for (int i=GetSiteMin(); i<GetSiteMax(); i++)    {
        logratesum += log(rate[i]);
        ratesum += rate[i];
    }
}

double GammaRateProcess::LogRatePrior() {
    return LogAlphaPrior() + LogSiteRatesPrior();
}

double GammaRateProcess::LogSiteRatesPrior()    {
    UpdateRateHyperSuffStat();
    return RateHyperSuffStatLogProb();
}

double GammaRateProcess::RateHyperSuffStatLogProb()    {
    return GetNsite() * (alpha*log(alpha) - rnd::GetRandom().logGamma(alpha)) + (alpha-1)*logratesum - alpha*ratesum;
}

double GammaRateProcess::LogAlphaPrior()	{
	double a = meanalpha * meanalpha / varalpha;
	double b = meanalpha / varalpha;
	double ret = a*log(b) - rnd::GetRandom().logGamma(a) + (a-1)*log(alpha) - b*alpha;
}

double GammaRateProcess::RateSuffStatLogProb()	{
	double total = 0;
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
        if (GetSiteRateSuffStatCount(i))    {
			total += GetSiteRateSuffStatCount(i) * log(rate[i]);
		}
		total -= GetSiteRateSuffStatBeta(i) * rate[i];
	}
	return total;
}

double GammaRateProcess::MoveRateParams(double tuning, int nrep)	{

	double naccepted = 0;
	for (int rep=0; rep<nrep; rep++)	{
		naccepted += MoveAlpha(tuning);
	}
	return naccepted / nrep;
}

double GammaRateProcess::MoveAlpha(double tuning)	{

	double bkalpha = alpha;
	double deltalogprob = -LogAlphaPrior() - RateHyperSuffStatLogProb();
	double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
	double e = exp(m);
    alpha *= e;
	deltalogprob += m + LogAlphaPrior() + RateHyperSuffStatLogProb();
	int accepted = (log(rnd::GetRandom().Uniform()) < deltalogprob);
	if (alpha < alphamin)	{
		accepted = 0;
	}
	if (!accepted)	{
        alpha /= e;
	}
	return ((double) accepted);
}

void GammaRateProcess::GlobalCollectSiteRates()  {

    if (GetNprocs() > 1)    {
        MPI_Status stat;
        MESSAGE signal = COLLECT_SRATE;
        MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
        for(int i=1; i<GetNprocs(); i++) {
            MPI_Recv(rate+GetProcSiteMin(i),GetProcSiteNumber(i),MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
        }
    }
}

void GammaRateProcess::SlaveCollectSiteRates()  {
	MPI_Send(rate+GetSiteMin(),GetSiteNumber(),MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

