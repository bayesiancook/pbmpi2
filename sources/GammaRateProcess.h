
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef CGAMRATE_H
#define CGAMRATE_H

#include "RateProcess.h"

class GammaRateProcess : public virtual RateProcess {

	public:

	GammaRateProcess() : rate(0), fixalpha(false), meanalpha(1), varalpha(1), alphamin(0.0) {}

	virtual ~GammaRateProcess() {}

	double GetAlpha() {return alpha;}

	int GetNrate(int site)	{
        return 1;
    }

	void SetFixAlpha(bool in)	{
		fixalpha = in;
	}

	bool FixAlpha()	{
		return fixalpha;
	}

    void BackupRate();
    void RestoreRate();

	double GetRate(int site, int cat = 0)	{
        if (cat)    {
            cerr << "error in GammaRateProcess: GetRate with non null category\n";
            exit(1);
        }
		return rate[site];
	}

	double GetRateWeight(int site, int cat)	{
        return 1.0;
    }

	void ActivateSumOverRateAllocations() {
		condflag = false;
	}

	void InactivateSumOverRateAllocations()	{
		condflag = true;
	}

	void SiteActivateSumOverRateAllocation(int site) {
		condflag = false;
	}

	void SiteInactivateSumOverRateAllocation(int site)	{
		condflag = true;
	}

	double GetPriorMeanRate()	{
        return 1.0;
    }

	virtual double Move()   {
        if (GetNprocs() > 1)    {
            MPIMove();
        }
        else    {
            NonMPIMove();
        }
        return 1.0;
    }

    virtual double MPIMove()    {
		GlobalUpdateSiteRateSuffStat();
        GlobalResampleSiteRates();
        GlobalUpdateRateHyperSuffStat();
		MoveRateParams(1,100);
		MoveRateParams(0.3,100);
		MoveRateParams(0.1,100);
        return 1.0;
	}

	double NonMPIMove() {
		UpdateSiteRateSuffStat();
        ResampleSiteRates();
        UpdateRateHyperSuffStat();
        /*
		MoveRateParams(1,100);
		MoveRateParams(0.3,100);
		MoveRateParams(0.1,100);
        */
		return 1.0;
	}

	// uses sufficient stats
    double MoveRateParams(double tuning, int nrep);
	double MoveAlpha(double tuning);

	void ToStream(ostream& os);
	void FromStream(istream& is);

	protected:

	virtual void Create();
	virtual void Delete();

	virtual void SampleRate();
	virtual void PriorSampleRate();

    void SampleSiteRates();
    void ResampleSiteRates();
    void GlobalResampleSiteRates();

    void GlobalUpdateRateHyperSuffStat();
    void SlaveUpdateRateHyperSuffStat();
    void UpdateRateHyperSuffStat();

    void GlobalCollectSiteRates();
    void SlaveCollectSiteRates();

	virtual double LogRatePrior();
    double LogAlphaPrior();
    double LogSiteRatesPrior();
    double RateHyperSuffStatLogProb();
	double RateSuffStatLogProb();

	double* rate;
    double* bkrate;
	double alpha;
	double bkalpha;
	
    double ratesum;
    double logratesum;

	bool fixalpha;

	public:

	double meanalpha;
	double varalpha;

	double alphamin;
};

#endif

