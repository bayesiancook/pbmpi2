
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef DGAMRATE_H
#define DGAMRATE_H

#include "RateProcess.h"

class DGamRateProcess : public virtual RateProcess {

	public:

	DGamRateProcess() : Ncat(0), rate(0) {}
	virtual ~DGamRateProcess() {}

	double GetAlpha() {return alpha;}

	int GetNrate(int site)	{
		return Ncat;
		/*
		if (SumOverRateAllocations())	{
			return Ncat;
		}
		return 1;
		*/
	}

	virtual int GetNrate() {return GetNcat();}

	int GetNcat() {return Ncat;}

	double GetRate(int site, int cat = 0)	{
		// cat should be == 0
		if (SumOverRateAllocations())	{
			return rate[cat];
		}
		return rate[alloc[site]];
	}

	double GetRateWeight(int site, int cat)	{
		if (SumOverRateAllocations())	{
			return 1.0/Ncat;
		}
		return 1.0;
	}

	void ActivateSumOverRateAllocations() {
		condflag = false;
	}

	void InactivateSumOverRateAllocations(int* ratealloc) {
		for (int i=0; i<GetNsite(); i++)	{
			if (ActiveSite(i))	{
				alloc[i] = ratealloc[i];
			}
		}
		condflag = true;
	}

	void SiteActivateSumOverRateAllocation(int site) {
		condflag = false;
	}

	void SiteInactivateSumOverRateAllocation(int site, int ratealloc) {
		alloc[site] = ratealloc;
		condflag = true;
	}

	double GetPriorMeanRate()	{
		double total = 0;
		for (int k=0; k<GetNcat(); k++)	{
			total += rate[k];
		}
		return total / GetNcat();
	}

	double Move(double tuning = 1, int nrep = 1)	{
		GlobalUpdateSiteRateSuffStat();
		chronorate.Start();
		return MoveAlpha(tuning, nrep);
		chronorate.Stop();
	}

	double NonMPIMove(double tuning = 1, int nrep = 1)	{
		UpdateSiteRateSuffStat();
		NonMPIMoveAlpha(tuning, nrep);
		return 1;
	}

	// uses suffisicent stats
	double MoveAlpha(double tuning, int nrep);
	double NonMPIMoveAlpha(double tuning, int nrep);
	
	void SetAlpha(double inalpha)	{
		alpha = inalpha;
		UpdateDiscreteCategories();
	}

	void ToStream(ostream& os);
	void FromStream(istream& is);

	protected:

	void SetNcat(int inncat)	{
		Ncat = inncat;
	}

	void Create();
	void Delete();

	void SampleRate();
	void PriorSampleRate();
	double LogRatePrior();


	void GlobalUpdateRateSuffStat();
	virtual void SlaveUpdateRateSuffStat();
	void UpdateRateSuffStat();
	double RateSuffStatLogProb();

	void UpdateDiscreteCategories();
	
	double* rate;
	double alpha;
	int* alloc;
	int* ratesuffstatcount;
	double* ratesuffstatbeta;

	int Ncat;
};

#endif

