
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef RATE_H
#define RATE_H

#include <iostream>
using namespace std;

#include "Chrono.h"
#include "MPIModule.h"

class RateProcess : public virtual MPIModule {

	public:

	RateProcess() : condflag(false) {}
	virtual ~RateProcess() {}

	virtual int GetNrate() {return 1;}

	virtual int GetNrate(int site)	{
		cerr << "in RateProcess::GetNrate\n";
		exit(1);
		return 0;
	}

	virtual double GetRate(int site, int cat = 0)	{
		cerr << "in RateProcess::GetRate\n";
		exit(1);
		return 0;
	}

	virtual double GetRateWeight(int site, int cat)	{
		cerr << "in RateProcess::GetRateWeight\n";
		exit(1);
		return 0;
	}

	double GetMeanRate();
	virtual double GetPriorMeanRate()	{
		cerr << "in RateProcess::GetPriorMeanRate\n";
		exit(1);
	}

	virtual void BackupRate() {
		cerr << "in RateProcess::BackupRate\n";
		exit(1);
	}

	virtual void RestoreRate() {
		cerr << "in RateProcess::RestoreRate\n";
		exit(1);
	}

	virtual double GetAlpha() {return 1;}

	virtual void SiteActivateSumOverRateAllocation(int site)	{
		cerr << "in RateProcess::SiteActivateSumOverRateAlloc\n";
		exit(1);
	}
	virtual void SiteInactivateSumOverRateAllocation(int site, int ratealloc)	{
		cerr << "in RateProcess::SiteINactivateSumOverRateAlloc\n";
		exit(1);
	}
	virtual void ActivateSumOverRateAllocations()	{
		cerr << "in RateProcess::ActivateSumOverRateAlloc\n";
		exit(1);
	}
	virtual void InactivateSumOverRateAllocations(int* ratealloc)	{
		cerr << "in RateProcess::InactivateSumOverRateAlloc\n";
		exit(1);
	}
	bool SumOverRateAllocations() {return ! condflag;}

	virtual double LogRatePrior()	{
		cerr << "in RateProcess::LogRatePrior\n";
		exit(1);
		return 0;
	}
	virtual void SampleRate()	{
		cerr << "in RateProcess::SampleRate\n";
		exit(1);
	}
	virtual void PriorSampleRate() {};

	virtual void ToStream(ostream& os)	{
		cerr << "in RateProcess::ToStream\n";
		exit(1);
	}
	virtual void FromStream(istream& is)	{
		cerr << "in RateProcess::FromStream\n";
		exit(1);
	}

	virtual void UpdateSiteRateSuffStat()	{
		cerr << "in RateProcess::UpdateSiteRateSuffStat\n";
		exit(1);
	}

	protected:

	// abstract classes will be implemented in phyloprocess
	virtual void GlobalUpdateSiteRateSuffStat()	{
		cerr << "in RateProcess::GlobalUpdateSiteRateSuffStat\n";
		exit(1);
	}
	virtual void SlaveUpdateSiteRateSuffStat()	{
		cerr << "in RateProcess::SlaveUpdateSiteRateSuffStat\n";
		exit(1);
	}

	virtual double GetSiteRateSuffStatBeta(int site)	{
		cerr << "in RateProcess::GetSiteRateSuffStatBeta\n";
		exit(1);
	}
	virtual int GetSiteRateSuffStatCount(int site)	{
		cerr << "in RateProcess::GetSiteRateSuffStatCount\n";
		exit(1);
	}

	void Create() {}
	void Delete() {}

	bool condflag;

	Chrono chronorate;
};

#endif

