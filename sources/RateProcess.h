
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

	virtual int GetNrate(int site) = 0;
	virtual double GetRate(int site, int cat = 0) = 0;
	virtual double GetRateWeight(int site, int cat) = 0;
	double GetMeanRate();
	virtual double GetPriorMeanRate() = 0;
	virtual double GetAlpha() {return 1;}

	virtual void SiteActivateSumOverRateAllocation(int site) = 0;
	virtual void SiteInactivateSumOverRateAllocation(int site, int ratealloc) = 0;
	virtual void ActivateSumOverRateAllocations() = 0;
	virtual void InactivateSumOverRateAllocations(int* ratealloc) = 0;
	bool SumOverRateAllocations() {return ! condflag;}

	virtual double LogRatePrior() = 0;
	virtual void SampleRate() = 0;
	virtual void PriorSampleRate() {};

	virtual void ToStream(ostream& os) = 0;
	virtual void FromStream(istream& is) = 0;

	protected:

	// abstract classes will be implemented in phyloprocess
	virtual void GlobalUpdateSiteRateSuffStat() = 0;
	virtual void SlaveUpdateSiteRateSuffStat() = 0;

	virtual void UpdateSiteRateSuffStat() = 0;
	virtual double GetSiteRateSuffStatBeta(int site) = 0;
	virtual int GetSiteRateSuffStatCount(int site) = 0;

	void Create() {}
	void Delete() {}

	bool condflag;

	Chrono chronorate;
};

#endif

