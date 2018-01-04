
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

	DGamRateProcess() : Ncat(0), rate(0), fixalpha(false), meanalpha(1), varalpha(1), alphamin(0.0), pinv(0.01), withpinv(false), fixpinv(false), meanpinv(0.5), invconcpinv(2) {}

	virtual ~DGamRateProcess() {}

	double GetAlpha() {return alpha;}
	double GetPinv() {return pinv;}

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

	void SetFixAlpha(bool in)	{
		fixalpha = in;
	}

	bool FixAlpha()	{
		return fixalpha;
	}

	void SetFixPinv(bool in)	{
		fixpinv = in;
	}

	bool FixPinv()	{
		return fixpinv;
	}

	virtual void BackupRate() {
		bkalpha = alpha;
		bkpinv = pinv;
	}

	virtual void RestoreRate() {
		SetRateParams(bkalpha,bkpinv);
	}

	double GetRate(int site, int cat = 0)	{
		// cat should be == 0
		if (SumOverRateAllocations())	{
			return rate[cat];
		}
		return rate[ratealloc[site]];
	}

	double GetRateWeight(int site, int cat)	{
		if (SumOverRateAllocations())	{
			if (withpinv)	{
				if (! cat)	{
					return pinv;
				}
				else	{
					return (1.0 - pinv) / (Ncat-1);
				}
			}
			else	{
				return 1.0/Ncat;
			}
		}
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
		double total = 0;
		for (int k=0; k<GetNcat(); k++)	{
			total += GetRateWeight(0,k) * rate[k];
		}
		return total;
	}

	virtual double Move(double tuning = 1, int nrep = 1)	{
		GlobalUpdateSiteRateSuffStat();
		chronorate.Start();
		return MoveRateParams(tuning, nrep);
		chronorate.Stop();
	}

	double NonMPIMove(double tuning = 1, int nrep = 1)	{
		UpdateSiteRateSuffStat();
		NonMPIMoveRateParams(tuning, nrep);
		return 1;
	}

	// uses suffisicent stats
	double MoveRateParams(double tuning, int nrep);
	double NonMPIMoveRateParams(double tuning, int nrep);

	double MoveAlpha(double tuning);
	double MovePinv(double tuning);
	
	void SetRateParams(double inalpha, double inpinv)	{
		alpha = inalpha;
		pinv = inpinv;
		UpdateDiscreteCategories();
	}

	void ToStream(ostream& os);
	void FromStream(istream& is);

	virtual void UpdateRateSuffStat();
	const double* GetRateSuffStatCount() {return ratesuffstatcount;}
	const double* GetRateSuffStatBeta() {return ratesuffstatbeta;}

	protected:

	void SetNcat(int inncat)	{
		if (withpinv)	{
			Ncat = inncat + 1;
		}
		else	{
			Ncat = inncat;
		}
	}

	virtual void Create();
	virtual void Delete();

	virtual void SampleRate();
	virtual void PriorSampleRate();
	virtual double LogRatePrior();

	void GlobalUpdateRateSuffStat();
	void SlaveUpdateRateSuffStat();
	double RateSuffStatLogProb();

	void UpdateDiscreteCategories();
	
	double* rate;
	double alpha;
	double bkalpha;
	
	double* ratesuffstatcount;
	double* ratesuffstatbeta;
	double Ninv;

	int Ncat;
	bool fixalpha;

	public:

	double meanalpha;
	double varalpha;

	double alphamin;

	double pinv;
	double bkpinv;
	bool withpinv;
	bool fixpinv;
	double meanpinv;
	double invconcpinv;
};

#endif

