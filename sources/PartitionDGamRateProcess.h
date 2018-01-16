
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef PARTDGAMRATE_H
#define PARTDGAMRATE_H

#include "RateProcess.h"

class PartitionDGamRateProcess : public virtual RateProcess {

	public:

	PartitionDGamRateProcess() : Ncat(0), rate(0), fixalpha(false), meanalpha(1), varalpha(1), alphamin(0.0), pinv(0), withpinv(false), fixpinv(false), meanpinv(0.5), invconcpinv(2) {}

	virtual ~PartitionDGamRateProcess() {}

	double GetAlpha(int part) {return alpha[part];}
	double GetPinv(int part) {return pinv[part];}

	double GetMeanAlpha()	{
		double mean = 0;
		for (int part=0; part<Npart; part++)	{
			mean += alpha[part];
		}
		mean /= Npart;
		return mean;
	}

	double GetVarAlpha()	{
		double mean = 0;
		double var = 0;
		for (int part=0; part<Npart; part++)	{
			mean += alpha[part];
			var += alpha[part]*alpha[part];
		}
		mean /= Npart;
		var /= Npart;
		var -= mean*mean;
		return var;
	}

	double GetMeanPinv()	{
		double mean = 0;
		for (int part=0; part<Npart; part++)	{
			mean += pinv[part];
		}
		mean /= Npart;
		return mean;
	}

	double GetVarPinv()	{
		double mean = 0;
		double var = 0;
		for (int part=0; part<Npart; part++)	{
			mean += pinv[part];
			var += pinv[part]*pinv[part];
		}
		mean /= Npart;
		var /= Npart;
		var -= mean*mean;
		return var;
	}

	int GetNrate(int site)	{
		return Ncat;
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
		for (int part=0; part<Npart; part++)	{
			bkalpha[part] = alpha[part];
			bkpinv[part] = pinv[part];
		}
	}

	virtual void RestoreRate() {
		SetRateParams(bkalpha,bkpinv);
	}

	double GetRate(int site, int cat = 0)	{
		// cat should be == 0
		int part = partalloc[site];
		if (SumOverRateAllocations())	{
			return rate[part][cat];
		}
		return rate[part][ratealloc[site]];
	}

	double GetRateWeight(int site, int cat)	{
		int part = partalloc[site];
		if (SumOverRateAllocations())	{
			if (withpinv)	{
				if (! cat)	{
					return pinv[part];
				}
				else	{
					return (1.0 - pinv[part]) / (Ncat-1);
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

	double GetPriorMeanRate(int part)	{
		double total = 0;
		for (int k=0; k<GetNcat(); k++)	{
			total += GetRateWeight(0,k) * rate[part][k];
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
	double MoveAlpha(int part, double tuning);
	double MovePinv(double tuning);
	double MovePinv(int part, double tuning);
	
	void SetRateParams(int part, double inalpha, double inpinv)	{
		alpha[part] = inalpha;
		pinv[part] = inpinv;
		UpdateDiscreteCategories(part);
	}

	void SetRateParams(double* inalpha, double* inpinv)	{
		for (int part=0; part<Npart; part++)	{
			alpha[part] = inalpha[part];
			pinv[part] = inpinv[part];
		}
		UpdateDiscreteCategories();
	}

	void ToStream(ostream& os);
	void FromStream(istream& is);

	virtual void UpdateRateSuffStat();
	const double* GetRateSuffStatCount(int part) {return ratesuffstatcount[part];}
	const double* GetRateSuffStatBeta(int part) {return ratesuffstatbeta[part];}

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
	double LogRatePrior(int part);

	void GlobalUpdateRateSuffStat();
	void SlaveUpdateRateSuffStat();

	double RateSuffStatLogProb();
	double RateSuffStatLogProb(int part);

	void UpdateDiscreteCategories();
	void UpdateDiscreteCategories(int part);
	
	double** rate;
	double* alpha;
	double* bkalpha;
	
	double** ratesuffstatcount;
	double** ratesuffstatbeta;
	double* allocratesuffstatcount;
	double* allocratesuffstatbeta;
	double* Ninv;

	int Ncat;
	bool fixalpha;

	public:

	double meanalpha;
	double varalpha;

	double alphamin;

	double* pinv;
	double* bkpinv;
	bool withpinv;
	bool fixpinv;
	double meanpinv;
	double invconcpinv;
};

#endif

