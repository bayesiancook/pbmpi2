
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef MIXTUREPROFILE_H
#define MIXTUREPROFILE_H

#include <cmath>
#include "ProfileProcess.h"

// general superclass for all finite process mixtures on site-specific profiles
class MixtureProfileProcess: public virtual ProfileProcess	{

	public:

	MixtureProfileProcess() : profile(0), nmodemax(5000), mtryalloc(0) {}
	virtual ~MixtureProfileProcess(){}

	double* GetProfile(int site)	{
		if (alloc[site] == -1)	{
			cerr << "error in MixtureProfileProcess::GetProfile(int site): alloc is -1\n";
			cerr << "site : " << site << '\n';
			exit(1);
		}
		return profile[alloc[site]];
	}

	// number of components of the mixture
	int GetNcomponent() { return Ncomponent;}

	virtual int GetNDisplayedComponent()	{
		return Ncomponent;
	}

	int GetNOccupiedComponent() {
		int n = 0;
		for (int k=0; k<Ncomponent; k++)	{
			if (occupancy[k])	{
				n++;
			}
		}
		return n;
	}

	// for allocation purpose
	virtual int GetNmodeMax() {
		/*
		if (Npart)	{
			return Npart;
		}
		*/
		return GetNsite() > nmodemax ? nmodemax : GetNsite();
	}
	virtual void SetNmodeMax(int n) {nmodemax = n;}

	void GlobalChooseMultipleTryAlloc();
	void SlaveChooseMultipleTryAlloc();
	void ChooseMultipleTryAlloc();

	void GlobalActivateSumOverComponents();
	void SlaveActivateSumOverComponents();
	void ActivateSumOverComponents();

	double GetStatEnt();

	protected:

	//------
	// priors
	//------

	double LogProfilePrior();
	virtual double LogStatPrior();
	virtual double LogStatPrior(int cat);

	virtual double LogAllocPrior()	{
		return 0;
	}

	//------
	// sampling from prior
	//------

	virtual void SampleProfile();
	virtual void PriorSampleProfile();
	virtual void SampleStat();
	virtual void SampleStat(int cat);

	virtual void SampleAlloc() = 0;

	double ResampleEmptyProfiles();

	//------
	// log likelihoods based on sufficient statistics (substitution mappings)
	//------

	double ProfileSuffStatLogProb();
	// the component suff stat log prob is yet to be implemented in subclasses
	virtual double ProfileSuffStatLogProb(int cat) = 0;

	// suffstat lnL of site <site> when allocated to component <cat>
	virtual double LogStatProb(int site, int cat) = 0;

	//------
	// moves
	//------

	// the following MPI move
	// assumes that master and all slaves are in sync concerning siteprofilesuffstats
	// (which will be the case upon calling GlobalUpdateSiteProfileSuffStat)
	virtual double GlobalMoveProfile(double tuning = 1, int n = 1, int nrep = 1);
	virtual void SlaveMoveProfile();
	double MoveProfile(double tuning = 1, int n = 1, int nrep = 1);
	double MoveProfile(int cat, double tuning, int n, int nrep);
	virtual double GlobalSMCAddSites();

	//------
	// create and delete
	//------

	virtual void Create();
	virtual void Delete();

	virtual void CreateComponent(int k) = 0;
	virtual void DeleteComponent(int k) = 0;
	virtual void UpdateComponent(int k) = 0;

	virtual void UpdateComponents()	{
		for (int k=0; k<GetNcomponent(); k++)	{
			UpdateComponent(k);
		}
	}

	virtual void AddSite(int site, int cat)	{
		alloc[site] = cat;
		occupancy[cat]++;
	}

	virtual void RemoveSite(int site, int cat)	{
		if (cat != -1)	{
			alloc[site] = -1;
			occupancy[cat]--;
		}
	}

	virtual double GetWeight(int cat)	{
		cerr << "in MixtureProfileProcess::GetWeight\n";
		exit(1);
		return 1;
	}

	virtual void SwapComponents(int cat1, int cat2);

	void UpdateOccupancyNumbers();

	double GetMeanStationaryEntropy() {return GetStatEnt();}
	double GetSiteStationaryEntropy(int site) {return GetStatEnt(alloc[site]);}

	// summary statistic: mean entropy over all profiles
	double GetStatEnt(int k);

	void RenormalizeProfiles();

	virtual double GetAllocEntropy()	{
		double total = 0;
		UpdateOccupancyNumbers();
		double totnsite = 0;
		for (int k=0; k<GetNcomponent(); k++)	{
			totnsite += occupancy[k];
		}
		for (int k=0; k<GetNcomponent(); k++)	{
			double tmp = ((double) occupancy[k]) / totnsite;
			if (tmp)	{
				total -= tmp * log(tmp);
			}
		}
		return total;
	}

	double** profile;
	double* allocprofile;
	int* alloc;
	int* occupancy;
	int Ncomponent;
	double* logstatprior;
	double* profilesuffstatlogprob;

	// sumovercomponents
	int** mtryalloc;
	double** mtryweight;

	int nmodemax;

	// 0 : flexible
	// 1 : rigid

	Chrono totchrono;
	Chrono profilechrono;
	Chrono incchrono;
};

#endif

