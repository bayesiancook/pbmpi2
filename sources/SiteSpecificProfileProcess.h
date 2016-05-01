
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef SSPROFILE_H
#define SSPROFILE_H

#include <cmath>
#include "ProfileProcess.h"

class SiteSpecificProfileProcess: public virtual ProfileProcess	{

	public:

	SiteSpecificProfileProcess() : profile(0) {}
	virtual ~SiteSpecificProfileProcess(){}

	double* GetProfile(int site)	{
		return profile[site];
	}

	protected:

	//------
	// priors
	//------

	double LogProfilePrior();
	virtual double LogStatPrior();
	virtual double LogStatPrior(int site);

	//------
	// sampling from prior
	//------

	void SampleProfile();
	virtual void SampleStat();
	virtual void SampleStat(int site);

	//------
	// log likelihoods based on sufficient statistics (substitution mappings)
	//------

	double ProfileSuffStatLogProb();
	// the site-specific suff stat log prob is yet to be implemented in subclasses
	virtual double ProfileSuffStatLogProb(int site) = 0;

	//------
	// moves
	//------

	double Move(double tuning, int n, int nrep);
	double MPIMove(double tuning, int n, int nrep);
	double NonMPIMove(double tuning, int n, int nrep);

	// the following MPI move
	// assumes that master and all slaves are in sync concerning siteprofilesuffstats
	// (which will be the case upon calling GlobalUpdateSiteProfileSuffStat)
	virtual double GlobalMoveProfile(double tuning = 1, int n = 1, int nrep = 1);
	virtual void SlaveMoveProfile();
	double MoveProfile(double tuning = 1, int n = 1, int nrep = 1);
	virtual double MoveProfile(int site, double tuning, int n, int nrep);

	//------
	// create and delete
	//------

	// called at the beginning and end of the run (see PhyloProcess)
	virtual void Create();

	virtual void Delete();

	// create, delete, update, all variables associated to each site (e.g. substitution matrix)
	virtual void CreateSite(int i) = 0;
	virtual void DeleteSite(int i) = 0;
	virtual void UpdateSite(int i) = 0;

	virtual void UpdateSites()	{
		for (int k=0; k<GetNsite(); k++)	{
			if (ActiveSite(k))	{
				UpdateSite(k);
			}
		}
	}

	void RenormalizeProfiles();

	double GetMeanStationaryEntropy() {return GetStatEnt();}
	double GetSiteStationaryEntropy(int site) {return GetStatEnt(site);}

	// summary statistic: mean entropy over all profiles
	double GetStatEnt();
	double GetStatEnt(int k);

	double** profile;
	double* allocprofile;
	double* logstatprior;
	double* profilesuffstatlogprob;

	Chrono totchrono;
	Chrono profilechrono;
};

#endif

