
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef PROFILE_H
#define PROFILE_H

#include "Chrono.h"
#include "StateSpace.h"
#include "SequenceAlignment.h"

#include <iostream>
#include <map>
#include "MPIModule.h"

using namespace std;

// this class takes care of all aspects related to substitution processes (not including branch lengths and site-specific relative rates)
// this includes global parameters (such as global relative exchangeabilities for CATGTR, or mutational parameters for Codon Mutation Selection models)
// and hyperparameters of the mixtures of the mixture models
//
// this code allows for variations across sites only of the profiles
// all other parameters should be homogeneous along the sequence

const double stateps = 1e-100;

class ProfileProcess : public virtual MPIModule {

	public:

	ProfileProcess() : dim(0), activesuffstat(false), statinfcount(0), totstatcount(0), profilepriortype(0), proposemode(0), allocmode(0) {}
	virtual ~ProfileProcess() {}

	//------
	// basic accessors
	//------

	// dimension of the profile(s)
	int GetDim() {return dim;}

	virtual int GetNstate() = 0;

	virtual StateSpace* GetStateSpace() = 0;

	// profile associated with given site
	virtual double* GetProfile(int site) = 0;


	//------
	// priors
	//------

	// total prior associated to all aspects of the substitution processes across sites
	virtual double LogProfilePrior() = 0;

	// prior on hyperparameters
	virtual double LogHyperPrior() = 0;

	// prior on profiles, given hyperparameters
	virtual double LogStatPrior() = 0;

	// basic tool for getting the log prob associated to a profile
	// can be a frequency profile (Dirichlet prior) or a log-frequency profile (normal prior)
	virtual double LogFrequencyStatPrior(double* prof) = 0;

	virtual double GetMinStat(double* profile, int site) = 0;


	//------
	// sampling from prior
	//------

	// sample all parameters (including hyperparameters)
	// strictly from prior ...
	virtual void PriorSampleProfile() {
		// .. although .. 
		SampleProfile();
	}

	// with some quirks to make draw more reasonable for MCMC
	virtual void SampleProfile() = 0;

	// sampling hyperparameters
	virtual void SampleHyper() = 0;

	// sampling profiles, given hyperparameters
	virtual void SampleStat() = 0;

	// basic tool for sampling a profile, depending on the prior
	// can be a frequency profile (Dirichlet prior) or a log-frequency profile (normal prior)
	virtual void SampleFrequencyStat(double* prof) = 0;

	//------
	// log likelihoods based on sufficient statistics (substitution mappings)
	//------

	virtual double ProfileSuffStatLogProb() = 0;

	//------
	// moves
	//------

	// implemented in specialized phyloprocess classes
	// will collect the sufficient statistics necessary for updating substitution parameters
	// separately for each site
	// these statistics will be accessed through specific functions declared and implemented in subclasses
	virtual void UpdateSiteProfileSuffStat() = 0;
	virtual void GlobalUpdateSiteProfileSuffStat() = 0;
	virtual void SlaveUpdateSiteProfileSuffStat() = 0;

	virtual void GlobalUpdateModeProfileSuffStat() {}
	virtual void SlaveUpdateModeProfileSuffStat() {}
	virtual void UpdateModeProfileSuffStat() {}

	virtual void GlobalUpdateParameters() = 0;
	virtual void SlaveUpdateParameters() = 0;

	// MCMC MOVES

	// generic Move function
	virtual double Move(double tuning = 1, int n = 1, int nrep = 1) = 0;

	// move on hyperparameters
	virtual double MoveHyper(double tuning, int nrep) = 0;

	// auxiliary: for Metropolis Moves on profiles
	virtual double ProfileProposeMove(double* profile, double tuning, int n, int K, int cat, double statmin);

	// move all parameters (except component-specific or site-specific frequency profiles) 
	virtual double GlobalParametersMove() = 0;

	// STREAMS

	virtual void ToStream(ostream& os) = 0;
	virtual void FromStream(istream& is) = 0;


	// more specialized methods

	// summary statistics
	virtual double GetMeanStationaryEntropy() = 0;
	virtual double GetSiteStationaryEntropy(int site) = 0;

	double GetStatInfCount() {
		double tmp = 0;
		if (totstatcount)	{
			tmp = ((double) statinfcount) / totstatcount;
		}
		statinfcount = 0;
		totstatcount = 0;
		return tmp;
	}

	protected:

	virtual double GetNormalizationFactor()	{cerr << "should not be here\n"; exit(1); return 1;}

	// called at the beginning of the run only
	// but potentially called several times due to multiple inheritance
	// thus, should internally control that no object is created twice
	// (by checking that pointers are null before creating)

	void SetDim(int indim)	{
		dim = indim;
	}

	virtual SequenceAlignment* GetData() = 0;

	virtual void Create() {}
	virtual void Delete() {}

	virtual double* GetEmpiricalFreq() = 0;

	virtual double GlobalSMCAddSites();
	virtual double SMCAddSites() {}

	virtual void SampleSiteMapping(int site) = 0;
	virtual double SiteLogLikelihood(int site) = 0;

	int dim;
	bool activesuffstat;

	Chrono chronodp;
	Chrono chronostat;
	Chrono chronorr;

	int statinfcount;
	int totstatcount;

	int profilepriortype;

	int proposemode;
	double profacc;
	double proftry;
	int allocmode;
	double rracc;
	double rrtry;
};


#endif
