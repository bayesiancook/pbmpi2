
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef SUBS_H
#define SUBS_H

#include "DGamRateProcess.h"
#include "ProfileProcess.h"
#include "BranchSitePath.h"
#include "Chrono.h"
#include <algorithm>

// ----
// Substitution Process is the class gathering nearly all CPU-intensive methods of the program
//
// A substitution process object is aware of the "across site" dimension, but not of the "across branch" dimension (which is taken care of by PhyloProcess)
//
// most important accessors:
//
// int GetNstate(int i) : number of states for the substitution process at site i
//
// int GetNrate(int i) : number of rate categories at site i
// double GetRate(int i) : gives the rate at site i (when sites are explicitly allocated to one of the rate categories)
// 
// const double* GetProfile(int i) : equilibrium frequency profile at site i

// how these rates and profiles are distributed across sites is not managed by SubstitutionProcess itself,
// but by its two super classes RateProcess and ProfileProcess

// in addition, SubstitutionProcess implements a whole series of CPU intenstive methods for 
// 	- propagating conditional likelihood vectors
//	- sampling rates at each site proportionally to the conditional likelihoods
//	- sampling states (amino-acids) at internal nodes of the tree
// 	- sampling substitution mappings along branches

// since SubstitutionProcess is not aware of the overarching tree structure
// it does not deal with the overall structure of those computations 
// which are in fact managed by the PhyloProcess subclass

// In other words. PhyloProcess is responsible for the higher-level algorithmic part
// while SubstitutionProcess makes the "across-site" micro-managment part of the computation

// Thus, for all these methods, PhyloProcess organizes the computations
// and repeteadly asks SubstitutionProcess to perform one pass over all sites for one of those computational items
// each time, providing SubstitutionProcess with the necessary information
// (e.g. length of the branch, state of the process at both ends for each site, etc.)

// since phylogenomics is about using many aligned positions (large value of GetNsite())
// about everything that needs to be parallelized is here
// the only exception is "

class SubstitutionProcess : public virtual RateProcess, public virtual ProfileProcess {

	public:

	SubstitutionProcess() : condsitelogL(0), sitelogL(0), infprobcount(0), suboverflowcount(0) {}
	virtual ~SubstitutionProcess() {}

	// those two accessors are abstract: what they return depends on whether we use recoded Poisson or GTR substitution processes
	virtual int GetNstate(int site) {return GetNstate();}
	virtual int GetNstate() {return GetDim();}
	virtual const double* GetStationary(int site)	{
		cerr << "in SubstitutionProcess::GetStationary\n";
		exit(1);
		return 0;
	}

	int GetInfProbCount() {return infprobcount;}

	virtual void ActivateZip() {}
	virtual void InactivateZip() {}

	protected:

	void Create();
	void Delete();

	// basic modules for creating deleting arrays of conditional likelihoods
	// used by PhyloProcess
	double*** CreateConditionalLikelihoodVector();
	void DeleteConditionalLikelihoodVector(double*** condl);

	double* CreateProbVector()	{
		return new double[GetSiteMax() - GetSiteMin()];
	}

	void CreateCondSiteLogL();
	void DeleteCondSiteLogL();

	// ------------------
	// various computatonial accessory methods
	// used by PhyloProcess

	// CPU : level 1
	// if aux==0, assumes likelihoods have been computed
	void DrawAllocations(double*** aux);

	// in the following
	// bool condalloc = true means that we want to make the computation, for each site,
	// only for the category specified for that site by double* ratealloc

	// CPU : level 1
	void Reset(double*** condl, bool condalloc, int phase = 1);
	void Multiply(double*** from, double*** to, bool condalloc);
	void Add(double*** from, double*** to, bool condalloc);
	void MultiplyByStationaries(double*** from, bool condalloc);
	void Offset(double*** condl, bool condalloc);
	virtual void Initialize(double*** condl, const int* leafstates, bool condalloc);

    int GetParsimonyBranchScore(int nchildren, double*** t, int* sitescores);
    void ChooseParsimonyStates(const int* upstates, int* downstates, double*** t);
    virtual void ParsimonyProfiles(const int* upstates, const int* downstates, double** profile);

	// CPU : level 2
	double ComputeLikelihood(double*** aux, bool condalloc);

	// CPU : level 3
	// implemented in GTR or POisson Substitution process
	virtual void Propagate(double*** from, double*** to, double time, bool condalloc)	{
		cerr << "in SubstitutionProcess::Propagate\n";
		exit(1);
	}

	// CPU : level 1
	// implemented in GTR or POisson Substitution process
	// here, assumes that each site is under the rate category defined by double* ratealloc
	virtual void ChooseStates(double*** aux, int* states);
	

	// site-wise versions of the same functions

	// CPU : level 1
	// if aux==0, assumes likelihoods have been computed
	void SiteDrawAllocations(int site, double** aux);

	// CPU : level 1
	void SiteReset(int site, double** condl, bool condalloc = false);
	void SiteMultiply(int site, double** from, double** to, bool condalloc = false);
	void SiteMultiplyByStationaries(int site, double** from, bool condalloc = false);
	void SiteOffset(int site, double** condl, bool condalloc = false);
	virtual void SiteInitialize(int site, double** condl, const int leafstate, bool condalloc = false);

	// CPU : level 2
	double SiteComputeLikelihood(int site, double** aux, bool condalloc = false);

	// CPU : level 3
	// implemented in GTR or POisson Substitution process
	virtual void SitePropagate(int site, double** from, double** to, double time, bool condalloc = false)	{
		cerr << "in SubstitutionProcess::SitePropagate\n";
		exit(1);
	}

	// CPU : level 1
	// implemented in GTR or POisson Substitution process
	// here, assumes that each site is under the rate category defined by double* ratealloc
	virtual int SiteChooseState(int site, double** aux);
	
	// CPU : level 3
	// implemented in GTR or POisson Substitution process
	void SamplePaths(BranchSitePath** patharray, int* stateup, int* statedown, double time);
	void SampleRootPaths(BranchSitePath** patharray, int* rootstate);

	virtual BranchSitePath* SampleSitePath(int site, int stateup, int statedown, double time)	{
		cerr << "in SubstitutionProcess::SampleSitePath\n";
		exit(1);
		return 0;
	}

	virtual BranchSitePath* SampleRootSitePath(int site, int rootstate)	{
		cerr << "in SubstitutionProcess::SampleRootSitePath\n";
		exit(1);
		return 0;
	}

	// -------------------------

	double** condsitelogL;
	double* sitelogL;
	double logL;

	int infprobcount;
	int suboverflowcount;
};

#endif
