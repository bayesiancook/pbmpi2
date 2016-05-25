
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef SBDPPROFILE_H
#define SBDPPROFILE_H

#include <cmath>
#include "DPProfileProcess.h"

// general superclass for all finite process mixtures on site-specific profiles
class SBDPProfileProcess: public virtual DPProfileProcess	{

	using MixtureProfileProcess::LogStatPrior;

	public:

	SBDPProfileProcess() : DPProfileProcess(), V(0), maxweighterror(0), InitIncremental(0) {}
	virtual ~SBDPProfileProcess(){}

	virtual double Move(double tuning = 1, int nmix = 5, int nrep = 1, int nalloc = 1);

	int GetNDisplayedComponent()	{
		return GetNOccupiedComponent();
	}

	protected:

	double GetMaxWeightError() {return maxweighterror;}
	void ResetMaxWeightError() {maxweighterror = 0;}

	virtual void Create();
	virtual void Delete();

	virtual double MPIMove(double tuning = 1, int nmix = 1, int nrep = 1, int nalloc = 1);
	virtual double NonMPIMove(double tuning = 1, int nmix = 1, int nrep = 1, int nalloc = 1);

	virtual double GlobalSMCAddSites();
	virtual double SMCAddSites();

	double GetWeightEnt()	{
		double tot = 0;
		for (int k=0; k<GetNcomponent(); k++)	{
			if (weight[k] > 1e-8)	{
				tot -= weight[k] * log(weight[k]);
			}
		}
		return tot;
	}

	int GetLastOccupiedComponent()	{
		int kmax = 0;
		for (int i=0; i<GetNsite(); i++)	{
			if (ActiveSite(i))	{
				if (kmax < alloc[i])	{
					kmax = alloc[i];
				}
			}
		}
		return kmax;
	}

	int GetNCutoff(double cutoff)	{
		int n = (int) (GetNOccupiedComponent() * (1 - cutoff));
		int k = GetLastOccupiedComponent();
		int tot = occupancy[k];
		while (k && (tot < n))	{
			k--;
			if (k)	{
				tot += occupancy[k];
			}
		}
		return k;
	}
		
	virtual void SwapComponents(int cat1, int cat2);

	virtual double GetWeight(int cat)	{
		return weight[cat];
	}

	// double LogAllocPrior();
	double LogIntegratedAllocProb();
	double MoveKappa(double tuning, int nrep);

	void SampleAlloc();
	void IncrementalSampleAlloc();

	void SampleWeights();
	void ResampleWeights();

	/*
	virtual void SampleHyper()	{
		DPProfileProcess::SampleHyper();
		SampleWeights();
	}

	virtual void PriorSampleHyper()	{
		DPProfileProcess::PriorSampleHyper();
		SampleWeights();
	}
	*/

	// void ResampleLastWeight();
	double MoveOccupiedCompAlloc(int nrep = 1);
	double MoveAdjacentCompAlloc(int nrep = 1);

	double LogStatPrior();

	virtual double MixMove(int nrep, int nallocrep, double epsilon, int nprofilerep);
	virtual double GlobalMixMove(int nrep, int nallocrep, double epsilon, int nprofilerep);
	virtual void SlaveMixMove();

	double IncrementalDPMove(int nrep, double epsilon);

	double* V;
	double* weight;

	double maxweighterror;

	int InitIncremental;
};

#endif

