
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef POISSONMIXTUREPROFILE_H
#define POISSONMIXTUREPROFILE_H

#include "PoissonProfileProcess.h"
#include "DirichletMixtureProfileProcess.h"

// superclass for Poisson (F81) implementations
class PoissonMixtureProfileProcess: public virtual PoissonProfileProcess, public virtual DirichletMixtureProfileProcess	{

	public:

	PoissonMixtureProfileProcess() : profilesuffstatcount(0) {}
	virtual ~PoissonMixtureProfileProcess() {}

	protected:

	virtual void Create();
	virtual void Delete();

	virtual void CreateComponent(int k) {
		occupancy[k] = 0;
		SampleStat(k);
	}

	virtual void UpdateZip(int site) = 0;

	virtual void DeleteComponent(int k) {
	}
	virtual void UpdateComponent(int k) {
		for (int i=0; i<GetNsite(); i++)	{
			if (ActiveSite(i))	{
				if (alloc[i] == k)	{
					UpdateZip(i);
				}
			}
		}
	}

	// posterior
	// collects sufficient statistics across sites, pools them componentwise
	void UpdateModeProfileSuffStat();
	void GlobalUpdateModeProfileSuffStat();
	void SlaveUpdateModeProfileSuffStat();

	// virtual void CreateComponent(int k)	{}

	// suffstat lnL of all sites allocated to component cat
	double ProfileSuffStatLogProb(int cat);

	// difference between:
	// suffstat lnL of all sites allocated to component cat when site <site> is among them, and
	// suffstat lnL of all sites allocated to component cat when site <site> is not among them
	double DiffLogSampling(int cat, int site);
	virtual double LogStatProb(int site, int cat);
	double LogStatIntPrior(int cat);
	double LogStatIntPrior();
	double MoveDirWeights(double tuning, int nrep);

	double MoveProfile();
	double MoveProfile(int cat);

	void SwapComponents(int cat1, int cat2);
	void AddSite(int site, int cat);
	void RemoveSite(int site, int cat);

	double GetNormRate(int k)	{

		double tot = 0;
		for (int i=0; i<GetDim(); i++)	{
			for (int j=i+1; j<GetDim(); j++)	{
				tot += profile[k][i] * profile[k][j];
			}
		}
		return 2*tot;
	}

	virtual double GetNormalizationFactor()	{
		UpdateOccupancyNumbers();
		double norm = 0;
		int tot = 0;
		for (int k=0; k<GetNcomponent(); k++)	{
			if (occupancy[k])	{
				double tmp = GetNormRate(k);
				norm += (occupancy[k] + 1) * tmp;
				tot += occupancy[k] + 1;
			}
		}
		if (! tot)	{
			return 1;
		}
		norm /= tot;
		return norm;
	}

	// private:
	double** profilesuffstatcount;
	double* allocprofilesuffstatcount;
	double** tmpprofilesuffstatcount;
	double* alloctmpprofilesuffstatcount;
};

#endif

