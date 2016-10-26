
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef EXPCONPARTGTRMIXTUREPROFILE_H
#define EXPCONPARTGTRMIXTUREPROFILE_H

#include "PartitionGTRMixtureProfileProcess.h"
#include "ExpoConjugatePartitionGTRProfileProcess.h"
#include "DirichletProfileProcess.h"

// Exponential conjugate GTR models
class ExpoConjugatePartitionGTRMixtureProfileProcess : public virtual PartitionGTRMixtureProfileProcess, public virtual ExpoConjugatePartitionGTRProfileProcess, public virtual DirichletProfileProcess {

	public:

	ExpoConjugatePartitionGTRMixtureProfileProcess() : profilesuffstatcount(0), profilesuffstatbeta(0) {}
	virtual ~ExpoConjugatePartitionGTRMixtureProfileProcess() {}

	protected:

	// matrices are not used during Sufficient-statistics based updates
	// all matrices are deleted upon 'Collapsing' (pruning->suffstat)
	// and re-created upon 'Unfolding' (suffstat->pruning) state-vectors
	// using CreateMatrices() and DeleteMatrices()
	virtual void CreateComponent(int k) {
		occupancy[k] = 0;
		if (activesuffstat)	{
			for (int l=0; l<GetDim(); l++)	{
				profilesuffstatcount[k][l] = 0;
				profilesuffstatbeta[k][l] = 0;
			}
		}
		else	{
		}
		SampleStat(k);
	}
	virtual void DeleteComponent(int k) {
	}
	virtual void UpdateComponent(int k) {
		// update suff stats ?
	};

	double ProfileSuffStatLogProb(int cat);
	void SwapComponents(int cat1, int cat2);
	virtual double LogStatProb(int site, int cat);

	virtual double LogProxy(int site, int cat)	{
		cerr << "in ExpoConjugateGTRMixtureProfileProcess::LogProxy\n";
		exit(1);
	}

	virtual void Create();
	virtual void Delete();

	// collects site-specific suffstats and pools them componentwise
	virtual void UpdateModeProfileSuffStat();
	virtual void GlobalUpdateModeProfileSuffStat();
	virtual void SlaveUpdateModeProfileSuffStat();

	// component-specific sufficient statistics
	int** profilesuffstatcount;
	double** profilesuffstatbeta;
	int* allocprofilesuffstatcount;
	double* allocprofilesuffstatbeta;
	int** tmpprofilesuffstatcount;
	double** tmpprofilesuffstatbeta;
	int* alloctmpprofilesuffstatcount;
	double* alloctmpprofilesuffstatbeta;

	double PoissonDiffLogSampling(int cat, int site);
};

#endif

