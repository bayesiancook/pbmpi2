
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef DIRPROFILE_H
#define DIRPROFILE_H

#include "ProfileProcess.h"

class DirichletProfileProcess : public virtual ProfileProcess {

	public:

	DirichletProfileProcess() : dirweight(0), statweight(0) {}
	virtual ~DirichletProfileProcess() {}

	double GetMeanDirWeight();
	double GetCenterStatEnt();

	virtual void SampleFrequencyStat(double* prof);
	virtual double LogFrequencyStatPrior(double* prof);

	// prior of hyperparameters
	virtual double LogHyperPrior();
	virtual void SampleHyper();
	virtual void PriorSampleHyper();

	void SetEmpiricalPriorStatMix();

	// move on hyperparameters
	virtual double MoveHyper(double tuning, int nrep);

	double GetPostStatProb(double* prof, double* post);

	virtual void GetPostCount(int* count)	{
		cerr << "in DirichletProfileProcess::GetPostCount, generic version\n";
		exit(1);
	}

	double MoveStatWeight();
	double MoveDirWeights(int cat, double tuning, int nrep);

	protected:

	virtual void Create();
	virtual void Delete();

    // int Nstatcomp;
	double* statweight;
	double** dirweight;
};


#endif
