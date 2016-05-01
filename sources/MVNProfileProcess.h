
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef MVNPROFILE_H
#define MVNPROFILE_H

#include "ProfileProcess.h"
#include "CovMatrix.h"

class MVNProfileProcess : public virtual ProfileProcess {

	public:

	MVNProfileProcess() : kappa(0), df(0), covmatrix(0), mean(0) {}
	virtual ~MVNProfileProcess() {}

	virtual void SampleFrequencyStat(double* logprof);
	virtual double LogFrequencyStatPrior(double* logprof);

	virtual void UpdateFrequencyStat(double* logprof, double* prof);

	// prior of hyperparameters
	virtual double LogHyperPrior();
	virtual void SampleHyper();

	// move on hyperparameters
	virtual double MoveHyper(double tuning, int nrep)	{
		ResampleCovMatrix();
		MoveKappa(tuning,nrep);
	}

	double MoveKappa(double tuning, int nrep);

	// auxiliary: for Metropolis Moves on profiles
	double LogProfileProposeMove(double* logprofile, double tuning, int n);

	// conjugate inverse-wishart / normal Gibbs sampling of the covariance matrix
	// will be implemented in site-specific or component-specific sub-classes
	virtual void ResampleCovMatrix() = 0;

	protected:

	virtual void Create();
	virtual void Delete();

	double* kappa;
	int df;
	InverseWishartMatrix* covmatrix;
	double* mean;
};


#endif
