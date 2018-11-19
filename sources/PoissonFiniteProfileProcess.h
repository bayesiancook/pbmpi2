
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef POISSONFINITEPROFILE_H
#define POISSONFINITEPROFILE_H

#include "PoissonMixtureProfileProcess.h"
#include "FiniteProfileProcess.h"

// superclass for Poisson (F81) implementations
class PoissonFiniteProfileProcess: public virtual PoissonMixtureProfileProcess, public virtual FiniteProfileProcess	{

	public:

	PoissonFiniteProfileProcess() {}
	virtual ~PoissonFiniteProfileProcess() {}

    void EM_UpdateProfiles();

	virtual double Move(double tuning = 1, int n = 1, int nrep = 1);
	virtual double MPIMove(double tuning = 1, int n = 1, int nrep = 1);
	virtual double NonMPIMove(double tuning = 1, int n = 1, int nrep = 1);

	virtual void ToStream(ostream& os);
	virtual void FromStream(istream& is);

	virtual double GlobalParametersMove()	{
		return 1;
	}

	protected:

	virtual void Create()	{
		FiniteProfileProcess::Create();
		PoissonMixtureProfileProcess::Create();
	}

	virtual void Delete()	{
		PoissonMixtureProfileProcess::Delete();
		FiniteProfileProcess::Delete();
	}

	double IncrementalFiniteMove(int nrep);

	double GlobalIncrementalFiniteMove(int nrep);
	double SlaveIncrementalFiniteMove();


};

#endif

