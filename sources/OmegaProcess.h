
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef OMEGA_H
#define OMEGA_H

#include <cstdlib>
#include <iostream>
#include <cmath>
using namespace std;

#include "MPIModule.h"

class OmegaProcess : public virtual MPIModule	{

	public:

	OmegaProcess() : fixomega(0) {}
	virtual ~OmegaProcess() {}

	protected:

	// omega
	virtual double LogOmegaPrior()	{
		cerr << "in OmegaProcess::LogOmegaPrior\n";
		exit(1);
	}

	virtual double OmegaSuffStatLogProb()	{
		cerr << "in OmegaProcess::OmegaSuffStatLogProb\n";
		exit(1);
	}

	virtual void SampleOmega()	{
		cerr << "in OmegaProcess::SampleOmega\n";
		exit(1);
	}

	virtual double MoveOmega(double tuning)	{
		cerr << "in OmegaProcess::MoveOmega\n";
		exit(1);
	}

	virtual void UpdateOmega()	{
		cerr << "in OmegaProcess::UpdateOmega\n";
		exit(1);
	}

	int fixomega;
};


class SingleOmegaProcess : public virtual OmegaProcess	{

	public:

	SingleOmegaProcess() : omega(0), omegaprior(0)	{}
	virtual ~SingleOmegaProcess() {}

	double GetOmega()	{
		return *omega;
	}

	// omega
	virtual double LogOmegaPrior();
	virtual void SampleOmega();
	double MoveOmega(double tuning); 

	protected:

	virtual void Create()	{
		if (! omega)	{
			omega = new double;
		}
	}

	virtual void Delete()	{
		if (omega)	{
			delete omega;
			omega = 0;
		}
	}

	double* omega;
	int omegaprior;
};

#endif



