
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

	OmegaProcess() : fixomega(0), omegaprior(0) {}
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
	int omegaprior;
};


class SingleOmegaProcess : public virtual OmegaProcess	{

	public:

	SingleOmegaProcess() : omega(0)	{}
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
};

class MultipleOmegaProcess : public virtual OmegaProcess	{

	public: 
	MultipleOmegaProcess() : omega(0) {}

	virtual ~MultipleOmegaProcess() {}

	int GetNomega() {return Nomega;}

	double GetOmega(int k)	{
		return omega[k];
	}

	double GetMeanOmega()	{
		double tot = 0;
		for (int i=0; i<GetNsite(); i++)	{
			tot += GetOmega(i);
		}
		return tot / GetNsite();
	}

	double* GetOmegaPointer(int k)	{
		return &omega[k];
	}

	double GetSiteOmega(int site)	{
		return omega[omegaalloc[site]];
	}

	int GetOmegaSiteAlloc(int site)	{
		return omegaalloc[site];
	}

	// still to be implemented
	virtual void SampleOmega()	{
	}

	// still to be implemented
	virtual double LogOmegaPrior()	{
		return 0;
	}

	virtual double OmegaSuffStatLogProb()	{

		double total = 0;
		for (int l=0; l<Nomega; l++)	{
			total += OmegaSuffStatLogProb(l);
		}
		return total;
	}

	// still to be implemented
	virtual double OmegaSuffStatLogProb(int l)	{
		return 0;
	}

	protected:

	virtual void Create()	{
		if (! omega)	{
			omega = new double[Nomega];
			omegaalloc = new int[GetNsite()];
			SampleOmega();
		}
	}

	virtual void Delete()	{
		if (omega)	{
			delete[] omega;
			delete[] omegaalloc;
			omega = 0;
		}
	}

	double* omega;
	int* omegaalloc;
	int Nomega;
};

#endif



