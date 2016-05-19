
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
#include "Random.h"

class OmegaProcess : public virtual MPIModule	{

	public:

	OmegaProcess() : fixomega(0), omegaprior(0) {}
	virtual ~OmegaProcess() {}

	protected:

	// omega

	virtual void UpdateOmegaSuffStat()	{
		cerr << "in OmegaProcess::UpdateOmegaSuffStat\n";
		exit(1);
	}

	virtual void GlobalUpdateOmegaSuffStat()	{
		cerr << "in OmegaProcess::GlobalUpdateOmegaSuffStat\n";
		exit(1);
	}

	virtual void SlaveUpdateOmegaSuffStat()	{
		cerr << "in OmegaProcess::SlaveUpdateOmegaSuffStat\n";
		exit(1);
	}	

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

	virtual double GetSiteOmega(int site)	{
		cerr << "in OmegaProcess::GetSiteOmega\n";
		exit(1);
	}

	virtual void UpdateSiteOmegaSuffStat()	{
		cerr << "in OmegaProcess::UpdateSiteOmegaSuffStat\n";
	}

	int fixomega;
	int omegaprior;
	double* siteomegasuffstatbeta;
	int* siteomegasuffstatcount;
};


class SingleOmegaProcess : public virtual OmegaProcess	{

	public:

	SingleOmegaProcess() : omega(0)	{}
	virtual ~SingleOmegaProcess() {}

	double GetOmega()	{
		return *omega;
	}

	double GetSiteOmega(int site)	{
		GetOmega();
	}
	
	virtual double OmegaSuffStatLogProb()	{
		//return ProfileSuffStatLogProb();
		return omegasuffstatcount * log(*omega) - omegasuffstatbeta * *omega;
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
	int omegasuffstatcount;
	double omegasuffstatbeta;
	
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

	virtual void SampleOmega();
	virtual void SampleOmegaWeights();
	virtual double LogOmegaPrior();
	double MoveOmega(double tuning); 
	double MoveOmegaValues(double tuning); 
	void ResampleOmegaWeights();
	double GlobalOmegaIncrementalFiniteMove(int nrep);
	double SlaveOmegaIncrementalFiniteMove();
	void UpdateOmegaSuffStat();
	void GlobalUpdateOmegaSuffStat();
	void SlaveUpdateOmegaSuffStat();
	virtual double OmegaSuffStatLogProb()	{

		double total = 0;
		for (int l=0; l<Nomega; l++)	{
			total += OmegaSuffStatLogProb(l);
		}
		return total;
	}

	virtual double OmegaSuffStatLogProb(int l)	{
		return omegasuffstatcount[l] * log(omega[l]) - omegasuffstatbeta[l] * *omega;
	}

	double OmegaLogStatProb(int site, int omegaalloc)	{
		return siteomegasuffstatcount[site] * log(omega[omegaalloc]) - siteomegasuffstatbeta[site] * omega[omegaalloc];
	}


	protected:

	virtual void Create()	{
		if (! omega)	{
			omega = new double[Nomega];
			omegaalloc = new int[GetNsite()];
			omegasuffstatbeta = new double[Nomega];
			omegasuffstatcount = new int[Nomega];
			siteomegasuffstatbeta = new double[GetNsite()];
			siteomegasuffstatcount = new int[GetNsite()];
			//siteomegasuffstatbeta = new double*[GetNsite()];
			//siteomegasuffstatcount = new int*[GetNsite()];
			omegaweight = new double[Nomega];
			//for (int i=0; i<GetNsite(); i++)	{
			//	siteomegasuffstatbeta[i] = new double[Nomega];
			//	siteomegasuffstatcount[i] = new int[Nomega];
			//}
			SampleOmega();
		}
	}

	virtual void Delete()	{
		if (omega)	{
			delete[] omega;
			delete[] omegaalloc;
			delete[] omegasuffstatbeta;
			delete[] omegasuffstatcount;
			delete[] siteomegasuffstatbeta;
			delete[] siteomegasuffstatcount;
			delete[] omegaweight;
			omega = 0;
		}
	}

	double* omega;
	int* omegaalloc;
	int Nomega;
	double* omegasuffstatbeta;
	int* omegasuffstatcount;
	double* omegaweight;
	double omegaweightalpha;
};

#endif



