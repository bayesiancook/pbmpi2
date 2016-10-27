
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef MIXOMEGA_H
#define MIXOMEGA_H

#include "OmegaProcess.h"

class MixtureOmegaProcess : public virtual OmegaProcess	{

	public: 
	MixtureOmegaProcess() : omega(0), omegaprior(1) {}

	virtual ~MixtureOmegaProcess() {}

	// implementing the interface declared by OmegaProcess
	double GetSiteOmega(int site)	{
		//return omega[omegaalloc[site]];
		return GetOmega(omegaalloc[site]);
	}

	double* GetSiteOmegaPtr(int site)	{
		return &omega[omegaalloc[site]];
	}

	// more specific methods defined by mixtures of omega
	int GetNomega() {return Nomega;}

	double GetOmega(int k)	{
		return omega[k];
	}

	double* GetOmegaPointer(int k)	{
		return &omega[k];
	}

	int GetOmegaSiteAlloc(int site)	{
		return omegaalloc[site];
	}

	// generic: valid for all mixtures, finite and infinite

	// choice between alternative priors
	// definition of prior hyperparameters
	// what is shared by all mixtures, what is specific to sbdp or finite?

	virtual double LogOmegaPrior();

	virtual double LogOmegaPrior(int component)	{
		cerr << "implement MixtureProfile::LogOmegaPrior\n";
		cerr << "for a single component\n";
		exit(1);
	}

	virtual void SampleOmegaAlloc();

	// should depend on omegaprior choice
	double LogOmegaHyperPrior()	{
		return -omegaalpha - omegabeta;
	}

	virtual void SampleOmegaHyper()	{
		omegaalpha = 1.0;
		omegabeta = 1.0;
	}

	// overall sampling of the omega part of the model
	virtual void SampleOmega() = 0;

	virtual void SampleOmegas()	{
		if (omegaprior == 0)	{ //ratio of exponential random variables
			double r1, r2;
			for (int i=0; i<GetNomega(); i++)	{
				r1 = rnd::GetRandom().Uniform();
				r2 = rnd::GetRandom().Uniform();
				omega[i] = r1/r2;
			}
			
		}
		else if (omegaprior == 1)	{ // gamma
			omegaalpha = omegabeta = 1.0;
			for (int i=0; i<GetNomega(); i++)	{
				omega[i] = rnd::GetRandom().Gamma(omegaalpha,omegabeta);
			}
		}
		else	{
			cerr << "omega prior not correctly set\n";
			exit(1);
		}
	}

	// -----------------
	// -----------------
	// MOVES

	// this one is implemented in MixtureOmegaProcess.cpp
	// just calls MPI or NonMPIMoveOmega
	virtual double MoveOmega(double tuning);

	// those ones are implemented in Finite (Multiple) and SBDP
	virtual double MPIMoveOmega(double tuning, int nrep) = 0;
	virtual double NonMPIMoveOmega(double tuning, int nrep) = 0;

	virtual void SampleOmegaWeights() = 0;
	virtual void ResampleOmegaWeights() = 0;

	// reallocation move: same for all types of mixtures
	double GlobalOmegaIncrementalFiniteMove(int nrep);
	double SlaveOmegaIncrementalFiniteMove();

	// moving omega values
	double MoveOmegas(double tuning, int nrep);
	double MoveCompOmega(double tuning, int component);

	double MoveOmegaHyper(double tuning, int nrep);
	double MoveOmegaAlpha(double tuning);
	double MoveOmegaBeta(double tuning);

	double MoveOmegaValuesAndHyperParameters(double tuning, int nrep)	{

		// should decide which procedure depending on prior on omegas
		// conjugate gibbs: only under gamma prior
		// (1) pure MH
		// resample omegas
		MoveOmegas(tuning, nrep);
		// resample hyperparameters
		MoveOmegaHyper(tuning,nrep);
		MoveOmegaHyper(0.1*tuning,nrep);

		// (2) Gibbs resampling the omegas + MH on hyperparams
		// ResampleOmegas();
		// move hyperparam
		// MoveOmegaHyper(tuning,nhyperrep);
		// MoveOmegaHyper(0.1*tuning,nhyperrep);

		// (3) only if correlations between alpha,beta and the omega_k's are strong
		// MoveOmegaHyperIntegrated(tuning,nhyperrep);
		// MoveOmegaHyperIntegrated(0.1*tuning,nhyperrep);
		// ResampleOmegas();
	}

	// -----------------
	// -----------------
	// -----------------

	void UpdateOmegaSuffStat();
	void GlobalUpdateOmegaSuffStat();
	void SlaveUpdateOmegaSuffStat();

	// total log prob summed over components
	virtual double OmegaSuffStatLogProb()	{
		double total = 0;
		for (int l=0; l<Nomega; l++)	{
			total += OmegaSuffStatLogProb(l);
		}
		return total;
	}

	// total log prob for a given component
	virtual double OmegaSuffStatLogProb(int l)	{
		return compomegasuffstatcount[l] * log(omega[l]) - compomegasuffstatbeta[l] * omega[l];
	}

	// used in realloc move
	double OmegaLogStatProb(int site, int omegaalloc)	{
		return siteomegasuffstatcount[site] * log(omega[omegaalloc]) - siteomegasuffstatbeta[site] * omega[omegaalloc];
	}

	protected:

	virtual void Create()	{
		if (! omega)	{
			OmegaProcess::Create();
			omega = new double[Nomega];
			omegaalloc = new int[GetNsite()];
			compomegasuffstatbeta = new double[Nomega];
			compomegasuffstatcount = new int[Nomega];
			omegaweight = new double[Nomega];
			SampleOmega();
		}
	}

	virtual void Delete()	{
		if (omega)	{
			delete[] omega;
			delete[] omegaalloc;
			delete[] compomegasuffstatbeta;
			delete[] compomegasuffstatcount;
			delete[] omegaweight;
			omega = 0;
			OmegaProcess::Delete();
		}
	}

	double* omega;
	int* omegaalloc;
	int Nomega;
	double* compomegasuffstatbeta;
	int* compomegasuffstatcount;
	double* omegaweight;

	int omegaprior;
	double omegaalpha;
	double omegabeta;
};

#endif
