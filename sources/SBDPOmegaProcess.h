/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef SBDPOMEGA_H
#define SBDPOMEGA_H

#include "MixtureOmegaProcess.h"

class SBDPOmegaProcess : public virtual MixtureOmegaProcess	{

	public:
	SBDPOmegaProcess() : omegakappa(1), omegamovekappa(true), omegakappaprior(0), omegaV(0), omegamaxweighterror(0) {}
	virtual ~SBDPOmegaProcess() {}

	protected:

	double GetOmegaKapp()	{
		return omegakappa;
	}
	
	int GetNDisplayedOmegaComponent()	{
		return GetNOccupiedOmegaComponent();
	}

	int GetNOccupiedOmegaComponent()	{
		int n = 0;
		for (int k=0; k<GetNomegaMax(); k++)	{
			if (omegaoccupancy[k])	{
				n++;
			}
		}
		return n;
	}
		
	virtual int GetNomegaMax() {
		return GetNsite() > nomegamax ? nomegamax : GetNsite();
	}
	virtual void SetNomegaMax(int n) {nomegamax = 1000;}


	// overall sampling of the omega part of the model
	void SampleOmega();

	// MOVES

	double MPIMoveOmega(double tuning, int nrep);
	double NonMPIMoveOmega(double tuning, int nrep);

	virtual void SampleOmegaWeights();
	void ResampleOmegaWeights();

	// double MoveOccupiedCompAlloc(int nrep);
	// double MoveAdjacentCompAlloc(int nrep);
	
	// double MoveKappa...

	// are those ones really useful?
	/*
	double MixMoveOmega(int nmix, double tuning, int nsitenrep, int nhyperrep);
	double GlobalMixMoveOmega(int nmix, double tuning, int nsitenrep, int nhyperrep);
	void SlaveMixMoveOmega();
	*/

	protected:

	virtual void Create()	{
		if (! omegaV)	{
			MixtureOmegaProcess::Create();
			// create all data structures specific to sbdp
		}
	}

	virtual void Delete()	{
		if (omegaV)	{
			// specific deletes
			MixtureOmegaProcess::Delete();
		}
	}
	
	int nomegamax;
	double* omegaV;
	int* omegaoccupancy;
	double omegakappa;
	bool omegamovekappa;
	int omegakappaprior;

	double omegamaxweighterror;
};

#endif

