
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef FinitePROFILE_H
#define FinitePROFILE_H

#include <cmath>
#include "DirichletMixtureProfileProcess.h"

// general superclass for all finite process mixtures on site-specific profiles
class FiniteProfileProcess: public virtual DirichletMixtureProfileProcess	{

	public:

	FiniteProfileProcess() : weight(0), fixncomp(false), empmix(0), Nfixcomp(0), statfix(0), empweight(0), statcons(0), Ncons(0), conscutoff(0)  {
		weightalpha = 1;
		mixtype = "None";
		Ncomponent = -1;
	}

	virtual ~FiniteProfileProcess(){}

	// uniform prior on component number
	double GetLogNPrior() {return 0;}

	void SetFixedNcomponent(bool in = true)	{
		fixncomp = in;
	}

	protected:

	virtual double GlobalSMCAddSites();
	virtual double SMCAddSites();

	void ReadStatFix(string name);
	void SetStatFix();
	void BroadcastStatFix();
	void SlaveGetStatFix();

	virtual double Move(double tuning = 1, int n = 1, int nrep = 1);
	virtual double MPIMove(double tuning = 1, int n = 1, int nrep = 1);
	virtual double NonMPIMove(double tuning = 1, int n = 1, int nrep = 1);

	virtual double MoveHyper(double tuning, int nrep); // added virtual
	virtual double MoveWeightAlpha(double tuning, int nrep); // added virtual
	virtual double MoveStatFixAlpha(double tuning, int nrep); // added virtual
	double MoveNcomponent(int nrep);
	virtual void ResampleWeights(); // added virtual

	// the following MPI move
	// assumes that master and all slaves are in sync concerning siteprofilesuffstats
	// (which will be the case upon calling GlobalUpdateSiteProfileSuffStat)
	double GlobalIncrementalFiniteMove(int nrep);
	double SlaveIncrementalFiniteMove();
	double IncrementalFiniteMove(int nrep);

	// called at the beginning and end of the run (see PhyloProcess)

	void SetNcomponent(int inncomp)	{
		Ncomponent = inncomp;
	}

	virtual void Create();
	virtual void Delete();

	// multinomial 
	void SampleAlloc();
	void SampleStat();
	virtual void SampleStat(int cat);
	virtual void SampleHyper(); // added virtual
	virtual void SampleWeights(); // added virtual

	virtual void PriorSampleProfile();
	virtual void PriorSampleHyper();
	virtual void PriorSampleWeights();

	virtual double LogHyperPrior();

	virtual double LogStatPrior(int cat);
	virtual double LogStatPriorConstrained(int cat);
	// void ReadConstraints(string filename);
	void SetStatCons();

	double LogWeightPrior();
	double LogStatAlphaPrior();

	double GetWeightedStationaryEntropy()	{
		double total = 0;
		for (int k=0; k<Ncomponent; k++)	{
			total += weight[k] * GetStatEnt(k);
		}
		return total;
	}

	double* weight;
	double weightalpha;

	double statfixalpha;

	bool fixncomp;
	int empmix;
	
	int Nfixcomp;
	double** statfix;
	double* empweight;

	int** statcons;
	int Ncons;
	double* conscutoff;

	string mixtype;
};

#endif

