
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef DPPROFILE_H
#define DPPROFILE_H

#include <cmath>
#include "DirichletMixtureProfileProcess.h"

// general superclass for all finite process mixtures on site-specific profiles
class DPProfileProcess: public virtual DirichletMixtureProfileProcess	{

	public:

	DPProfileProcess() : kappa(1), movekappa(true), kappaprior(0), Nadd(30), Ninc(3) {}
	virtual ~DPProfileProcess(){}

	double GetKappa()	{
		return kappa;
	}

	protected:

	virtual double Move(double tuning = 1, int n = 1, int nrep = 1);

	virtual double IncrementalDPMove(int nrep, double epsilon);
	double MoveHyper(double tuning, int nrep);
	virtual double MoveKappa(double tuning, int nrep);

	// the following MPI move
	// assumes that master and all slaves are in sync concerning siteprofilesuffstats
	// (which will be the case upon calling GlobalUpdateSiteProfileSuffStat)
	// double GlobalIncrementalDPMove(int nrep);
	// double SlaveIncrementalDPMove();
	double IncrementalDPMove(int nrep);

	// in DP incremental move: Nadd additional components are proposed
	int Nadd;
	int Ninc;

	// static allocation of many component-specific variables
	// such as: profiles, occupancy number
	// basically everything except substitution matrices

	// multinomial 
	virtual double LogProxy(int site, int cat);
	virtual void SampleAlloc();
	void SampleHyper();

	virtual void PriorSampleProfile() 	{
		cerr << "in DPProfileProcess::PriorSampleProfile\n";
		exit(1);
	}

	// kappa has an exponential prior of mean 10
	double LogHyperPrior();
	virtual double LogAllocPrior();

	double kappa;
	bool movekappa;
	int kappaprior;
	// 0 : exponential of mean 20
	// 1 : jeffreys prior 
};

#endif

