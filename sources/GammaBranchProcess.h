
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef GAMMABRANCH_H
#define GAMMABRANCH_H

#include "TaxonSet.h"
#include "BranchProcess.h"

class GammaBranchProcess : public virtual BranchProcess	{

	public:

	GammaBranchProcess() : branchmean(0), branchrelvar(0), betaprior(0), hierarchicallengthprior(0) {}
	virtual ~GammaBranchProcess() {}

	virtual void Create();
	virtual void Delete();

	double LogBranchLengthPrior(const Branch* branch);
	virtual double LogLengthHyperPrior();

	void SampleBranchLength(const Branch* branch);
	virtual void PriorSampleLengthHyperParameters();
	virtual void SampleLengthHyperParameters();

	void ToStreamWithLengths(ostream& os, const Link* from);

	// conjugate sampling

	virtual double BranchProcessMove(double tuning = 1, int nrep=1)	{
		return Move(tuning,nrep);
	}

	virtual double Move(double tuning = 1, int nrep=1);
	virtual double NonMPIMove(double tuning = 1, int nrep=1);
	virtual double MoveLengthHyperParameters(double tuning, int nrep);

	double MPIMoveBranchLengths();
	double NonMPIMoveBranchLengths();

	double GetMeanLengthRelVar()	{

		if (! hierarchicallengthprior)	{
			cerr << "error: in GammaBranchProcess::GetMeanLengthRelVar, with non hierarchical prior\n";
			exit(1);
		}
		double tot = 0;
		for (int j=1; j<GetNbranch(); j++)	{
			tot += branchrelvar[j];
		}
		tot /= (GetNbranch() - 1);
		return tot;
	}

	void ToStream(ostream& os);
	void FromStream(istream& is);

	void SetHyperParameters(double* inmean, double* inrelvar)	{
		for (int j=1; j<GetNbranch(); j++)	{
			branchmean[j] = inmean[j];
			branchrelvar[j] = inrelvar[j];
		}
	}

	double branchalpha;
	double branchbeta;
	double* branchmean;
	double* branchrelvar;

	int betaprior;
	int hierarchicallengthprior;
};

#endif

