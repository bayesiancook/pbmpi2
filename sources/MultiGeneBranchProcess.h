/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef MULTIGENEBRANCHPROCESS_H
#define MULTIGENEBRANCHPROCESS_H

#include "MultiGeneMPIModule.h"
#include "GammaBranchProcess.h"

class MultiGeneBranchProcess : public virtual GammaBranchProcess, public virtual MultiGeneMPIModule	{

	public:

	MultiGeneBranchProcess() : geneblarray(0), mappsuffstat(1) {}
	virtual ~MultiGeneBranchProcess() {}

	virtual void Create();
	virtual void Delete();

	virtual void ToStream(ostream& os);
	virtual void FromStream(istream& is);

	void SetGlobalBranchLengths(int in)	{
		globalbl = in;
		if (! in)	{
			hierarchicallengthprior = 1;
		}
	}

	int GlobalBranchLengths()	{
		return globalbl;
	}

	double GetMeanTotalLength()	{
		if (globalbl)	{
			return GetTotalLength();
		}
		return GlobalGetMeanTotalLength();
	}

	double GlobalGetMeanTotalLength();
	void SlaveGetMeanTotalLength();

	void GlobalCollectGeneBranchLengths();
	void SlaveCollectGeneBranchLengths();

	void GlobalCollectGeneLengthMappingSuffStat();
	void SlaveCollectGeneLengthMappingSuffStat();

	void SampleLengthHyperParameters();
	void PriorSampleLengthHyperParameters();

	double LogLengthPrior();
	double LogLengthHyperPrior();
	double LogLengthHyperPrior(int j);
	double LogLengthHyperHyperPrior();
	double LogGeneLengthSuffStatPrior();
	double LogGeneLengthSuffStatPrior(int j);
	double LogGeneLengthMarginalSuffStatProb();
	double LogGeneLengthMarginalSuffStatProb(int j);
	void ComputeGeneLengthSuffStat();

	double Move(double tuning = 1, int nrep = 1);
	double MoveGeneLengthHyperParameters(double tuning = 1, int nrep = 1);
	double MoveGeneLengthHyperHyperParameters(double tuning = 1, int nrep = 1);

	void SlaveDetach(int,int);
	void SlaveAttach(int,int,int,int);
	void SlaveKnit();

	GammaBranchProcess* GetBranchProcess(int gene)	{

		GammaBranchProcess* tmp = dynamic_cast<GammaBranchProcess*>(process[gene]);
		if (!tmp)	{
			cerr << "error in GetBranchProcess\n";
			exit(1);
		}
		return tmp;
	}

	int globalbl;
	int mappsuffstat;

	double** geneblarray;
	double** tmpgeneblarray;
	double* allocgeneblarray;
	double* alloctmpgeneblarray;

	double** geneblcount;
	double** tmpgeneblcount;
	double* allocgeneblcount;
	double* alloctmpgeneblcount;

	double** geneblbeta;
	double** tmpgeneblbeta;
	double* allocgeneblbeta;
	double* alloctmpgeneblbeta;

	double* totlength;
	double* totloglength;

	double meanbranchmean;
	double relvarbranchmean;
	double meanbranchrelvar;
	double relvarbranchrelvar;
};

#endif
