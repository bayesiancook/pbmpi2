
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef MULTIGENEPHYLOPROCESS_H
#define MULTIGENEPHYLOPROCESS_H

#include "PhyloProcess.h"
#include "MultiGeneRateProcess.h"
#include "MultiGeneProfileProcess.h"
#include "MultiGeneBranchProcess.h"

class MultiGenePhyloProcess : public virtual PhyloProcess, public virtual MultiGeneRateProcess, public virtual MultiGeneBranchProcess, public virtual MultiGeneProfileProcess	{

	public:

	MultiGenePhyloProcess() {
	}

	virtual ~MultiGenePhyloProcess() {}

	virtual void Create();
	virtual void Delete();

	virtual void ToStream(ostream& os);
	virtual void FromStream(istream& is);
	void SlaveToStream();
	void SlaveFromStream();

	void AllocateAlignments(string datafile);

	virtual void New(int unfold = 1);
	virtual void Open(istream& is, int unfold = 1);
	void SlavePostNew();
	void SlavePostOpen();

	virtual int GetTotalNsite()	{
		return GlobalNsite;
	}

	protected:

	virtual void GlobalWriteSiteRankToStream(ostream& os);
	virtual void GlobalReadSiteRankFromStream(istream& is);
	virtual void SlaveWriteSiteRankToStream();
	virtual void SlaveReadSiteRankFromStream();

    virtual int SpecialSlaveExecute(MESSAGE);

	virtual void GlobalWriteSiteLogLToStream(ostream& os);
	virtual void SlaveWriteSiteLogLToStream();

	virtual void SlaveUpdateSiteRateSuffStat()	{
		MultiGeneRateProcess::SlaveUpdateSiteRateSuffStat();
	}

	virtual void GlobalUpdateSiteProfileSuffStat()	{
		MultiGeneProfileProcess::GlobalUpdateSiteProfileSuffStat();
	}

	virtual void SlaveUpdateSiteProfileSuffStat()	{
		MultiGeneProfileProcess::SlaveUpdateSiteProfileSuffStat();
	}

	virtual void UpdateBranchLengthSuffStat();

	void GlobalSample();
	void SlaveSample();

	virtual void GlobalGeneMove();
	virtual void SlaveGeneMove();

	virtual void GlobalReshuffleSites();
	virtual void SlaveReshuffleSites();

	// re-implement slave functions
	// should dispatch job over all genes,
	// collect results and send back to master

	void SetTree(string treefile)	{
		if (treefile == "None")	{
			tree = new Tree(GetData()->GetTaxonSet());
			if (GetMyid() == 0)	{
				tree->MakeRandomTree();
				GlobalBroadcastTree();
			}
			else	{
				PhyloProcess::SlaveBroadcastTree();
			}
		}
		else	{
			tree = new Tree(treefile);
		}
		tree->RegisterWith(GetData()->GetTaxonSet());
	}


	// void SlaveResetTree();
	void SlaveBroadcastTree();

	void SlaveSetBFFrac();
	void SlaveSetSISFrac();

	void SlaveUnfold();
	void SlaveCollapse();

	void SlaveCollectLogLikelihood();
	void SlaveCollectFullLogLikelihood();
	virtual void ReadPB(int argc, char* argv[]);
	virtual void ReadFullLogL(string name, int burnin, int every, int until);

	void SlaveActivateSumOverRateAllocations();
	void SlaveInactivateSumOverRateAllocations();

	void SlaveActivateZip();
	void SlaveInactivateZip();

	virtual void SlaveComputeTopoBFLogLikelihoodRatio();

	void SlaveUpdateConditionalLikelihoods();
	void SlaveComputeNodeLikelihood(int,int);
	virtual void SlaveGetFullLogLikelihood();

	void SlaveReset(int n, bool v);
	void SlaveMultiply(int n, int m, bool v);
	void SlaveMultiplyByStationaries(int n, bool v);
	void SlaveInitialize(int n, int m, bool v);
	void SlavePropagate(int n, int m, bool v, double t);

	void SlaveProposeMove(int n, double x);
	void SlaveRestoreBranch(int n);

	void SlaveRoot(int n);

	void SlaveBackupTree();
	void SlaveRestoreTree();
	void SlaveSwapTree();

	// in Gibbs.cpp
	void SlaveGibbsSPRScan(int idown, int iup);

	// in NNI.cpp
	void SlavePropagateOverABranch(int);
	virtual void LocalTryNNI(int l, int n, int* br, double* m, double* loglikelihood, int mimick);
	virtual void LocalFinalizeNNI(int n, int* br, int choice, int mimick);

	// in SMC.cpp
	void SlaveSetMinMax();

	// virtual void ReadSiteLogL(string name, int burnin, int every, int until);
	// virtual void PostPred(int ppredtype, string name, int burnin, int every, int until, int rateprior, int profileprior, int rootprior);
    
    void SlaveActivateCPO();
    void SlaveInactivateCPO();
    void SlaveCollectPseudoMarginalLogLikelihood();
};

#endif
