
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

/*
class MultiGeneMPIModule	{

	public:

	MultiGeneMPIModule() : Ngene(0)	{}
	virtual ~MultiGeneMPIModule() {}

	void AllocateAlignments(string datafile, string treefile);

	virtual StateSpace* GetStateSpace() {return statespace;}
	virtual int GetNstate() {return Nstate;}

	void SetMPI(int inmyid, int innprocs)	{
		myid = inmyid;
		nprocs = innprocs;
	}

	protected:

	int Nstate;
	StateSpace* statespace;

	int Ngene;
	int* genealloc;
	int* genesize;
	string* genename;

	int GlobalNsite;
	int* globalnsite;

	int myid;
	int nprocs;

};
*/

class MultiGenePhyloProcess : public virtual PhyloProcess	{

	public:

	MultiGenePhyloProcess() : genelnL(0) {}
	virtual ~MultiGenePhyloProcess() {}

	virtual void Create();
	virtual void Delete();

	void AllocateAlignments(string datafile, string treefile);

	virtual void New();

	protected:

        virtual int SpecialSlaveExecute(MESSAGE);


	void GlobalSample();
	void SlaveSample();

	virtual void GlobalGeneMove();
	virtual void SlaveGeneMove();

	// re-implement slave functions
	// should dispatch job over all genes,
	// collect results and send back to master

	void SlaveUpdateSiteRateSuffStat();

	void SlaveBroadcastTree();

	void SlaveUnfold();
	void SlaveCollapse();

	void SlaveActivateSumOverRateAllocations();
	void SlaveInactivateSumOverRateAllocations();

	void SlaveActivateZip();
	void SlaveInactivateZip();

	void SlaveUpdateConditionalLikelihoods();
	void SlaveComputeNodeLikelihood(int,int);
	void SlaveGetFullLogLikelihood();

	void SlaveReset(int n, bool v);
	void SlaveMultiply(int n, int m, bool v);
	void SlaveMultiplyByStationaries(int n, bool v);
	void SlaveInitialize(int n, int m, bool v);
	void SlavePropagate(int n, int m, bool v, double t);

	void SlaveProposeMove(int n, double x);
	void SlaveRestore(int n);

	void SlaveRoot(int n);

	void SlaveBackupTree();
	void SlaveRestoreTree();
	void SlaveSwapTree();

	// in Gibbs.cpp
	void SlaveGibbsSPRScan(int idown, int iup);

	// in NNI.cpp
	/*
	void SlaveKnit();
	void SlaveNNI(int,int);
	void SlavePropagateOverABranch(int);
	*/

	// in SMC.cpp
	void SlaveSetMinMax();

	// should be defined at the level of the BranchProcess
	void SlaveDetach(int,int);
	void SlaveAttach(int,int,int,int);
	void SlaveDetach1(int,int);
	void SlaveAttach1(int,int,int,int);
	void SlaveDetach2(int,int);
	void SlaveAttach2(int,int,int,int);
	void SlaveSwapRoot();


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
		CloneTree();
		tree2->RegisterWith(GetData()->GetTaxonSet());
	}


	PhyloProcess** process;
	double* genelnL;
	double* tmpgenelnL;

	int Ngene;
	int* genealloc;
	int* genesize;
	string* genename;

	int GlobalNsite;
	int* globalnsite;


};

#endif
