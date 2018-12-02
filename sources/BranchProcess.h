
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef BRANCH_H
#define BRANCH_H

#include "Tree.h"
#include "Chrono.h"
#include "MPIModule.h"

class BranchProcess : public NewickTree, public virtual MPIModule {

	public:

	BranchProcess() : tree(0), blarray(0), fixbl(0), BranchNcat(2), branchalloc(0) {}
	virtual ~BranchProcess() {}

	int GetBranchNcat()	{
		return BranchNcat;
	}

	NewickTree* GetLengthTree() {return this;}
	Tree* GetTree() {return tree;}
	Link* GetRoot() {
		return tree->GetRoot();
	}

	const Link* GetRoot() const {
		return tree->GetRoot();
	}

	string GetBranchName(const Link* link) const;
	string GetNodeName(const Link* link) const;

	void SetFixBL(bool in)	{
		fixbl = in;
	}

	bool FixBL()	{
		return fixbl;
	}

	double GetMinLength(const Link* from, double l)	{

		double min = -1;
		const Link* link = from;
		do	{
			if (! link->isRoot())	{
				double tmp = GetLength(link->GetBranch());
				if (tmp != l)	{
					if ((min == -1) || (min > tmp))	{
						min = tmp;
					}
				}
			}
			link = link->Next();
		}
		while (link != from);
		if (! min)	{
			cerr << "error in BranchProcess::GetMinLength\n";
			exit(1);
		}
		return min;
	}

	double GetLength(const Branch* branch) {return branch ? blarray[branch->GetIndex()] : 0;}
	double GetLength(const Branch* branch) const {return branch ? blarray[branch->GetIndex()] : 0;}
	void SetLength(const Branch* branch, double inlength) {
		if (! branch)	{
			cerr << "error in branch process: null branch\n";
			exit(1);
		}
		if (! blarray)	{
			cerr << "array not created\n";
			exit(1);
		}
		blarray[branch->GetIndex()] = inlength;
	}

	void SetBranchLengths(const double* inblarray)	{
		for (int j=0; j<GetNbranch(); j++)	{
			blarray[j] = inblarray[j];
		}
	}

	const double* GetBranchLengths()	{
		return blarray;
	}

	void BackupLength();
	void RestoreLength();
	// void SwapLength();

	// Move functions

	virtual double BranchProcessMove(double tuning = 1, int nrep=1) = 0;

	double ProposeMove(const Branch* branch, double tuning);
	void MoveBranch(const Branch* branch, double factor);
	void RestoreBranch(const Branch* branch);

	virtual double GetNormFactor()	{
		cerr << "in BranchProcess::GetNormFactor\n";
		exit(1);
		return 0;
	}

	// return number of branches
	int MoveAllBranches(double factor);
	int RecursiveMoveAllBranches(const Link* from, double e);

	double GetAllocTotalLength(int cat)	{
		double tot = 0;
		for (int j=1; j<GetNbranch(); j++)	{
			if (branchalloc[j] == cat)	{
				tot += blarray[j];
			}
		}
		return tot;
	}

	double GetTotalLength()	{
		return RecursiveTotalLength(GetRoot());
	}

	double GetRenormTotalLength()	{
		return RecursiveTotalLength(GetRoot()) * GetNormFactor();
	}

	void RenormalizeBranchLengths()	{
		double tmp = GetNormFactor();
		RecursiveNormalizeBranchLengths(GetRoot(),tmp);
	}

	void DenormalizeBranchLengths()	{
		double tmp = GetNormFactor();
		RecursiveNormalizeBranchLengths(GetRoot(),1.0 / tmp);
	}

	void RecursiveNormalizeBranchLengths(const Link* from, double factor);
	
	double RecursiveTotalLength(const Link* from);

	void SetLengthsFromNames()	{
		if (blarray)	{
			RecursiveSetLengthsFromNames(GetRoot());
		}
		else	{
			cerr << "set lengths from names called without blarray\n";
		}
	}
	
	void SetNamesFromLengths()	{
		if (blarray)	{
			RecursiveSetNamesFromLengths(GetRoot());
		}
		else	{
			cerr << "set names from length called without blarray\n";
		}
	}

	void RecursiveSetLengthsFromNames(const Link* from);
	void RecursiveSetNamesFromLengths(const Link* from);

	// sampling all params (lengths and hyperparams) from prior
	void PriorSampleLength()	{
		if (! fixbl)	{
			PriorSampleLengthHyperParameters();
			RecursiveSampleLength(GetRoot());
		}
	}

	// sampling all params (lengths and hyperparams) for starting mcmc
	void SampleLength()	{
		if (! fixbl)	{
			SampleLengthHyperParameters();
			RecursiveSampleLength(GetRoot());
		}
	}

	// recursive sampling of all branches 
	void RecursiveSampleLength(const Link* from);

	// sampling one branch length (to be implemented in specialized classes)
	virtual void SampleBranchLength(const Branch* branch)	{
		cerr << "in BranchProcess::SampleBranchLength\n";
		exit(1);
	}

	// sampling hyperparameters from prior
	virtual void PriorSampleLengthHyperParameters()	{
		cerr << "in branchProcess::PriorSampleLengthHyperParameters\n";
		exit(1);
	}

	// sampling hyperparameters for starting mcmc
	virtual void SampleLengthHyperParameters()	{
		cerr << "in branchProcess::SampleLengthHyperParameters\n";
		exit(1);
	}

	virtual double LogBranchProcessPrior()	{
		return LogLengthPrior() + LogLengthHyperPrior();
	}

	// prior on branch lengths, given hyperparameters
	virtual double LogLengthPrior()	{
		return RecursiveLogLengthPrior(GetRoot());
	}

	double RecursiveLogLengthPrior(const Link* from);

	// prior on the length of one specific branch
	virtual double LogBranchLengthPrior(const Branch* branch)	{
		cerr << "in BranchProcess::LogBranchLengthPrior\n";
		exit(1);
	}

	// prior on hyperparameters
	virtual double LogLengthHyperPrior()	{
		cerr << "in BranchProcess::LogHyperPrior\n";
		exit(1);
		return 0;
	}

	virtual void GlobalUpdateBranchLengthSuffStat()	{
		cerr << "in BranchProcess::GlobalUpdateBranchLengthSuffStat\n";
		exit(1);
	}

	virtual void SlaveUpdateBranchLengthSuffStat()	{
		cerr << "in BranchProcess::SlaveUpdateBranchLengthSuffStat\n";
		exit(1);
	}

	virtual void UpdateBranchLengthSuffStat()	{
		cerr << "in BranchProcess::UpdateBranchLengthSuffStat\n";
		exit(1);
	}


	virtual double GetBranchLengthSuffStatBeta(int index)	{
		cerr << "in BranchProcess::GetBranchLengthSuffStatBeta\n";
		exit(1);
		return 0;
	}

	virtual double GetBranchLengthSuffStatCount(int index)	{
		cerr << "in BranchProcess::GetBranchLengthSuffStatCount\n";
		exit(1);
		return 0;
	}

	void ResetBranchAlloc();

	void SetBranchAlloc(int branchindex, int alloc)	{
		branchalloc[branchindex] = alloc;
	}

	virtual double GetBranchScaling(int index)	{
		cerr << "error: in BranchProcess::GetBranchScaling\n";
		exit(1);
		return 1.0;
	}

	virtual void SetBranchScaling(double scale, int index) {
		cerr << "error: in BranchProcess::SetBranchScaling\n";
		exit(1);
	}

	virtual void RescaleBranchPrior(double factor, int index)	{
		SetBranchScaling(factor*GetBranchScaling(index),index);
	}
	/*
	void GlobalSetBranchAlloc(int branchindex, int alloc);
	void SlaveSetBranchAlloc();

	void GlobalRescaleBranchPrior(double factor, int index);
	void SlaveRescaleBranchPrior();
	*/

	/*
	virtual double GetBranchScalingFactor(int index)	{
		return 1.0;
	}
	*/
	
	// SetBranchAlloc(string taxon1, strong taxon2, int alloc);

	// implements a map<const Branch*, double>

	double LengthSuffStatLogProb();
    void EM_UpdateBranchLengths();

	virtual void ToStream(ostream& os)	{
		cerr << "in BranchProcess::ToStream\n";
		exit(1);
	}
	virtual void FromStream(istream& is)	{
		cerr << "in BranchProcess::FromStream\n";
		exit(1);
	}

	virtual const TaxonSet* GetTaxonSet() const = 0;

	// protected:

	int GetNbranch()	{
		return tree->GetNbranch();
	}

	int GetNnode()	{
		return tree->GetNnode();
	}

	int GetNlink()	{
		return tree->GetNlink();
	}

	void SetTree(Tree* intree)	{
		cerr << "in BranchProcess::SetTree\n";
		exit(1);
		tree = intree;
	}

	virtual void Create() {
		if (! blarray)	{
			if (! tree)	{
				cerr << "error: BranchProcess::tree has not been initialized\n";
				exit(1);
			}
			blarray = new double[GetNbranch()];
			bkarray = new double[GetNbranch()];
			bk2array = new double[GetNbranch()];
			blarray[0] = 0;
			SetLengthsFromNames();
			branchalloc = new int[GetNbranch()];
			for (int j=0; j<GetNbranch(); j++)	{
				branchalloc[j] = 0;
			}
		}
	}

	virtual void Delete() {
		if (blarray)	{
			delete[] blarray;
			blarray = 0;
			delete[] bkarray;
			delete[] bk2array;
			delete[] branchalloc;
		}
	}

	// translation tables : from pointers of type Link* Branch* and Node* to their index and vice versa
	// this translation is built when the Tree::RegisterWithTaxonSet method is called (in the model, in PB.cpp)
	Link* GetLink(int linkindex)	{
		if (! linkindex)	{
			//return GetRoot();
			return 0;
		}
		return GetTree()->GetLink(linkindex);
	}

	Link* GetLinkForGibbs(int linkindex)	{
		if (! linkindex)	{
			return GetRoot();
		}
		return GetTree()->GetLink(linkindex);
	}

	const Branch* GetBranch(int branchindex)	{
		return GetTree()->GetBranch(branchindex);
	}

	const Node* GetNode(int nodeindex)	{
		return GetTree()->GetNode(nodeindex);
	}


	int GetLinkIndex(const Link* link)	{
		return link ? link->GetIndex() : 0;
	}

	int GetBranchIndex(const Branch* branch)	{
		if (! branch)	{
			return 0;
			/*
			cerr << "error in get branch index\n";
			exit(1);
			*/
		}
		return branch->GetIndex();
	}

	int GetNodeIndex(const Node* node)	{
		return node->GetIndex();
	}

	Link* GlobalDetach(Link* down, Link* up);
	void GlobalAttach(Link* down, Link* up, Link* fromdown, Link* fromup);
	virtual void SlaveDetach(int,int);
	virtual void SlaveAttach(int,int,int,int);
	virtual void LocalDetach(int,int);
	virtual void LocalAttach(int,int,int,int);

	void GlobalKnit(Link*);
	virtual void SlaveKnit();
	void LocalKnit(int arg);

	void GetWeights(Link* from, map<pair<Link*,Link*>,double>& weights, double lambda);
	double WeightedDrawSubTree(double lambda, Link*& down, Link*& up);
	double GetSubTreeWeight(double lambda, Link* down, Link* up);

	Tree* tree;
	double* blarray;
	double* bkarray;
	double* bk2array;

	Chrono chronolength;

	int fixbl;
	int BranchNcat;
	int* branchalloc;

};

#endif

