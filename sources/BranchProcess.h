
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

	BranchProcess() : tree(0), blarray(0), swaproot(0), fixbl(0) {}
	virtual ~BranchProcess() {}

	NewickTree* GetLengthTree() {return this;}
	Tree* GetTree() {return tree;}
	Tree* GetTree2() {return tree2;}
	Link* GetRoot() {
		if (swaproot)	{
			return tree2->GetRoot();
		}
		return tree->GetRoot();
	}

	const Link* GetRoot() const {
		if (swaproot)	{
			return tree2->GetRoot();
		}
		return tree->GetRoot();
	}

	Link* GetRoot2() {return tree2->GetRoot();}
	const Link* GetRoot2() const {return tree2->GetRoot();}

	/*
	void GlobalSwapRoot();
	virtual void SlaveSwapRoot()	{
		SwapRoot();
	}

	void SwapRoot()	{
		swaproot = 1 - swaproot;
	}
	*/

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

	virtual int GetBranchLengthSuffStatCount(int index)	{
		cerr << "in BranchProcess::GetBranchLengthSuffStatCount\n";
		exit(1);
		return 0;
	}

	// implements a map<const Branch*, double>

	double LengthSuffStatLogProb();

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

	void CloneTree()	{
		tree2 = new Tree(tree);
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
		}
	}
	virtual void Delete() {
		delete[] blarray;
		delete[] bkarray;
		delete[] bk2array;
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

	Link* GetLink2(int linkindex)	{
		if (! linkindex)	{
			//return GetRoot();
			return 0;
		}
		return GetTree2()->GetLink(linkindex);
	}

	Link* GetLinkForGibbs(int linkindex)	{
		if (! linkindex)	{
			return GetRoot();
		}
		return GetTree()->GetLink(linkindex);
	}

	Link* GetLinkForGibbs2(int linkindex)	{
		if (! linkindex)	{
			return GetRoot2();
		}
		return GetTree2()->GetLink(linkindex);
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

	Link* GetCloneLink(const Link* link)	{
		return GetLinkForGibbs2(link->GetIndex());
	}

	Link* GlobalDetach(Link* down, Link* up);
	void GlobalAttach(Link* down, Link* up, Link* fromdown, Link* fromup);
	virtual void SlaveDetach(int,int);
	virtual void SlaveAttach(int,int,int,int);
	virtual void LocalDetach(int,int);
	virtual void LocalAttach(int,int,int,int);

	Link* GlobalDetach1(Link* down, Link* up);
	void GlobalAttach1(Link* down, Link* up, Link* fromdown, Link* fromup);
	virtual void SlaveDetach1(int,int);
	virtual void SlaveAttach1(int,int,int,int);
	virtual void LocalDetach1(int,int);
	virtual void LocalAttach1(int,int,int,int);

	Link* GlobalDetach2(Link* down, Link* up);
	void GlobalAttach2(Link* down, Link* up, Link* fromdown, Link* fromup);
	virtual void SlaveDetach2(int,int);
	virtual void SlaveAttach2(int,int,int,int);
	virtual void LocalDetach2(int,int);
	virtual void LocalAttach2(int,int,int,int);

	void GlobalKnit(Link*);
	virtual void SlaveKnit();
	void LocalKnit(int arg);

	void GetWeights(Link* from, map<pair<Link*,Link*>,double>& weights, double lambda);
	double WeightedDrawSubTree(double lambda, Link*& down, Link*& up);
	double GetSubTreeWeight(double lambda, Link* down, Link* up);

	Tree* tree;
	Tree* tree2;
	double* blarray;
	double* bkarray;
	double* bk2array;

	Chrono chronolength;

	int swaproot;

	int fixbl;

};

#endif

