
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "StringStreamUtils.h"

#include "Random.h"
#include "PhyloProcess.h"
#include <string>

#include "Parallel.h"
extern MPI_Datatype Propagate_arg;

#include "TexTab.h"


//-------------------------------------------------RG
//-------------------------------------------------
//------------------- NNI Moves -------------------
//-------------------------------------------------
//-------------------------------------------------


double PhyloProcess::GibbsNNI(double tuning, int type){

	if (!type)	{
		tuning=0;
	}
	GlobalRootAtRandom();
	GlobalUpdateConditionalLikelihoods();
	int success =0;
	int moves =0;

	int anumber = rnd::GetRandom().Uniform()*6+1;
	if (anumber > 1)	{
		GlobalKnit(GetRoot());
	}
	if (anumber > 3)	{
		GlobalKnit(GetRoot());
	}
	RecursiveGibbsNNI(GetRoot()->Next()->Out(),tuning,type,success,moves);
	for (int i=anumber%2+1; i; i--)	{	
		GlobalKnit(GetRoot());
	}
	RecursiveGibbsNNI(GetRoot()->Next()->Out(),tuning,type,success,moves);
	for (int i=anumber%2+1; i; i--)	{
		GlobalKnit(GetRoot());
	}
	RecursiveGibbsNNI(GetRoot()->Next()->Out(),tuning,type,success,moves);
	
	return success/(double)moves;
}

// The Local Gibbs sampling over the 3 topologies over the branches is done recursively.
void PhyloProcess::RecursiveGibbsNNI(Link* from, double tuning, int type, int& success, int& moves){

	// Forward
	if ((! from->isLeaf()) and (! from->isRoot()))	{
		if (rnd::GetRandom().Uniform() > 0.5)	{
			GlobalKnit(from);
		}
		for (int i=2; i; i--)	{
			GlobalKnit(from);
			if (! from->Next()->Out()->isLeaf())	{
				success += GlobalNNI(from,tuning,type);
				moves++;
			}
		}
	}
	// Recursion
	if (rnd::GetRandom().Uniform() > 0.5)	{
		GlobalKnit(from);
	}
	for (Link* link=from->Next(); link!=from; link=link->Next())	{
		GlobalPropagateOverABranch(link->Out());
		RecursiveGibbsNNI(link->Out(),tuning,type,success,moves);
	}

	if ((! from->isLeaf()) and (! from->isRoot()))	{
		if (rnd::GetRandom().Uniform() > 0.5)	{
			GlobalKnit(from);
		}
		// Backward
		for (int i=2; i; i--){
			GlobalKnit(from);
			if (! from->Next()->Out()->isLeaf())	{
				success += GlobalNNI(from,tuning,0);
				moves++;
			}
		}
		GlobalPropagateOverABranch(from->Out());
	}
}

// GLOBAL NNI
int PhyloProcess::GlobalNNI(Link* from, double tuning, int type)	{

	// MPI send NNI signal
	MESSAGE signal = NNI;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	ofstream tos((name + ".topo").c_str(),ios_base::app);

	/*
	ostringstream s1, s2, s3, s4;
	GetTree()->ToStreamStandardForm(s1);
	GetTree()->NNIturn(from);	
	GetTree()->ToStreamStandardForm(s2);
	GetTree()->NNIturn(from);	
	GetTree()->ToStreamStandardForm(s3);
	GetTree()->NNIturn(from);	
	GetTree()->ToStreamStandardForm(s4);
	if (s1.str() != s4.str())	{
		cerr << "error in NNI: nni turns do not fall back onto original tree\n";
		exit(1);
	}
	*/

	int n =0;
	if (type)	{
		n = 1 + (int) 5 * rnd::GetRandom().Uniform();
	}
	int args[] = {GetLinkIndex(from),n};

	// MPI send NNI spot
	MPI_Bcast(args,2,MPI_INT,0,MPI_COMM_WORLD);

	// MPI send randomly jittered branch lengths
	Link** branches;
	double logDiffPriorAndHastings = 0.0;
	if (type)	{
		branches = new Link*[n];
		logDiffPriorAndHastings = GlobalSendRandomBranches(from, tuning, branches, n);
	}
	// Initialise the loglikelihood vector
	double* loglikelihood = new double[3];
	loglikelihood[0]=logL;
	loglikelihood[1]=logDiffPriorAndHastings;
	loglikelihood[2]=logDiffPriorAndHastings;

	// MPI receive the new loglikelihood
	double* vec = new double[2];
	MPI_Status stat;
	for(int i=1; i<GetNprocs(); ++i) {
		MPI_Recv(vec,2,MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
		loglikelihood[1]+=vec[0];
		loglikelihood[2]+=vec[1];
	}
	delete[] vec;

	// Sample a configuration
	int choice = rnd::GetRandom().DrawFromLogDiscreteDistribution(loglikelihood, 3);

	// NNI broadcast choice
	MPI_Bcast(&choice,1,MPI_INT,0,MPI_COMM_WORLD);
	bool success = (choice != 0);

	// Update the logL of the model
	logL=loglikelihood[choice];
	if (success)	{
		logL-=logDiffPriorAndHastings;
	}
	delete[] loglikelihood;

	// Restore the branch length if the initial configuration was chosen
	if (type)	{
		if (! success)	{
			for(int i=0; i<n; ++i){
				Restore(branches[i]->GetBranch());
			}
		}
		delete[] branches;
	}

	if (success)	{
		for(int i = choice;i>0;i--){
			GetTree()->NNIturn(from);
			GetTree2()->NNIturn(GetCloneLink(from));
		}
	}

	/*
	ostringstream s5;
	GetTree()->ToStreamStandardForm(s5);

	if (s5.str() == s1.str())	{
		tos << "nni\t" << s1.str() << '\t' << s2.str() << '\t' << "reject\n";
		tos << "nni\t" << s1.str() << '\t' << s3.str() << '\t' << "reject\n";
	}
	if (s5.str() == s2.str())	{
		tos << "nni\t" << s1.str() << '\t' << s3.str() << '\t' << "reject\n";
		tos << "nni\t" << s1.str() << '\t' << s2.str() << '\t' << "accept\n";
	}
	if (s5.str() == s3.str())	{
		tos << "nni\t" << s1.str() << '\t' << s2.str() << '\t' << "reject\n";
		tos << "nni\t" << s1.str() << '\t' << s3.str() << '\t' << "accept\n";
	}
	*/
	
	return success;
}

// This function fill the table of link with n link pointing on the branches that will move
// It also send the movements to each of the slaves.
double PhyloProcess::GlobalSendRandomBranches(Link* from, double tuning, Link**& branches, int n){

	// Here we draw the n branches
	int* br = new int[n];
	rnd::GetRandom().DrawFromUrn (br , n, 5);

	double* m = new double[n];
	double logDiffPriorAndHastings = 0;
	for(int i = 0; i<n; i++){
		switch(br[i]){
			case 0:	branches[i] = from->Next()->Out()->Next()->Next();
				break;
			case 1: branches[i] = from->Next()->Out()->Next();
				break;
			case 2: branches[i] = from->Next()->Next();
				break;
			case 3: branches[i] = from;
				break;
			case 4: branches[i] = from->Next();
				break;
			default:cerr << "Error in PhyloProcess::SendRandomBranches " << br[i] << ' ';exit(1);
		}
		m[i] = tuning * (rnd::GetRandom().Uniform() - 0.5);
		logDiffPriorAndHastings += m[i];
		Branch* b = branches[i]->GetBranch();
		logDiffPriorAndHastings -= LogBranchLengthPrior(b);
		// logDiffPriorAndHastings += LogBranchLengthPrior(b);
		MoveBranch(b,m[i]);
		logDiffPriorAndHastings += LogBranchLengthPrior(b);
		// logDiffPriorAndHastings -= LogBranchLengthPrior(b);
		br[i]=GetLinkIndex(branches[i]);
	}
	MPI_Bcast(br,n,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(m,n,MPI_DOUBLE,0,MPI_COMM_WORLD);

	delete[] m;
	delete[] br;
	return logDiffPriorAndHastings;
}

void PhyloProcess::SlaveNNI()	{

	// MPI receive arguments
	int arg[2];
	MPI_Bcast(arg,2,MPI_INT,0,MPI_COMM_WORLD);

	int l = arg[0];
	int n = arg[1];
	int* br;
	double* m;
	if (n)	{
		br = new int[n];
		m = new double[n];

		// MPI receive new branch lengths
		MPI_Bcast(br,n,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(m,n,MPI_DOUBLE,0,MPI_COMM_WORLD);
	}

	// local computations 1
	double* loglikelihood = new double[2];
	loglikelihood[0] = 0;
	loglikelihood[1] = 0;

	// this can be overloaded by MultiGenePhyloProcess
	// distribute over genes and add up log likelihoods
	LocalTryNNI(l,n,br,m,loglikelihood,0);

	// MPI Send likelihood components
	MPI_Send(loglikelihood,2,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);

	// MPI Receive choice
	int choice;
	MPI_Bcast(&choice,1,MPI_INT,0,MPI_COMM_WORLD);

	// local computations 2
	// this can be overloaded by MultiGenePhyloProcess
	// just dispatch order across genes
	LocalFinalizeNNI(n,br,choice,0);

	if (n)	{
		delete[] br;
		delete[] m;
	}
	delete[] loglikelihood;
}

void PhyloProcess::LocalTryNNI(int l, int n, int* br, double* m, double* loglikelihood, int mimick)	{

	Link* from = GetLinkForGibbs(l);
	Link* up = from->Next();
	bknnifrom = from;
	bknniup = up;

	if (n)	{
		for(int i=0; i<n; ++i){
			Link* link = GetLinkForGibbs(br[i]);
			MoveBranch(link->GetBranch(),m[i]);
			if ((link!=up) && (! mimick))	{
				PropagateOverABranch(link);
			}
		}
		if (! mimick)	{
			PropagateOverABranch(up);
		}
	}
	
	GetTree()->NNIturn(from);
	GetTree2()->NNIturn(GetCloneLink(from));
	if (! mimick)	{
		PropagateOverABranch(from->Next());
		loglikelihood[0] += ComputeNodeLikelihood(from,-1);
	}

	GetTree()->NNIturn(from);
	GetTree2()->NNIturn(GetCloneLink(from));
	if (! mimick)	{
		PropagateOverABranch(from->Next());
		loglikelihood[1] += ComputeNodeLikelihood(from,-1);
	}
}

void PhyloProcess::LocalFinalizeNNI(int n, int* br, int choice, int mimick)	{

	// recover parameters
	Link* from = bknnifrom;
	Link* up = bknniup;

	if (n)	{
		if (choice == 0)	{
			for(int i=0; i<n; ++i){
				Link* link = GetLinkForGibbs(br[i]);
				Restore(link->GetBranch());
				if ((link!=up) && (!mimick))	{
					PropagateOverABranch(link);
				}
			}
		}
	}

	if (choice == 1)	{
		GetTree()->NNIturn(from);	
		GetTree2()->NNIturn(GetCloneLink(from));
	}
	if (choice != 2)	{
		GetTree()->NNIturn(from);
		GetTree2()->NNIturn(GetCloneLink(from));
		if (! mimick)	{
			PropagateOverABranch(up);
		}
	}
}

void PhyloProcess::GlobalPropagateOverABranch(Link* from)	{
	if (! from->isLeaf())	{
		MESSAGE signal = BRANCHPROPAGATE;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
		int args[] = {GetLinkIndex(from)};
		MPI_Bcast(args,1,MPI_INT,0,MPI_COMM_WORLD);
	}
}

void PhyloProcess::SlavePropagateOverABranch(int l)	{

	const Link* from = GetLinkForGibbs(l);
	PropagateOverABranch(from);
}

// Update the conditional Likelihoods associated to from
// Assume that everything else is updated
void PhyloProcess::PropagateOverABranch(const Link* from){


	if (from->Out()->isLeaf())	{
		Initialize(condlmap[0],GetData(from->Out()),false);
	}
	else{
		Reset(condlmap[0],false);
		for (const Link* link=from->Out()->Next(); link!=from->Out(); link=link->Next())	{
			if (!link->isRoot())	{
				Multiply(GetConditionalLikelihoodVector(link),condlmap[0],false);
			}
		}
		Offset(condlmap[0],false);
	}
	Propagate(condlmap[0],GetConditionalLikelihoodVector(from),GetLength(from->GetBranch()),false);
}

