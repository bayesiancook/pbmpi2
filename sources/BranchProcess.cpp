
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#include "BranchProcess.h"
#include "Random.h"
#include "Parallel.h"

#include <sstream>
#include <iomanip>

string BranchProcess::GetNodeName(const Link* link)  const	{
	if (link->isLeaf())	{
		return link->GetNode()->GetName();
	}
	else	{
		return "";
	}
}

string BranchProcess::GetBranchName(const Link* link) const	{
	if (link->isRoot())	{
		return "";
	}
	else	{
		ostringstream s;
		s << GetLength(link->GetBranch());
		return s.str();
	}
}	

void BranchProcess::BackupLength()	{

	for (int i=0; i<GetNbranch(); i++)	{
		bk2array[i] = blarray[i];
	}
}

void BranchProcess::RestoreLength()	{

	for (int i=0; i<GetNbranch(); i++)	{
		blarray[i] = bk2array[i];
	}
}

void BranchProcess::ResetBranchAlloc()	{

	for (int i=0; i<GetNbranch(); i++)	{
		branchalloc[i] = 0;
	}
}

/*
void BranchProcess::SwapLength()	{

	for (int i=0; i<GetNbranch(); i++)	{
		double tmp = bkarray[i];
		blarray[i] = bk2array[i];
		bk2array[i] = tmp;
	}
}
*/

/*
void BranchProcess::SetBranchAlloc(string taxon1, string taxon2, int alloc)	{

	Link* link = tree->GetLCA(taxon1,taxon2);
	if (!link)	{
		cerr << "error in GammaBranchProcess::SetBranchAlloc: did not find LCA of " << taxon1 << " and " << taxon2 << '\n';
		exit(1);
	}
	if (link->isRoot())	{
		cerr << "error in GammaBranchProcess::SetBranchAlloc: LCA of " << taxon1 << " and " << taxon2 << " is root\n";
		exit(1);
	}
	branchalloc[link->GetBranch()->GetIndex()] = alloc;
}
*/

/*
void BranchProcess::GlobalSetBranchAlloc(int branchindex, int alloc)	{

	MESSAGE signal = SETBRANCHALLOC;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&branchindex,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&alloc,1,MPI_INT,0,MPI_COMM_WORLD);
	branchalloc[branchindex] = alloc;
}


void BranchProcess::SlaveSetBranchAlloc()	{

	int branchindex,alloc;
	MPI_Bcast(&branchindex,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&alloc,1,MPI_INT,0,MPI_COMM_WORLD);
	branchalloc[branchindex] = alloc;
}
*/

/*
void BranchProcess::GlobalRescaleBranchPrior(double factor, int alloc)	{

	MESSAGE signal = RESCALEBRANCH;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&factor,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&alloc,1,MPI_INT,0,MPI_COMM_WORLD);
	RescaleBranchPrior(factor,alloc);
}


void BranchProcess::SlaveRescaleBranchPrior()	{

	double factor;
	int alloc;
	MPI_Bcast(&factor,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&alloc,1,MPI_INT,0,MPI_COMM_WORLD);
	RescaleBranchPrior(factor,alloc);
}
*/


void BranchProcess::RestoreBranch(const Branch* branch)	{
	if (! branch)	{
		cerr << "error in branchprocess::Movebranch: null branch\n";
		exit(1);
	}
	int index = branch->GetIndex();
	blarray[index] = bkarray[index];
}

void BranchProcess::MoveBranch(const Branch* branch, double m)	{

	if (! branch)	{
		cerr << "error in branchprocess::Movebranch: null branch\n";
		exit(1);
	}
	int index = branch->GetIndex();
	bkarray[index] = blarray[index];
	blarray[index] *= exp(m);
}

double BranchProcess::ProposeMove(const Branch* branch, double tuning)	{

	if (! branch)	{
		cerr << "error in branchprocess::Movebranch: null branch\n";
		exit(1);
	}
	int index = branch->GetIndex();
	bkarray[index] = blarray[index];
	double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
	blarray[index] *= exp(m);
	return m;
}

int BranchProcess::MoveAllBranches(double e)	{
	return RecursiveMoveAllBranches(GetRoot(),e);
}

int BranchProcess::RecursiveMoveAllBranches(const Link* from, double e)	{

	int n = 0;
	if (!from->isRoot())	{
		int index = from->GetBranch()->GetIndex();
		blarray[index] *= e;
		n++;
	}
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		n += RecursiveMoveAllBranches(link->Out(),e);
	}
	return n;
}

double BranchProcess::RecursiveLogLengthPrior(const Link* from)	{

	double total = 0;
	if (! from->isRoot())	{
		total += LogBranchLengthPrior(from->GetBranch());
	}
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		total += RecursiveLogLengthPrior(link->Out());
	}
	return total;
}

void BranchProcess::RecursiveSampleLength(const Link* from)	{

	if (! from->isRoot())	{
		SampleBranchLength(from->GetBranch());
	}
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		RecursiveSampleLength(link->Out());
	}
}

void BranchProcess::RecursiveNormalizeBranchLengths(const Link* from, double factor)	{

	if (! from->isRoot())	{
		SetLength(from->GetBranch(),GetLength(from->GetBranch())*factor);
	}
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		RecursiveNormalizeBranchLengths(link->Out(),factor);
	}
}

double BranchProcess::RecursiveTotalLength(const Link* from)	{
	double total = 0;
	if (! from->isRoot())	{
		total += GetLength(from->GetBranch());
	}
	else	{
		if (GetLength(from->GetBranch()))	{
			cerr << "error: non null branch length for root\n";
			cerr << GetLength(from->GetBranch()) << '\n';
			exit(1);
		}
	}
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		total += RecursiveTotalLength(link->Out());
	}
	return total;
}

void BranchProcess::RecursiveSetLengthsFromNames(const Link* from)	{
	if (! from->isRoot())	{
		double l = atof(from->GetBranch()->GetName().c_str());
		if (l <= 0)	{
			l = 0.1;
			/*
			cerr << "error in BranchProcess::SetLengthsFromFile : " << l << '\n';
			exit(1);
			*/
		}
		SetLength(from->GetBranch(),l);
	}
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		RecursiveSetLengthsFromNames(link->Out());
	}
}

void BranchProcess::RecursiveSetNamesFromLengths(const Link* from)	{

	if (! from->isRoot())	{
		ostringstream s;
		s.precision(15);
		s << GetLength(from->GetBranch());
		from->GetBranch()->SetName(s.str());
	}
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		RecursiveSetNamesFromLengths(link->Out());
	}
}

double BranchProcess::LengthSuffStatLogProb()	{

	double total = 0;
	for (int i=1; i<GetNbranch(); i++)	{
		total += GetBranchLengthSuffStatCount(i) * log(blarray[i]);
		total -= GetBranchLengthSuffStatBeta(i) * blarray[i];
	}
	return total;
}

void BranchProcess::EM_UpdateBranchLengths()    {
    for (int j=1; j<GetNbranch(); j++)  {
        blarray[j] = GetBranchLengthSuffStatCount(j) / GetBranchLengthSuffStatBeta(j);
    }
}

void BranchProcess::GlobalKnit(Link* from)	{

	MESSAGE signal = KNIT;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	int args[] = {GetLinkIndex(from)};
	MPI_Bcast(args,1,MPI_INT,0,MPI_COMM_WORLD);
	from->Knit();
}

void BranchProcess::SlaveKnit()	{
	int arg;
	MPI_Bcast(&arg,1,MPI_INT,0,MPI_COMM_WORLD);
	LocalKnit(arg);
}

void BranchProcess::LocalKnit(int arg)	{
	GetLinkForGibbs(arg)->Knit();
}

Link* BranchProcess::GlobalDetach(Link* down, Link* up)	{

	if (GetNprocs() > 1)	{
		// MPI
		// master and all slaved should call 
		// GetTree()->Detach(down,up,fromdown,fromup);
		// but message passing will again  use link to index, then index to link, translations.
		MESSAGE signal = DETACH;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
		int args[] = {GetLinkIndex(down),GetLinkIndex(up)};
		// int args[] = {down->GetIndex(),up->GetIndex(),fromdown->GetIndex(),fromup->GetIndex()};
		MPI_Bcast(args,2,MPI_INT,0,MPI_COMM_WORLD);
	}
	return GetTree()->Detach(down,up);
}

void BranchProcess::GlobalAttach(Link* down, Link* up, Link* fromdown, Link* fromup)	{

	if (GetNprocs() > 1)	{
		// MPI
		// same thing as for detach
		MESSAGE signal = ATTACH;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
		int args[] = {GetLinkIndex(down),GetLinkIndex(up),GetLinkIndex(fromdown),GetLinkIndex(fromup)};
		// int args[] = {down->GetIndex(),up->GetIndex(),fromdown->GetIndex(),fromup->GetIndex()};
		MPI_Bcast(args,4,MPI_INT,0,MPI_COMM_WORLD);
	}
	GetTree()->Attach(down,up,fromdown,fromup);
}


void BranchProcess::SlaveDetach(int n,int m) {
	LocalDetach(n,m);
}

void BranchProcess::LocalDetach(int n,int m) {
	Link* down = GetLinkForGibbs(n);
	Link* up = GetLinkForGibbs(m);
	GetTree()->Detach(down,up);
}

void BranchProcess::SlaveAttach(int n,int m,int p,int q) {
	LocalAttach(n,m,p,q);
}

void BranchProcess::LocalAttach(int n,int m,int p,int q) {
	Link* down = GetLinkForGibbs(n);
	Link* up = GetLinkForGibbs(m);
	Link* fromdown = GetLinkForGibbs(p);
	Link* fromup = GetLinkForGibbs(q);
	GetTree()->Attach(down,up,fromdown,fromup);
}

void BranchProcess::GetWeights(Link* from, map<pair<Link*,Link*>,double>& weights, double lambda)	{

	for (Link* link=from->Next(); link!=from; link=link->Next())	{
		double w = exp(-lambda * GetMinLength(link,GetLength(link->GetBranch())));
		weights[pair<Link*,Link*>(from,link->Out())] = w;
		GetWeights(link->Out(),weights,lambda);
	}
}

double BranchProcess::WeightedDrawSubTree(double lambda, Link*& down, Link*& up)	{

	map<pair<Link*,Link*>,double> weights;
	for (const Link* link=GetRoot()->Next(); link!=GetRoot(); link=link->Next())	{
		GetWeights(link->Out(),weights,lambda);
	}

	double totweight = 0;
	for (map<pair<Link*,Link*>,double>::iterator i=weights.begin(); i!=weights.end(); i++)	{
		totweight += i->second;
	}
	double u = totweight * rnd::GetRandom().Uniform();
	map<pair<Link*,Link*>, double>::iterator i = weights.begin();
	double cumul = i->second;
	while ((i!=weights.end()) && (cumul < u))	{
		i++;
		cumul += i->second;
	}
	if (i == weights.end())	{
		cerr << "error in BranchProcess::WeightedDrawSubTree: overflow\n";
		exit(1);
	}
	
	down = i->first.second;
	up = i->first.first;

	if ((! down) || (down->isRoot()))	{
		cerr << "error in tree::drawsubtree: down\n";
		cerr << down << '\n';
		exit(1);
	}
	if ((! up) || (up->isLeaf()))	{
		cerr << "error in tree::drawsubtree: up\n";
		cerr << up << '\n';
		exit(1);
	}
	return i->second/totweight;
}

double BranchProcess::GetSubTreeWeight(double lambda, Link* down, Link* up)	{

	lambda = 0;
	map<pair<Link*,Link*>,double> weights;
	for (const Link* link=GetRoot()->Next(); link!=GetRoot(); link=link->Next())	{
		GetWeights(link->Out(),weights,lambda);
	}

	double totweight = 0;
	for (map<pair<Link*,Link*>,double>::iterator i=weights.begin(); i!=weights.end(); i++)	{
		totweight += i->second;
	}
	return weights[pair<Link*,Link*>(up,down)] / totweight;
}

