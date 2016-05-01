
#include "StringStreamUtils.h"

#include "Random.h"
#include "PhyloProcess.h"
#include <string>

#include "Parallel.h"
extern MPI_Datatype Propagate_arg;

#include "TexTab.h"

double PhyloProcess::GibbsSPR(int nrep)	{

	double naccepted = 0;

	if (GetNprocs() > 1)	{
		for (int rep=0; rep<nrep; rep++)	{
			naccepted += MPIGibbsSPR();
		}
	}
	else	{
		for (int rep=0; rep<nrep; rep++)	{
			naccepted += NonMPIGibbsSPR();
		}
	}

	return naccepted / nrep;
}

double PhyloProcess::GibbsMHSPR(double lambda, int nrep)	{

	double naccepted = 0;

	if (GetNprocs() > 1)	{
		for (int rep=0; rep<nrep; rep++)	{
			naccepted += MPIGibbsMHSPR(lambda);
		}
	}
	else	{
		for (int rep=0; rep<nrep; rep++)	{
			naccepted += NonMPIGibbsMHSPR(lambda);
		}
	}

	return naccepted / nrep;
}

double PhyloProcess::TemperedGibbsSPR(double lambda, double mu, int nstep, int nrep)	{

	double naccepted = 0;

	if (GetNprocs() > 1)	{
		for (int rep=0; rep<nrep; rep++)	{
			naccepted += MPITemperedGibbsSPR(lambda,mu,nstep);
		}
	}
	else	{
		for (int rep=0; rep<nrep; rep++)	{
			naccepted += NonMPITemperedGibbsSPR(lambda,mu,nstep);
		}
	}

	return naccepted / nrep;
}

double PhyloProcess::BPPSPR(int nrep)	{

	double naccepted = 0;

	if (GetNprocs() > 1)	{
		for (int rep=0; rep<nrep; rep++)	{
			naccepted += MPIBPPSPR();
		}
	}
	else	{
		for (int rep=0; rep<nrep; rep++)	{
			naccepted += NonMPIBPPSPR();
		}
	}

	return naccepted / nrep;
}

double PhyloProcess::TemperedBPPSPR(int nrep, int nstep)	{

	double naccepted = 0;

	if (GetNprocs() > 1)	{
		for (int rep=0; rep<nrep; rep++)	{
			naccepted += MPITemperedBPPSPR(nstep);
		}
	}
	else	{
		for (int rep=0; rep<nrep; rep++)	{
			naccepted += NonMPITemperedBPPSPR(nstep);
		}
	}

	return naccepted / nrep;
}

double PhyloProcess::SpecialSPR(int nrep)	{

	double naccepted = 0;

	if (GetNprocs() > 1)	{
		for (int rep=0; rep<nrep; rep++)	{
			naccepted += MPISpecialSPR();
		}
	}
	else	{
		for (int rep=0; rep<nrep; rep++)	{
			naccepted += NonMPISpecialSPR();
		}
	}

	return naccepted / nrep;
}

double PhyloProcess::TemperedSpecialSPR(int nstep, int nrep)	{

	double naccepted = 0;

	if (GetNprocs() > 1)	{
		for (int rep=0; rep<nrep; rep++)	{
			naccepted += MPITemperedSpecialSPR(nstep);
		}
	}
	else	{
		for (int rep=0; rep<nrep; rep++)	{
			naccepted += NonMPITemperedSpecialSPR(nstep);
		}
	}

	return naccepted / nrep;
}


// MPI and Non MPI versions

int PhyloProcess::NonMPITemperedGibbsSPR(double lambda, double mu, int nstep)	{

	cerr << "in PhyloProcess::NonMPITemperedGibbsSPR\n";
	exit(1);
	return 0;
}

int PhyloProcess::MPITemperedGibbsSPR(double lambda, double mu, int nstep)	{

	Link* up = 0;
	Link* down = 0;

	ofstream tos((name + ".topo").c_str(),ios_base::app);
	tos << "tspr\t";
	GetTree()->ToStreamStandardForm(tos);
	tos << '\t';

	GlobalRootAtRandom();

	Backup();

	// choose subtree to be pruned and regrafted
	// returns probability of choosing that subtree
	double q1 = 1.0;
	if (lambda)	{
		if (BPP)	{
			q1 = BPP->ProposeFragileSPR(GetTree(),down,up);
		}
		else	{
			cerr << "weighted draw subtree\n";
			exit(1);
			q1 = WeightedDrawSubTree(lambda,down,up);
		}
	}
	else	{
		GetTree()->DrawSubTree(down,up);
	}

	if (down->isRoot())	{
		cerr << "down is root\n";
		return 0;
	}
	if (up->isRoot())	{
		cerr << "up is root\n";
		return 0;
	}

	Link* fromdown = GlobalDetach(down,up);
	
	GlobalUpdateConditionalLikelihoods();
	
	GlobalGibbsSPRScan(down,up,loglarray);
	map<pair<Link*,Link*>, double> loglmap;
	int n = 0;
	RecursiveGibbsFillMap(GetRoot(),GetRoot(),loglmap,loglarray,n);

	// check for numerical problems
	for (map<pair<Link*,Link*>,double>::iterator i=loglmap.begin(); i!=loglmap.end(); i++)	{
		if (isnan(i->second))	{
			cerr << "nan log prob in gibbs\n";
			exit(1);
		}
		if (isinf(i->second))	{
			cerr << "inf log prob in gibbs\n";
			exit(1);
		}
	}

	// choose among all target positions different from current position
	double max1 = 0;
	int notfound = 1;
	Link* fromup = 0;
	double prologp1 = 0;
	for (map<pair<Link*,Link*>,double>::iterator i=loglmap.begin(); i!=loglmap.end(); i++)	{
		// grep values of current position
		if (i->first.first == fromdown)	{
			fromup = i->first.second;
			prologp1 = i->second;
		}
		// find max, excluding current position
		if ((i->first.first != fromdown) && ((notfound) || (max1 < i->second)))	{
			max1 = i->second;
			notfound = 0;
		}
	}

	// compute total prob mass
	double total1 = 0;
	int nchoice = 0;
	for (map<pair<Link*,Link*>,double>::iterator i=loglmap.begin(); i!=loglmap.end(); i++)	{
		if (i->first.first != fromdown)	{
			nchoice++;
			double tmp = exp(i->second - max1);
			total1 += tmp;
		}
		if (isinf(total1))	{
			cerr << "error in gibbs: inf\n";
			cerr << "total1\n";
			exit(1);
		}
		if (isnan(total1))	{
			cerr << "error in gibbs: nan\n";
		}
	}
	if (! nchoice)	{
		cerr << "error in gibbs: no possible choice\n";
		exit(1);
	}

	// randomly choose final position, duly excluding current position
	double u = total1 * rnd::GetRandom().Uniform();
	map<pair<Link*,Link*>, double>::iterator i = loglmap.begin();
	double cumul = 0;
	if (i->first.first != fromdown)	{
		cumul += exp(i->second-max1);
	}
	while ((i!=loglmap.end()) && (cumul < u))	{
		i++;
		if (i == loglmap.end())	{
			cerr << "error in gibbs spr: overflow\n";
			exit(1);
		}
		if (i->first.first != fromdown)	{
			cumul += exp(i->second-max1);
		}
	}
	
	// store values of target position
	Link* todown = i->first.first;
	Link* toup = i->first.second;
	double logp2 = i->second;

	double deltalogp = 0;

	double reldifffwd = 2 * fabs(logp2 - prologp1) / fabs(logp2 + prologp1);
	double ptempfwd = exp(-reldifffwd / mu);
	int maketempmove = mu ? (rnd::GetRandom().Uniform() < ptempfwd) : 0;

	tsprtot++;
	if (maketempmove)	{

		tsprtmp++;
		// do the tempered move between the two topologies
		GlobalAttach(down,up,fromdown,fromup);

		if (! sumratealloc)	{
			GlobalActivateSumOverRateAllocations();

			sumratealloc = 1;

			GlobalUpdateConditionalLikelihoods();
			deltalogp = GlobalTemperedTreeMoveLogProb(nstep,down,up,fromdown,fromup,todown,toup);

			GlobalInactivateSumOverRateAllocations();

			sumratealloc = 0;
		}
		else	{

			GlobalUpdateConditionalLikelihoods();
			deltalogp = GlobalTemperedTreeMoveLogProb(nstep,down,up,fromdown,fromup,todown,toup);

		}

		// reverse probability 
		GlobalDetach(down,up);
		GlobalUpdateConditionalLikelihoods();

		GlobalGibbsSPRScan(down,up,loglarray);
		loglmap.clear();
		n = 0;
		RecursiveGibbsFillMap(GetRoot(),GetRoot(),loglmap,loglarray,n);

		// check for numerical problems
		for (map<pair<Link*,Link*>,double>::iterator i=loglmap.begin(); i!=loglmap.end(); i++)	{
			if (isnan(i->second))	{
				cerr << "nan log prob in gibbs\n";
				exit(1);
			}
			if (isinf(i->second))	{
				cerr << "inf log prob in gibbs\n";
				exit(1);
			}
		}
	}

	else	{
		deltalogp = logp2 - prologp1;
	}

	// calculate probability of reverse move
	double max2 = 0;
	notfound = 1;
	double prologp2 = 0;
	for (map<pair<Link*,Link*>,double>::iterator i=loglmap.begin(); i!=loglmap.end(); i++)	{
		// grep values of current position
		if (i->first.first == todown)	{
			prologp2 = i->second;
		}
		if ((i->first.first != todown) && ((notfound) || (max2 < i->second)))	{
			max2 = i->second;
			notfound = 0;
		}
	}

	double logp1 = 0;
	double total2 = 0;
	for (map<pair<Link*,Link*>,double>::iterator i=loglmap.begin(); i!=loglmap.end(); i++)	{
		if (i->first.first != todown)	{
			double tmp = exp(i->second - max2);
			total2 += tmp;
		}
		if (i->first.first == fromdown)	{
			logp1 = i->second;
		}
		if (isinf(total2))	{
			cerr << "error in gibbs: inf\n";
		}
		if (isnan(total2))	{
			cerr << "error in gibbs: nan\n";
		}
	}

	// reverse probability of choosing subtree to be pruned and regrafted
	GlobalAttach(down,up,todown,toup);

	GetTree()->ToStreamStandardForm(tos);
	tos << '\t';

	double q2 = 1.0;
	if (lambda)	{
		if (BPP)	{
			q2 = BPP->GetFragileSPRProposalProb(GetTree(),todown,toup);
		}
		else	{
			q2 = GetSubTreeWeight(lambda,todown,toup);
		}
	}

	double reldiffrev = 2 * fabs(logp1 - prologp2) / fabs(logp1 + prologp2);
	double ptemprev = exp(-reldiffrev / mu);
	if (! maketempmove)	{
		if (fabs(ptemprev - ptempfwd) > 1e-6)	{
			cerr << "error: fwd rev move should have same probability of being tempered if tempered move was not chosen\n";
			exit(1);
		}
	}
	double loghtemp = 0;
	if (maketempmove)	{
		loghtemp = log(ptemprev) - log(ptempfwd);
	}

	// hastings log ratio
	double logpfwd = logp2 - max1 - log(total1);
	double logprev = logp1 - max2 - log(total2);
	double logh = logprev - logpfwd + log(q2/q1) + loghtemp;

	// MH log ratio
	double logratio = logh + deltalogp;

	int accepted = (log(rnd::GetRandom().Uniform()) < logratio);

	if (accepted)	{
		tos << "accept";
	}
	else	{
		tos << "reject";
	}

	if (maketempmove)	{
		tos << "\ttempered";
		tos << '\t' << logratio << '\t' << deltalogp << '\t' << logp2 - prologp1 << '\t' << logprev - logpfwd << '\t' << loghtemp;
		double pseudologratio = logprev - logpfwd + log(q2/q1) + logp2 - prologp1;
		tos << '\t' << pseudologratio;
		double l1 = (logratio < 0) ? logratio : 0;
		double l2 = (pseudologratio < 0) ? pseudologratio : 0;
		double p1 = exp(l1);
		double p2 = exp(l2);
		tsprtmpacc00 += (1-p1) * (1-p2);
		tsprtmpacc01 += (1-p1) * p2;
		tsprtmpacc10 += p1 * (1-p2);
		tsprtmpacc11 += p1 * p2;
		
		ofstream os((name + ".temperedmove").c_str(),ios_base::app);
		os << logratio << '\t' << logp2 - prologp1 << '\t' << prologp2 - prologp1 << '\t' << logp1 - prologp2 << '\t' << logp1 - logp2 << '\t' << accepted << '\n';
		os.close();
	}
	else	{
		tos << "\tdirect";
		tos << '\t' << logratio << '\t' << deltalogp << '\t' << logprev - logpfwd << '\t' << loghtemp;
	}
	tos << '\n';

	if (! accepted)	{
		GlobalDetach(down,up);
		GlobalAttach(down,up,fromdown,fromup);
		if (maketempmove)	{
			Restore();
			GlobalUpdateParameters();
		}
	}
	if (maketempmove &&  (! sumratealloc))	{
		GlobalActivateSumOverRateAllocations();
		sumratealloc = 1;
		GlobalUpdateConditionalLikelihoods();
		GlobalInactivateSumOverRateAllocations();
		sumratealloc = 0;
	}
	// GlobalUpdateConditionalLikelihoods();
	return accepted;
}

int PhyloProcess::NonMPIGibbsMHSPR(double lambda)	{
}

int PhyloProcess::MPIGibbsMHSPR(double lambda)	{

	Link* up = 0;
	Link* down = 0;

	ofstream tos((name + ".topo").c_str(),ios_base::app);
	tos << "mhspr\t";
	GetTree()->ToStreamStandardForm(tos);
	tos << '\t';

	GlobalRootAtRandom();
	double q1 = 1.0;
	if (lambda)	{
		if (BPP)	{
			q1 = BPP->ProposeFragileSPR(GetTree(),down,up);
		}
		else	{
			cerr << "weighted draw subtree\n";
			exit(1);
			q1 = WeightedDrawSubTree(lambda,down,up);
		}
	}
	else	{
		GetTree()->DrawSubTree(down,up);
	}

	if (down->isRoot())	{
		cerr << "down is root\n";
		return 0;
	}
	if (up->isRoot())	{
		cerr << "up is root\n";
		return 0;
	}
	int sizebefore = GetTree()->GetSize();
	int subtreesize = GetTree()->GetSize(down);

	Link* fromdown = GlobalDetach(down,up);
	
	int sizeafter = GetTree()->GetSize();

	if (sizebefore-subtreesize-sizeafter)	{
		cerr << "error in gibbs spr: non matching subtree sizes\n";
		exit(1);
	}

	GlobalUpdateConditionalLikelihoods();
	
	GlobalGibbsSPRScan(down,up,loglarray);
	map<pair<Link*,Link*>, double> loglmap;
	int n = 0;
	RecursiveGibbsFillMap(GetRoot(),GetRoot(),loglmap,loglarray,n);

	// check for numerical problems
	for (map<pair<Link*,Link*>,double>::iterator i=loglmap.begin(); i!=loglmap.end(); i++)	{
		if (isnan(i->second))	{
			cerr << "nan log prob in gibbs\n";
			exit(1);
		}
		if (isinf(i->second))	{
			cerr << "inf log prob in gibbs\n";
			exit(1);
		}
	}

	// choose among all target positions different from current position
	double max1 = 0;
	int notfound = 1;
	double logp1 = 0;
	Link* fromup = 0;
	for (map<pair<Link*,Link*>,double>::iterator i=loglmap.begin(); i!=loglmap.end(); i++)	{
		// grep values of current position
		if (i->first.first == fromdown)	{
			logp1 = i->second;
			fromup = i->first.second;
		}
		// find max, excluding current position
		if ((i->first.first != fromdown) && ((notfound) || (max1 < i->second)))	{
			max1 = i->second;
			notfound = 0;
		}
	}

	// compute total prob mass
	double total1 = 0;
	int nchoice = 0;
	for (map<pair<Link*,Link*>,double>::iterator i=loglmap.begin(); i!=loglmap.end(); i++)	{
		if (i->first.first != fromdown)	{
			nchoice++;
			double tmp = exp(i->second - max1);
			total1 += tmp;
		}
		if (isinf(total1))	{
			cerr << "error in gibbs: inf\n";
			cerr << "total1\n";
			exit(1);
		}
		if (isnan(total1))	{
			cerr << "error in gibbs: nan\n";
		}
	}
	if (! nchoice)	{
		cerr << "error in gibbs: no possible choice\n";
		exit(1);
	}

	// randomly choose final position, duly excluding current position
	double u = total1 * rnd::GetRandom().Uniform();
	map<pair<Link*,Link*>, double>::iterator i = loglmap.begin();
	double cumul = 0;
	if (i->first.first != fromdown)	{
		cumul += exp(i->second-max1);
	}
	while ((i!=loglmap.end()) && (cumul < u))	{
		i++;
		if (i == loglmap.end())	{
			cerr << "error in gibbs spr: overflow\n";
			exit(1);
		}
		if (i->first.first != fromdown)	{
			cumul += exp(i->second-max1);
		}
	}
	
	// store values of target position
	Link* todown = i->first.first;
	Link* toup = i->first.second;
	double logp2 = i->second;

	// calculate probability of reverse move
	double max2 = 0;
	notfound = 1;
	for (map<pair<Link*,Link*>,double>::iterator i=loglmap.begin(); i!=loglmap.end(); i++)	{
		if ((i->first.first != todown) && ((notfound) || (max2 < i->second)))	{
			max2 = i->second;
			notfound = 0;
		}
	}

	double total2 = 0;
	for (map<pair<Link*,Link*>,double>::iterator i=loglmap.begin(); i!=loglmap.end(); i++)	{
		if (i->first.first != todown)	{
			double tmp = exp(i->second - max2);
			total2 += tmp;
		}
		if (isinf(total2))	{
			cerr << "error in gibbs: inf\n";
		}
		if (isnan(total2))	{
			cerr << "error in gibbs: nan\n";
		}
	}

	GlobalAttach(down,up,todown,toup);
	GetTree()->ToStreamStandardForm(tos);
	tos << '\t';
	// reverse probability of choosing subtree to be pruned and regrafted
	double q2 = 1.0;
	if (lambda)	{
		if (BPP)	{
			q2 = BPP->GetFragileSPRProposalProb(GetTree(),todown,toup);
		}
		else	{
			q2 = GetSubTreeWeight(lambda,todown,toup);
		}
	}

	// hastings log ratio
	double logpfwd = logp2 - max1 - log(total1);
	double logprev = logp1 - max2 - log(total2);
	double logh = logprev - logpfwd + log(q2 / q1);

	// MH log ratio
	double logratio = logh + logp2 - logp1;

	int accepted = (log(rnd::GetRandom().Uniform()) < logratio);
	if (accepted)	{
		// GlobalAttach(down,up,todown,toup);
	}
	else	{
		GlobalDetach(down,up);
		GlobalAttach(down,up,fromdown,fromup);
	}
	if (accepted)	{
		tos << "accept\n";
	}
	else	{
		tos << "reject\n";
	}

	// GlobalUpdateConditionalLikelihoods();
	return accepted;
}

int PhyloProcess::MPIGibbsSPR()	{

	Link* up = 0;
	Link* down = 0;
	GlobalRootAtRandom();
	GetTree()->DrawSubTree(down,up);

	if (down->isRoot())	{
		cerr << "down is root\n";
		return 0;
	}
	if (up->isRoot())	{
		cerr << "up is root\n";
		return 0;
	}
	int sizebefore = GetTree()->GetSize();
	int subtreesize = GetTree()->GetSize(down);

	Link* fromdown = GlobalDetach(down,up);
	
	int sizeafter = GetTree()->GetSize();

	if (sizebefore-subtreesize-sizeafter)	{
		cerr << "error in gibbs spr: non matching subtree sizes\n";
		exit(1);
	}

	GlobalUpdateConditionalLikelihoods();
	
	GlobalGibbsSPRScan(down,up,loglarray);
	map<pair<Link*,Link*>, double> loglmap;
	int n = 0;
	RecursiveGibbsFillMap(GetRoot(),GetRoot(),loglmap,loglarray,n);

	double max = 0;
	int j = 0;
	for (map<pair<Link*,Link*>,double>::iterator i=loglmap.begin(); i!=loglmap.end(); i++)	{
		if (isnan(i->second))	{
			cerr << "nan log prob in gibbs\n";
			exit(1);
		}
		if (isinf(i->second))	{
			cerr << "inf log prob in gibbs\n";
			exit(1);
		}
		if ((i==loglmap.begin()) || (max < i->second))	{
			max = i->second;
		}
		j++;
	}
	if (j != 2*sizeafter-3)	{
		cerr << "error in gibbs: number of cases not matching 2p-3\n";
		cerr << j << '\t' << sizeafter << '\n';
		exit(1);
	}

	double total = 0;
	for (map<pair<Link*,Link*>,double>::iterator i=loglmap.begin(); i!=loglmap.end(); i++)	{
		total += exp(i->second - max);
		if (isinf(total))	{
			cerr << "error in gibbs: inf\n";
		}
		if (isnan(total))	{
			cerr << "error in gibbs: nan\n";
		}
	}

	double u = total * rnd::GetRandom().Uniform();
	map<pair<Link*,Link*>, double>::iterator i = loglmap.begin();
	double cumul = exp(i->second-max);
	while ((i!=loglmap.end()) && (cumul < u))	{
		i++;
		if (i == loglmap.end())	{
			cerr << "error in gibbs spr: overflow\n";
			exit(1);
		}
		cumul += exp(i->second-max);
	}
	
	int accepted = (fromdown != i->first.first);
	if (i->second - max < -40)	{
		cerr << "suspicious choice in gibbs : " << i->second - max << '\t' << exp(i->second - max) << '\n';
		exit(1);
	}
	GlobalAttach(down,up,i->first.first,i->first.second);
	// GlobalUpdateConditionalLikelihoods();

	return accepted;
}

double PhyloProcess::NonMPIBPPSPR()	{

	cerr << "in NonMPIBPPSPR\n";
	exit(1);
	return 0;

}

double PhyloProcess::MPIBPPSPR()	{

	GlobalRootAtRandom();
	// GlobalUpdateConditionalLikelihoods();

	Link* subtree = 0;
	Link* subtreeup = 0;
	Link* fromdown = 0;
	Link* fromup = 0;
	Link* down = 0;
	Link* up = 0;

	double logp1 = logL;

	double logh = BPP->ProposeFullGibbsSPR(tree,subtree,subtreeup,fromdown,fromup,down,up);

	GlobalDetach(subtree,subtreeup);
	GlobalAttach(subtree,subtreeup,down,up);

	GlobalUpdateConditionalLikelihoods();
	double logp2 = logL;

	double logratio = logp2 - logp1 + logh;
	
	int accepted = (log(rnd::GetRandom().Uniform()) < logratio);

	if (! accepted)	{
		GlobalDetach(subtree,subtreeup);
		GlobalAttach(subtree,subtreeup,fromdown,fromup);
		// GlobalUpdateConditionalLikelihoods();
	}
	return accepted;
}

double PhyloProcess::NonMPITemperedBPPSPR(int nstep)	{

	cerr << "in NonMPIBPPSPR\n";
	exit(1);
	return 0;

}

double PhyloProcess::MPITemperedBPPSPR(int nstep)	{

	GlobalRootAtRandom();
	Backup();

	Link* subtree = 0;
	Link* subtreeup = 0;
	Link* fromdown = 0;
	Link* fromup = 0;
	Link* todown = 0;
	Link* toup = 0;

	double logh = BPP->ProposeFullGibbsSPR(tree,subtree,subtreeup,fromdown,fromup,todown,toup);

	double deltalogprior = -LogBranchProcessPrior();

	GlobalUpdateConditionalLikelihoods();

	double deltalogp = GlobalTemperedTreeMoveLogProb(nstep,subtree,subtreeup,fromdown,fromup,todown,toup);

	deltalogprior += LogBranchProcessPrior();

	double logratio = logh + deltalogprior + deltalogp;
	
	int accepted = (log(rnd::GetRandom().Uniform()) < logratio);

	if (accepted)	{
		GlobalDetach(subtree,subtreeup);
		GlobalAttach(subtree,subtreeup,todown,toup);
	}
	else	{
		Restore();
		GlobalUpdateParameters();
	}

	// GlobalUpdateConditionalLikelihoods();

	return accepted;
}

double PhyloProcess::NonMPITemperedSpecialSPR(int nstep)	{

	cerr << "in PhyloProcess::NonMPITemperedSpecialSPR\n";
	exit(1);
	return 0;
}

double PhyloProcess::NonMPISpecialSPR()	{

	cerr << "in PhyloProcess::NonMPISpecialSPR\n";
	exit(1);
	return 0;
}

double PhyloProcess::MPISpecialSPR()	{

	// GlobalUpdateConditionalLikelihoods();
	// double deltalogp = -logL;
	double deltalogp = - GlobalGetFullLogLikelihood();

	string taxname = "Strongyloc";
	Link* down = GetTree()->GetLCA(taxname,taxname);
	Link* up = GetTree()->GetAncestor(down);
	
	string taxname2 = "Branchiost";
	Link* down2 = GetTree()->GetLCA(taxname2,taxname2);
	Link* up2 = GetTree()->GetAncestor(down);
	
	Link* fromdown = GlobalDetach(down,up);
	Link* fromup = GetTree()->GetAncestor(fromdown);

	Link* todown = 0;
	Link* toup = 0;
	if (down2 == fromdown)	{
		currenttopo = 0;
		todown = fromup;
		toup = GetTree()->GetAncestor(todown);
	}
	else	{
		currenttopo = 1;
		todown = down2;
		toup = GetTree()->GetAncestor(todown);
	}

	GlobalAttach(down,up,todown,toup);

	// GlobalUpdateConditionalLikelihoods();
	// deltalogp += logL;
	deltalogp += GlobalGetFullLogLikelihood();

	// MH log ratio
	double logratio = deltalogp;

	int accepted = (log(rnd::GetRandom().Uniform()) < logratio);

	if (! accepted)	{
		GlobalDetach(down,up);
		GlobalAttach(down,up,fromdown,fromup);
		Restore();
		GlobalGetFullLogLikelihood();
		// GlobalUpdateParameters();
	}
	GlobalUpdateConditionalLikelihoods();
	return ((double) accepted);
}

double PhyloProcess::MPITemperedSpecialSPR(int nstep)	{

	Backup();
	GlobalGetFullLogLikelihood();
	// GlobalUpdateConditionalLikelihoods();

	string taxname = "Strongyloc";
	Link* down = GetTree()->GetLCA(taxname,taxname);
	Link* up = GetTree()->GetAncestor(down);

	string taxname2 = "Branchiost";
	Link* down2 = GetTree()->GetLCA(taxname2,taxname2);
	Link* up2 = GetTree()->GetAncestor(down);
	
	Link* fromdown = GlobalDetach(down,up);
	Link* fromup = GetTree()->GetAncestor(fromdown);

	Link* todown = 0;
	Link* toup = 0;
	if (down2 == fromdown)	{
		currenttopo = 0;
		todown = fromup;
		toup = GetTree()->GetAncestor(todown);
	}
	else	{
		currenttopo = 1;
		todown = down2;
		toup = GetTree()->GetAncestor(todown);
	}

	GlobalAttach(down,up,fromdown,fromup);

	double deltalogp = GlobalTemperedTreeMoveLogProb(nstep,down,up,fromdown,fromup,todown,toup);

	int accepted = (log(rnd::GetRandom().Uniform()) < deltalogp);

	if (accepted)	{
		GlobalDetach(down,up);
		GlobalAttach(down,up,todown,toup);
	}
	else	{
		Restore();
		GlobalUpdateParameters();
	}
	// GlobalUpdateConditionalLikelihoods();
	GlobalGetFullLogLikelihood();

	return ((double) accepted);
}

// Helper functions

void PhyloProcess::RecursiveGibbsSPRScan(Link* from, Link* fromup, Link* down, Link* up, double* loglarray, int& n)	{

	if (! from->isRoot())	{
		GetTree()->Attach(down,up,from,fromup);
		double*** aux = condlmap[0];
		Reset(aux,false);
		for (const Link* link=up->Next(); link!=up; link=link->Next())	{
			if (link->isRoot())	{
				cerr << "ROOT\n";
				exit(1);
			}
			Multiply(GetConditionalLikelihoodVector(link),aux,false);
		}
		Propagate(aux,GetConditionalLikelihoodVector(up->Out()),GetLength(up->GetBranch()),false);
		double logl = ComputeNodeLikelihood(up->Out(),0);
		if (n >= GetNbranch())	{
			cerr << "branch overflow\n";
			exit(1);
		}
		loglarray[n] = logl;
		n++;
		Link* tmp1 = GetTree()->Detach(down,up);
	}
	Link* trailer = from;
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		RecursiveGibbsSPRScan(link->Out(),trailer,down,up,loglarray,n);
		trailer = trailer->Next();
	}
}

void PhyloProcess::RecursiveGibbsFillMap(Link* from, Link* fromup, map<pair<Link*,Link*>,double>& loglmap, double* loglarray, int& n)	{

	if (! from->isRoot())	{
		if (n >= GetNbranch())	{
			cerr << "branch overflow\n";
			exit(1);
		}
		loglmap[pair<Link*,Link*>(from,fromup)] = loglarray[n];
		n++;
	}
	Link* trailer = from;
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		RecursiveGibbsFillMap(link->Out(),trailer,loglmap,loglarray,n);
		trailer = trailer->Next();
	}
}

double PhyloProcess::NonMPIGibbsSPR()	{

	GetTree()->RootAtRandom();

	Link* up = 0;
	Link* down = 0;

	GetTree()->DrawSubTree(down,up);
	
	int sizebefore = GetTree()->GetSize();
	int subtreesize = GetTree()->GetSize(down);
	Link* fromdown = GetTree()->Detach(down,up);
	
	int sizeafter = GetTree()->GetSize();

	if (sizebefore-subtreesize-sizeafter)	{
		cerr << "error in gibbs spr: non matching subtree sizes\n";
		exit(1);
	}

	UpdateConditionalLikelihoods();
	
	map<pair<Link*,Link*>, double> loglmap;
	RecursiveNonMPIGibbsSPRScan(GetRoot(),GetRoot(),down,up,loglmap);

	double max=0;
	int j = 0;
	for (map<pair<Link*,Link*>,double>::iterator i=loglmap.begin(); i!=loglmap.end(); i++)	{
		if ((i==loglmap.begin()) || (max < i->second))	{
			max = i->second;
		}
		j++;
	}
	if (j != 2*sizeafter-3)	{
		cerr << "error in gibbs: number of cases not matching 2p-3\n";
		cerr << j << '\t' << sizeafter << '\n';
		exit(1);
	}

	double total = 0;
	for (map<pair<Link*,Link*>,double>::iterator i=loglmap.begin(); i!=loglmap.end(); i++)	{
		total += exp(i->second - max);
		if (isinf(total))	{
			cerr << "error in gibbs: inf\n";
		}
		if (isnan(total))	{
			cerr << "error in gibbs: nan\n";
		}
	}
	double u = total * rnd::GetRandom().Uniform();
	map<pair<Link*,Link*>, double>::iterator i = loglmap.begin();
	double cumul = exp(i->second-max);
	while ((i!=loglmap.end()) && (cumul < u))	{
		i++;
		if (i == loglmap.end())	{
			cerr << "error in gibbs spr: overflow\n";
			exit(1);
		}
		cumul += exp(i->second-max);
	}
	
	int accepted = (fromdown != i->first.first);
	if (i->second - max < -40)	{
		cerr << "suspicious choice in gibbs : " << i->second - max << '\t' << exp(i->second - max) << '\n';
		exit(1);
	}
	GetTree()->Attach(down,up,i->first.first,i->first.second);
	UpdateConditionalLikelihoods();
	return accepted;
}

void PhyloProcess::RecursiveNonMPIGibbsSPRScan(Link* from, Link* fromup, Link* down, Link* up, map<pair<Link*,Link*>,double>& loglmap)	{

	if (! from->isRoot())	{
		GetTree()->Attach(down,up,from,fromup);
		double*** aux = condlmap[0];
		Reset(aux,false);
		for (const Link* link=up->Next(); link!=up; link=link->Next())	{
			if (link->isRoot())	{
				cerr << "ROOT\n";
				exit(1);
			}
			Multiply(GetConditionalLikelihoodVector(link),aux,false);
		}
		Propagate(aux,GetConditionalLikelihoodVector(up->Out()),GetLength(up->GetBranch()),false);
		double logl = ComputeNodeLikelihood(up->Out(),0);
		loglmap[pair<Link*,Link*>(from,fromup)] = logl;
		Link* tmp1 = GetTree()->Detach(down,up);
	}
	Link* trailer = from;
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		RecursiveNonMPIGibbsSPRScan(link->Out(),trailer,down,up,loglmap);
		trailer = trailer->Next();
	}
}

void PhyloProcess::GlobalGibbsSPRScan(Link* down, Link* up, double* loglarray)  {

	int args[2];
	MPI_Status stat;
	MESSAGE signal = SCAN;
	args[0] = GetLinkIndex(down);
	args[1] = GetLinkIndex(up);

	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(args,2,MPI_INT,0,MPI_COMM_WORLD);

	double dvector[GetNbranch()];
	for(int i=0; i<GetNbranch(); i++) {
		loglarray[i] = 0.0;
	}
	for(int i=1; i<GetNprocs(); i++) {
		MPI_Recv(dvector,GetNbranch(),MPI_DOUBLE,MPI_ANY_SOURCE,TAG1,MPI_COMM_WORLD,&stat);
		for(int j=0; j<GetNbranch(); j++) {
			loglarray[j] += dvector[j];
		}
	}
}

void PhyloProcess::SlaveGibbsSPRScan(int idown, int iup)	{

	LocalGibbsSPRScan(idown,iup);
	MPI_Send(loglarray,GetNbranch(),MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

void PhyloProcess::LocalGibbsSPRScan(int idown, int iup)	{

	int n = 0;
	Link* down = GetLink(idown);
	Link* up = GetLink(iup);
	RecursiveGibbsSPRScan(GetRoot(),GetRoot(),down,up,loglarray,n);
}
