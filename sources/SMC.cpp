
#include "PhyloProcess.h"

void PhyloProcess::GlobalIncrementNsite(int indelta)	{

	MESSAGE signal = SETNSITE;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&indelta,1,MPI_INT,0,MPI_COMM_WORLD);
	IncrementNsite(indelta);
}

void PhyloProcess::GlobalResetNsite()	{

	MESSAGE signal = RESETNSITE;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	ResetNsite();
}

void PhyloProcess::NonMPISMCBurnin(string name, int deltansite, int shortcycle, int longcycle, int cutoffsize, int nrep)	{

	cerr << "in PhyloProcess::NonMPISMCBurnin\n";
	exit(1);
}

void PhyloProcess::MPISMCBurnin(string name, int deltansite, int shortcycle, int longcycle, int cutoffsize, int nrep)	{

	ofstream los((name + ".logweights").c_str());
	ofstream cos((name + ".burnintimes").c_str());
	ofstream os((name + ".burnintrace").c_str());
	TraceHeader(os);
	os.flush();

	ofstream tos((name + ".burnintreelist").c_str());

	GlobalCollapse();

	Sample();
	// PriorSample();
	GlobalUnfold();
	GlobalCollapse();

	GlobalResetNsite();
	GlobalUnfold();

	int i = 0;
	int cycle = shortcycle;

	double logweight = 0;

	int delta = deltansite * shortcycle;

	while (i < GetNsite())	{

		Chrono ch;
		ch.Start();

		GlobalCollapse();

		int bki = i;
		i += delta;
		if (i > GetNsite())	{
			i = GetNsite();
		}
			

		GlobalIncrementNsite(i - bki);
		double deltalogp = GlobalSMCAddSites();

		logweight += deltalogp;
		los << i << '\t' << deltalogp << '\t' << logweight << '\n';
		// cerr << i << '\t' << delta << '\t' << logweight << '\n';
		los.flush();

		GlobalRestrictedMoveCycle();
		GlobalUnfold();
		TopoMoveCycle(nrep,1.0);

		Trace(os);
		os.flush();

		SetNamesFromLengths();
		RenormalizeBranchLengths();
		GetTree()->ToStream(tos);
		DenormalizeBranchLengths();
		tos.flush();

		if (cutoffsize && (i > cutoffsize))	{
			delta = deltansite * longcycle;
			cutoffsize = 0;
		}

		ch.Stop();
		cos << i << '\t' << ch.GetTime() / 1000 << '\n';
		cos.flush();
	}
}

/*
void PhyloProcess::MPISMCBurnin(string name, int deltansite, int shortcycle, int longcycle, int cutoffsize, int nrep)	{

	ofstream los((name + ".logweights").c_str());
	ofstream cos((name + ".burnintimes").c_str());
	ofstream os((name + ".burnintrace").c_str());
	TraceHeader(os);
	os.flush();

	ofstream tos((name + ".burnintreelist").c_str());

	GlobalCollapse();

	Sample();
	// PriorSample();
	GlobalUnfold();
	GlobalCollapse();

	GlobalResetNsite();
	GlobalUnfold();

	int i = 0;
	int cycle = shortcycle;

	double logweight = 0;

	while (i < finalnsite)	{

		Chrono ch;
		ch.Start();

		GlobalCollapse();

		for (int k=0; k<cycle; k++)	{

			if (i < finalnsite)	{

				int bki = i;
				i += deltansite;
				if (i > finalnsite)	{
					i = finalnsite;
				}
					

				GlobalIncrementNsite(i - bki);
				double delta = GlobalSMCAddSites();
				logweight += delta;
				los << i << '\t' << delta << '\t' << logweight << '\n';
				// cerr << i << '\t' << delta << '\t' << logweight << '\n';
				los.flush();
				GlobalRestrictedMoveCycle();


				Trace(os);
				os.flush();

			}
		}

		GlobalUnfold();

		TopoMoveCycle(nrep);

		SetNamesFromLengths();
		RenormalizeBranchLengths();
		GetTree()->ToStream(tos);
		DenormalizeBranchLengths();
		tos.flush();

		if (cutoffsize && (i > cutoffsize))	{
			cycle = longcycle;
			cutoffsize = 0;
		}

		ch.Stop();
		cos << i << '\t' << ch.GetTime() / 1000 << '\n';
		cos.flush();
	}
}
*/

void PhyloProcess::CreateSiteConditionalLikelihoods()	{

	if (! sitecondlmap)	{
		sitecondlmap = new double**[GetNlink()];
		for (int j=0; j<GetNlink(); j++)	{
			sitecondlmap[j] = new double*[GetNrate()];
			for (int k=0; k<GetNrate(); k++)	{
				sitecondlmap[j][k] = new double[GetNstate()+1];
			}
		}
	}
}

void PhyloProcess::DeleteSiteConditionalLikelihoods()	{

	if (sitecondlmap)	{
		for (int j=0; j<GetNlink(); j++)	{
			for (int k=0; k<GetNrate(); k++)	{
				delete[] sitecondlmap[j][k];
			}
			delete[] sitecondlmap[j];
		}
		delete[] sitecondlmap;
		sitecondlmap = 0;
	}
}

double PhyloProcess::SiteLogLikelihood(int site)	{

	PrepareSiteLogLikelihood(site);
	SiteActivateSumOverRateAllocation(site);
	SitePostOrderPruning(site,GetRoot());
	SiteMultiplyByStationaries(site,sitecondlmap[0]);
	SiteComputeLikelihood(site,sitecondlmap[0]);
	return sitelogL[site];
}

void PhyloProcess::SitePostOrderPruning(int site, const Link* from)	{

	if (from->isLeaf())	{
		SiteInitialize(site,sitecondlmap[0],GetData(from)[site]);
	}
	else	{
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			SitePostOrderPruning(site,link->Out());
			SitePropagate(site,sitecondlmap[0],sitecondlmap[GetLinkIndex(link)],GetLength(link->GetBranch()));
		}
		SiteReset(site,sitecondlmap[0]);
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			SiteMultiply(site,sitecondlmap[GetLinkIndex(link)],sitecondlmap[0]);
		}
		SiteOffset(site,sitecondlmap[0]);
	}
	if (from->isRoot())	{
	}	
}

void PhyloProcess::SampleSiteMapping(int site)	{

	SiteDrawAllocations(site,0);
	SampleSiteNodeStates(site);
	SiteInactivateSumOverRateAllocation(site,ratealloc[site]);
	SampleSiteSubstitutionMapping(site,GetRoot());
}


void PhyloProcess::SampleSiteNodeStates(int site)	{
	SampleSiteNodeStates(site, GetRoot(),sitecondlmap[0]);
}


void PhyloProcess::SampleSiteNodeStates(int site, const Link* from, double** aux)	{
	
	if (from->isLeaf())	{
		SiteInitialize(site,aux,GetData(from)[site]);
	}
	else	{
		SiteReset(site,aux,true);
	}
	// make product of conditional likelihoods around node
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		SiteMultiply(site,sitecondlmap[GetLinkIndex(link)],aux,true);
	}
	if (!from->isRoot())	{
		SiteMultiply(site,sitecondlmap[GetLinkIndex(from)],aux,true);
	}
	SiteMultiplyByStationaries(site,aux,true);
	// let substitution process choose states based on this vector
	// this should collapse the vector into 1s and 0s
	GetStates(from->GetNode())[site] = SiteChooseState(site,aux);

	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		// propagate forward
		SitePropagate(site,aux,sitecondlmap[GetLinkIndex(link->Out())],GetLength(link->GetBranch()),true);
	}
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		SampleSiteNodeStates(site,link->Out(),aux);
	}
}

void PhyloProcess::GlobalActivateFastTopo()	{

	GlobalSetMinMax(0,fasttopofracmin);
	if (fastcondrate)	{
		GlobalInactivateSumOverRateAllocations();
	}
	GlobalUpdateConditionalLikelihoods();
}

void PhyloProcess::GlobalInactivateFastTopo()	{

	GlobalSetMinMax(0,1);
	if (fastcondrate)	{
		GlobalActivateSumOverRateAllocations();
	}
	GlobalUpdateConditionalLikelihoods();
}

double PhyloProcess::FastTopoMoveCycle(int nrep, double mu)	{

	fasttopotry++;

	GlobalReshuffleSites();

	GlobalBackupTree();

	GlobalUpdateConditionalLikelihoods();
	double logp1 = logL;

	GlobalActivateFastTopo();
	double subdeltalogp = -logL;
	// double deltalogprior = -LogLengthPrior();

	double preacc = 0;
	preacc = SPRMove(10);
	// preacc = NNIMove(5,0);
	/*
	if (rnd::GetRandom().Uniform() < 0.5)	{
		preacc = SPRMove(10);
	}
	else	{
		preacc = NNIMove(5,0);
	}
	*/
	// double preacc = SimpleTopoMoveCycle(nrep,tuning);

	GlobalUpdateConditionalLikelihoods();
	subdeltalogp += logL;
	// deltalogprior += LogLengthPrior();

	GlobalInactivateFastTopo();

	double logp2 = logL;

	if (fabs(subdeltalogp) < 1e-6)	{
		// topo has not changed
		// always accepted
		return 1;
	}

	int accepted = 0;

	double reldifffwd = 2 * fabs(logp2 - logp1) / fabs(logp2 + logp1);
	double ptempfwd = exp(-reldifffwd / mu);
	int maketempmove = (rnd::GetRandom().Uniform() < ptempfwd);

	double logratio = -subdeltalogp;

	if ((fasttoponstep > 1) && maketempmove)	{

		GlobalSwapTree();
		GlobalUpdateConditionalLikelihoods();
		double tempered =  GlobalTemperedTreeMoveLogProb(fasttoponstep);
		logratio += tempered;

		// reverse move:
		// likelihood under initial tree
		double logq1 = logL;
		// likelihood under final tree
		GlobalSwapTree();
		GlobalUpdateConditionalLikelihoods();
		double logq2 = logL;
		// GlobalSwapTree();
		
		// prob of reverse move
		double reldiffrev = 2 * fabs(logq2 - logq1) / fabs(logq2 + logq1);
		double ptemprev = exp(-reldiffrev / mu);
		double loghtemp = log(ptemprev) - log(ptempfwd);

		logratio += loghtemp;
		// logratio -= deltalogprior;

		accepted = (log(rnd::GetRandom().Uniform()) < logratio);

		ofstream os((name + ".temperedmove").c_str(),ios_base::app);
		os << logratio << '\t' << tempered << '\t' << logp2 - logp1 << '\t' << subdeltalogp << '\t' << loghtemp << '\t' << logq2 - logp1 << '\t' << logq2 - logq1 << '\t' << logp2 - logq1 << '\t' << accepted << '\n';
		os.close();
	}
	else	{
		logratio += logp2 - logp1;
		accepted = (log(rnd::GetRandom().Uniform()) < logratio);
	}

	if (! accepted)	{
		GlobalRestoreTree();
		GlobalUpdateConditionalLikelihoods();
	}

	fasttopochange++;
	if (accepted)	{
		fasttopoacc++;
	}

	if (fasttoponstep > 1)	{
		anntot++;
		if (maketempmove)	{
			anntmp++;
			// would have been accepted without tempering ?
			double pseudologratio = logp2 - logp1 - subdeltalogp;
			double l1 = (logratio < 0) ? logratio : 0;
			double l2 = (pseudologratio < 0) ? pseudologratio : 0;
			double p1 = exp(l1);
			double p2 = exp(l2);
			anntmpacc00 += (1-p1) * (1-p2);
			anntmpacc01 += (1-p1) * p2;
			anntmpacc10 += p1 * (1-p2);
			anntmpacc11 += p1 * p2;
		}
	}

	return ((double) accepted);
}

/*
double PhyloProcess::AnnealedTopoMoveCycle(int nrep, double tuning, double frac0, int nstep)	{

	int accepted = 0;

	GlobalReshuffleSites();

	GlobalBackupTree();

	double deltalogp = 0;

	if (nstep == 1)	{
		GlobalUpdateConditionalLikelihoods();
		deltalogp -= logL;
	}

	GlobalEnterFastLikelihood();
	GlobalSetMinMax(0,frac0);
	GlobalUpdateConditionalLikelihoods();
	double subdeltalogp = -logL;
	double deltalogprior = -LogLengthPrior();

	double preacc = 0;
	// preacc = SPRMove(10);
	// preacc = NNIMove(5,0);
	if (rnd::GetRandom().Uniform() < 0.5)	{
		preacc = SPRMove(10);
	}
	else	{
		preacc = NNIMove(5,0);
	}
	// double preacc = SimpleTopoMoveCycle(nrep,tuning);

	// for (int rep=0; rep<nrep; rep++)	{
	//	GlobalCollapse();
	//	GlobalRestrictedTemperedMove();
	//	GlobalUnfold();
	//	double preacc = SimpleTopoMoveCycle(1,tuning);
	//}
	//GlobalCollapse();
	//GlobalRestrictedTemperedMove();
	//GlobalUnfold();

	GlobalUpdateConditionalLikelihoods();
	subdeltalogp += logL;
	deltalogprior += LogLengthPrior();
	GlobalLeaveFastLikelihood();

	GlobalSetMinMax(0,1);

	if (nstep == 1)	{
		GlobalUpdateConditionalLikelihoods();
		deltalogp += logL;
	}
	else	{
		GlobalSwapTree();
		GlobalUpdateConditionalLikelihoods();
		deltalogp = GlobalTemperedTreeMoveLogProb(nstep);
		GlobalSwapTree();
		deltalogp -= deltalogprior;
	}

	double logratio = deltalogp - subdeltalogp;

	accepted = (log(rnd::GetRandom().Uniform()) < logratio);
	if (! accepted)	{
		GlobalRestoreTree();
	}
	GlobalUpdateConditionalLikelihoods();

	if (accepted)	{
		annealtopoacc++;
	}
	annealtopotry++;

	return ((double) accepted);
}
*/

/*
double PhyloProcess::AnnealedTopoMoveCycle(int nrep, double tuning, double frac0, int nstep)	{

	int accepted = 0;

	GlobalReshuffleSites();

	if (nstep == 1)	{

		GlobalBackupTree();
		GlobalUpdateConditionalLikelihoods();
		double logp1 = logL;

		GlobalSetMinMax(0,frac0);
		GlobalUpdateConditionalLikelihoods();
		double sublogp1 = logL;

		double preacc = SimpleTopoMoveCycle(nrep,tuning);
		// for (int rep=0; rep<nrep; rep++)	{
		//	GlobalCollapse();
		//	GlobalRestrictedTemperedMove();
		//	GlobalUnfold();
		//	double preacc = SimpleTopoMoveCycle(1,tuning);
		//}
		//GlobalCollapse();
		//GlobalRestrictedTemperedMove();
		//GlobalUnfold();

		GlobalUpdateConditionalLikelihoods();
		double sublogp2 = logL;
		
		GlobalSetMinMax(0,1);
		GlobalUpdateConditionalLikelihoods();
		double logp2 = logL;

		double logratio = (logp2 - logp1) - (sublogp2 - sublogp1);

		accepted = (log(rnd::GetRandom().Uniform()) < logratio);
		if (! accepted)	{
			GlobalRestoreTree();
			GlobalUpdateConditionalLikelihoods();
		}
	}
	else	{

		double stepsize = (1.0 - frac0) / nstep;

		double deltalogp = 0;

		// back up topology and branchlengths
		GlobalBackupTree();

		// reshuffle sites

		// mappings for all sites
		GlobalUpdateConditionalLikelihoods();
		double logp1 = logL;
		GlobalCollapse();

		// reduce dataset down to fracmin
		// along the process, add log likelihoods for each fraction, as they are being removed
		for (int step=nstep-1; step>=0; step--)	{

			double fracmin = frac0 + stepsize*step;
			double fracmax = fracmin + stepsize;
			GlobalSetMinMax(fracmin,fracmax);

			GlobalUnfold();
			deltalogp -= logL;

			GlobalSetMinMax(0,0);
			GlobalCollapse();
			GlobalSetMinMax(0,fracmin);

			if (step > 0)	{
				GlobalRestrictedTemperedMove();
			}
		}

		GlobalUnfold();
		double preacc = SimpleTopoMoveCycle(nrep,tuning);
		GlobalCollapse();

		// expand dataset up to 1
		// along the process, add log likelihoods for each fraction, as they are being added back
		for (int step=0; step<nstep; step++)	{

			if (step > 0)	{
				GlobalRestrictedTemperedMove();
			}

			GlobalSetMinMax(0,0);
			GlobalUnfold();

			double fracmin = frac0 + stepsize*step;
			double fracmax = fracmin + stepsize;

			GlobalSetMinMax(fracmin,fracmax);
			// GlobalUnfold();
			GlobalUpdateConditionalLikelihoods();
			deltalogp += logL;

			GlobalCollapse();
			
			GlobalSetMinMax(0,fracmax);

		}

		accepted = (log(rnd::GetRandom().Uniform()) < deltalogp);
		if (! accepted)	{
			// restore topology and branch lengths
			GlobalRestoreTree();
			// GlobalUpdateConditionalLikelihoods();
		}

		GlobalUnfold();
		double logp2 = logL;
		if (! accepted)	{
			if (fabs(logp2 - logp1) > 1e-6)	{
				cerr << "error in annealed topo move: incorrect restore\n";
				exit(1);
			}
		}
	}

	if (accepted)	{
		annealtopoacc++;
	}
	annealtopotry++;

	return ((double) accepted);
}
*/

double PhyloProcess::GlobalTemperedTreeMoveLogProb(int nstep, Link* down, Link* up, Link* fromdown, Link* fromup, Link* todown, Link* toup)	{

	// assumes the two trees have already been set up

	double deltalogp = 0;

	// under old topology (all sites)
	GlobalCollapse();

	for (int step=0; step<nstep; step++)	{

		double fracmin = ((double) step) / nstep;
		double fracmax = ((double) step+1) / nstep;

		GlobalSetMinMax(fracmin,fracmax);

		// under old topology (current fraction of sites)
		GlobalUnfold();

		if (sumovercomponents)	{
			deltalogp -= GlobalGetFullLogLikelihood();
		}
		else	{
			deltalogp -= logL;
		}

		// switch to new topology
		GlobalDetach(down,up);
		GlobalAttach(down,up,todown,toup);

		GlobalUpdateConditionalLikelihoods();

		if (sumovercomponents)	{
			deltalogp += GlobalGetFullLogLikelihood();
		}
		else	{
			deltalogp += logL;
		}

		// under new topology (current fraction of sites)
		GlobalCollapse();

		// switch back to old topology
		GlobalDetach(down,up);
		GlobalAttach(down,up,fromdown,fromup);

		GlobalSetMinMax(0,1);

		// for all sites
		if (step < nstep-1)	{
			GlobalRestrictedTemperedMove();
		}
	}

	// switch to new topology
	GlobalDetach(down,up);
	GlobalAttach(down,up,todown,toup);

	// under new topology (all sites)
	GlobalUnfold();
	return deltalogp;
}

double PhyloProcess::GlobalTemperedTreeMoveLogProb(int nstep)	{

	// assumes the two trees have already been set up

	// backup

	double deltalogp = 0;

	// under old topology (all sites)
	GlobalCollapse();

	for (int step=0; step<nstep; step++)	{

		double fracmin = ((double) step) / nstep;
		double fracmax = ((double) step+1) / nstep;

		GlobalSetMinMax(fracmin,fracmax);

		// under old topology (current fraction of sites)
		GlobalUnfold();

		if (sumovercomponents)	{
			deltalogp -= GlobalGetFullLogLikelihood();
		}
		else	{
			deltalogp -= logL;
		}

		GlobalSwapTree();

		GlobalUpdateConditionalLikelihoods();

		if (sumovercomponents)	{
			deltalogp += GlobalGetFullLogLikelihood();
		}
		else	{
			deltalogp += logL;
		}

		// under new topology (current fraction of sites)
		GlobalCollapse();

		GlobalSwapTree();

		GlobalSetMinMax(0,1);

		// for all sites
		if (step < nstep-1)	{
			GlobalRestrictedTemperedMove();
		}
	}

	GlobalSwapTree();

	// under new topology (all sites)
	GlobalUnfold();
	
	return deltalogp;
}

double PhyloProcess::GlobalRestrictedTemperedMove()	{

	BranchProcessMove(1.0,0);
	GlobalUpdateParameters();
}

void PhyloProcess::GlobalSetMinMax(double min, double max)	{

	if (GetNprocs() > 1)	{
		MESSAGE signal = MINMAX;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
		double minmax[2];
		minmax[0] = min;
		minmax[1] = max;
		MPI_Bcast(&minmax,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
	}
	SetMinMax(min,max);
}

void PhyloProcess::SlaveSetMinMax()	{

	double minmax[2];
	MPI_Bcast(minmax,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
	SetMinMax(minmax[0],minmax[1]);
}


