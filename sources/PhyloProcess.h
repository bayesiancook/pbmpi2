
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef PHYLOPROCESS_H
#define PHYLOPROCESS_H

#include "SequenceAlignment.h"
#include "CodonSequenceAlignment.h"
#include "ZippedSequenceAlignment.h"

#include "SubstitutionProcess.h"
#include "BranchProcess.h"

#include "Parallel.h"

#include "CCP.h"

#include <map>
#include <vector>

class PhyloProcess : public virtual SubstitutionProcess, public virtual BranchProcess {

	public:

	virtual void WaitLoop();

	virtual void SlaveExecute(MESSAGE);

        virtual void SlaveRoot(int);

	virtual void SlaveGibbsSPRScan(int,int);
	void LocalGibbsSPRScan(int,int);

	virtual void SlaveProposeMove(int,double);
	virtual void SlaveRestore(int);
	virtual void SlaveReset(int,bool);
	virtual void SlaveMultiplyByStationaries(int,bool);
	virtual void SlaveMultiply(int,int,bool);
	virtual void SlaveInitialize(int,int,bool);
	virtual void SlavePropagate(int,int,bool,double);

	// virtual void SlaveUpdate();

	// default constructor: pointers set to nil
	PhyloProcess() :  sitecondlmap(0), siteratesuffstatcount(0), siteratesuffstatbeta(0), branchlengthsuffstatcount(0), branchlengthsuffstatbeta(0), size(0), totaltime(0), currenttopo(0), sumovercomponents(0), data(0), fasttopo(0) {
		empfreq = 0;
		spracc = sprtry = 0;
		mhspracc = mhsprtry = 0;
		tspracc = tsprtry = tsprtmp = tsprtot = 0;
		tsprtmpacc00 = 0;
		tsprtmpacc01 = 0;
		tsprtmpacc10 = 0;
		tsprtmpacc11 = 0;
		anntot = anntmp = 0;
		anntmpacc00 = 0;
		anntmpacc01 = 0;
		anntmpacc10 = 0;
		anntmpacc11 = 0;
		nniacc = nnitry = 0;
		bppspracc = bppsprtry = 0;
		tbppspracc = tbppsprtry = 0;
		specacc = spectry = 0;
		tspecacc = tspectry = 0;
		profacc = proftry = 0;
		rracc = rrtry = 0;
		fasttopoacc = fasttopotry = fasttopochange = 0;
		ziptopoacc = ziptopotry = 0;
		nprelim = 0;
	}

	virtual ~PhyloProcess() {}

	// performs one full cycle of MCMC
	// returns average success rate
	virtual double Move(double tuning = 1.0) = 0;
	virtual double AugmentedMove(double tuning = 1.0) {}

	virtual double RestrictedMoveCycle(int nrep = 1, double tuning = 1.0) {
		cerr << "in default restricted move cycle\n";
		exit(1);
	}

	virtual double GlobalRestrictedMoveCycle(int nrep = 1, double tuning = 1.0) {
		cerr << "in default global restricted move cycle\n";
		exit(1);
	}

	virtual double GlobalExtendedMoveCycle(int nrep = 1, double tuning = 1.0) {
		return TopoMoveCycle(nrep,tuning);
	}

	virtual double TopoMoveCycle(int nrep, double tuning)	{
		if (fasttopo)	{
			return FastTopoMoveCycle(nrep,tuning);
		}
		return SimpleTopoMoveCycle(nrep,tuning);
	}

	virtual double SimpleTopoMoveCycle(int nrep, double tuning);

	double FastTopoMoveCycle(int nrep, double tuning);

	double MoveTopo();
	double SPRMove(int nrep);
	double NNIMove(int nrep, double tuning);

	void GlobalBackupTree();
	void GlobalRestoreTree();
	void GlobalSwapTree();
	virtual void SlaveBackupTree();
	virtual void SlaveRestoreTree();
	virtual void SlaveSwapTree();

	virtual void GlobalActivateFastTopo();
	virtual void GlobalInactivateFastTopo();

	/*
	virtual void GlobalEnterFastLikelihood() {}
	virtual void SlaveEnterFastLikelihood() {}
	virtual void GlobalLeaveFastLikelihood() {}
	virtual void SlaveLeaveFastLikelihood() {}
	*/

	virtual void PrepareSiteLogLikelihood(int site) {
		cerr << "in default PrepareSiteLogLikelihood\n";
		exit(1);
	}

	// weird but simpler for multiple inheritance reasons
	virtual void GlobalActivateZip();
	virtual void GlobalInactivateZip();

	virtual void SlaveActivateZip();
	virtual void SlaveInactivateZip();

	virtual double SiteLogLikelihood(int site);
	virtual void SampleSiteMapping(int site);

	void SMCBurnin(string name, int deltansite, int shortcycle, int longcycle, int cutoffsize, int nrep)	{

		if (GetNprocs() > 1)	{
			MPISMCBurnin(name,deltansite,shortcycle,longcycle,cutoffsize,nrep);
		}
		else	{
			NonMPISMCBurnin(name,deltansite,shortcycle,longcycle,cutoffsize,nrep);
		}
	}

	void MPISMCBurnin(string name, int deltansite, int shortcycle, int longcycle, int cutoffsize, int nrep);
	void NonMPISMCBurnin(string name, int deltansite, int shortcycle, int longcycle, int cutoffsize, int nrep);

	/*
	double  SMC(string name, int shortcycle = 1, int mediumcycle = 10, int longcycle = 100, int mediumsize = 1000, int maxsize = 10000)	{

		double ret = 0;
		if (GetNprocs() > 1)	{
			ret = MPISMC(name,shortcycle,mediumcycle,longcycle,mediumsize,maxsize);
		}
		else	{
			ret = NonMPISMC(name,shortcycle,mediumcycle,longcycle,mediumsize,maxsize);
		}
		return ret;
	}

	double  NonMPISMC(string name, int shortcycle = 1, int mediumcycle = 10, int longcycle = 100, int mediumsize = 1000, int maxsize = 10000);
	double  MPISMC(string name, int shortcycle = 1, int mediumcycle = 10, int longcycle = 100, int mediumsize = 1000, int maxsize = 10000);
	*/

	void GlobalIncrementNsite(int innsite);
	void GlobalResetNsite();

	// sample from prior
	virtual void Sample()	{
		SampleRate();
		SampleLength();
		SampleProfile();
	}

	virtual void PriorSample()	{
		PriorSampleRate();
		PriorSampleLength();
		PriorSampleProfile();
	}

	double GlobalTemperedTreeMoveLogProb(int nstep, Link* from, Link* up, Link* fromdown, Link* fromup, Link* todown, Link* toup);
	double GlobalTemperedTreeMoveLogProb(int nstep);
	double GlobalRestrictedTemperedMove();	

	void GlobalSetMinMax(double min, double max);
	virtual void SlaveSetMinMax();

	// print out the first line (header) of the trace file
	virtual void TraceHeader(ostream& os) = 0;

	// print out one line of trace (summary statistics such as logprob, treelength, totaltime, etc)
	virtual void Trace(ostream& os) = 0;

	virtual void Monitor(ostream& os)  {
		os << "matrix uni" << '\t' << SubMatrix::GetUniSubCount() << '\n';
		os << "inf prob  " << '\t' << GetInfProbCount() << '\n';
		os << "stat inf  " << '\t' << GetStatInfCount() << '\n';
		if (sprtry)	{
			os << "spr " << '\t' << (100 * spracc) / sprtry << '\n';
		}
		if (mhsprtry)	{
			os << "mhspr " << '\t' << (100 * mhspracc) / mhsprtry << '\n';
		}
		if (tsprtry)	{
			os << "tspr " << '\t' << (100 * tspracc) / tsprtry << '\n';
			os << "fraction of temperedmoves: " << (100 * tsprtmp) / tsprtot << '\n';
			if (tsprtmp)	{
				os << "10" << '\t' << (100 * tsprtmpacc10) / tsprtmp << '\n';
				os << "01" << '\t' << (100 * tsprtmpacc01) / tsprtmp << '\n';
				os << "11" << '\t' << (100 * tsprtmpacc11) / tsprtmp << '\n';
				os << "00" << '\t' << (100 * tsprtmpacc00) / tsprtmp << '\n';
			}
		}
		if (nnitry)	{
			os << "nni " << '\t' << (100 * nniacc) / nnitry << '\n';
		}
		if (spectry)	{
			os << "special spr : " << '\t' << (100 * specacc) / spectry << '\n';
		}
		if (tspectry)	{
			os << "special tempered spr : " << '\t' << (100 * tspecacc) / tspectry << '\n';
		}
		if (bppsprtry)	{
			os << "bppspr " << '\t' << (100 * bppspracc) / bppsprtry << '\n';
		}
		if (tbppsprtry)	{
			os << "tbppspr " << '\t' << (100 * tbppspracc) / tbppsprtry << '\n';
		}
		if (proftry)	{
			os << "profile moves " << '\t' << (100 * profacc) / proftry << '\n';
		}
		if (rrtry)	{
			os << "rr moves " << '\t' << (100 * rracc) / rrtry << '\n';
		}
		if (ziptopotry)	{
			os << "zip topo : " << '\t' << (100 * ziptopoacc) / ziptopotry << '\n';
		}
		if (fasttopotry)	{
			os << "fast topo moves: \n";
			os << "topo changed : " << '\t' << (100 * fasttopochange) / fasttopotry << '\n';
			if (fasttopochange)	{
				os << "accepted     : " << '\t' << (100 * fasttopoacc) / fasttopochange << '\n';
			}
			if (anntot)	{
				os << "tempered fraction: " << '\t' << (100 * anntmp) / anntot << '\n';
				if (anntmp)	{
					os << "10" << '\t' << (100 * anntmpacc10) / anntmp << '\n';
					os << "01" << '\t' << (100 * anntmpacc01) / anntmp << '\n';
					os << "11" << '\t' << (100 * anntmpacc11) / anntmp << '\n';
					os << "00" << '\t' << (100 * anntmpacc00) / anntmp << '\n';
				}
			}
		}
		if (! fixtopo)	{
			double totaltime = nnichrono.GetTime() + sprchrono.GetTime() + tsprchrono.GetTime();
			os << "nni  time : " << nnichrono.GetTime() / totaltime << '\n';
			os << "spr  time : " << sprchrono.GetTime() / totaltime << '\n';
			os << "tspr time : " << tsprchrono.GetTime() / totaltime << '\n';
		}
	}

	void SetParameters(string indatafile, string intreefile, int iniscodon, GeneticCodeType incodetype, int infixtopo, int inNSPR, int inNMHSPR, int inNTSPR, double intopolambda, double intopomu, int intoponstep, int inNNNI, int innspec, int inntspec, int inbpp, int innbpp, int inntbpp, int inbppnstep, string inbppname, double inbppcutoff, double inbppbeta, int inprofilepriortype, int indc, int infixbl, int insumovercomponents, int inproposemode, int inallocmode, int insumratealloc,int infasttopo, double infasttopofracmin, int infasttoponstep, int infastcondrate)	{

		datafile = indatafile;
		treefile = intreefile;
		iscodon = iniscodon;
		codetype = incodetype;
		fixtopo = infixtopo;
		NSPR = inNSPR;
		NMHSPR = inNMHSPR;
		NTSPR = inNTSPR;
		topolambda = intopolambda;
		topomu = intopomu;
		toponstep = intoponstep;
		NNNI = inNNNI;
		nspec = innspec;
		ntspec = inntspec;
		bpp = inbpp;
		nbpp = innbpp;
		ntbpp = inntbpp;
		bppnstep = inbppnstep;
		bppname = inbppname;
		bppcutoff = inbppcutoff;
		bppbeta = inbppbeta;
		
		BPP = 0;
		if (bpp == 1)	{
			cerr << "make new BPP\n";
			BPP = new UnrootedBPP(bppname,bppcutoff,bppbeta);
		}
		else if (bpp == 2)	{
			cerr << "make new CCP\n";
			BPP = new UnrootedCCP(bppname,bppcutoff,bppbeta);
		}
		else if (bpp == 3)	{
			cerr << "make new CCP\n";
			BPP = new UnrootedCCP(bppname,bppcutoff,bppbeta);
		}

		profilepriortype = inprofilepriortype;
		dc = indc;
		fixbl = infixbl;
		sumovercomponents = insumovercomponents;
		proposemode = inproposemode;
		allocmode = inallocmode;
		sumratealloc = insumratealloc;
		fasttopo = infasttopo;
		fasttopofracmin = infasttopofracmin;
		fasttoponstep = infasttoponstep;
		fastcondrate = infastcondrate;
	}

	void SetMPI(int inmyid, int innprocs)	{
		myid = inmyid;
		nprocs = innprocs;
	}

	virtual void ToStreamHeader(ostream& os)	{
		os << version << '\n';
		propchrono.ToStream(os);
		chronototal.ToStream(os);
		os << datafile << '\n';
		os << iscodon << '\n';
		os << codetype << '\n';
		os << fixtopo << '\n';
		os << NSPR << '\t' << NMHSPR << '\t' << NTSPR << '\n';
		os << topolambda << '\t' << topomu << '\t' << toponstep << '\n';
		os << NNNI << '\n';
		os << nspec << '\t' << ntspec << '\n';
		os << bpp << '\t' << nbpp << '\t' << ntbpp << '\t' << bppnstep << '\t' << bppname << '\t' << bppcutoff << '\t' << bppbeta << '\n';
		os << dc << '\n';
		os << fixbl << '\n';
		os << proposemode << '\n';
		os << allocmode << '\n';
		os << sumratealloc << '\n';
		os << fasttopo << '\t' << fasttopofracmin << '\t' << fasttoponstep << '\n';
		os << fastcondrate << '\n';
		os << sumovercomponents << '\n';
		SetNamesFromLengths();
		GetTree()->ToStream(os);
	}

	virtual void FromStreamHeader(istream& is)	{
		is >> version;
		if (atof(version.substr(0,3).c_str()) < 1.2)	{
			cerr << "error: version is too old : " << version << '\n';
			exit(1);
		}
		propchrono.FromStream(is);
		chronototal.FromStream(is);
		string indatafile;
		is >> datafile;
		is >> iscodon;
		is >> codetype;
		is >> fixtopo;
		is >> NSPR >> NMHSPR >> NTSPR;
		is >> topolambda >> topomu >> toponstep;
		is >> NNNI;
		is >> nspec >> ntspec;
		is >> bpp >> nbpp >> ntbpp >> bppnstep >> bppname >> bppcutoff >> bppbeta;
		is >> dc;
		is >> fixbl;
		is >> proposemode;
		is >> allocmode;
		is >> sumratealloc;
		is >> fasttopo >> fasttopofracmin >> fasttoponstep;
		is >> fastcondrate;
		is >> sumovercomponents;
		is >> treestring;
	}

	virtual void ToStream(ostream& os)	{
		cerr << "error: in phyloprocess::ToStream\n";
		exit(1);
	}

	virtual void FromStream(istream& is)	{
		cerr << "error: in phyloprocess::FromStream\n";
		exit(1);
	}

	const TaxonSet* GetTaxonSet() const {
		return GetData()->GetTaxonSet();
	}

	void GlobalUpdateConditionalLikelihoods();
	virtual void SlaveUpdateConditionalLikelihoods();

	double GlobalComputeNodeLikelihood(const Link* from, int auxindex = -1);
	virtual void SlaveComputeNodeLikelihood(int,int);
	double LocalComputeNodeLikelihood(int,int);

	// Feb 1st, 2016
	// special device set up for calculating, on the fly, the likelihood summed over profile allocations
	// valid only under specific settings
	virtual double GlobalGetFullLogLikelihood();
	virtual void SlaveGetFullLogLikelihood();
	virtual double GetFullLogLikelihood()	{
		return logL;
	}

	protected:

	virtual SequenceAlignment* GetData() {
		if (! data)	{
			cerr << "error in PhyloProcess::GetData: no data\n";
			exit(1);
		}
		return data;
	}

	virtual const SequenceAlignment* GetData() const {
		if (! data)	{
			cerr << "error in PhyloProcess::GetData: no data\n";
			exit(1);
		}
		return data;
	}

	virtual StateSpace* GetStateSpace() {return GetData()->GetStateSpace();}

	double GetMinStat(double* profile, int site)	{
		double min = 1;
		for (int k=0; k<GetDim(); k++)	{
			if (observedarray[site][k])	{
				if (min > profile[k])	{
					min = profile[k];
				}
			}
		}
		return min;
	}

	// returns total number of taxa in the analysis
	int GetNtaxa()	{
		return GetData()->GetNtaxa();
	}

	double* GetEmpiricalFreq()	{
		return empfreq;
	}

	// the following methods are particularly important for MPI
	// Create / Delete / Unfold and Collapse should probably be specialized
	// according to whether this is a slave or the master processus
	virtual void Create();
	virtual void Delete();

	public :

	virtual void Unfold();
	virtual void Collapse();

	void GlobalBroadcastTree();
	virtual void SlaveBroadcastTree();

	// MPI Global dispatcher functions
	// all methods with a "Global" prefix are the ones through which MPI parallelization is supposed to go
	// these methods should dispatch the computation over slaves, and possibly, collect the result of the computation
	// in the non MPI version, they just call their non MPI counterpart (or nearly so)

	void SetFixTopo(bool in)	{
		fixtopo = in;
	}

	bool FixTopo()	{
		return fixtopo;
	}
	
	public: 

	void QuickUpdate()	{

		MPI_Status stat;
		MESSAGE signal = BCAST_TREE;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
		GlobalBroadcastTree();
		
		GlobalCollapse();
		GlobalUnfold();
	}

	virtual void SetDataFromLeaves()	{
		for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
			if (ActiveSite(i))	{
				RecursiveSetDataFromLeaves(i,GetRoot());
			}
		}
	}

	void RecursiveSetDataFromLeaves(int site, const Link* from)	{
		if (from->isLeaf())	{
			if (GetData()->GetBKState(GetNodeIndex(from->GetNode()),site) != -1)	{
				GetData()->SetState(GetNodeIndex(from->GetNode()),site,nodestate[GetNodeIndex(from->GetNode())][site]);
			}
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveSetDataFromLeaves(site, link->Out());
		}
	}

	void GlobalUnclamp();
	void GlobalRestoreData();
	void GlobalSetDataFromLeaves();
	double GlobalGetMeanDiversity();
	void SlaveUnclamp();
	void SlaveRestoreData();
	void SlaveSetDataFromLeaves();
	void SlaveGetMeanDiversity();

	void GlobalSetNodeStates();
	void SlaveSetNodeStates();
	void WriteNodeStates(ostream& os, const Link* from);

	virtual void ReadPB(int argc, char* argv[]);
	virtual void Read(string name, int burnin, int every, int until);
	virtual void ReadSiteLogL(string name, int burnin, int every, int until);
	virtual void ReadCV(string testdatafile, string name, int burnin, int every, int until, int iscodon = 0, GeneticCodeType codetype = Universal);
	virtual void PostPred(int ppredtype, string name, int burnin, int every, int until);

	void ReadSiteRates(string name, int burnin, int every, int until);

	// The following methids are here to write the mappings.
	void ReadMap(string name, int burnin, int every, int until);
	void ReadPostPredMap(string name, int burnin, int every, int until);
	void GlobalWriteMappings(string name);
	virtual void SlaveWriteMappings();
	void WriteTreeMapping(ostream& os, const Link* from, int i);



	virtual void GlobalSetTestData();
	virtual void SlaveSetTestData();
	void SetTestSiteMinAndMax();
	virtual void SlaveComputeCVScore() {
		cerr << "error: in PhyloProcess::native slave compute cv score\n";
		exit(1);
	}
	
	virtual void SlaveComputeSiteLogL() {
		cerr << "error: in PhyloProcess::slave compute site logL\n";
		exit(1);
	}
	
	double GetLogProb()	{
		return GetLogPrior() + GetLogLikelihood();
	}

	double GetLogPrior()	{
		return 0;
	}

	double GetLogLikelihood()	{
		return logL;
	}

	virtual void GlobalUnfold();
	virtual void GlobalCollapse();

	virtual void SlaveUnfold();
	virtual void SlaveCollapse();

	void GlobalReset(const Link* from, bool condalloc = false);
	void GlobalMultiply(const Link* from, const Link* to, bool condalloc = false);
	void GlobalMultiplyByStationaries(const Link* from, bool condalloc = false);
	void GlobalInitialize(const Link* from, const Link* link, bool condalloc = false);

	void GlobalPropagate(const Link* from, const Link* to, double time, bool condalloc = false);
	double GlobalProposeMove(const Branch* branch, double tuning);
	void GlobalRestore(const Branch* branch);

	void GlobalRootAtRandom();

	void RecursiveGibbsNNI(Link* from, double tuning, int type, int& success, int& moves);
	double GibbsNNI(double tuning, int);
	int  GlobalNNI(Link*,double,int);
	void GlobalPropagateOverABranch(Link*);
	virtual void SlavePropagateOverABranch(int);
	void SlaveNNI(int,int);
	void PropagateOverABranch(const Link*);
	int SlaveSendNNILikelihood(Link*);
	double SendRandomBranches(Link*,double,Link**&, int);

	// MCMC on branch lengths
	double BranchLengthMove(double tuning);
	double NonMPIBranchLengthMove(double tuning);

	// MCMC on topology
	virtual void CreateSuffStat();
	virtual void DeleteSuffStat();

	// if  auxindex == -1 then create auxiliary array, otherwise use condlmap[auxindex] as the auxiliary array
	double ComputeNodeLikelihood(const Link* from, int auxindex);
	// double ComputeNodeLikelihood(const Link* from, int auxindex = -1);

	// assumes that rate allocations have already been defined
	// and that conditional likelihoods are updated
	// those conditional likelihoods will be corrupted
	void SampleNodeStates();
	void SampleNodeStates(const Link* from, double*** aux);

	void SampleSiteNodeStates(int site);
	void SampleSiteNodeStates(int site, const Link* from, double** aux);

	// assumes that states at nodes have been sampled (using ResampleState())
	void SampleSubstitutionMappings(const Link* from);
	void SampleSiteSubstitutionMapping(int site, const Link* from);

	// conditional likelihood propagations
	void PostOrderPruning(const Link* from, double*** aux);
	void PreOrderPruning(const Link* from, double*** aux);

	void SitePostOrderPruning(int site, const Link* from);

	double RecursiveBranchLengthMove(const Link* from, double tuning, int& n);
	double RecursiveNonMPIBranchLengthMove(const Link* from, double tuning, int& n);
	double LocalBranchLengthMove(const Link* from, double tuning);
	double LocalNonMPIBranchLengthMove(const Link* from, double tuning);

	double SpecialSPR(int nrep);
	double NonMPISpecialSPR();
	double MPISpecialSPR();

	double TemperedSpecialSPR(int nrep, int nstep);
	double NonMPITemperedSpecialSPR(int nstep);
	double MPITemperedSpecialSPR(int nstep);

	double BPPSPR(int nrep);
	double NonMPIBPPSPR();
	double MPIBPPSPR();

	double TemperedBPPSPR(int nrep, int nstep);
	double NonMPITemperedBPPSPR(int nstep);
	double MPITemperedBPPSPR(int nstep);

	double TemperedGibbsSPR(double lambda, double mu, int nstep, int nrep);
	int MPITemperedGibbsSPR(double lambda, double mu, int nstep);
	int NonMPITemperedGibbsSPR(double lambda, double mu, int nstep);

	double GibbsSPR(int nrep);
	int MPIGibbsSPR();
	double NonMPIGibbsSPR();

	double GibbsMHSPR(double lambda, int nrep);
	int MPIGibbsMHSPR(double lambda);
	int NonMPIGibbsMHSPR(double lambda);

	void GlobalGibbsSPRScan(Link* down, Link* up, double* loglarray);
	void RecursiveGibbsSPRScan(Link* from, Link* fromup, Link* down, Link* up, double* loglarray, int& n);
	void RecursiveGibbsFillMap(Link* from, Link* fromup, map<pair<Link*,Link*>,double>& loglmap, double* loglarray, int& n);
	void RecursiveNonMPIGibbsSPRScan(Link* from, Link* fromup, Link* down, Link* up, map<pair<Link*,Link*>,double>& loglmap);

	// various protected accessors
	// used for computation and maintenance from within PhyloProcess classes
	const int* GetData(int index)	{
		return GetData()->GetState(index);
	}

	const int* GetData(const Link* from)	{
		if (! from->isLeaf())	{
			cerr << "error in PhyloProcess::GetData\n";
			exit(1);
		}
		return GetData()->GetState(from->GetNode()->GetIndex());
	}

	int* GetStates(const Node* node)	{
		return nodestate[GetNodeIndex(node)];
	}

	void CreateSiteConditionalLikelihoods();
	void DeleteSiteConditionalLikelihoods();

	void CreateConditionalLikelihoods();
	void DeleteConditionalLikelihoods();
	void UpdateConditionalLikelihoods();

	double*** GetConditionalLikelihoodVector(const Link* link)	{
		return condlmap[GetLinkIndex(link)];
	}

	void CreateMappings();
	void DeleteMappings();
	void FullDeleteMappings();

	void CreateNodeStates();
	void DeleteNodeStates();
	
	// sufficient statistics for rates and branch lengths (do not depend on the model)
	int GetSiteRateSuffStatCount(int site) {return siteratesuffstatcount[site];}
	double GetSiteRateSuffStatBeta(int site) {return siteratesuffstatbeta[site];}

	int GetBranchLengthSuffStatCount(int index) {return branchlengthsuffstatcount[index];}
	double GetBranchLengthSuffStatBeta(int index) {return branchlengthsuffstatbeta[index];}

	const int* GetBranchLengthSuffStatCount() {return branchlengthsuffstatcount;}
	const double* GetBranchLengthSuffStatBeta() {return branchlengthsuffstatbeta;}

	void GlobalUpdateSiteRateSuffStat();
	void GlobalUpdateBranchLengthSuffStat();

	virtual void SlaveUpdateSiteRateSuffStat();
	void SlaveUpdateBranchLengthSuffStat();

	void GlobalGetMeanSiteRate();
	void SlaveSendMeanSiteRate();

	void GlobalActivateSumOverRateAllocations();
	void GlobalInactivateSumOverRateAllocations();
	virtual void SlaveActivateSumOverRateAllocations();
	virtual void SlaveInactivateSumOverRateAllocations();

	virtual int CountMapping();
	virtual int CountMapping(int site);
	virtual int GlobalCountMapping();
	void SlaveCountMapping();

	void ReadData(string indatafile)	{
		datafile = indatafile;
		SequenceAlignment* plaindata;
		if (iscodon)	{
			SequenceAlignment* tempdata = new FileSequenceAlignment(datafile);
			data = new CodonSequenceAlignment(tempdata,true,codetype);
		}
		else	{
			data = new FileSequenceAlignment(datafile);
		}
		if (dc)	{
			data->DeleteConstantSites();
		}

		MakeObservedArray();
	}

	virtual void MakeObservedArray()	{

		observedarray = new int*[data->GetNsite()];
		for (int i=0; i<data->GetNsite(); i++)	{
			observedarray[i] = new int[data->GetNstate()];
			for (int k=0; k<data->GetNstate(); k++)	{
				observedarray[i][k] = 0;
			}
		}
		for (int i=0; i<data->GetNsite(); i++)	{
			for (int j=0; j<GetNtaxa(); j++)	{
				if (data->GetState(j,i))	{
					observedarray[i][data->GetState(j,i)] = 1;
				}
			}
		}	
	}

	virtual void SetTree(string treefile)	{
		if (treefile == "None")	{
			tree = new Tree(GetData()->GetTaxonSet());
			if (GetMyid() == 0)	{
				tree->MakeRandomTree();
				GlobalBroadcastTree();
			}
			else	{
				SlaveBroadcastTree();
			}
		}
		else	{
			tree = new Tree(treefile);
		}
		tree->RegisterWith(GetData()->GetTaxonSet());
		CloneTree();
		tree2->RegisterWith(GetData()->GetTaxonSet());
	}

	virtual void SetProfileDim()	{
		SetDim(GetData()->GetNstate());
	}

	virtual void New();

	void Open(istream& is);

	virtual double GetObservedCompositionalHeterogeneity()	{
		return GetData()->CompositionalHeterogeneity(0);
	}

	virtual double GetCompositionalHeterogeneity()	{
		return GetData()->CompositionalHeterogeneity(0);
	}

	double*** sitecondlmap;
	double**** condlmap;
	BranchSitePath*** submap;
	int** nodestate;

	// sufficient statistics for rates and branch lengths (do not depend on the model)
	int* siteratesuffstatcount;
	double* siteratesuffstatbeta;
	int* branchlengthsuffstatcount;
	double* branchlengthsuffstatbeta;

	string datafile;
	string treefile;
	string treestring;

	double* empfreq;

	Chrono propchrono;
	Chrono chronototal;
	Chrono nnichrono;
	Chrono sprchrono;
	Chrono tsprchrono;

	double* loglarray;

	int size;
	void IncSize()	{size++;}
	int GetSize() {return size;}
	void SetSize(int insize) {size = insize;}

	double GetNormFactor() {return GetNormalizationFactor();}

	double totaltime;

	SequenceAlignment* testdata;
	int testnsite;
	int testsitemin;
	int testsitemax;

	int fixtopo;
	int NSPR;
	int NMHSPR;
	int NTSPR;
	double topolambda;
	double topomu;
	int toponstep;
	int NNNI;
	int bpp;
	int nbpp;
	int ntbpp;
	int bppnstep;
	string bppname;
	double bppcutoff;
	double bppbeta;
	UnrootedBPP* BPP;
	int nspec;
	int ntspec;
	int dc;
	int iscodon;
	GeneticCodeType codetype;

	int** observedarray;

	double spracc;
	double sprtry;
	double mhspracc;
	double mhsprtry;
	double tspracc;
	double tsprtry;
	double tsprtmp;
	double tsprtmpacc00;
	double tsprtmpacc01;
	double tsprtmpacc10;
	double tsprtmpacc11;
	double tsprtot;
	double anntot;
	double anntmp;
	double anntmpacc00;
	double anntmpacc01;
	double anntmpacc10;
	double anntmpacc11;
	double nniacc;
	double nnitry;
	double bppsprtry;
	double bppspracc;
	double tbppsprtry;
	double tbppspracc;
	double specacc;
	double spectry;
	double tspecacc;
	double tspectry;
	double fasttopoacc;
	double fasttopotry;
	double fasttopochange;
	int fastcondrate;
	double ziptopotry;
	double ziptopoacc;

	int nprelim;
	int currenttopo;
	int sumovercomponents;
	int sumratealloc;

	int fasttopo;
	double fasttopofracmin;
	int fasttoponstep;

	SequenceAlignment* data;
};



#endif
