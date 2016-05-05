
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef MULTIGENERASCATGTRSBDPGAMMAPHYLOPROCESS_H
#define MULTIGENERASCATGTRSBDPGAMMAPHYLOPROCESS_H

#include "MultiGenePhyloProcess.h"
#include "RASCATGTRSBDPGammaPhyloProcess.h"
#include "ExpoConjugateGTRProfileProcess.h"

class MultiGeneExpoConjugateGTRProfileProcess : public virtual ExpoConjugateGTRProfileProcess, public virtual MultiGeneProfileProcess	{

	public:

	MultiGeneExpoConjugateGTRProfileProcess()	{}
	virtual ~MultiGeneExpoConjugateGTRProfileProcess() {}

	// re-implement: should gather gene-specific suff stats
	virtual void UpdateRRSuffStat();

	ExpoConjugateGTRProfileProcess* GetGTRProfileProcess(int gene)	{
		ExpoConjugateGTRProfileProcess* tmp = dynamic_cast<ExpoConjugateGTRProfileProcess*>(process[gene]);
		if (! tmp)	{
			cerr << "error in GetGTRProfileProcess\n";
			exit(1);
		}
		return tmp;
	}

};


class MultiGeneRASCATGTRSBDPGammaPhyloProcess : public virtual MultiGenePhyloProcess, public virtual MultiGeneExpoConjugateGTRProfileProcess, public virtual RASCATGTRSBDPGammaPhyloProcess	{

	public:

        virtual void SlaveExecute(MESSAGE);

	// broadcast only the global parameters (not the gene-specific ones)
	virtual void GlobalUpdateParameters();
	// receive global parameters and then send them to genes
	virtual void SlaveUpdateParameters();

	virtual void SlaveUpdateSiteRateSuffStat()	{
		MultiGenePhyloProcess::SlaveUpdateSiteRateSuffStat();
	}

	virtual void GlobalUpdateSiteProfileSuffStat()	{
		MultiGenePhyloProcess::GlobalUpdateSiteProfileSuffStat();
	}

	virtual void SlaveUpdateSiteProfileSuffStat()	{
		MultiGenePhyloProcess::SlaveUpdateSiteProfileSuffStat();
	}

	virtual void UpdateBranchLengthSuffStat()	{
		MultiGenePhyloProcess::UpdateBranchLengthSuffStat();
	}

	virtual void UpdateRateSuffStat()	{
		MultiGenePhyloProcess::UpdateRateSuffStat();
	}

	virtual void UpdateRRSuffStat()	{
		MultiGeneExpoConjugateGTRProfileProcess::UpdateRRSuffStat();
	}

	double GlobalGetMeanNcomponent();
	double GlobalGetMeanStatEnt();
	double GlobalGetMeanStatAlpha();
	double GlobalGetMeanKappa();

	void SlaveGetMeanNcomponent();
	void SlaveGetMeanStatEnt();
	void SlaveGetMeanStatAlpha();
	void SlaveGetMeanKappa();

	virtual void Create();
	virtual void Delete();

	RASCATGTRSBDPGammaPhyloProcess* GetProcess(int gene)	{
		RASCATGTRSBDPGammaPhyloProcess* proc = dynamic_cast<RASCATGTRSBDPGammaPhyloProcess*> (process[gene]);
		if ((! proc) && (process[gene]))	{
			cerr << "error in phyloprocess conversion\n";
			exit(1);
		}
		return proc;
	}

	MultiGeneRASCATGTRSBDPGammaPhyloProcess() {}

	MultiGeneRASCATGTRSBDPGammaPhyloProcess(int nratecat, string inrrtype, int inkappaprior, int inglobalalpha, int inglobalbl)	{

		Ncat = nratecat;
		rrtype = inrrtype;
		kappaprior = inkappaprior;
		globalalpha = inglobalalpha;
		globalbl = inglobalbl;
		if (! globalbl)	{
			hierarchicallengthprior = 1;
		}
	}

	MultiGeneRASCATGTRSBDPGammaPhyloProcess(istream& is, int inmyid, int innprocs)	{

		// generic
		FromStreamHeader(is);

		// should be different
		SetMPI(inmyid,innprocs);

		// specific
		is >> Ncat;
		is >> kappaprior;
		is >> rrtype;
		is >> globalalpha;
		is >> globalbl;
		if (! globalbl)	{
			hierarchicallengthprior = 1;
		}

		Open(is);

	}

	~MultiGeneRASCATGTRSBDPGammaPhyloProcess() {
		Delete();
	}

	void ToStreamHeader(ostream& os)	{
		PhyloProcess::ToStreamHeader(os);
		os << Ncat << '\n';
		os << kappaprior << '\n';
		os << rrtype << '\n';
		os << globalalpha << '\n';
		os << globalbl << '\n';
	}

	void TraceHeader(ostream& os)	{
		os << "#iter\ttime\ttopo\tloglik";
		os << "\tlength";
		if (! GlobalBranchLengths())	{
			os << "\tlengthvar";
		}
		os << "\talpha";
		if (! GlobalAlpha())	{
			os << "\tvaralpha";
		}
		os << "\tNmode\tstatent\tstatalpha\tkappa";
		if (! fixrr)	{
			os << "\trrent\trrmean";
		}
		os << '\n'; 
	}

	virtual double GetNormalizationFactor()	{
		return 1.0;
	}

	void Trace(ostream& os)	{

		os << GetSize();
		if (chronototal.GetTime())	{
			os << '\t' << chronototal.GetTime() / 1000;
			os << '\t' << ((int) (propchrono.GetTime() / chronototal.GetTime() * 100));
			chronototal.Reset();
			propchrono.Reset();
		}
		else	{
			os << '\t' << 0;
			os << '\t' << 0;
		}

		os << '\t' << GetLogLikelihood();
		os << '\t' << GetMeanTotalLength();
		if (! GlobalBranchLengths())	{
			os << '\t' << GetMeanLengthRelVar();
		}
		os << '\t' << GetMeanAlpha();
		if (! GlobalAlpha())	{
			os << '\t' << varalpha;
		}
		os << '\t' << GlobalGetMeanNcomponent();
		os << '\t' << GlobalGetMeanStatEnt();
		os << '\t' << GlobalGetMeanStatAlpha();
		os << '\t' << GlobalGetMeanKappa();
		if (! fixrr)	{
			os << '\t' << GetRREntropy();
			os << '\t' << GetRRMean();
		}
		os << '\n';

	}

	virtual void Monitor(ostream& os)  {
		PhyloProcess::Monitor(os);
	}

	double Move(double tuning = 1.0);

	virtual double GlobalRestrictedMoveCycle(int nrep = 1, double tuning = 1.0)	{

		for (int rep=0; rep<nrep; rep++)	{

			GlobalUpdateParameters();
			GammaBranchProcess::Move(tuning,10);

			GlobalUpdateParameters();
			DGamRateProcess::Move(0.3*tuning,10);
			DGamRateProcess::Move(0.03*tuning,10);
			//ExpoConjugateGTRSBDPProfileProcess::Move(1,1,1);
			GlobalUpdateParameters();
			MoveRR();
			GlobalUpdateParameters();
			GlobalGeneMove();
		}
		return 1;
	}

	virtual void PrepareSiteLogLikelihood(int site) {
		cerr << "in multi gene: prepare site log likelihood\n";
		exit(1);
		/*
		int cat = ExpoConjugateGTRSBDPProfileProcess::alloc[site];
		if (! matrixarray[cat])	{
			CreateMatrix(cat);
		}
		UpdateMatrix(cat);
		*/
	}

	void ToStream(ostream& os)	{
		GammaBranchProcess::ToStream(os);
		DGamRateProcess::ToStream(os);
		// ExpoConjugateGTRProfileProcess::ToStream(os);
		// ExpoConjugateGTRSBDPProfileProcess::ToStream(os);
	}

	void FromStream(istream& is)	{
		GammaBranchProcess::FromStream(is);
		DGamRateProcess::FromStream(is);
		// ExpoConjugateGTRProfileProcess::FromStream(is);
		// ExpoConjugateGTRSBDPProfileProcess::FromStream(is);
		GlobalUpdateParameters();
	}

	/*
	virtual void ReadPB(int argc, char* argv[]);
	void ReadNocc(string name, int burnin, int every, int until);
	void ReadTestProfile(string name, int nrep, double tuning, int burnin, int every, int until);
	void ReadRelRates(string name, int burnin, int every, int until);
	void SlaveComputeCVScore();
	void SlaveComputeSiteLogL();
	*/
};

#endif

