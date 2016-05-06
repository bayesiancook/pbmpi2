
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef MULTIGENERASCATSBDPGAMMAPHYLOPROCESS_H
#define MULTIGENERASCATSBDPGAMMAPHYLOPROCESS_H

#include "MultiGenePhyloProcess.h"
#include "RASCATSBDPGammaPhyloProcess.h"

class MultiGeneRASCATSBDPGammaPhyloProcess : public virtual MultiGenePhyloProcess, public virtual RASCATSBDPGammaPhyloProcess	{

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

	RASCATSBDPGammaPhyloProcess* GetProcess(int gene)	{
		RASCATSBDPGammaPhyloProcess* proc = dynamic_cast<RASCATSBDPGammaPhyloProcess*> (process[gene]);
		if ((! proc) && (process[gene]))	{
			cerr << "error in phyloprocess conversion\n";
			exit(1);
		}
		return proc;
	}

	MultiGeneRASCATSBDPGammaPhyloProcess() {}

	MultiGeneRASCATSBDPGammaPhyloProcess(int nratecat, int inkappaprior, int inglobalalpha, int inglobalbl)	{

		Ncat = nratecat;
		kappaprior = inkappaprior;
		globalalpha = inglobalalpha;
		globalbl = inglobalbl;
		if (! globalbl)	{
			hierarchicallengthprior = 1;
		}
	}

	MultiGeneRASCATSBDPGammaPhyloProcess(istream& is, int inmyid, int innprocs)	{

		// generic
		FromStreamHeader(is);

		// should be different
		SetMPI(inmyid,innprocs);

		// specific
		is >> Ncat;
		is >> kappaprior;
		is >> globalalpha;
		is >> globalbl;
		if (! globalbl)	{
			hierarchicallengthprior = 1;
		}

		Open(is);

	}

	~MultiGeneRASCATSBDPGammaPhyloProcess() {
		Delete();
	}

	void ToStreamHeader(ostream& os)	{
		PhyloProcess::ToStreamHeader(os);
		os << Ncat << '\n';
		os << kappaprior << '\n';
		os << globalalpha << '\n';
		os << globalbl << '\n';
	}

	void TraceHeader(ostream& os)	{
		os << "#iter\ttime\ttopo\tloglik";
		os << "\tlength";
		if (! GlobalBranchLengths())	{
			os << "\tlengthvar";
			os << "\tmean\trelvar\tmean\tvrelvar";
		}
		os << "\talpha";
		if (! GlobalAlpha())	{
			os << "\tvaralpha";
		}
		os << "\tNmode\tstatent\tstatalpha\tkappa";
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
			os << '\t' << meanbranchmean << '\t' << relvarbranchmean << '\t' << meanbranchrelvar << '\t' << relvarbranchrelvar;
		}
		os << '\t' << GetMeanAlpha();
		if (! GlobalAlpha())	{
			os << '\t' << varalpha;
		}
		os << '\t' << GlobalGetMeanNcomponent();
		os << '\t' << GlobalGetMeanStatEnt();
		os << '\t' << GlobalGetMeanStatAlpha();
		os << '\t' << GlobalGetMeanKappa();
		os << '\n';
	}

	virtual void Monitor(ostream& os)  {
		PhyloProcess::Monitor(os);
	}

	double Move(double tuning = 1.0);

	virtual double GlobalRestrictedMoveCycle(int nrep = 1, double tuning = 1.0)	{

		cerr << "in multi gene: globalrestriced move\n";
		exit(1);
		for (int rep=0; rep<nrep; rep++)	{

			GlobalUpdateParameters();
			GammaBranchProcess::Move(tuning,10);

			GlobalUpdateParameters();
			DGamRateProcess::Move(0.3*tuning,10);
			DGamRateProcess::Move(0.03*tuning,10);
			GlobalUpdateParameters();

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

};

#endif
