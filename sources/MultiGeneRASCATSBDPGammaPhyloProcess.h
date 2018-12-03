
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
#include "MultiGeneSBDPProfileProcess.h"
#include "RASCATSBDPGammaPhyloProcess.h"

class MultiGeneRASCATSBDPGammaPhyloProcess : public virtual MultiGenePhyloProcess, public virtual MultiGeneSBDPProfileProcess, public virtual RASCATSBDPGammaPhyloProcess	{

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

	virtual double GlobalGetFullLogLikelihood()	{
		return PhyloProcess::GlobalGetFullLogLikelihood();
	}

	virtual void SlaveGetFullLogLikelihood()	{
		MultiGenePhyloProcess::SlaveGetFullLogLikelihood();
	}

	virtual double GetFullLogLikelihood()	{
		cerr << "error: in multi gene GeFullLogLikelihood\n";
		exit(1);
	}

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

	MultiGeneRASCATSBDPGammaPhyloProcess(int nratecat, int inkappaprior, int innmodemax, int inglobalalpha, int inglobalbl, int inmappsuffstat)	{

		Ncat = nratecat;
		kappaprior = inkappaprior;
		nmodemax = innmodemax;
		globalalpha = inglobalalpha;
		globalbl = inglobalbl;
		mappsuffstat = inmappsuffstat;
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
		is >> nmodemax;
		is >> globalalpha;
		is >> globalbl;
		if (! globalbl)	{
			hierarchicallengthprior = 1;
			is >> mappsuffstat;
		}

		Open(is,0);

	}

	~MultiGeneRASCATSBDPGammaPhyloProcess() {
		Delete();
	}

	void ToStreamHeader(ostream& os)	{
		PhyloProcess::ToStreamHeader(os);
		os << Ncat << '\n';
		os << kappaprior << '\n';
		os << nmodemax << '\n';
		os << globalalpha << '\n';
		os << globalbl << '\n';
		if (! globalbl)	{
			os << mappsuffstat << '\n';
		}
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
		if (kappaprior == 2)	{
			os << "\tkappamean\trelvar";
		}
		os << '\n'; 
	}

	virtual double GetNormalizationFactor()	{
		return 1.0;
	}

	double GetLogLikelihood()	{
		return logL;
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
		if (kappaprior == 2)	{
			os << '\t' << kappamean;
			os << '\t' << kapparelvar;
		}
		os << '\n';
	}

	virtual void Monitor(ostream& os)  {
        if (monitorlogl)    {
            GlobalWriteSiteLogLToStream(os);
        }
        else    {
		/*
            GlobalCollectKappas();
            for (int gene=0; gene<Ngene; gene++)    {
                os << nmodearray[gene] << '\t';
            }
            for (int gene=0; gene<Ngene; gene++)    {
                os << statalphaarray[gene] << '\t';
            }
            for (int gene=0; gene<Ngene; gene++)    {
                os << kappaarray[gene] << '\t';
            }
            os << '\n';
	    */
	PhyloProcess::Monitor(os);
        }
	}

	double Move(double tuning = 1.0);

	virtual double GlobalRestrictedTemperedMove();	

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
			cerr << "error in prepare site log likelihood: matrix is not allocated\n";
			exit(1);
		}
		UpdateMatrix(cat);
		*/
	}

	void ToStream(ostream& os)	{
		MultiGenePhyloProcess::ToStream(os);
	}

	void FromStream(istream& is)	{
		MultiGenePhyloProcess::FromStream(is);
	}

	virtual void ReadPB(int argc, char* argv[]) {
        MultiGenePhyloProcess::ReadPB(argc,argv);
    }


	int nmodemax;
};

#endif

