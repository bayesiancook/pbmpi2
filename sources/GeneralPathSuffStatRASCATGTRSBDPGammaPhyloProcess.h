
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef GENPATHSSRASCATGTRSBDP_H
#define GENPATHSSRASCATGTRSBDP_H

#include "GTRSBDPProfileProcess.h"
#include "GeneralPathSuffStatMatrixMixtureProfileProcess.h"
#include "GeneralPathSuffStatMatrixPhyloProcess.h"
#include "GeneralPathSuffStatMatrixSubstitutionProcess.h"
#include "DGamRateProcess.h"
#include "GammaBranchProcess.h"

class GeneralPathSuffStatGTRSBDPProfileProcess : public virtual GTRSBDPProfileProcess, public virtual GeneralPathSuffStatMatrixMixtureProfileProcess	{

	public:

	GeneralPathSuffStatGTRSBDPProfileProcess() {}
	virtual ~GeneralPathSuffStatGTRSBDPProfileProcess() {}

	protected:

	void Create()	{
		GTRSBDPProfileProcess::Create();
		GeneralPathSuffStatMatrixMixtureProfileProcess::Create();
	}

	virtual void Delete()	{
		GeneralPathSuffStatMatrixMixtureProfileProcess::Delete();
		GTRSBDPProfileProcess::Delete();
	}
};

class GeneralPathSuffStatRASCATGTRSBDPSubstitutionProcess : public virtual GeneralPathSuffStatMatrixSubstitutionProcess, public virtual DGamRateProcess, public virtual GeneralPathSuffStatGTRSBDPProfileProcess {

	public:

	GeneralPathSuffStatRASCATGTRSBDPSubstitutionProcess() {}
	virtual ~GeneralPathSuffStatRASCATGTRSBDPSubstitutionProcess() {}

	virtual void Create()	{
		GeneralPathSuffStatMatrixSubstitutionProcess::Create();
		DGamRateProcess::Create();
		GeneralPathSuffStatGTRSBDPProfileProcess::Create();
	}

	virtual void Delete()	{
		GeneralPathSuffStatGTRSBDPProfileProcess::Delete();
		DGamRateProcess::Delete();
		GeneralPathSuffStatMatrixSubstitutionProcess::Delete();
	}

};

class GeneralPathSuffStatRASCATGTRSBDPGammaPhyloProcess : public virtual GeneralPathSuffStatMatrixPhyloProcess, public virtual GeneralPathSuffStatRASCATGTRSBDPSubstitutionProcess, public virtual GammaBranchProcess	{

	public:

	GeneralPathSuffStatRASCATGTRSBDPGammaPhyloProcess() {}

	GeneralPathSuffStatRASCATGTRSBDPGammaPhyloProcess(int nratecat, string inrrtype, int inkappaprior)	{

		Ncat = nratecat;
		rrtype = inrrtype;
		kappaprior = inkappaprior;
	}

	GeneralPathSuffStatRASCATGTRSBDPGammaPhyloProcess(istream& is, int inmyid, int innprocs)	{

		// generic
		FromStreamHeader(is);
		SetMPI(inmyid,innprocs);

		// specific
		is >> Ncat;
		is >> kappaprior;
		is >> rrtype;

		Open(is);

	}

	~GeneralPathSuffStatRASCATGTRSBDPGammaPhyloProcess() {
		Delete();
	}

	void ToStreamHeader(ostream& os)	{
		PhyloProcess::ToStreamHeader(os);
		os << Ncat << '\n';
		os << kappaprior << '\n';
		os << rrtype << '\n';
	}
	virtual void SlaveUpdateParameters();
	virtual void GlobalUpdateParameters();

	virtual void SlaveExecute(MESSAGE);
	virtual void UpdateRRSuffStat() {}
	virtual void GlobalUpdateRRSuffStat() {}
	virtual void SlaveUpdateRRSuffStat() {}

	void TraceHeader(ostream& os)	{
		os << "#iter\ttime\ttopo\tsuffstat\tloglik\tlength\talpha\tNmode\tstatent\tstatalpha\tkappa";
		if (! fixrr)	{
			os << "\trrent\trrmean";
		}
		os << '\n'; 
	}

	void Trace(ostream& os)	{

		UpdateOccupancyNumbers();

		os << GetSize();
		if (chronototal.GetTime())	{
			os << '\t' << chronototal.GetTime() / 1000;
			os << '\t' << ((int) (propchrono.GetTime() / chronototal.GetTime() * 100));
			os << '\t' << ((int) (sschrono.GetTime() / chronototal.GetTime() * 100));
			chronototal.Reset();
			propchrono.Reset();
			sschrono.Reset();
		}
		else	{
			os << '\t' << 0;
			os << '\t' << 0;
		}

		os << '\t' << GetLogLikelihood();
		os << '\t' << GetRenormTotalLength();
		os << '\t' << GetAlpha();
		os << '\t' << GetNDisplayedComponent();
		os << '\t' << GetStatEnt();
		os << '\t' << GetMeanDirWeight();
		os << '\t' << kappa;
		if (! fixrr)	{
			os << '\t' << GetRREntropy();
			os << '\t' << GetRRMean();
		}

		os << '\n';
	}

	virtual double Move(double tuning = 1.0)	{

		chronototal.Start();

		propchrono.Start();
		BranchLengthMove(tuning);
		BranchLengthMove(0.1 * tuning);
		if (! fixtopo)	{
			MoveTopo();
		}
		propchrono.Stop();

		GlobalCollapse();

		GammaBranchProcess::Move(tuning,10);

		GlobalUpdateParameters();
		DGamRateProcess::Move(0.3*tuning,10);
		DGamRateProcess::Move(0.03*tuning,10);

		GeneralPathSuffStatGTRSBDPProfileProcess::Move(1,1,10);

		if (! fixrr)	{
			LengthRelRateMove(1,10);
			LengthRelRateMove(0.1,10);
			LengthRelRateMove(0.01,10);
		}

		GlobalUnfold();

		chronototal.Stop();

		return 1;
	}

	void ToStream(ostream& os)	{
		GammaBranchProcess::ToStream(os);
		DGamRateProcess::ToStream(os);
		GeneralPathSuffStatGTRSBDPProfileProcess::ToStream(os);
	}

	void FromStream(istream& is)	{
		GammaBranchProcess::FromStream(is);
		DGamRateProcess::FromStream(is);
		GeneralPathSuffStatGTRSBDPProfileProcess::FromStream(is);
	}


	protected:

	virtual void Create()	{
		GeneralPathSuffStatMatrixPhyloProcess::Create();
		GeneralPathSuffStatRASCATGTRSBDPSubstitutionProcess::Create();
		GammaBranchProcess::Create();
	}
		
	virtual void Delete()	{
		GammaBranchProcess::Delete();
		GeneralPathSuffStatRASCATGTRSBDPSubstitutionProcess::Delete();
		GeneralPathSuffStatMatrixPhyloProcess::Delete();
	}


	double LengthRelRateMove(double tuning, int nrep)	{

		double naccept = 0;
		for (int rep=0; rep<nrep; rep++)	{
			double deltalogratio = - LogRRPrior() - LogLengthPrior();
			double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
			double e = exp(m);
			int nbranch = MoveAllBranches(e);
			for (int i=0; i<GetNrr(); i++)	{
				rr[i] /= e;
			}
			deltalogratio += LogRRPrior() + LogLengthPrior();
			deltalogratio += (nbranch-GetNrr()) * m;

			int accepted = (log(rnd::GetRandom().Uniform()) < deltalogratio);

			if (accepted)	{
				naccept++;
			}
			else	{
				MoveAllBranches(1.0/e);
				for (int i=0; i<GetNrr(); i++)	{
					rr[i] *= e;
				}
			}	
		}
		return naccept / nrep;
	}
};

#endif

