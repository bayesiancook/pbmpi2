
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef AASSRASCATSBDP_H
#define AASSRASCATSBDP_H

#include "GTRSBDPProfileProcess.h"
#include "AASubSelMixtureProfileProcess.h"
#include "GeneralPathSuffStatMatrixMixtureProfileProcess.h"
#include "GeneralPathSuffStatMatrixPhyloProcess.h"
#include "GeneralPathSuffStatMatrixSubstitutionProcess.h"
#include "DGamRateProcess.h"
#include "GammaBranchProcess.h"

class AASubSelSBDPProfileProcess : public virtual GTRSBDPProfileProcess, public virtual AASubSelMixtureProfileProcess, public virtual GeneralPathSuffStatMatrixMixtureProfileProcess	{

	public:

	AASubSelSBDPProfileProcess() {}
	virtual ~AASubSelSBDPProfileProcess() {}

	protected:

	virtual void SwapComponents(int cat1, int cat2)	{
		SBDPProfileProcess::SwapComponents(cat1,cat2);
	}

	void Create()	{
		GTRSBDPProfileProcess::Create();
		AASubSelMixtureProfileProcess::Create();
		GeneralPathSuffStatMatrixMixtureProfileProcess::Create();
	}

	virtual void Delete()	{
		GeneralPathSuffStatMatrixMixtureProfileProcess::Delete();
		AASubSelMixtureProfileProcess::Delete();
		GTRSBDPProfileProcess::Delete();
	}


};

class AASubSelRASCATSBDPSubstitutionProcess : public virtual GeneralPathSuffStatMatrixSubstitutionProcess, public virtual DGamRateProcess, public virtual AASubSelSBDPProfileProcess {

	public:

	AASubSelRASCATSBDPSubstitutionProcess() {}
	virtual ~AASubSelRASCATSBDPSubstitutionProcess() {}

	virtual void Create()	{
		GeneralPathSuffStatMatrixSubstitutionProcess::Create();
		DGamRateProcess::Create();
		AASubSelSBDPProfileProcess::Create();
	}

	virtual void Delete()	{
		AASubSelSBDPProfileProcess::Delete();
		DGamRateProcess::Delete();
		GeneralPathSuffStatMatrixSubstitutionProcess::Delete();
	}

};

class AASubSelRASCATSBDPGammaPhyloProcess : public virtual GeneralPathSuffStatMatrixPhyloProcess, public virtual AASubSelRASCATSBDPSubstitutionProcess, public virtual GammaBranchProcess	{

	public:

	AASubSelRASCATSBDPGammaPhyloProcess(int nratecat, int inkappaprior)	{

		Ncat = nratecat;
		kappaprior = inkappaprior;
	}

	AASubSelRASCATSBDPGammaPhyloProcess(istream& is, int inmyid, int innprocs)	{

		// generic
		FromStreamHeader(is);
		SetMPI(inmyid,innprocs);

		// specific
		is >> Ncat;
		is >> kappaprior;

		Open(is);

	}

	virtual ~AASubSelRASCATSBDPGammaPhyloProcess() {
		Delete();
	}

	void ToStreamHeader(ostream& os)	{
		PhyloProcess::ToStreamHeader(os);
		os << Ncat << '\n';
		os << kappaprior << '\n';
	}

	virtual void ReadPB(int argc, char* argv[]);
	void ReadNocc(string name, int burnin, int every, int until);
	void ReadTestProfile(string name, int nrep, double tuning, int burnin, int every, int until);
	void ReadRelRates(string name, int burnin, int every, int until);
	void SlaveComputeCVScore();
	void SlaveComputeSiteLogL();

	void SlaveUpdateParameters();
	void GlobalUpdateParameters();

	virtual void SlaveExecute(MESSAGE);
	void UpdateRRSuffStat() {}
	void GlobalUpdateRRSuffStat() {}
	void SlaveUpdateRRSuffStat() {}

	void TraceHeader(ostream& os)	{
		os << "#time\ttimeperccyle\ttopo\tlnL\tlength\talpha\tNmode\tstatent\tstatalpha";
		if (! fixrr)	{
			os << "\trrent\trrmean";
		}
		os << '\n'; 
	}

	void Trace(ostream& os)	{

		UpdateOccupancyNumbers();

		os << ((int) (chronototal.GetTime() / 1000));
		if (chronototal.GetTime())	{
			os << '\t' << ((double) ((int) (chronototal.GetTime() / (1 + GetSize())))) / 1000;
			os << '\t' << ((int) (propchrono.GetTime() / chronototal.GetTime() * 100));
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
		if (! fixrr)	{
			os << '\t' << GetRREntropy();
			os << '\t' << GetRRMean();
		}
		os << '\n';
	}

	double Move(double tuning = 1.0)	{

		chronototal.Start();

		propchrono.Start();
		if (! FixBL())	{
			BranchLengthMove(tuning);
			BranchLengthMove(0.1 * tuning);
		}

		if (! fixtopo)	{
			MoveTopo();
		}
		propchrono.Stop();

		for (int rep=0; rep<5; rep++)	{
			GlobalCollapse();
			AugmentedMove(tuning);
			GlobalUnfold();
		}

		chronototal.Stop();

		return 1;
	}

	double AugmentedMove(double tuning = 1)	{

		if (! FixBL())	{
			GammaBranchProcess::Move(tuning,10);
			GlobalUpdateParameters();
		}


		if (! FixAlpha())	{
			DGamRateProcess::Move(tuning,10);
			DGamRateProcess::Move(0.3*tuning,10);
			DGamRateProcess::Move(0.03*tuning,10);
		}

		GlobalUpdateParameters();
		AASubSelSBDPProfileProcess::Move(1,1,2);
		GlobalUpdateParameters();

		if (! FixRR()){
			LengthRelRateMove(1,10);
			LengthRelRateMove(0.1,10);
			LengthRelRateMove(0.01,10);
		}
	}

	void ToStream(ostream& os)	{
		GammaBranchProcess::ToStream(os);
		DGamRateProcess::ToStream(os);
		AASubSelSBDPProfileProcess::ToStream(os);
	}

	void FromStream(istream& is)	{
		GammaBranchProcess::FromStream(is);
		DGamRateProcess::FromStream(is);
		AASubSelSBDPProfileProcess::FromStream(is);
	}


	protected:

	virtual void Create()	{
		GeneralPathSuffStatMatrixPhyloProcess::Create();
		AASubSelRASCATSBDPSubstitutionProcess::Create();
		GammaBranchProcess::Create();
	}
		
	virtual void Delete()	{
		GammaBranchProcess::Delete();
		AASubSelRASCATSBDPSubstitutionProcess::Delete();
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

