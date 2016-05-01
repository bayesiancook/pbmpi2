
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef GENPATHSSRASCATGTR_H
#define GENPATHSSRASCATGTR_H

#include "GTRFiniteProfileProcess.h"
#include "GeneralPathSuffStatMatrixMixtureProfileProcess.h"
#include "GeneralPathSuffStatMatrixPhyloProcess.h"
#include "GeneralPathSuffStatMatrixSubstitutionProcess.h"
#include "DGamRateProcess.h"
#include "GammaBranchProcess.h"

class GeneralPathSuffStatGTRFiniteProfileProcess : public virtual GTRFiniteProfileProcess, public virtual GeneralPathSuffStatMatrixMixtureProfileProcess	{

	public:

	GeneralPathSuffStatGTRFiniteProfileProcess() {}
	virtual ~GeneralPathSuffStatGTRFiniteProfileProcess() {}

	protected:

	virtual void Create()	{
		GTRFiniteProfileProcess::Create();
		GeneralPathSuffStatMatrixMixtureProfileProcess::Create();
	}

	virtual void Delete()	{
		GeneralPathSuffStatMatrixMixtureProfileProcess::Delete();
		GTRFiniteProfileProcess::Delete();
	}
};

class GeneralPathSuffStatRASCATGTRFiniteSubstitutionProcess : public virtual GeneralPathSuffStatMatrixSubstitutionProcess, public virtual DGamRateProcess, public virtual GeneralPathSuffStatGTRFiniteProfileProcess {

	public:

	GeneralPathSuffStatRASCATGTRFiniteSubstitutionProcess() {}
	virtual ~GeneralPathSuffStatRASCATGTRFiniteSubstitutionProcess() {}

	virtual void Create()	{
		GeneralPathSuffStatMatrixSubstitutionProcess::Create();
		DGamRateProcess::Create();
		GeneralPathSuffStatGTRFiniteProfileProcess::Create();
	}

	virtual void Delete()	{
		GeneralPathSuffStatGTRFiniteProfileProcess::Delete();
		DGamRateProcess::Delete();
		GeneralPathSuffStatMatrixSubstitutionProcess::Delete();
	}

};

class GeneralPathSuffStatRASCATGTRFiniteGammaPhyloProcess : public virtual GeneralPathSuffStatMatrixPhyloProcess, public virtual GeneralPathSuffStatRASCATGTRFiniteSubstitutionProcess, public virtual GammaBranchProcess	{

	public:

	GeneralPathSuffStatRASCATGTRFiniteGammaPhyloProcess() {}

	GeneralPathSuffStatRASCATGTRFiniteGammaPhyloProcess(int nratecat, int ncat, int infixncomp, int inempmix, string inmixtype, string inrrtype)	{

		Ncat = nratecat;
		Ncomponent = ncat;
		fixncomp = infixncomp;
		empmix = inempmix;
		mixtype = inmixtype;
		rrtype = inrrtype;
	}

	GeneralPathSuffStatRASCATGTRFiniteGammaPhyloProcess(istream& is, int inmyid, int innprocs)	{

		// generic
		FromStreamHeader(is);
		SetMPI(inmyid,innprocs);

		// specific
		is >> Ncat;
		is >> fixncomp;
		is >> empmix;
		is >> mixtype;
		is >> rrtype;

		Open(is);
	}

	void ToStreamHeader(ostream& os)	{
		PhyloProcess::ToStreamHeader(os);
		os << Ncat << '\n';
		os << fixncomp << '\t' << empmix << '\t' << mixtype << '\n';
		os << rrtype << '\n';
	}

	~GeneralPathSuffStatRASCATGTRFiniteGammaPhyloProcess() {
		Delete();
	}

	virtual void SlaveUpdateParameters();
	virtual void GlobalUpdateParameters();

	virtual void SlaveExecute(MESSAGE);
	virtual void UpdateRRSuffStat() {}
	virtual void GlobalUpdateRRSuffStat() {}
	virtual void SlaveUpdateRRSuffStat() {}

	void TraceHeader(ostream& os)	{
		os << "#time\ttimeperccyle\ttopo\tlnL\tlength\talpha\tNmode\tstatent\tstatalpha";
		if (! fixrr)	{
			os << "\trrent\trrmean";
		}
		os << '\n'; 
	}

	void Trace(ostream& os)	{

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
		BranchLengthMove(tuning);
		BranchLengthMove(0.1 * tuning);
		if (!fixtopo)	{
			MoveTopo();
		}
		propchrono.Stop();

		GlobalCollapse();

		GammaBranchProcess::Move(tuning,10);

		GlobalUpdateParameters();
		DGamRateProcess::Move(0.3*tuning,10);
		DGamRateProcess::Move(0.03*tuning,10);

		GeneralPathSuffStatGTRFiniteProfileProcess::Move(1,1,5);

		if (! fixrr)	{
			LengthRelRateMove(1,10);
			LengthRelRateMove(0.1,10);
			LengthRelRateMove(0.01,10);
		}

		GlobalUnfold();

		chronototal.Stop();

		return 1;
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

	void ToStream(ostream& os)	{

		os << datafile << '\n';
		os << GetNcat() << '\n';
		GetTree()->ToStream(os);

		GammaBranchProcess::ToStream(os);
		DGamRateProcess::ToStream(os);
		GeneralPathSuffStatGTRFiniteProfileProcess::ToStream(os);
	}

	void FromStream(istream& is)	{
		GammaBranchProcess::FromStream(is);
		DGamRateProcess::FromStream(is);
		GeneralPathSuffStatGTRFiniteProfileProcess::FromStream(is);
	}


	protected:

	virtual void Create()	{
		GeneralPathSuffStatMatrixPhyloProcess::Create();
		GeneralPathSuffStatRASCATGTRFiniteSubstitutionProcess::Create();
		GammaBranchProcess::Create();
	}
		
	virtual void Delete()	{
		GammaBranchProcess::Delete();
		GeneralPathSuffStatRASCATGTRFiniteSubstitutionProcess::Delete();
		GeneralPathSuffStatMatrixPhyloProcess::Delete();
	}
};

#endif

