
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef RASCATSBDP_H
#define RASCATSBDP_H

#include "RASCATGammaPhyloProcess.h"
#include "PoissonSBDPProfileProcess.h"
#include "CodonSequenceAlignment.h"

class RASCATSBDPSubstitutionProcess : public virtual RASCATSubstitutionProcess, public virtual PoissonSBDPProfileProcess {

	public:

	RASCATSBDPSubstitutionProcess() {}
	virtual ~RASCATSBDPSubstitutionProcess() {}

	protected:

	virtual void Create()	{
		RASCATSubstitutionProcess::Create();
		PoissonSBDPProfileProcess::Create();
	}

	virtual void Delete()	{
		PoissonSBDPProfileProcess::Delete();
		RASCATSubstitutionProcess::Delete();
	}

};

class RASCATSBDPGammaPhyloProcess : public virtual RASCATGammaPhyloProcess, public virtual RASCATSBDPSubstitutionProcess {

	public:

	RASCATSBDPGammaPhyloProcess() {}

	RASCATSBDPGammaPhyloProcess(int nratecat, int inwithpinv, int inkappaprior) : RASCATGammaPhyloProcess(nratecat, inwithpinv, inkappaprior)	{}

	RASCATSBDPGammaPhyloProcess(istream& is, int inmyid, int innprocs)	{

		// generic
		FromStreamHeader(is);
		SetMPI(inmyid,innprocs);

		// specific
		is >> Ncat;
		// is >> withpinv;
		is >> kappaprior;

		Open(is);
	}

	virtual void GlobalUpdateParameters();
	virtual void SlaveUpdateParameters();
	void SlaveComputeCVScore();
	void SlaveComputeSiteLogL();

	virtual double GlobalGetFullLogLikelihood();
	virtual void SlaveGetFullLogLikelihood();
	virtual double GetFullLogLikelihood();

	void ReadStatMin(string name, int burnin, int every, int until);

	void FromStream(istream& is)	{
		GammaBranchProcess::FromStream(is);
		DGamRateProcess::FromStream(is);
		PoissonSBDPProfileProcess::FromStream(is);
		GlobalUpdateParameters();
	}

	virtual double Move(double tuning = 1.0)	{

		chronototal.Start();
		propchrono.Start();

		if (! FixBL())	{
			if ((topobf != 3) && ((topobf != 1) || (size < bfburnin)))	{
				BranchLengthMove(tuning);
				BranchLengthMove(0.1 * tuning);
			}
		}

		if (! FixTopo())	{
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

	virtual double AugmentedMove(double tuning = 1.0)	{

		// important to start with that one
		// if marginal suff stat move is done before that in a multi gene context
		if (! FixBL())	{
			GammaBranchProcess::Move(tuning,10);
			GammaBranchProcess::Move(0.1*tuning,10);
			GlobalUpdateParameters();
		}

		if (! FixAlpha())	{
			DGamRateProcess::Move(tuning,10);
			DGamRateProcess::Move(0.3*tuning,10);
			DGamRateProcess::Move(0.03*tuning,10);
		}

		PoissonSBDPProfileProcess::Move(1,5,1,1);
		GlobalUpdateParameters();
		/*
		if (iscodon)	{
			PoissonSBDPProfileProcess::Move(0.1,1,15);
			PoissonSBDPProfileProcess::Move(0.01,1,15);
		}
		*/
	}

	virtual double GlobalRestrictedTemperedMove();

	virtual double GlobalRestrictedMoveCycle(int nrep = 1, double tuning = 1.0)	{

		for (int rep=0; rep<nrep; rep++)	{

			GlobalUpdateParameters();
			GammaBranchProcess::Move(tuning,10);

			GlobalUpdateParameters();
			DGamRateProcess::Move(0.3*tuning,10);
			DGamRateProcess::Move(0.03*tuning,10);

			GlobalUpdateParameters();
			PoissonSBDPProfileProcess::Move(1,1,1);
			GlobalUpdateParameters();
		}
		return 1;
	}

	virtual void Monitor(ostream& os)  {
		PhyloProcess::Monitor(os);
		os << "weight " << '\t' << GetMaxWeightError() << '\n';
		ResetMaxWeightError();
	}

	virtual void Create()	{
		RASCATSBDPSubstitutionProcess::Create();
		RASCATGammaPhyloProcess::Create();
	}
		
	virtual void Delete()	{
		RASCATGammaPhyloProcess::Delete();
		RASCATSBDPSubstitutionProcess::Delete();
	}

	void SlaveExecute(MESSAGE signal);
};

#endif

