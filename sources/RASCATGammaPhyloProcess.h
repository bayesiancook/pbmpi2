
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef RASCAT_H
#define RASCAT_H

#include "PoissonSubstitutionProcess.h"
#include "PoissonPhyloProcess.h"
#include "DGamRateProcess.h"
#include "PoissonDPProfileProcess.h"
#include "GammaBranchProcess.h"

class RASCATSubstitutionProcess : public virtual PoissonSubstitutionProcess, public virtual DGamRateProcess, public virtual PoissonDPProfileProcess {

	using PoissonSubstitutionProcess::UpdateZip;

	public:

	RASCATSubstitutionProcess() {}
	virtual ~RASCATSubstitutionProcess() {}

	protected:

	virtual void UpdateZip(int site)	{
		PoissonSubstitutionProcess::UpdateZip(site);
	}

	virtual void Create()	{
		PoissonSubstitutionProcess::Create();
		DGamRateProcess::Create();
		PoissonDPProfileProcess::Create();
	}

	virtual void Delete()	{
		PoissonDPProfileProcess::Delete();
		DGamRateProcess::Delete();
		PoissonSubstitutionProcess::Delete();
	}

};

class RASCATGammaPhyloProcess : public virtual PoissonPhyloProcess, public virtual RASCATSubstitutionProcess, public virtual GammaBranchProcess	{

	public:

        virtual void SlaveExecute(MESSAGE);
	virtual void GlobalUpdateParameters();
	virtual void SlaveUpdateParameters();

	RASCATGammaPhyloProcess() {}

	RASCATGammaPhyloProcess(int nratecat, int inwithpinv, int inkappaprior)	{

		Ncat = nratecat;
		withpinv = inwithpinv;
		kappaprior = inkappaprior;
	}

	RASCATGammaPhyloProcess(istream& is, int inmyid, int innprocs)	{

		// generic
		FromStreamHeader(is);
		SetMPI(inmyid,innprocs);

		// specific
		is >> Ncat;
		is >> withpinv;
		is >> kappaprior;

		Open(is);
	}

	void ToStreamHeader(ostream& os)	{
		PhyloProcess::ToStreamHeader(os);
		os << Ncat << '\n';
		os << kappaprior << '\n';
	}

	~RASCATGammaPhyloProcess() {
		Delete();
	}

	void TraceHeader(ostream& os)	{
		os << "#iter\ttime\ttopo\tloglik\tlength\talpha";
		// if (withpinv)	{
			os << "\tpinv";
		// }
		os << "\tNmode\tstatent\tstatalpha";
		os << "\tkappa";
		os << '\n'; 
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

		os << '\t' << GetLogLikelihood() << '\t' << GetRenormTotalLength() << '\t' << GetAlpha();
		// if (withpinv)	{
			os << '\t' << GetPinv();
		// }
		os << '\t' << GetNOccupiedComponent() << '\t' << GetStatEnt();
		os << '\t' << GetMeanDirWeight();
		os << '\t' << kappa;
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

		PoissonDPProfileProcess::Move(1,1,5);

		GlobalUnfold();
		chronototal.Stop();

		return 1;
	
	}

	virtual void ReadPB(int argc, char* argv[]);
	virtual void ReadStatMin(string name, int burnin, int every, int until);
	virtual void ReadProfileDistribution(string name, int burnin, int every, int until, int ndisc, double cialpha, int nsample);
	void ReadMeanDirWeight(string name, int burnin, int every, int until);

	void ToStream(ostream& os)	{
		GammaBranchProcess::ToStream(os);
		DGamRateProcess::ToStream(os);
		PoissonDPProfileProcess::ToStream(os);
	}

	void FromStream(istream& is)	{
		GammaBranchProcess::FromStream(is);
		DGamRateProcess::FromStream(is);
		PoissonDPProfileProcess::FromStream(is);
		GlobalUpdateParameters();
	}


	virtual void Create()	{
		PoissonPhyloProcess::Create();
		RASCATSubstitutionProcess::Create();
		GammaBranchProcess::Create();
	}
		
	virtual void Delete()	{
		GammaBranchProcess::Delete();
		RASCATSubstitutionProcess::Delete();
		PoissonPhyloProcess::Delete();
	}
};

#endif

