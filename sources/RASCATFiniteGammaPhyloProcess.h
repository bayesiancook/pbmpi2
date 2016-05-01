
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef RASCATFINITE_H
#define RASCATFINITE_H

#include "PoissonSubstitutionProcess.h"
#include "PoissonPhyloProcess.h"
#include "DGamRateProcess.h"
#include "PoissonFiniteProfileProcess.h"
#include "GammaBranchProcess.h"

class RASCATFiniteSubstitutionProcess : public virtual PoissonSubstitutionProcess, public virtual DGamRateProcess, public virtual PoissonFiniteProfileProcess {

	using PoissonSubstitutionProcess::UpdateZip;

	public:

	RASCATFiniteSubstitutionProcess() {}
	virtual ~RASCATFiniteSubstitutionProcess() {}

	protected:

	virtual void UpdateZip(int site)	{
		PoissonSubstitutionProcess::UpdateZip(site);
	}

	virtual void Create()	{
		PoissonSubstitutionProcess::Create();
		DGamRateProcess::Create();
		PoissonFiniteProfileProcess::Create();
	}

	virtual void Delete()	{
		PoissonFiniteProfileProcess::Delete();
		DGamRateProcess::Delete();
		PoissonSubstitutionProcess::Delete();
	}

};

class RASCATFiniteGammaPhyloProcess : public virtual PoissonPhyloProcess, public virtual RASCATFiniteSubstitutionProcess, public virtual GammaBranchProcess	{

	public:

        virtual void SlaveExecute(MESSAGE);
	virtual void GlobalUpdateParameters();
	virtual void SlaveUpdateParameters();

	RASCATFiniteGammaPhyloProcess() {}

	RASCATFiniteGammaPhyloProcess(int nratecat, int ncat, int infixncomp, int inempmix, string inmixtype)	{

		Ncat = nratecat;
		Ncomponent = ncat;
		fixncomp = infixncomp;
		empmix = inempmix;
		mixtype = inmixtype;
	}

	RASCATFiniteGammaPhyloProcess(istream& is, int inmyid, int innprocs)	{

		// generic
		FromStreamHeader(is);
		SetMPI(inmyid,innprocs);

		// specific
		is >> Ncat;
		is >> fixncomp;
		is >> empmix;
		is >> mixtype;

		Open(is);
	}

	void ToStreamHeader(ostream& os)	{
		PhyloProcess::ToStreamHeader(os);
		os << Ncat << '\n';
		os << fixncomp << '\t' << empmix << '\t' << mixtype << '\n';
	}

	~RASCATFiniteGammaPhyloProcess() {
		Delete();
	}

	void TraceHeader(ostream& os)	{
		os << "#iter\ttime\ttopo\tloglik\tlength\talpha\tNmode\tstatent\tstatalpha\ttopo";
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
		os << '\t' << GetNOccupiedComponent() << '\t' << GetStatEnt();
		os << '\t' << GetMeanDirWeight();
		os << '\t' << branchalpha << '\t' << branchbeta;
		os << '\t' << currenttopo;
		os << '\n';
	}

	virtual double Move(double tuning = 1.0)	{

		chronototal.Start();
		propchrono.Start();
		BranchLengthMove(tuning);
		BranchLengthMove(0.1 * tuning);
		if (! fixtopo)	{
			TopoMoveCycle(1,tuning);
		}
		propchrono.Stop();

		GlobalCollapse();

		GammaBranchProcess::Move(tuning,10);

		// this one is important 
		GlobalUpdateParameters();
		DGamRateProcess::Move(0.3*tuning,10);
		DGamRateProcess::Move(0.03*tuning,10);
		// RASCATSubstitutionProcess::MoveRate(tuning);

		// this one is not useful
		// because uniformized process:
		// conditional on discrete substitution mapping
		// profiles do not depend on branch lengths and site rates
		// GlobalUpdateParameters();

		PoissonFiniteProfileProcess::Move(1,1,5);

		GlobalUnfold();
		chronototal.Stop();

		// Trace(cerr);

		return 1;
	
	}

	virtual double GlobalGetFullLogLikelihood();
	virtual void SlaveGetFullLogLikelihood();
	virtual double GetFullLogLikelihood();

	virtual double GlobalRestrictedMoveCycle(int nrep = 1, double tuning = 1.0)	{

		for (int rep=0; rep<nrep; rep++)	{

			GlobalUpdateParameters();
			GammaBranchProcess::Move(tuning,10);

			GlobalUpdateParameters();
			DGamRateProcess::Move(0.3*tuning,10);
			DGamRateProcess::Move(0.03*tuning,10);

			// if (! empmix)	{
			PoissonFiniteProfileProcess::Move(1,1,1);
			// }
		}
		return 1;
	}

	virtual void ReadPB(int argc, char* argv[]);
	void SlaveComputeCVScore();
	void SlaveComputeSiteLogL();

	void ToStream(ostream& os)	{
		GammaBranchProcess::ToStream(os);
		DGamRateProcess::ToStream(os);
		PoissonFiniteProfileProcess::ToStream(os);
	}

	void FromStream(istream& is)	{
		GammaBranchProcess::FromStream(is);
		DGamRateProcess::FromStream(is);
		PoissonFiniteProfileProcess::FromStream(is);
		GlobalUpdateParameters();
	}


	virtual void Create()	{
		PoissonPhyloProcess::Create();
		RASCATFiniteSubstitutionProcess::Create();
		GammaBranchProcess::Create();
	}
		
	virtual void Delete()	{
		GammaBranchProcess::Delete();
		RASCATFiniteSubstitutionProcess::Delete();
		PoissonPhyloProcess::Delete();
	}
};

#endif

