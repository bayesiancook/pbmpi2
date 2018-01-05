
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef IIDDIRGAM_H
#define IIDDIRGAM_H

#include "PoissonSubstitutionProcess.h"
#include "PoissonPhyloProcess.h"
#include "GammaRateProcess.h"
#include "PoissonSiteSpecificProfileProcess.h"
#include "GammaBranchProcess.h"

class IIDDirichletIIDGammaSubstitutionProcess : public virtual PoissonSubstitutionProcess, public virtual GammaRateProcess, public virtual PoissonSiteSpecificProfileProcess {

	using PoissonSubstitutionProcess::UpdateZip;

	public:

	IIDDirichletIIDGammaSubstitutionProcess() {}
	virtual ~IIDDirichletIIDGammaSubstitutionProcess() {}

	protected:

	virtual void UpdateZip(int site)	{
		PoissonSubstitutionProcess::UpdateZip(site);
	}

	virtual void Create()	{
		PoissonSubstitutionProcess::Create();
		GammaRateProcess::Create();
		PoissonSiteSpecificProfileProcess::Create();
	}

	virtual void Delete()	{
		PoissonSiteSpecificProfileProcess::Delete();
		GammaRateProcess::Delete();
		PoissonSubstitutionProcess::Delete();
	}

};

class IIDDirichletIIDGammaPhyloProcess : public virtual PoissonPhyloProcess, public virtual IIDDirichletIIDGammaSubstitutionProcess, public virtual GammaBranchProcess	{

	public:

    virtual void SlaveExecute(MESSAGE);
	virtual void GlobalUpdateParameters();
	virtual void SlaveUpdateParameters();

	IIDDirichletIIDGammaPhyloProcess() {}

	IIDDirichletIIDGammaPhyloProcess(istream& is, int inmyid, int innprocs)	{

		// generic
		FromStreamHeader(is);
		SetMPI(inmyid,innprocs);
		Open(is);
	}

	void ToStreamHeader(ostream& os)	{
		PhyloProcess::ToStreamHeader(os);
	}

	~IIDDirichletIIDGammaPhyloProcess() {
		Delete();
	}

	void TraceHeader(ostream& os)	{
		os << "#iter\ttime\ttopo\tloglik\tlength\trate\talpha";
		os << "\tstatent\tstatalpha";
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

		os << '\t' << GetLogLikelihood() << '\t' << GetRenormTotalLength() << '\t' << GetMeanRate() << '\t' << GetAlpha();
		os << '\t' << GetStatEnt();
		os << '\t' << GetMeanDirWeight();
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

		for (int rep=0; rep<5; rep++)	{
            GlobalCollapse();

            GammaBranchProcess::Move(tuning,10);
            GlobalUpdateParameters();

            GammaRateProcess::Move();
            GlobalUpdateParameters();

            PoissonSiteSpecificProfileProcess::Move();
            GlobalUpdateParameters();

            GlobalUnfold();
        }

        GlobalCollectSiteRates();
        GlobalCollectSiteProfiles();

        chronototal.Stop();
		return 1;
	
	}

    virtual void VarBayes();
    double GetVarLogMarginalLikelihood();
    void UpdateVarLengths();
    void UpdateVarRates();
    void UpdateVarProfiles();

    void UpdateSite(int i)  {
        UpdateZip(i);
    }

	// virtual void ReadPB(int argc, char* argv[]);

	void ToStream(ostream& os)	{
		GammaBranchProcess::ToStream(os);
		GammaRateProcess::ToStream(os);
		PoissonSiteSpecificProfileProcess::ToStream(os);
	}

	void FromStream(istream& is)	{
		GammaBranchProcess::FromStream(is);
		GammaRateProcess::FromStream(is);
		PoissonSiteSpecificProfileProcess::FromStream(is);
		GlobalUpdateParameters();
	}


	virtual void Create()	{
		PoissonPhyloProcess::Create();
		IIDDirichletIIDGammaSubstitutionProcess::Create();
		GammaBranchProcess::Create();
	}
		
	virtual void Delete()	{
		GammaBranchProcess::Delete();
		IIDDirichletIIDGammaSubstitutionProcess::Delete();
		PoissonPhyloProcess::Delete();
	}
};

#endif

