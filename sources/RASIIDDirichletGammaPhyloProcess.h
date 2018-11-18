
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef IIDDIR_H
#define IIDDIR_H

#include "PoissonSubstitutionProcess.h"
#include "PoissonPhyloProcess.h"
#include "DGamRateProcess.h"
#include "PoissonSiteSpecificProfileProcess.h"
#include "GammaBranchProcess.h"

class RASIIDDirichletSubstitutionProcess : public virtual PoissonSubstitutionProcess, public virtual DGamRateProcess, public virtual PoissonSiteSpecificProfileProcess {

	using PoissonSubstitutionProcess::UpdateZip;

	public:

	RASIIDDirichletSubstitutionProcess() {}
	virtual ~RASIIDDirichletSubstitutionProcess() {}

	protected:

	virtual void UpdateZip(int site)	{
		PoissonSubstitutionProcess::UpdateZip(site);
	}

	virtual void Create()	{
		PoissonSubstitutionProcess::Create();
		DGamRateProcess::Create();
		PoissonSiteSpecificProfileProcess::Create();
	}

	virtual void Delete()	{
		PoissonSiteSpecificProfileProcess::Delete();
		DGamRateProcess::Delete();
		PoissonSubstitutionProcess::Delete();
	}

};

class RASIIDDirichletGammaPhyloProcess : public virtual PoissonPhyloProcess, public virtual RASIIDDirichletSubstitutionProcess, public virtual GammaBranchProcess	{

	public:

    virtual void SlaveExecute(MESSAGE);
	virtual void GlobalUpdateParameters();
	virtual void SlaveUpdateParameters();

	RASIIDDirichletGammaPhyloProcess() {}

	RASIIDDirichletGammaPhyloProcess(int inNcat, string insitefreq, int infixprofile, double inpseudocount) {
        Ncat = inNcat;
        sitefreq = insitefreq;
        fixprofile = infixprofile;
        pseudocount = inpseudocount;
        initprofile = 0;
    }

	RASIIDDirichletGammaPhyloProcess(istream& is, int inmyid, int innprocs)	{

		// generic
		FromStreamHeader(is);
		SetMPI(inmyid,innprocs);
        is >> Ncat;
        is >> sitefreq >> fixprofile >> pseudocount;
		Open(is);
	}

	void ToStreamHeader(ostream& os)	{
		PhyloProcess::ToStreamHeader(os);
        os << Ncat << '\n';
        os << sitefreq << '\t' << fixprofile << '\t' << pseudocount << '\n';
	}

	~RASIIDDirichletGammaPhyloProcess() {
		Delete();
	}

	void TraceHeader(ostream& os)	{
		os << "#iter\ttime\ttopo\tloglik\tlength\trate\talpha";
        if (! fixprofile)   {
            os << "\tstatent\tstatalpha";
        }
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

		os << '\t' << GetLogLikelihood() << '\t' << GetTotalLength() << '\t' << GetAlpha();
        if (! fixprofile)   {
            os << '\t' << GetStatEnt();
            os << '\t' << GetMeanDirWeight();
        }
		os << '\n';
	}

	virtual double Move(double tuning = 1.0)	{

		chronototal.Start();
		propchrono.Start();
		// BranchLengthMove(tuning);
		// BranchLengthMove(0.1 * tuning);
		if (! fixtopo)	{
			MoveTopo();
		}
		propchrono.Stop();

		for (int rep=0; rep<5; rep++)	{
            GlobalCollapse();

            GammaBranchProcess::Move(tuning,10);
            GlobalUpdateParameters();

            DGamRateProcess::Move(tuning,10);
            DGamRateProcess::Move(0.3*tuning,10);
            DGamRateProcess::Move(0.03*tuning,10);
            GlobalUpdateParameters();

            if (! fixprofile)   {
                PoissonSiteSpecificProfileProcess::Move();
                GlobalUpdateParameters();
            }

            GlobalUnfold();
        }

        if (! fixprofile)   {
            GlobalCollectSiteProfiles();
        }

        chronototal.Stop();
		return 1;
	
	}

    void UpdateSite(int i)  {
        UpdateZip(i);
    }

	// virtual void ReadPB(int argc, char* argv[]);

	void ToStream(ostream& os)	{
		GammaBranchProcess::ToStream(os);
		DGamRateProcess::ToStream(os);
        // if (! fixprofile)   {
            PoissonSiteSpecificProfileProcess::ToStream(os);
        // }
	}

	void FromStream(istream& is)	{
		GammaBranchProcess::FromStream(is);
		DGamRateProcess::FromStream(is);
        // if (! fixprofile)   {
            PoissonSiteSpecificProfileProcess::FromStream(is);
        // }
		GlobalUpdateParameters();
	}


	virtual void Create()	{
		PoissonPhyloProcess::Create();
		RASIIDDirichletSubstitutionProcess::Create();
		GammaBranchProcess::Create();
        if (sitefreq != "free") {
            SetSiteProfiles();
        }
    }

	virtual void Delete()	{
		GammaBranchProcess::Delete();
		RASIIDDirichletSubstitutionProcess::Delete();
		PoissonPhyloProcess::Delete();
	}

    void EM(double cutoff, int nrep);
    double EMUpdateMeanSuffStat();
    void EM_UpdateProfiles(double pseudocount);
	double** modesitelogL;
	double** modesitepostprob;
	int initprofile;
};

#endif

