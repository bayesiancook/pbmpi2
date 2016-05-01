
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef AACODONMUTSELMVNSSPHYLO_H
#define AACODONMUTSELMVNSSPHYLO_H

#include "AACodonMutSelMVNSiteSpecificSubstitutionProcess.h"
#include "GeneralPathSuffStatMatrixPhyloProcess.h"
#include "GammaBranchProcess.h"

class AACodonMutSelMVNSiteSpecificPhyloProcess : public virtual AACodonMutSelMVNSiteSpecificSubstitutionProcess, public virtual GeneralPathSuffStatMatrixPhyloProcess, public virtual GammaBranchProcess	{

	public:

	AACodonMutSelMVNSiteSpecificPhyloProcess(int infixcodonprofile, int infixomega, int inomegaprior)	{

		fixcodonprofile = infixcodonprofile;
		fixomega = infixomega;
		omegaprior = inomegaprior;
	}

	AACodonMutSelMVNSiteSpecificPhyloProcess(istream& is, int inmyid, int innprocs)	{

		// generic
		FromStreamHeader(is);
		SetMPI(inmyid,innprocs);

		// specific
		is >> fixcodonprofile;
		is >> fixomega;
		is >> omegaprior;

		Open(is);
	}

	virtual ~AACodonMutSelMVNSiteSpecificPhyloProcess()	{
		Delete();
	}

	void ToStreamHeader(ostream& os)	{
		PhyloProcess::ToStreamHeader(os);
		os << fixcodonprofile << '\n';
		os << fixomega << '\n';
		os << omegaprior << '\n';
	}

        virtual void SlaveExecute(MESSAGE);
	void SlaveUpdateParameters();
	void GlobalUpdateParameters();

	void TraceHeader(ostream& os)	{
		os << "#iter\ttime\tpruning\tlnL\tlength\tcodonent\tomega\tstatent\tnucsA\tnucsC\tnucsG\tnucsT\n";
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

		os << '\t' <<  GetLogLikelihood();
		os << '\t' << GetTotalLength();
		os << '\t' << GetCodonProfileEntropy();
		os << '\t' << GetOmega();
		os << '\t' << GetStatEnt();
		os << '\t' << GetNucStat(0) << '\t' << GetNucStat(1) << '\t' << GetNucStat(2) << '\t' << GetNucStat(3);
		os << '\n';
	}

	void ToStream(ostream& os)	{
		GammaBranchProcess::ToStream(os);
		AACodonMutSelMVNSiteSpecificProfileProcess::ToStream(os);
	}

	void FromStream(istream& is)	{
		GammaBranchProcess::FromStream(is);
		AACodonMutSelMVNSiteSpecificProfileProcess::FromStream(is);
		GlobalUpdateParameters();
	}

	virtual void ReadPB(int argc, char* argv[]);

	double Move(double tuning = 1.0)	{
		chronototal.Start();
		propchrono.Start();
		if (! fixbl)	{
			BranchLengthMove(tuning);
			BranchLengthMove(0.1 * tuning);
		}
		if (! fixtopo)	{
			MoveTopo();
		}
		propchrono.Stop();

		chronosuffstat.Start();

		chronocollapse.Start();
		GlobalCollapse();
		chronocollapse.Stop();
		if (! fixbl)	{
			GammaBranchProcess::Move(0.1 * tuning,10);
			GammaBranchProcess::Move(tuning,10);
		}

		GlobalUpdateParameters();
		AACodonMutSelMVNSiteSpecificProfileProcess::Move(tuning,1,10);
		chronosuffstat.Stop();

		chronounfold.Start();
		GlobalUnfold();
		chronounfold.Stop();

		chronototal.Stop();
		return 1;
	}


	protected:

	virtual void SetProfileDim()	{
		SetDim(20);
	}

	virtual void Create()	{
		AACodonMutSelMVNSiteSpecificSubstitutionProcess::Create();
		GeneralPathSuffStatMatrixPhyloProcess::Create();
		GammaBranchProcess::Create();
	}

	virtual void Delete()	{
		GeneralPathSuffStatMatrixPhyloProcess::Delete();
		AACodonMutSelMVNSiteSpecificSubstitutionProcess::Delete();
		GammaBranchProcess::Delete();
	}

	int fixcodonprofile;
	int fixomega;

	Chrono chronopruning;
	Chrono chronosuffstat;
	Chrono chronototal;
	Chrono chronocollapse;
	Chrono chronounfold;

};

#endif

