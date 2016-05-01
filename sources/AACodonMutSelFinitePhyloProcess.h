
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef AACODONMUTSELFINITEPHYLO_H
#define AACODONMUTSELFINITEPHYLO_H

#include "AACodonMutSelFiniteSubstitutionProcess.h"
#include "GeneralPathSuffStatMatrixPhyloProcess.h"
#include "GammaBranchProcess.h"

class AACodonMutSelFinitePhyloProcess : public virtual AACodonMutSelFiniteSubstitutionProcess, public virtual GeneralPathSuffStatMatrixPhyloProcess, public virtual GammaBranchProcess	{

	public:

	AACodonMutSelFinitePhyloProcess(int ncat, int infixncomp, int inempmix, string inmixtype, int infixcodonprofile, int infixomega, int inomegaprior)	{

		fixcodonprofile = infixcodonprofile;
		fixomega = infixomega;
		omegaprior = inomegaprior;

		Ncomponent = ncat;
		fixncomp = infixncomp;
		empmix = inempmix;
		mixtype = inmixtype;
	}

	AACodonMutSelFinitePhyloProcess(istream& is, int inmyid, int innprocs)	{

		// generic
		FromStreamHeader(is);
		SetMPI(inmyid,innprocs);

		// specific
		is >> fixcodonprofile;
		is >> fixomega;
		is >> omegaprior;
		is >> fixncomp;
		is >> empmix;
		is >> mixtype;

		Open(is);
	}

	virtual ~AACodonMutSelFinitePhyloProcess()	{
		Delete();
	}

	void ToStreamHeader(ostream& os)	{
		PhyloProcess::ToStreamHeader(os);
		os << fixcodonprofile << '\n';
		os << fixomega << '\n';
		os << omegaprior << '\n';
		os << fixncomp << '\n';
		os << empmix << '\n';
		os << mixtype << '\n';
	}

        virtual void SlaveExecute(MESSAGE);
	void SlaveComputeCVScore();
	void SlaveUpdateParameters();
	void GlobalUpdateParameters();

	void TraceHeader(ostream& os)	{
		os << "#iter\ttime\tpruning\tlnL\tlength\tcodonent\tomega\tNmode\tstatent\tstatalpha\tnucsA\tnucsC\tnucsG\tnucsT\tnucrrAC\tnucrrAG\tnucrrAT\tnucrrCG\tnucrrCT\tnucrrGT";
		os << "\n";
	}

	void Trace(ostream& os)	{
		UpdateOccupancyNumbers();

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
		os << '\t' << GetNDisplayedComponent();
		os << '\t' << GetStatEnt();
		os << '\t' << GetMeanDirWeight();
		os << '\t' << GetNucStat(0) << '\t' << GetNucStat(1) << '\t' << GetNucStat(2) << '\t' << GetNucStat(3);
		os << '\t' << GetNucRR(0) << '\t' << GetNucRR(1) << '\t' << GetNucRR(2) << '\t' << GetNucRR(3) << '\t' << GetNucRR(4) << '\t' << GetNucRR(5);

		os << '\n';
	}

	void ToStream(ostream& os)	{
		GammaBranchProcess::ToStream(os);
		AACodonMutSelFiniteProfileProcess::ToStream(os);
	}

	void FromStream(istream& is)	{
		GammaBranchProcess::FromStream(is);
		AACodonMutSelFiniteProfileProcess::FromStream(is);
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
		AACodonMutSelFiniteProfileProcess::Move(tuning,1,10);
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
		AACodonMutSelFiniteSubstitutionProcess::Create();
		GeneralPathSuffStatMatrixPhyloProcess::Create();
		GammaBranchProcess::Create();
	}

	virtual void Delete()	{
		GeneralPathSuffStatMatrixPhyloProcess::Delete();
		AACodonMutSelFiniteSubstitutionProcess::Delete();
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

