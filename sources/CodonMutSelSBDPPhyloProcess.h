
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef CODONMUTSELSBDPPHYLO_H
#define CODONMUTSELSBDPPHYLO_H

#include "CodonMutSelSBDPSubstitutionProcess.h"
#include "GeneralPathSuffStatMatrixPhyloProcess.h"
#include "GammaBranchProcess.h"

class CodonMutSelSBDPPhyloProcess : public virtual CodonMutSelSBDPSubstitutionProcess, public virtual GeneralPathSuffStatMatrixPhyloProcess, public virtual GammaBranchProcess	{

	public:

	CodonMutSelSBDPPhyloProcess(int inkappaprior)	{

		kappaprior = inkappaprior;

	}

	CodonMutSelSBDPPhyloProcess(istream& is, int inmyid, int innprocs)	{

		// generic
		FromStreamHeader(is);
		SetMPI(inmyid,innprocs);

		// specific
		is >> kappaprior;

		Open(is);

	}

	virtual ~CodonMutSelSBDPPhyloProcess()	{
		Delete();
	}

	void ToStreamHeader(ostream& os)	{
		PhyloProcess::ToStreamHeader(os);
		os << kappaprior << '\n';
	}

        virtual void SlaveExecute(MESSAGE);
	void SlaveUpdateParameters();
	void GlobalUpdateParameters();

	void TraceHeader(ostream& os)	{
		os << "#time\ttimeperccyle\tpruning\tlnL\tlength\talpha\tNmode\tstatent\tstatalpha\tnucsA\tnucsC\tnucsG\tnucsT\tnucrrAC\tnucrrAG\tnucrrAT\tnucrrCG\tnucrrCT\tnucrrGT";
		os << '\n'; 
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
		os << '\t' << GetAlpha();
		os << '\t' << GetNDisplayedComponent();
		os << '\t' << GetStatEnt();
		os << '\t' << GetMeanDirWeight();
		os << '\t' << GetNucStat(0) << '\t' << GetNucStat(1) << '\t' << GetNucStat(2) << '\t' << GetNucStat(3);
		os << '\t' << GetNucRR(0) << '\t' << GetNucRR(1) << '\t' << GetNucRR(2) << '\t' << GetNucRR(3) << '\t' << GetNucRR(4) << '\t' << GetNucRR(5);
		os << '\n';
	}

	void ToStream(ostream& os)	{
		GammaBranchProcess::ToStream(os);
		CodonMutSelSBDPProfileProcess::ToStream(os);
	}
	void FromStream(istream& is)	{
		GammaBranchProcess::FromStream(is);
		CodonMutSelSBDPProfileProcess::FromStream(is);
		GlobalUpdateParameters();
	}

	virtual void ReadPB(int argc, char* argv[]);
	void Read(string name, int burnin, int every, int until);

	double Move(double tuning = 1.0)	{
		chronototal.Start();
		propchrono.Start();
		if (! fixbl)	{
			BranchLengthMove(0.1 * tuning);
			BranchLengthMove(tuning);
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
		CodonMutSelSBDPProfileProcess::Move(tuning,1,10);
		chronosuffstat.Stop();

		chronounfold.Start();
		GlobalUnfold();
		chronounfold.Stop();

		chronototal.Stop();
		return 1;
	}

	protected:

	virtual void Create()	{
		CodonMutSelSBDPSubstitutionProcess::Create();
		GeneralPathSuffStatMatrixPhyloProcess::Create();
		GammaBranchProcess::Create();
	}

	virtual void Delete()	{
		GeneralPathSuffStatMatrixPhyloProcess::Delete();
		CodonMutSelSBDPSubstitutionProcess::Delete();
		GammaBranchProcess::Delete();
	}

	Chrono chronopruning;
	Chrono chronosuffstat;
	Chrono chronototal;
	Chrono chronocollapse;
	Chrono chronounfold;

};

#endif

