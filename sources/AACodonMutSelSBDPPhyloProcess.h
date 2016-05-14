
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef AACODONMUTSELSBDPPHYLO_H
#define AACODONMUTSELSBDPPHYLO_H

#include "AACodonMutSelSBDPSubstitutionProcess.h"
#include "GeneralPathSuffStatMatrixPhyloProcess.h"
#include "GammaBranchProcess.h"

class AACodonMutSelSBDPPhyloProcess : public virtual AACodonMutSelSBDPSubstitutionProcess, public virtual GeneralPathSuffStatMatrixPhyloProcess, public virtual GammaBranchProcess	{

	public:

	AACodonMutSelSBDPPhyloProcess(int infixcodonprofile, int infixomega, int inomegaprior, int inkappaprior)	{

		fixcodonprofile = infixcodonprofile;
		fixomega = infixomega;
		omegaprior = inomegaprior;

		kappaprior = inkappaprior;
	}

	AACodonMutSelSBDPPhyloProcess(istream& is, int inmyid, int innprocs)	{

		// generic
		FromStreamHeader(is);
		SetMPI(inmyid,innprocs);

		// specific
		is >> fixcodonprofile;
		is >> fixomega;
		is >> omegaprior;
		is >> kappaprior;

		Open(is);
	}

	virtual ~AACodonMutSelSBDPPhyloProcess()	{
		Delete();
	}

	void ToStreamHeader(ostream& os)	{
		PhyloProcess::ToStreamHeader(os);
		os << fixcodonprofile << '\n';
		os << fixomega << '\n';
		os << omegaprior << '\n';
		os << kappaprior << '\n';
	}

        virtual void SlaveExecute(MESSAGE);
	void SlaveUpdateParameters();
	void GlobalUpdateParameters();
	
	void SlaveComputeCVScore();

	void TraceHeader(ostream& os)	{
		os << "#iter\ttime\tpruning\tlnL\tlength\tcodonent\tomega\tNmode\tstatent\tstatalpha\tnucsA\tnucsC\tnucsG\tnucsT\tnucrrAC\tnucrrAG\tnucrrAT\tnucrrCG\tnucrrCT\tnucrrGT";
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
			//os << '\t' << ((double) ((int) (chronototal.GetTime() / (GetSize())))) / 1000;
			//os << '\t' << ((int) (100 * chronopruning.GetTime() /chronototal.GetTime()));
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

	virtual void Monitor(ostream& os)  {
		PhyloProcess::Monitor(os);
		os << "weight " << '\t' << GetMaxWeightError() << '\n';
		ResetMaxWeightError();
	}

	void ToStream(ostream& os)	{
		GammaBranchProcess::ToStream(os);
		AACodonMutSelSBDPProfileProcess::ToStream(os);
	}
	void FromStream(istream& is)	{
		GammaBranchProcess::FromStream(is);
		AACodonMutSelSBDPProfileProcess::FromStream(is);
		GlobalUpdateParameters();
	}

	virtual void ReadPB(int argc, char* argv[]);
	void Read(string name, int burnin, int every, int until);
	void ReadMapStats(string name, int burnin, int every, int until);
	int CountNonSynMapping(int i);
	int CountNonSynMapping();
	int GlobalNonSynMapping();
	virtual void SlaveNonSynMapping();

	// primary scheduler

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
		AACodonMutSelSBDPProfileProcess::Move(tuning,1,15);
		chronosuffstat.Stop();

		chronounfold.Start();
		GlobalUnfold();
		chronounfold.Stop();

		chronototal.Stop();
		return 1;
	}


	protected:

	virtual void Create()	{
		AACodonMutSelSBDPSubstitutionProcess::Create();
		GeneralPathSuffStatMatrixPhyloProcess::Create();
		GammaBranchProcess::Create();
	}

	virtual void Delete()	{
		GammaBranchProcess::Delete();
		GeneralPathSuffStatMatrixPhyloProcess::Delete();
		AACodonMutSelSBDPSubstitutionProcess::Delete();
	}

	virtual void SetProfileDim()	{
		SetDim(20);
	}

	Chrono chronopruning;
	Chrono chronosuffstat;
	Chrono chronototal;
	Chrono chronocollapse;
	Chrono chronounfold;

};



#endif

