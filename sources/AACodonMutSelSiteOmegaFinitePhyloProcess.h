
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef AACODONMUTSELSITEOMFINITEPHYLO_H
#define AACODONMUTSELSITEOMFINITEPHYLO_H

#include "AACodonMutSelSiteOmegaFiniteSubstitutionProcess.h"
#include "GeneralPathSuffStatMatrixPhyloProcess.h"
#include "GammaBranchProcess.h"

class AACodonMutSelSiteOmegaFinitePhyloProcess : public virtual AACodonMutSelSiteOmegaFiniteSubstitutionProcess, public virtual GeneralPathSuffStatMatrixPhyloProcess, public virtual GammaBranchProcess	{

	public:

	AACodonMutSelSiteOmegaFinitePhyloProcess(int ncat, int infixncomp, int inempmix, string inmixtype, int infixcodonprofile, int infixomega, int inomegaprior)	{

		fixcodonprofile = infixcodonprofile;
		fixomega = infixomega;
		omegaprior = inomegaprior;

		Ncomponent = ncat;
		fixncomp = infixncomp;
		empmix = inempmix;
		mixtype = inmixtype;
	}

	AACodonMutSelSiteOmegaFinitePhyloProcess(istream& is, int inmyid, int innprocs)	{

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

	virtual ~AACodonMutSelSiteOmegaFinitePhyloProcess()	{
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
	void SlaveUpdateParameters();
	void GlobalUpdateParameters();
	
	void SlaveComputeCVScore();

	void TraceHeader(ostream& os)	{
		os << "#iter\ttime\tpruning\tlnL\tlength\tcodonent\tmeanomega\tvaromega\tpos\tNmode\tstatent\tstatalpha\tnucsA\tnucsC\tnucsG\tnucsT\tnucrrAC\tnucrrAG\tnucrrAT\tnucrrCG\tnucrrCT\tnucrrGT";
		for (int i=0; i<GetNsite(); i++)	{
			os << "\t" << i+1;
		}
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
		os << '\t' << GetMeanOmega();
		os << '\t' << GetRelVarOmega();
		os << '\t' << GetProportionOmegaGreaterThan(1.0);
		os << '\t' << GetNDisplayedComponent();
		os << '\t' << GetStatEnt();
		os << '\t' << GetMeanDirWeight();
		os << '\t' << GetNucStat(0) << '\t' << GetNucStat(1) << '\t' << GetNucStat(2) << '\t' << GetNucStat(3);
		os << '\t' << GetNucRR(0) << '\t' << GetNucRR(1) << '\t' << GetNucRR(2) << '\t' << GetNucRR(3) << '\t' << GetNucRR(4) << '\t' << GetNucRR(5);
		for (int i=0; i<GetNsite(); i++)	{
			os << '\t' << GetSiteOmega(i);
		}
		os << '\n';
	}

	void ToStream(ostream& os)	{
		GammaBranchProcess::ToStream(os);
		AACodonMutSelSiteOmegaFiniteProfileProcess::ToStream(os);
	}
	void FromStream(istream& is)	{
		GammaBranchProcess::FromStream(is);
		AACodonMutSelSiteOmegaFiniteProfileProcess::FromStream(is);
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

	/*
	using MixtureProfileProcess::ProfileSuffStatLogProb;
	using MixtureProfileProcess::CountProfileSuffStatLogProb;
	using MixtureProfileProcess::BetaProfileSuffStatLogProb;
	*/

	double Move(double tuning = 1.0)	{
		/*
		GlobalCollapse();
		GlobalUpdateSiteProfileSuffStat();
		GlobalUpdateModeProfileSuffStat();
		GlobalUpdateOmegaSuffStat();
		cerr << "suff stat log prob : " << CountProfileSuffStatLogProb() + CountOmegaSuffStatLogProb() << '\t' << BetaProfileSuffStatLogProb() << '\t' << ProfileSuffStatLogProb() << '\n';
		exit(1);
		*/

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
		AACodonMutSelSiteOmegaFiniteProfileProcess::Move(tuning,1,15);
		chronosuffstat.Stop();

		chronounfold.Start();
		GlobalUnfold();
		chronounfold.Stop();

		chronototal.Stop();
		return 1;
	}

	protected:

	virtual void Create()	{
		AACodonMutSelSiteOmegaFiniteSubstitutionProcess::Create();
		GeneralPathSuffStatMatrixPhyloProcess::Create();
		GammaBranchProcess::Create();
		// if (GetMyid())	{
			CreateSiteMatrices();
		// }
	}

	virtual void Delete()	{
		// if (GetMyid())	{
			DeleteSiteMatrices();
		// }
		GammaBranchProcess::Delete();
		GeneralPathSuffStatMatrixPhyloProcess::Delete();
		AACodonMutSelSiteOmegaFiniteSubstitutionProcess::Delete();
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

