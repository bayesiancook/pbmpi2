
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef RASCATGTRSBDP_H
#define RASCATGTRSBDP_H

#include "GTRSBDPProfileProcess.h"
#include "ExpoConjugateGTRMixtureProfileProcess.h"
#include "ExpoConjugateGTRPhyloProcess.h"
#include "DGamRateProcess.h"
#include "GammaBranchProcess.h"
#include "CodonSequenceAlignment.h"


// Exponential conjugate GTR models
class ExpoConjugateGTRSBDPProfileProcess : public virtual GTRSBDPProfileProcess, public virtual ExpoConjugateGTRMixtureProfileProcess {

	public:

	ExpoConjugateGTRSBDPProfileProcess() {}
	virtual ~ExpoConjugateGTRSBDPProfileProcess() {}

	protected:

	virtual void SwapComponents(int cat1, int cat2)	{
		GTRSBDPProfileProcess::SwapComponents(cat1,cat2);
	}

	virtual void Create()	{
		GTRSBDPProfileProcess::Create();
		ExpoConjugateGTRMixtureProfileProcess::Create();
	}

	virtual void Delete()	{
		ExpoConjugateGTRMixtureProfileProcess::Delete();
		GTRSBDPProfileProcess::Delete();
	}
};

class RASCATGTRSBDPSubstitutionProcess : public virtual ExpoConjugateGTRSubstitutionProcess, public virtual DGamRateProcess, public virtual ExpoConjugateGTRSBDPProfileProcess {

	public:

	RASCATGTRSBDPSubstitutionProcess() {}
	virtual ~RASCATGTRSBDPSubstitutionProcess() {}

	virtual void Create()	{
		ExpoConjugateGTRSubstitutionProcess::Create();
		DGamRateProcess::Create();
		ExpoConjugateGTRSBDPProfileProcess::Create();
	}

	virtual void Delete()	{
		ExpoConjugateGTRSBDPProfileProcess::Delete();
		DGamRateProcess::Delete();
		ExpoConjugateGTRSubstitutionProcess::Delete();
	}
};


class RASCATGTRSBDPGammaPhyloProcess : public virtual ExpoConjugateGTRPhyloProcess, public virtual RASCATGTRSBDPSubstitutionProcess, public virtual GammaBranchProcess	{

	public:

        virtual void SlaveExecute(MESSAGE);
	virtual void GlobalUpdateParameters();
	virtual void SlaveUpdateParameters();

	RASCATGTRSBDPGammaPhyloProcess() {}

	RASCATGTRSBDPGammaPhyloProcess(int nratecat, string inrrtype, int inkappaprior)	{

		Ncat = nratecat;
		rrtype = inrrtype;
		kappaprior = inkappaprior;
	}

	RASCATGTRSBDPGammaPhyloProcess(istream& is, int inmyid, int innprocs)	{

		// generic
		FromStreamHeader(is);
		SetMPI(inmyid,innprocs);

		// specific
		is >> Ncat;
		is >> kappaprior;
		is >> rrtype;

		Open(is);

	}

	~RASCATGTRSBDPGammaPhyloProcess() {
		Delete();
	}

	void ToStreamHeader(ostream& os)	{
		PhyloProcess::ToStreamHeader(os);
		os << Ncat << '\n';
		os << kappaprior << '\n';
		os << rrtype << '\n';
	}

	void TraceHeader(ostream& os)	{
		os << "#iter\ttime\ttopo\tloglik\tlength\talpha\tNmode\tstatent\tstatalpha\tkappa";
		if (! fixrr)	{
			os << "\trrent\trrmean";
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
		}
		else	{
			os << '\t' << 0;
			os << '\t' << 0;
		}

		os << '\t' << GetLogLikelihood();
		os << '\t' << GetRenormTotalLength();
		os << '\t' << GetAlpha();
		os << '\t' << GetNDisplayedComponent();
		os << '\t' << GetStatEnt();
		os << '\t' << GetMeanDirWeight();
		os << '\t' << kappa;
		if (! fixrr)	{
			os << '\t' << GetRREntropy();
			os << '\t' << GetRRMean();
		}
		os << '\n';

	}

	virtual void Monitor(ostream& os)  {
		PhyloProcess::Monitor(os);
		os << "weight " << '\t' << GetMaxWeightError() << '\n';
		ResetMaxWeightError();
	}

	/*
	double MoveProfiles()	{
		ExpoConjugateGTRSBDPProfileProcess::Move(1,1,10);
	}
	*/

	double Move(double tuning = 1.0)	{

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

		ExpoConjugateGTRSBDPProfileProcess::Move(1,1,10);
		if (iscodon){
			ExpoConjugateGTRSBDPProfileProcess::Move(0.1,1,15);
			ExpoConjugateGTRSBDPProfileProcess::Move(0.01,1,15);
		}

		if (! fixrr){
			LengthRelRateMove(1,10);
			LengthRelRateMove(0.1,10);
			LengthRelRateMove(0.01,10);
		}

		GlobalUnfold();

		chronototal.Stop();

		return 1;
	}

	virtual double GlobalRestrictedMoveCycle(int nrep = 1, double tuning = 1.0)	{

		for (int rep=0; rep<nrep; rep++)	{

			GlobalUpdateParameters();
			GammaBranchProcess::Move(tuning,10);

			GlobalUpdateParameters();
			DGamRateProcess::Move(0.3*tuning,10);
			DGamRateProcess::Move(0.03*tuning,10);
			ExpoConjugateGTRSBDPProfileProcess::Move(1,1,1);
		}
		return 1;
	}

	virtual void PrepareSiteLogLikelihood(int site) {
		int cat = ExpoConjugateGTRSBDPProfileProcess::alloc[site];
		if (! matrixarray[cat])	{
			CreateMatrix(cat);
		}
		UpdateMatrix(cat);
	}

	void ToStream(ostream& os)	{
		GammaBranchProcess::ToStream(os);
		DGamRateProcess::ToStream(os);
		ExpoConjugateGTRSBDPProfileProcess::ToStream(os);
	}

	void FromStream(istream& is)	{
		GammaBranchProcess::FromStream(is);
		DGamRateProcess::FromStream(is);
		ExpoConjugateGTRSBDPProfileProcess::FromStream(is);
		GlobalUpdateParameters();
	}

	virtual void ReadPB(int argc, char* argv[]);
	void ReadNocc(string name, int burnin, int every, int until);
	void ReadTestProfile(string name, int nrep, double tuning, int burnin, int every, int until);
	void ReadRelRates(string name, int burnin, int every, int until);
	void SlaveComputeCVScore();
	void SlaveComputeSiteLogL();

	protected:

	virtual void Create()	{
		ExpoConjugateGTRPhyloProcess::Create();
		RASCATGTRSBDPSubstitutionProcess::Create();
		GammaBranchProcess::Create();
	}
		
	virtual void Delete()	{
		GammaBranchProcess::Delete();
		RASCATGTRSBDPSubstitutionProcess::Delete();
		ExpoConjugateGTRPhyloProcess::Delete();
	}

};

#endif

