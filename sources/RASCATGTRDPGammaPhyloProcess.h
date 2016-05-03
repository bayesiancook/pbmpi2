
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef RASCATGTRDP_H
#define RASCATGTRDP_H

#include "GTRDPProfileProcess.h"
#include "ExpoConjugateGTRMixtureProfileProcess.h"
#include "ExpoConjugateGTRPhyloProcess.h"
#include "DGamRateProcess.h"
#include "GammaBranchProcess.h"

// Exponential conjugate GTR models
class ExpoConjugateGTRDPProfileProcess : public virtual GTRDPProfileProcess, public virtual ExpoConjugateGTRMixtureProfileProcess {

	public:

	ExpoConjugateGTRDPProfileProcess() {}
	virtual ~ExpoConjugateGTRDPProfileProcess() {}

	protected:

	virtual void Create()	{
		GTRDPProfileProcess::Create();
		ExpoConjugateGTRMixtureProfileProcess::Create();
	}

	virtual void Delete()	{
		ExpoConjugateGTRMixtureProfileProcess::Delete();
		GTRDPProfileProcess::Delete();
	}
};
// this is the final class implementing the CATGTR phyloprocess

class RASCATGTRDPSubstitutionProcess : public virtual ExpoConjugateGTRSubstitutionProcess, public virtual DGamRateProcess, public virtual ExpoConjugateGTRDPProfileProcess {

	public:

	RASCATGTRDPSubstitutionProcess() {}
	virtual ~RASCATGTRDPSubstitutionProcess() {}

	virtual void Create()	{
		ExpoConjugateGTRSubstitutionProcess::Create();
		DGamRateProcess::Create();
		ExpoConjugateGTRDPProfileProcess::Create();
	}

	virtual void Delete()	{
		ExpoConjugateGTRDPProfileProcess::Delete();
		DGamRateProcess::Delete();
		ExpoConjugateGTRSubstitutionProcess::Delete();
	}
};


class RASCATGTRDPGammaPhyloProcess : public virtual ExpoConjugateGTRPhyloProcess, public virtual RASCATGTRDPSubstitutionProcess, public virtual GammaBranchProcess	{

	public:

        virtual void SlaveExecute(MESSAGE);
	virtual void GlobalUpdateParameters();
	virtual void SlaveUpdateParameters();


	RASCATGTRDPGammaPhyloProcess(int nratecat, string inrrtype, int inkappaprior)	{

		Ncat = nratecat;
		kappaprior = inkappaprior;
		rrtype = inrrtype;
	}

	RASCATGTRDPGammaPhyloProcess(istream& is, int inmyid, int innprocs)	{

		// generic
		FromStreamHeader(is);

		// specific
		is >> Ncat;
		is >> kappaprior;
		is >> rrtype;

		Open(is);
	}

	void ToStreamHeader(ostream& os)	{
		PhyloProcess::ToStreamHeader(os);
		os << Ncat << '\n';
		os << kappaprior << '\n';
		os << rrtype << '\n';
	}

	~RASCATGTRDPGammaPhyloProcess() {
		Delete();
	}

	void TraceHeader(ostream& os)	{
		os << "#time\ttime\ttopo\tloglik\tlength\talpha\tNmode\tstatent\tstatalpha";
		if (! fixrr)	{
			os << "\trrent\trrmean";
		}
		os << "\tkappa\tallocent";
		os << '\n'; 
	}

	void Trace(ostream& os)	{

		UpdateOccupancyNumbers();

		os << ((int) (chronototal.GetTime() / 1000));
		if (chronototal.GetTime())	{
			os << '\t' << ((double) ((int) (chronototal.GetTime() / (1 + GetSize())))) / 1000;
			os << '\t' << ((int) (propchrono.GetTime() / chronototal.GetTime() * 100));
			// os << '\t' << ((int) (chronosuffstat.GetTime() / chronototal.GetTime() * 100));
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
		if (! fixrr)	{
			os << '\t' << GetRREntropy();
			os << '\t' << GetRRMean();
		}
		// os << '\t' << kappa << '\t' << GetAllocEntropy();
		os << '\n';
	}

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

		ExpoConjugateGTRDPProfileProcess::Move(1,1,10);

		if (! fixrr)	{
			LengthRelRateMove(1,10);
			LengthRelRateMove(0.1,10);
			LengthRelRateMove(0.01,10);
		}

		GlobalUnfold();

		chronototal.Stop();

		return 1;
	}

	void ToStream(ostream& os)	{
		GammaBranchProcess::ToStream(os);
		DGamRateProcess::ToStream(os);
		ExpoConjugateGTRDPProfileProcess::ToStream(os);
	}

	void FromStream(istream& is)	{
		GammaBranchProcess::FromStream(is);
		DGamRateProcess::FromStream(is);
		ExpoConjugateGTRDPProfileProcess::FromStream(is);
		GlobalUpdateParameters();
	}

	protected:

	virtual void Create()	{
		ExpoConjugateGTRPhyloProcess::Create();
		RASCATGTRDPSubstitutionProcess::Create();
		GammaBranchProcess::Create();
	}
		
	virtual void Delete()	{
		GammaBranchProcess::Delete();
		RASCATGTRDPSubstitutionProcess::Delete();
		ExpoConjugateGTRPhyloProcess::Delete();
	}

	int fixrr;
};

#endif

