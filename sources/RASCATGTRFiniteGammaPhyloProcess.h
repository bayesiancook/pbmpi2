
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef RASCATGTR_H
#define RASCATGTR_H

#include "GTRFiniteProfileProcess.h"
#include "ExpoConjugateGTRMixtureProfileProcess.h"
#include "ExpoConjugateGTRPhyloProcess.h"
#include "DGamRateProcess.h"
#include "GammaBranchProcess.h"

// this is the final class implementing the CATGTR phyloprocess

class ExpoConjugateGTRFiniteProfileProcess : public virtual GTRFiniteProfileProcess, public virtual ExpoConjugateGTRMixtureProfileProcess {

	public:

	ExpoConjugateGTRFiniteProfileProcess() {}
	virtual ~ExpoConjugateGTRFiniteProfileProcess() {}

	protected:

	virtual void Create()	{
		GTRFiniteProfileProcess::Create();
		ExpoConjugateGTRMixtureProfileProcess::Create();
	}

	virtual void Delete()	{
		ExpoConjugateGTRMixtureProfileProcess::Delete();
		GTRFiniteProfileProcess::Delete();
	}
};

class RASCATGTRFiniteSubstitutionProcess : public virtual ExpoConjugateGTRSubstitutionProcess, public virtual DGamRateProcess, public virtual ExpoConjugateGTRFiniteProfileProcess {

	public:

	RASCATGTRFiniteSubstitutionProcess() {}
	virtual ~RASCATGTRFiniteSubstitutionProcess() {}

	virtual void Create()	{
		ExpoConjugateGTRSubstitutionProcess::Create();
		DGamRateProcess::Create();
		ExpoConjugateGTRFiniteProfileProcess::Create();
	}

	virtual void Delete()	{
		ExpoConjugateGTRFiniteProfileProcess::Delete();
		DGamRateProcess::Delete();
		ExpoConjugateGTRSubstitutionProcess::Delete();
	}
};

class RASCATGTRFiniteGammaPhyloProcess : public virtual ExpoConjugateGTRPhyloProcess, public virtual RASCATGTRFiniteSubstitutionProcess, public virtual GammaBranchProcess	{

	public:

        virtual void SlaveExecute(MESSAGE);

	virtual void GlobalUpdateParameters();

	RASCATGTRFiniteGammaPhyloProcess() {}

	RASCATGTRFiniteGammaPhyloProcess(int nratecat, int ncat, int infixncomp, int inempmix, string inmixtype, string inrrtype)	{

		Ncat = nratecat;
		Ncomponent = ncat;
		fixncomp = infixncomp;
		empmix = inempmix;
		mixtype = inmixtype;
		rrtype = inrrtype;
	}

	RASCATGTRFiniteGammaPhyloProcess(istream& is, int inmyid, int innprocs)	{

		// generic
		FromStreamHeader(is);
		SetMPI(inmyid,innprocs);

		// specific
		is >> Ncat;
		is >> fixncomp;
		is >> empmix;
		is >> mixtype;
		is >> rrtype;

		Open(is);
	}

	void ToStreamHeader(ostream& os)	{
		PhyloProcess::ToStreamHeader(os);
		os << Ncat << '\n';
		os << fixncomp << '\t' << empmix << '\t' << mixtype << '\n';
		os << rrtype << '\n';
	}

	~RASCATGTRFiniteGammaPhyloProcess() {
		Delete();
	}

	virtual void SlaveUpdateParameters();

	void TraceHeader(ostream& os)	{
		os << "#iter\ttime\ttopo\tloglik\tlength\talpha\tNmode\tstatent\tstatalpha";
		if (! fixrr)	{
			os << "\trrent\trrmean";
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

		os << '\t' << GetLogLikelihood();
		os << '\t' << GetRenormTotalLength();
		os << '\t' << GetAlpha();
		os << '\t' << GetNDisplayedComponent();
		os << '\t' << GetStatEnt();
		if (empmix == 2)	{
			os << '\t' << statfixalpha;
			// os << '\t' << weightalpha;
		}
		else	{
			os << '\t' << GetMeanDirWeight();
		}
		if (! fixrr)	{
			os << '\t' << GetRREntropy();
			os << '\t' << GetRRMean();
		}

		os << '\n';
	}

	double Move(double tuning = 1.0)	{

		chronototal.Start();
		propchrono.Start();

		/*
		double x = rnd::GetRandom().Uniform();
		if (x < 0.1)	{
			BranchLengthMove(tuning);
		}
		else if (x < 0.2)	{
			BranchLengthMove(0.1 * tuning);
		}
		*/
		if (! fixtopo)	{
			if (fasttopo)	{
				FastTopoMoveCycle(1,topomu);
				/*
				double x = rnd::GetRandom().Uniform();
				if (x < 0.25)	{
					SPRMove(2);
				}
				else if (x < 0.5)	{
					NNIMove(1,0.1);
				}
				else {
					FastTopoMoveCycle(1,topomu);
				}
				*/
			}
			else	{
				SimpleTopoMoveCycle(1,tuning);
			}
		}
		GlobalUpdateConditionalLikelihoods();

		propchrono.Stop();

		GlobalCollapse();

		GammaBranchProcess::Move(tuning,10);

		GlobalUpdateParameters();
		DGamRateProcess::Move(0.3*tuning,10);
		DGamRateProcess::Move(0.03*tuning,10);

		ExpoConjugateGTRFiniteProfileProcess::Move(1,1,10);

		if (! fixrr)	{
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
			ExpoConjugateGTRFiniteProfileProcess::Move(1,1,1);
		}
		return 1;
	}

	virtual void PrepareSiteLogLikelihood(int site) {
		int cat = ExpoConjugateGTRFiniteProfileProcess::alloc[site];
		if (! matrixarray[cat])	{
			cerr << "error in prepare site log likelihood: matrix is not allocated\n";
			exit(1);
			// CreateMatrix(cat);
		}
		UpdateMatrix(cat);
	}

	void ToStream(ostream& os)	{
		GammaBranchProcess::ToStream(os);
		DGamRateProcess::ToStream(os);
		ExpoConjugateGTRFiniteProfileProcess::ToStream(os);
	}

	void FromStream(istream& is)	{
		GammaBranchProcess::FromStream(is);
		DGamRateProcess::FromStream(is);
		ExpoConjugateGTRFiniteProfileProcess::FromStream(is);
		// GlobalUpdateParameters(); ?
	}

	virtual void ReadPB(int argc, char* argv[]);
	void ReadModeProfiles(string name, int burnin, int every, int until);
	void SlaveComputeCVScore();
	void SlaveComputeSiteLogL();

	protected:

	virtual void Create()	{
		ExpoConjugateGTRPhyloProcess::Create();
		RASCATGTRFiniteSubstitutionProcess::Create();
		GammaBranchProcess::Create();
	}
		
	virtual void Delete()	{
		GammaBranchProcess::Delete();
		RASCATGTRFiniteSubstitutionProcess::Delete();
		ExpoConjugateGTRPhyloProcess::Delete();
	}

};

#endif

