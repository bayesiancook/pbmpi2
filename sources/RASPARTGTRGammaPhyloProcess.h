
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef RASPARTGTR_H
#define RASPARTGTR_H

#include "PartitionGTRFiniteProfileProcess.h"
#include "ExpoConjugatePartitionGTRMixtureProfileProcess.h"
#include "ExpoConjugatePartitionGTRPhyloProcess.h"
#include "PartitionDGamRateProcess.h"
#include "GammaBranchProcess.h"

class ExpoConjugatePartitionGTRFiniteProfileProcess : public virtual PartitionGTRFiniteProfileProcess, public virtual ExpoConjugatePartitionGTRMixtureProfileProcess {

	public:

	ExpoConjugatePartitionGTRFiniteProfileProcess() {}
	virtual ~ExpoConjugatePartitionGTRFiniteProfileProcess() {}

	protected:

	virtual void Create()	{
		PartitionGTRFiniteProfileProcess::Create();
		ExpoConjugatePartitionGTRMixtureProfileProcess::Create();
	}

	virtual void Delete()	{
		ExpoConjugatePartitionGTRMixtureProfileProcess::Delete();
		PartitionGTRFiniteProfileProcess::Delete();
	}
};

class RASPARTGTRSubstitutionProcess : public virtual ExpoConjugatePartitionGTRSubstitutionProcess, public virtual PartitionDGamRateProcess, public virtual ExpoConjugatePartitionGTRFiniteProfileProcess {

	public:

	RASPARTGTRSubstitutionProcess() {}
	virtual ~RASPARTGTRSubstitutionProcess() {}

	virtual void Create()	{
		ExpoConjugatePartitionGTRSubstitutionProcess::Create();
		PartitionDGamRateProcess::Create();
		ExpoConjugatePartitionGTRFiniteProfileProcess::Create();
	}

	virtual void Delete()	{
		ExpoConjugatePartitionGTRFiniteProfileProcess::Delete();
		PartitionDGamRateProcess::Delete();
		ExpoConjugatePartitionGTRSubstitutionProcess::Delete();
	}
};

class RASPARTGTRGammaPhyloProcess : public virtual ExpoConjugatePartitionGTRPhyloProcess, public virtual RASPARTGTRSubstitutionProcess, public virtual GammaBranchProcess	{

	public:

        virtual void SlaveExecute(MESSAGE);

	virtual void GlobalUpdateParameters();

	RASPARTGTRGammaPhyloProcess() {}

	RASPARTGTRGammaPhyloProcess(int nratecat, int inwithpinv, string inrrtype)	{

		Ncat = nratecat;
		withpinv = inwithpinv;
		rrtype = inrrtype;
	}

	RASPARTGTRGammaPhyloProcess(istream& is, int inmyid, int innprocs)	{

		// generic
		FromStreamHeader(is);
		SetMPI(inmyid,innprocs);

		// specific
		is >> Ncat;
		// is >> withpinv;
		is >> rrtype;

		Open(is);
	}

	void ToStreamHeader(ostream& os)	{
		PhyloProcess::ToStreamHeader(os);
		os << Ncat << '\n';
		os << rrtype << '\n';
	}

	~RASPARTGTRGammaPhyloProcess() {
		Delete();
	}

	virtual void SlaveUpdateParameters();

	void TraceHeader(ostream& os)	{
		os << "#iter\ttime\ttopo\tloglik\tlength\tmeanalpha\tvaralpha";
		if (withpinv)	{
			os << "\tmeanpinv\tvarpinv";
		}
		os << "\tNmode\tstatent\tstatalpha";
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
		os << '\t' << GetMeanAlpha();
		os << '\t' << GetVarAlpha();
		if (withpinv)	{
			os << '\t' << GetMeanPinv();
			os << '\t' << GetVarPinv();
		}
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

		if (! FixBL())	{
			if ((topobf != 3) && ((topobf != 1) || (size < bfburnin)))	{
				BranchLengthMove(tuning);
				BranchLengthMove(0.1 * tuning);
			}
		}

		if (! FixTopo())	{
			MoveTopo();
		}

		propchrono.Stop();

		// for (int rep=0; rep<5; rep++)	{
			GlobalCollapse();
			AugmentedMove(tuning);
			GlobalUnfold();
		// }

		chronototal.Stop();

		return 1;
	}

	double AugmentedMove(double tuning = 1)	{

		if (! FixBL())	{
			GammaBranchProcess::Move(tuning,10);
			GammaBranchProcess::Move(0.1*tuning,10);
			GlobalUpdateParameters();
		}


		if (! FixAlpha())	{
			PartitionDGamRateProcess::Move(tuning,10);
			PartitionDGamRateProcess::Move(0.3*tuning,10);
			PartitionDGamRateProcess::Move(0.03*tuning,10);
		}

		GlobalUpdateParameters();
		ExpoConjugatePartitionGTRFiniteProfileProcess::Move(1,1,10);
		GlobalUpdateParameters();

		if (! FixRR()){
			LengthRelRateMove(1,10);
			LengthRelRateMove(0.1,10);
			LengthRelRateMove(0.01,10);
		}
		GlobalUpdateParameters();
	}


	virtual double GlobalRestrictedMoveCycle(int nrep = 1, double tuning = 1.0)	{

		for (int rep=0; rep<nrep; rep++)	{

			GlobalUpdateParameters();
			GammaBranchProcess::Move(tuning,10);

			GlobalUpdateParameters();
			PartitionDGamRateProcess::Move(0.3*tuning,10);
			PartitionDGamRateProcess::Move(0.03*tuning,10);
			ExpoConjugatePartitionGTRFiniteProfileProcess::Move(1,1,1);
		}
		return 1;
	}

	virtual void PrepareSiteLogLikelihood(int site) {
		int cat = ExpoConjugatePartitionGTRFiniteProfileProcess::alloc[site];
		if (! matrixarray[cat])	{
			cerr << "error in prepare site log likelihood: matrix is not allocated\n";
			exit(1);
			// CreateMatrix(cat);
		}
		UpdateMatrix(cat);
	}

	void ToStream(ostream& os)	{
		GammaBranchProcess::ToStream(os);
		PartitionDGamRateProcess::ToStream(os);
		ExpoConjugatePartitionGTRFiniteProfileProcess::ToStream(os);
	}

	void FromStream(istream& is)	{
		GammaBranchProcess::FromStream(is);
		PartitionDGamRateProcess::FromStream(is);
		ExpoConjugatePartitionGTRFiniteProfileProcess::FromStream(is);
		// GlobalUpdateParameters(); ?
	}

	virtual void ReadPB(int argc, char* argv[]);
	void ReadModeProfiles(string name, int burnin, int every, int until);
	void SlaveComputeCVScore();
	void SlaveComputeSiteLogL();

	protected:

	virtual void Create()	{
		ExpoConjugatePartitionGTRPhyloProcess::Create();
		RASPARTGTRSubstitutionProcess::Create();
		GammaBranchProcess::Create();
	}
		
	virtual void Delete()	{
		GammaBranchProcess::Delete();
		RASPARTGTRSubstitutionProcess::Delete();
		ExpoConjugatePartitionGTRPhyloProcess::Delete();
	}

};

#endif

