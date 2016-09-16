
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef ZIPRASCATGTRSBDPPHYLO_H
#define ZIPRASCATGTRSBDPPHYLO_H

#include "RASCATGTRSBDPGammaPhyloProcess.h"
#include "ZipExpoConjugateGTRMixtureProfileProcess.h"
#include "ZipExpoConjugateGTRPhyloProcess.h"

class ZipExpoConjugateGTRSBDPProfileProcess : public virtual ExpoConjugateGTRSBDPProfileProcess, public virtual ZipExpoConjugateGTRMixtureProfileProcess	{

	public:

	ZipExpoConjugateGTRSBDPProfileProcess() {}
	virtual ~ZipExpoConjugateGTRSBDPProfileProcess() {}

	protected:

	virtual void Create()	{
		ZipExpoConjugateGTRMixtureProfileProcess::Create();
		ExpoConjugateGTRSBDPProfileProcess::Create();
	}

	virtual void Delete()	{
		ExpoConjugateGTRSBDPProfileProcess::Delete();
		ZipExpoConjugateGTRMixtureProfileProcess::Delete();
	}

};

class ZipRASCATGTRSBDPSubstitutionProcess : public virtual RASCATGTRSBDPSubstitutionProcess, public virtual ZipExpoConjugateGTRSubstitutionProcess, public virtual ZipExpoConjugateGTRSBDPProfileProcess	{

	public:

	ZipRASCATGTRSBDPSubstitutionProcess() {}
	virtual ~ZipRASCATGTRSBDPSubstitutionProcess() {}

	protected:

	virtual void Create()	{
		RASCATGTRSBDPSubstitutionProcess::Create();
		ZipExpoConjugateGTRSBDPProfileProcess::Create();
		ZipExpoConjugateGTRSubstitutionProcess::Create();
	}

	virtual void Delete()	{
		ZipExpoConjugateGTRSubstitutionProcess::Delete();
		ZipExpoConjugateGTRSBDPProfileProcess::Delete();
		RASCATGTRSBDPSubstitutionProcess::Delete();
	}
};

class ZipRASCATGTRSBDPGammaPhyloProcess : public virtual RASCATGTRSBDPGammaPhyloProcess, public virtual ZipRASCATGTRSBDPSubstitutionProcess, public virtual ZipExpoConjugateGTRPhyloProcess	{

	public:

	ZipRASCATGTRSBDPGammaPhyloProcess(int nratecat, string inrrtype, int inkappaprior)	{

		Ncat = nratecat;
		rrtype = inrrtype;
		kappaprior = inkappaprior;
	}

	ZipRASCATGTRSBDPGammaPhyloProcess(istream& is, int inmyid, int innprocs)	{

		// generic
		FromStreamHeader(is);
		SetMPI(inmyid,innprocs);

		// specific
		is >> Ncat;
		is >> kappaprior;
		is >> rrtype;

		Open(is);

	}

	~ZipRASCATGTRSBDPGammaPhyloProcess() {
		Delete();
	}

	void ToStreamHeader(ostream& os)	{
		PhyloProcess::ToStreamHeader(os);
		os << Ncat << '\n';
		os << kappaprior << '\n';
		os << rrtype << '\n';
	}

	virtual int GetNstate() {
		return ZipExpoConjugateGTRSubstitutionProcess::GetNstate();
	}

	virtual int GetNstate(int i) {
		return ZipExpoConjugateGTRSubstitutionProcess::GetNstate(i);
	}

	/*
	virtual double Move(double tuning = 1.0)	{

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

		ZipExpoConjugateGTRSBDPProfileProcess::Move(1,1,10);

		if (! fixrr)	{
			LengthRelRateMove(1,10);
			LengthRelRateMove(0.1,10);
			LengthRelRateMove(0.01,10);
		}


		GlobalUnfold();

		chronototal.Stop();

		return 1;
	}
	*/

	protected:

	// virtual void SlaveExecute(MESSAGE signal);

	virtual void Create()	{
		ZipExpoConjugateGTRPhyloProcess::Create();
		ZipRASCATGTRSBDPSubstitutionProcess::Create();
		RASCATGTRSBDPGammaPhyloProcess::Create();
	}
		
	virtual void Delete()	{
		RASCATGTRSBDPGammaPhyloProcess::Delete();
		ZipRASCATGTRSBDPSubstitutionProcess::Delete();
		ZipExpoConjugateGTRPhyloProcess::Delete();
	}
};

#endif

