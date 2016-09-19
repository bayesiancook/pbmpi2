
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef ZIPRASCATGTRFinitePHYLO_H
#define ZIPRASCATGTRFinitePHYLO_H

#include "RASCATGTRFiniteGammaPhyloProcess.h"
#include "ZipExpoConjugateGTRMixtureProfileProcess.h"
#include "ZipExpoConjugateGTRPhyloProcess.h"

class ZipExpoConjugateGTRFiniteProfileProcess : public virtual ExpoConjugateGTRFiniteProfileProcess, public virtual ZipExpoConjugateGTRMixtureProfileProcess	{

	public:

	ZipExpoConjugateGTRFiniteProfileProcess() {}
	virtual ~ZipExpoConjugateGTRFiniteProfileProcess() {}

	protected:

	virtual void Create()	{
		ZipExpoConjugateGTRMixtureProfileProcess::Create();
		ExpoConjugateGTRFiniteProfileProcess::Create();
	}

	virtual void Delete()	{
		ExpoConjugateGTRFiniteProfileProcess::Delete();
		ZipExpoConjugateGTRMixtureProfileProcess::Delete();
	}

};

class ZipRASCATGTRFiniteSubstitutionProcess : 
public virtual RASCATGTRFiniteSubstitutionProcess, 
public virtual ZipExpoConjugateGTRSubstitutionProcess, 
public virtual ZipExpoConjugateGTRFiniteProfileProcess	{

	public:

	ZipRASCATGTRFiniteSubstitutionProcess() {}
	virtual ~ZipRASCATGTRFiniteSubstitutionProcess() {}

	protected:

	virtual void Create()	{
		RASCATGTRFiniteSubstitutionProcess::Create();
		ZipExpoConjugateGTRFiniteProfileProcess::Create();
		ZipExpoConjugateGTRSubstitutionProcess::Create();
	}

	virtual void Delete()	{
		ZipExpoConjugateGTRSubstitutionProcess::Delete();
		ZipExpoConjugateGTRFiniteProfileProcess::Delete();
		RASCATGTRFiniteSubstitutionProcess::Delete();
	}
};

class ZipRASCATGTRFiniteGammaPhyloProcess : public virtual RASCATGTRFiniteGammaPhyloProcess, public virtual ZipRASCATGTRFiniteSubstitutionProcess, public virtual ZipExpoConjugateGTRPhyloProcess	{

	public:

	ZipRASCATGTRFiniteGammaPhyloProcess(int nratecat, int inwithpinv, int ncat, int infixncomp, int inempmix, string inmixtype, string inrrtype)	{

		Ncat = nratecat;
		withpinv = inwithpinv;
		Ncomponent = ncat;
		fixncomp = infixncomp;
		empmix = inempmix;
		mixtype = inmixtype;
		rrtype = inrrtype;
	}

	ZipRASCATGTRFiniteGammaPhyloProcess(istream& is, int inmyid, int innprocs)	{

		// generic
		FromStreamHeader(is);
		SetMPI(inmyid,innprocs);

		// specific
		is >> Ncat;
		is >> withpinv;
		is >> fixncomp;
		is >> empmix;
		is >> mixtype;
		is >> rrtype;

		Open(is);

	}

	~ZipRASCATGTRFiniteGammaPhyloProcess() {
		Delete();
	}

	void ToStreamHeader(ostream& os)	{
		PhyloProcess::ToStreamHeader(os);
		os << Ncat << '\n';
		os << fixncomp << '\t' << empmix << '\t' << mixtype << '\n';
		os << rrtype << '\n';
	}

	virtual int GetNstate() {
		return ZipExpoConjugateGTRSubstitutionProcess::GetNstate();
	}

	virtual int GetNstate(int i) {
		return ZipExpoConjugateGTRSubstitutionProcess::GetNstate(i);
	}

	/*
	double Move(double tuning = 1.0)	{

		chronototal.Start();
		propchrono.Start();

		ZipTopoMoveCycle(1,tuning);

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
	*/

	protected:

	// virtual void SlaveExecute(MESSAGE signal);

	virtual void Create()	{
		ZipExpoConjugateGTRPhyloProcess::Create();
		ZipRASCATGTRFiniteSubstitutionProcess::Create();
		RASCATGTRFiniteGammaPhyloProcess::Create();
	}
		
	virtual void Delete()	{
		RASCATGTRFiniteGammaPhyloProcess::Delete();
		ZipRASCATGTRFiniteSubstitutionProcess::Delete();
		ZipExpoConjugateGTRPhyloProcess::Delete();
	}
};

#endif

