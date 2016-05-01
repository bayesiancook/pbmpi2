
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef GPSSZIPRASCATGTRFINITEPHYLO_H
#define GPSSZIPRASCATGTRFINITEPHYLO_H

#include "GeneralPathSuffStatRASCATGTRFiniteGammaPhyloProcess.h"
#include "ZipGeneralPathSuffStatMatrixMixtureProfileProcess.h"
#include "ZipMatrixPhyloProcess.h"

class ZipGeneralPathSuffStatGTRFiniteProfileProcess : public virtual GeneralPathSuffStatGTRFiniteProfileProcess, public virtual ZipGeneralPathSuffStatMatrixMixtureProfileProcess	{

	public:

	ZipGeneralPathSuffStatGTRFiniteProfileProcess() {}
	virtual ~ZipGeneralPathSuffStatGTRFiniteProfileProcess() {}

	protected:

	virtual void Create()	{
		GeneralPathSuffStatGTRFiniteProfileProcess::Create();
		ZipGeneralPathSuffStatMatrixMixtureProfileProcess::Create();
	}

	virtual void Delete()	{
		ZipGeneralPathSuffStatMatrixMixtureProfileProcess::Delete();
		GeneralPathSuffStatGTRFiniteProfileProcess::Delete();
	}

};

class ZipGeneralPathSuffStatRASCATGTRFiniteSubstitutionProcess : public virtual GeneralPathSuffStatRASCATGTRFiniteSubstitutionProcess, public virtual ZipMatrixSubstitutionProcess, public virtual ZipGeneralPathSuffStatGTRFiniteProfileProcess	{

	public:

	ZipGeneralPathSuffStatRASCATGTRFiniteSubstitutionProcess() {}
	virtual ~ZipGeneralPathSuffStatRASCATGTRFiniteSubstitutionProcess() {}

	int GetNstate(int site) {return ZipMatrixSubstitutionProcess::GetNstate(site);}

	protected:

	virtual void Create()	{
		ZipMatrixSubstitutionProcess::Create();
		GeneralPathSuffStatRASCATGTRFiniteSubstitutionProcess::Create();
		ZipGeneralPathSuffStatGTRFiniteProfileProcess::Create();
	}

	virtual void Delete()	{
		ZipGeneralPathSuffStatGTRFiniteProfileProcess::Delete();
		GeneralPathSuffStatRASCATGTRFiniteSubstitutionProcess::Delete();
		ZipMatrixSubstitutionProcess::Delete();
	}
};

class ZipGeneralPathSuffStatRASCATGTRFiniteGammaPhyloProcess : public virtual GeneralPathSuffStatRASCATGTRFiniteGammaPhyloProcess, public virtual ZipGeneralPathSuffStatRASCATGTRFiniteSubstitutionProcess, public virtual ZipMatrixPhyloProcess	{

	public:

	ZipGeneralPathSuffStatRASCATGTRFiniteGammaPhyloProcess(int nratecat, int ncat, int infixncomp, int inempmix, string inmixtype, string inrrtype)	{

		Ncat = nratecat;
		Ncomponent = ncat;
		fixncomp = infixncomp;
		empmix = inempmix;
		mixtype = inmixtype;
		rrtype = inrrtype;
	}

	ZipGeneralPathSuffStatRASCATGTRFiniteGammaPhyloProcess(istream& is, int inmyid, int innprocs)	{

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

	~ZipGeneralPathSuffStatRASCATGTRFiniteGammaPhyloProcess() {
		Delete();
	}

	protected:

	virtual void SlaveExecute(MESSAGE signal);

	virtual void Create()	{
		ZipMatrixPhyloProcess::Create();
		ZipGeneralPathSuffStatRASCATGTRFiniteSubstitutionProcess::Create();
		GeneralPathSuffStatRASCATGTRFiniteGammaPhyloProcess::Create();
	}
		
	virtual void Delete()	{
		GeneralPathSuffStatRASCATGTRFiniteGammaPhyloProcess::Delete();
		ZipGeneralPathSuffStatRASCATGTRFiniteSubstitutionProcess::Delete();
		ZipMatrixPhyloProcess::Delete();
	}

};

#endif

