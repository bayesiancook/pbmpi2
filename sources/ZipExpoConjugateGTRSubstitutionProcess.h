
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef ZIPEXPCONGTRSUB_H
#define ZIPEXPCONGTRSUB_H

#include "ExpoConjugateGTRSubstitutionProcess.h"
#include "ZipExpoConjugateGTRProfileProcess.h"
#include "ZipMatrixSubstitutionProcess.h"

class ZipExpoConjugateGTRSubstitutionProcess : public virtual ExpoConjugateGTRSubstitutionProcess, public virtual ZipExpoConjugateGTRProfileProcess, public virtual ZipMatrixSubstitutionProcess	{

	public:

	ZipExpoConjugateGTRSubstitutionProcess() : zipmode(0) {}
	virtual ~ZipExpoConjugateGTRSubstitutionProcess() {}

	// need to redefine ? 
	virtual int GetNstate() {return GetDim();}
	virtual int GetNstate(int i) {
		if (zipmode)	{
			return ZipMatrixSubstitutionProcess::GetNstate(i);
		}
		return GetNstate();
	}

	virtual SubMatrix* GetMatrix(int site)	{
		if (zipmode)	{
			return ZipMatrixSubstitutionProcess::GetMatrix(site);
		}
		return GetOriginalMatrix(site);
	}

	void ActivateZip()	{
		zipmode = 1;
	}

	void InactivateZip()	{
		zipmode = 0;
	}

	protected:

	void Create()	{
		ZipExpoConjugateGTRProfileProcess::Create();
		ZipMatrixSubstitutionProcess::Create();
		ExpoConjugateGTRSubstitutionProcess::Create();
	}

	void Delete() {
		ExpoConjugateGTRSubstitutionProcess::Delete();
		ZipExpoConjugateGTRProfileProcess::Delete();
		ZipMatrixSubstitutionProcess::Delete();
	}

	int zipmode;
};

#endif

