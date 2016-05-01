
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef ZIPEXPOCONJGTRMIXTUREPROFILE_H
#define ZIPEXPOCONJGTRMIXTUREPROFILE_H

#include "ZipExpoConjugateGTRProfileProcess.h"
#include "ExpoConjugateGTRMixtureProfileProcess.h"
#include "ZipMatrixMixtureProfileProcess.h"

// general subclass for all Matrix Mixture mixtures using generic sufficient statistics
class ZipExpoConjugateGTRMixtureProfileProcess : public virtual ZipExpoConjugateGTRProfileProcess, public virtual ZipMatrixMixtureProfileProcess, public virtual ExpoConjugateGTRMixtureProfileProcess  {

	public:

	ZipExpoConjugateGTRMixtureProfileProcess() {}
	virtual ~ZipExpoConjugateGTRMixtureProfileProcess() {}

	protected:

	virtual void Create()	{
		ZipMatrixMixtureProfileProcess::Create();
		ExpoConjugateGTRMixtureProfileProcess::Create();
		ZipExpoConjugateGTRProfileProcess::Create();
	}

	virtual void Delete()	{
		ZipExpoConjugateGTRProfileProcess::Delete();
		ExpoConjugateGTRMixtureProfileProcess::Delete();
		ZipMatrixMixtureProfileProcess::Delete();
	}

};

#endif

