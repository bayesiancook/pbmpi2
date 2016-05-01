
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef GENPATHSSMATSSPROFILE_H
#define GENPATHSSMATSSPROFILE_H

#include "GeneralPathSuffStatMatrixProfileProcess.h"
#include "MatrixSiteSpecificProfileProcess.h"

// general subclass for all Matrix SiteSpecific mixtures using generic sufficient statistics
class GeneralPathSuffStatMatrixSiteSpecificProfileProcess : public virtual MatrixSiteSpecificProfileProcess, public virtual GeneralPathSuffStatMatrixProfileProcess  {

	public:

	GeneralPathSuffStatMatrixSiteSpecificProfileProcess() {}
	virtual ~GeneralPathSuffStatMatrixSiteSpecificProfileProcess() {}

	protected:

	virtual void Create()	{
		MatrixSiteSpecificProfileProcess::Create();
		GeneralPathSuffStatMatrixProfileProcess::Create();
	}

	virtual void Delete()	{
		GeneralPathSuffStatMatrixProfileProcess::Delete();
		MatrixSiteSpecificProfileProcess::Delete();
	}

	virtual void CreateSite(int k)	{
		SampleStat(k);
		CreateMatrix(k);
		UpdateMatrix(k);
	}

	virtual void DeleteSite(int k)	{
		DeleteMatrix(k);
	}

	virtual double ProfileSuffStatLogProb(int cat);

};

#endif

