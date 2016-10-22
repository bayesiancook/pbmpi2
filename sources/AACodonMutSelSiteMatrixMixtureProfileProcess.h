
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef AACODONMUTSELSITEMATMIX_H
#define AACODONMUTSELSITEMATMIX_H

#include "AACodonMutSelProfileProcess.h"
#include "GeneralPathSuffStatMatrixProfileProcess.h"
#include "SiteMatrixMixtureProfileProcess.h"

class AACodonMutSelSiteMatrixMixtureProfileProcess : public virtual AACodonMutSelProfileProcess, public virtual SiteMatrixMixtureProfileProcess, public virtual GeneralPathSuffStatMatrixProfileProcess	{

	public:

	AACodonMutSelSiteMatrixMixtureProfileProcess() {}
	virtual ~AACodonMutSelSiteMatrixMixtureProfileProcess() {}

	virtual void Create() {}
	virtual void Delete() {}

	SubMatrix* GetMatrix(int site)	{
		return matrixarray[site];
	}

	virtual void UpdateModeProfileSuffStat() {}
	virtual void GlobalUpdateModeProfileSuffStat() {}
	virtual void SlaveUpdateModeProfileSuffStat() {}

	virtual double ProfileSuffStatLogProb(int cat) {}
	virtual double ProfileSuffStatLogProb() {}

	virtual double LogStatProb(int site, int cat) {}

};

#endif

