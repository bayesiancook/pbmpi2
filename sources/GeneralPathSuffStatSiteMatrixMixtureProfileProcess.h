
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef GENPATHSSSITEMATMIXTUREPROFILE_H
#define GENPATHSSSITEMATMIXTUREPROFILE_H


#include "GeneralPathSuffStatMatrixMixtureProfileProcess.h"
#include "SiteMatrixMixtureProfileProcess.h"

class GeneralPathSuffStatSiteMatrixMixtureProfileProcess : public virtual GeneralPathSuffStatMatrixMixtureProfileProcess, public virtual SiteMatrixMixtureProfileProcess	{

	public:

	GeneralPathSuffStatSiteMatrixMixtureProfileProcess() {}
	virtual ~GeneralPathSuffStatSiteMatrixMixtureProfileProcess() {}

	protected:

	virtual void Create()	{
		SiteMatrixMixtureProfileProcess::Create();
		GeneralPathSuffStatMatrixMixtureProfileProcess::Create();
	}

	virtual void Delete()	{
		GeneralPathSuffStatMatrixMixtureProfileProcess::Delete();
		SiteMatrixMixtureProfileProcess::Delete();
	}

	virtual void UpdateModeProfileSuffStat()	{
		// this one should be reimplemented
		cerr << "error: in GPSSSiteMatMixProfileProcess::UpdateModeProfileSuffStat\n";
		exit(1);
	}

	virtual double LogStatProb(int site, int cat)	{
		// this one should be reimplemented
		cerr << "error: in GPSSSiteMatMixProfileProcess::LogStatProb\n";
		exit(1);
	}
};

#endif
