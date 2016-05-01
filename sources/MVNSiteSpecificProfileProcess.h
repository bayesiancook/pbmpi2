
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef MVNSSPROFILE_H
#define MVNSSPROFILE_H

#include "SiteSpecificProfileProcess.h"
#include "MVNProfileProcess.h"

class MVNSiteSpecificProfileProcess : public virtual MVNProfileProcess, public virtual SiteSpecificProfileProcess {

	public:

	MVNSiteSpecificProfileProcess() : logprofile(0), alloclogprofile(0)  {}
	virtual ~MVNSiteSpecificProfileProcess() {}

	void ResampleCovMatrix();

	protected:

	virtual void Create();
	virtual void Delete();

	virtual void UpdateSite(int k)	{
		UpdateFrequencyStat(logprofile[k],profile[k]);
	}

	virtual void SampleStat(int site);

	virtual double MoveProfile(int site, double tuning, int n, int nrep);

	double** logprofile;
	double* alloclogprofile;
};


#endif
