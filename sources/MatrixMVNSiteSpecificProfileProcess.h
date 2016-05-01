
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef MATMVNSSPROFILE_H
#define MATMVNSSPROFILE_H

#include "MVNSiteSpecificProfileProcess.h"
#include "MatrixSiteSpecificProfileProcess.h"

class MatrixMVNSiteSpecificProfileProcess : public virtual MVNSiteSpecificProfileProcess, public virtual MatrixSiteSpecificProfileProcess {

	public:

	MatrixMVNSiteSpecificProfileProcess() {}
	virtual ~MatrixMVNSiteSpecificProfileProcess() {}

	protected:

	virtual void UpdateSite(int k)	{
		MVNSiteSpecificProfileProcess::UpdateSite(k);
		MatrixSiteSpecificProfileProcess::UpdateSite(k);
	}
	
	virtual void Create()	{
		MVNSiteSpecificProfileProcess::Create();
		MatrixSiteSpecificProfileProcess::Create();
	}

	virtual void Delete()	{
		MatrixSiteSpecificProfileProcess::Delete();
		MVNSiteSpecificProfileProcess::Delete();
	}
};

#endif
