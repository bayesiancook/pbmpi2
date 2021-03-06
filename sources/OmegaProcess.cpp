
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#include "Parallel.h"
#include "OmegaProcess.h"
#include "GeneralPathSuffStatMatrixProfileProcess.h"
#include "Random.h"

void OmegaProcess::GlobalUpdateSiteOmegaSuffStat()	{
	if (GetNprocs() > 1)	{
		// MPI2
		// should ask the slaves to call their UpdateRateSuffStat
		// and then gather the statistics;
		MPI_Status stat;
		MESSAGE signal = UPDATE_SITEOMEGA;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	}
	else	{
		UpdateSiteOmegaSuffStat();
	}
}

void OmegaProcess::SlaveUpdateSiteOmegaSuffStat()	{
	UpdateSiteOmegaSuffStat();
}
			

