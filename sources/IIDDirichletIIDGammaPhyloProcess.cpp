
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "StringStreamUtils.h"
#include "IIDDirichletIIDGammaPhyloProcess.h"
#include "Parallel.h"
#include <string>

void IIDDirichletIIDGammaPhyloProcess::GlobalUpdateParameters()	{

	if (GetNprocs() > 1)	{

        RenormalizeProfiles();

        int nd = 1 + GetNsite() + GetNbranch() + GetNsite()*GetDim() + GetDim();
        double* dvector = new double[nd];

        MESSAGE signal = PARAMETER_DIFFUSION;
        MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

        // First we assemble the vector of doubles for distribution
        int index = 0;
        dvector[index] = GetAlpha();
        index++;
        for (int i=0; i<GetNsite(); i++)    {
            dvector[index] = rate[i];
            index++;
        }

        for(int i=0; i<GetNbranch(); i++) {
            dvector[index] = blarray[i];
            index++;
        }
        
        for(int i=0; i<GetNsite(); ++i) {
            for(int j=0; j<GetDim(); ++j) {
                dvector[index] = profile[i][j];
                index++;
            }
        }
        for (int i=0; i<GetDim(); i++)	{
            dvector[index] = dirweight[i];
            index++;
        }

        // Now send out the doubles and ints over the wire...
        MPI_Bcast(dvector,nd,MPI_DOUBLE,0,MPI_COMM_WORLD);

        delete[] dvector;

	}
	else	{
		UpdateZip();
	}
}

void IIDDirichletIIDGammaPhyloProcess::SlaveUpdateParameters()	{

	int nd = 1 + GetNsite() + GetNbranch() + GetNsite()*GetDim() + GetDim();
    double* dvector = new double[nd];

	MPI_Bcast(dvector,nd,MPI_DOUBLE,0,MPI_COMM_WORLD);

	int index = 0;
    alpha = dvector[index];
    index++;

    for (int i=0; i<GetNsite(); i++)    {
        rate[i] = dvector[index];
        index++;
    }

	for(int i=0; i<GetNbranch(); ++i) {
		blarray[i] = dvector[index];
		index++;
	}

	for(int i=0; i<GetNsite(); ++i) {
		for(int j=0; j<GetDim(); ++j) {
			profile[i][j] = dvector[index];
			index++;
		}
	}
	for (int i=0; i<GetDim(); i++)	{
		dirweight[i] = dvector[index];
		index++;
	}

	delete[] dvector;

	UpdateZip();
}

void IIDDirichletIIDGammaPhyloProcess::SlaveExecute(MESSAGE signal)	{

	switch(signal) {

	case UPDATE_HYPERRATE:
		SlaveUpdateRateHyperSuffStat();
		break;
	case UPDATE_HYPERPROFILE:
		SlaveUpdateProfileHyperSuffStat();
		break;
    case MOVE_SRATE:
        ResampleSiteRates();
        break;
    case MOVE_SPROFILE:
        ResampleSiteProfiles();
        break;
    case COLLECT_SRATE:
        SlaveCollectSiteRates();
        break;
    case COLLECT_SPROFILE:
        SlaveCollectSiteProfiles();
        break;
	default:
		PhyloProcess::SlaveExecute(signal);
	}
}

