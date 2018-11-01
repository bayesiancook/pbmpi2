
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#include "RASCATGTRDPGammaPhyloProcess.h"
#include "Parallel.h"
#include <string>

void RASCATGTRDPGammaPhyloProcess::GlobalUpdateParameters()	{

	if (GetNprocs() > 1)	{

        int nd = 2 + GetNbranch() + GetNrr() + GetNmodeMax()*GetDim() + Nstatcomp*(GetDim()+1) + 1;
        int ni = 1 + GetNsite();
        int* ivector = new int[ni];
        double* dvector = new double[nd]; 

        MESSAGE signal = PARAMETER_DIFFUSION;
        MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

        // GlobalBroadcastTree();

        // First we assemble the vector of doubles for distribution
        int index = 0;
        dvector[index] = GetAlpha();
        index++;
        dvector[index] = GetPinv();
        index++;
        for(int i=0; i<GetNbranch(); i++)    {
            dvector[index] = blarray[i];
            index++;
        }
        
        for (int i=0; i<GetNrr(); i++)  {
            dvector[index] = rr[i];
            index++;
        }

        for (int i=0; i<GetNmodeMax(); i++) {
            for (int j=0; j<GetDim(); j++)  {
                dvector[index] = profile[i][j];
                index++;
            }
        }

        for (int k=0; k<Nstatcomp; k++) {
            dvector[index] = statweight[k];
            index++;
            for (int i=0; i<GetDim(); i++)	{
                dvector[index] = dirweight[k][i];
                index++;
            }
        }

        dvector[index] = kappa;
        index++;

        if (index != nd)    {
            cerr << "error in global update params: non matching dim\n";
            exit(1);
        }

        // Now the vector of ints
        ivector[0] = GetNcomponent();
        for(int i=0; i<GetNsite(); i++) {
            ivector[1+i] = DPProfileProcess::alloc[i];
        }

        // Now send out the doubles and ints over the wire...
        MPI_Bcast(ivector,ni,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(dvector,nd,MPI_DOUBLE,0,MPI_COMM_WORLD);

        delete[] dvector;
        delete[] ivector;
    }

    else    {
        UpdateMatrices();
    }

}


void RASCATGTRDPGammaPhyloProcess::SlaveUpdateParameters()	{

	// SlaveBroadcastTree();

	int nd = 2 + GetNbranch() + GetNrr() + GetNmodeMax()*GetDim() + Nstatcomp*(GetDim()+1) + 1;
	int ni = 1 + GetNsite();
	int* ivector = new int[ni];
	double* dvector = new double[nd];

	MPI_Bcast(ivector,ni,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(dvector,nd,MPI_DOUBLE,0,MPI_COMM_WORLD);

	int index = 0;
	SetRateParams(dvector[index],dvector[index+1]);
	index+=2;

    for (int i=0; i<GetNbranch(); i++)   {
		blarray[i] = dvector[index];
		index++;
	}

    for (int i=0; i<GetNrr(); i++)  {
		rr[i] = dvector[index];
		index++;
	}

    for (int i=0; i<GetNmodeMax(); i++) {
        for (int j=0; j<GetDim(); j++)  {
			profile[i][j] = dvector[index];
			index++;
		}
	}

    for (int k=0; k<Nstatcomp; k++) {
        statweight[k] = dvector[index];
        index++;
        for (int i=0; i<GetDim(); i++)	{
            dirweight[k][i] = dvector[index];
            index++;
        }
    }
	kappa = dvector[index];
	index++;
    
    if (index != nd)    {
        cerr << "error in slave update params: non matching dim\n";
        exit(1);
    }

	Ncomponent = ivector[0];
	for(int i=0; i<GetNsite(); i++) {
		DPProfileProcess::alloc[i] = ivector[1+i];
	}

	delete[] dvector;
	delete[] ivector;
	// some upate here ?
}

void RASCATGTRDPGammaPhyloProcess::SlaveExecute(MESSAGE signal)	{

	switch(signal) {

	case UPDATE_RATE:
		SlaveUpdateRateSuffStat();
		break;
	case UPDATE_RRATE:
		SlaveUpdateRRSuffStat();
		break;
	case PROFILE_MOVE:
		SlaveMoveProfile();
		break;
	default:
		PhyloProcess::SlaveExecute(signal);
	}
}
