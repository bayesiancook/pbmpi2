
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
#include "AACodonMutSelMVNSiteSpecificPhyloProcess.h"
#include "Parallel.h"
#include <string.h>


// MPI: these two functions are responsible for broadcasting/receiving the current state of the parameter vector
// are model dependent
// should be implemented in .cpp file
void AACodonMutSelMVNSiteSpecificPhyloProcess::SlaveUpdateParameters()	{

	int i,j,L1,L2,nd,nbranch = GetNbranch(),nnucrr = GetNnucrr(),nnucstat = 4;
	L1 = GetNsite();
	L2 = GetDim();
	int nstate = GetData()->GetNstate();
	nd = 2+ nbranch + nnucrr + nnucstat + L2 + L2*(L2+1)/2 + L1*L2 + nstate + 1; // check if these last terms are correct in this context...
	double* dvector = new double[nd];
	MPI_Bcast(dvector,nd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	int index = 0;
	dvector[index] = branchalpha;
	index++;
	dvector[index] = branchbeta;
	index++;
	for(i=0; i<nbranch; ++i) {
		blarray[i] = dvector[index];
		index++;
	}
	for(i=0; i<nnucrr; ++i) {
		nucrr[i] = dvector[index];
		index++;
	}
	for(i=0; i<nnucstat; ++i) {
		nucstat[i] = dvector[index];
		index++;
	}
	for (i=0; i<nstate; i++)	{
		codonprofile[i] = dvector[index];
		index++;
	}
	for(i=0; i<L1; ++i) {
		for(j=0; j<L2; ++j) {
			logprofile[i][j] = dvector[index];
			index++;
		}
		// UpdateSite(i);
	}
	for (int i=0; i<GetDim(); i++)	{
		kappa[i] = dvector[index];
		index++;
	}
	for (int i=0; i<GetDim(); i++)	{
		for (int j=i; j<GetDim(); j++)	{
			(*covmatrix)[i][j] = dvector[index];
			if (i != j)	{
				(*covmatrix)[j][i] = (*covmatrix)[i][j];
			}
			index++;
		}
	}
			
	*omega = dvector[index];
	index++;
	delete[] dvector;

	UpdateSites();
	CreateMatrices();
	UpdateMatrices();
}


void AACodonMutSelMVNSiteSpecificPhyloProcess::SlaveExecute(MESSAGE signal)	{

	switch(signal) {

		case PROFILE_MOVE:
			SlaveMoveProfile();
			break;
		default:
			PhyloProcess::SlaveExecute(signal);
	}
}

void AACodonMutSelMVNSiteSpecificPhyloProcess::GlobalUpdateParameters() {

	if (GetNprocs() > 1)	{
	// MPI2
	// should send the slaves the relevant information
	// about model parameters
	// for this model, should broadcast
	// (but should first call PutBranchLengthsIntoArray())
	// 
	// upon receiving this information
	// slave should 
	// store it in the local copies of the variables
	// and then call
	// SetBranchLengthsFromArray()
	int i,j,nnucrr,nnucstat,nbranch = GetNbranch(),nd,L1,L2;
	nnucrr = GetNnucrr();
	nnucstat = 4;	
	L1 = GetNsite();
	L2 = GetDim();
	int nstate = GetData()->GetNstate();
	nd = 2+ nbranch + nnucrr + nnucstat + L2 + L2*(L2+1)/2 + L1*L2 + nstate + 1; // check if these last terms are correct in this context...
	double dvector[nd]; 
	MESSAGE signal = PARAMETER_DIFFUSION;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	
	int index = 0;
	index++;
	dvector[index] = branchalpha;
	index++;
	dvector[index] = branchbeta;
	// First we assemble the vector of doubles for distribution
	for(i=0; i<nbranch; ++i) {
		dvector[index] = blarray[i];
		index++;
	}

	for(i=0; i<nnucrr; ++i) {
		dvector[index] = nucrr[i];
		index++;
	}
	for(i=0; i<nnucstat; ++i) {
		dvector[index] = nucstat[i];
		index++;
	}

	for (i=0; i<nstate; i++)	{
		dvector[index] = codonprofile[i];
		index++;
	}

	for(i=0; i<L1; ++i) {
		for(j=0; j<L2; ++j) {
			dvector[index] = logprofile[i][j];
			index++;
		}
	}
	for (int i=0; i<GetDim(); i++)	{
		dvector[index] = kappa[i];
		index++;
	}
	for (int i=0; i<GetDim(); i++)	{
		for (int j=i; j<GetDim(); j++)	{
			dvector[index] = (*covmatrix)[i][j];
			index++;
		}
	}

	dvector[index] = *omega;
	index++;

	// Now send out the doubles and ints over the wire...
	MPI_Bcast(dvector,nd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	}
}


void AACodonMutSelMVNSiteSpecificPhyloProcess::ReadPB(int argc, char* argv[])	{


	// Needs updating!

	string name = "";

	int burnin = 0;
	int every = 1;
	int until = -1;
	int ppred = 0;
	// 1 : plain ppred (outputs simulated data)
	// 2 : diversity statistic
	// 3 : compositional statistic

	try	{

		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];
			if ( (s == "-x") || (s == "-extract") )	{
				i++;
				if (i == argc) throw(0);
				burnin = atoi(argv[i]);
				i++;
				if (i == argc) throw(0);
				s = argv[i];
				if (IsInt(s))	{
					every = atoi(argv[i]);
					i++;
					if (i == argc) throw(0);
					string tmp = argv[i];
					if (IsInt(tmp))	{
						until = atoi(argv[i]);
					}
					else	{
						i--;
					}
				}
				else {
					i--;
				}
			}
			else	{
				if (i != (argc -1))	{
					throw(0);
				}
				name = argv[i];
			}
			i++;
		}
	}
	catch(...)	{
		cerr << "error in command\n";
		cerr << '\n';
		exit(1);
	}

	if (until == -1)	{
		until = GetSize();
	}

	Read(name,burnin,every,until);
}
