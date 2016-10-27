/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "SBDPOmegaProcess.h"
#include "Parallel.h"


void SBDPOmegaProcess::SampleOmega()	{
	/*
	for (int k=0; k<GetNomega(); k++)	{
		omega[k] = 1.0;
	}
	*/
	SampleOmegas();
	// some other things to either sample or initialize
	SampleOmegaWeights();
	SampleOmegaAlloc();
}

double SBDPOmegaProcess::MPIMoveOmega(double tuning, int nrep)	{

	double total=0;

	for (int rep=0; rep<nrep; rep++)	{

		total += GlobalOmegaIncrementalFiniteMove(1);
		/*
		MoveOccupiedCompAlloc(5);
		MoveAdjacentCompAlloc(5);
		*/

		GlobalUpdateOmegaSuffStat();

		total += MoveOmegaValuesAndHyperParameters(tuning,nrep);

		// Move kappa etc..

		ResampleOmegaWeights();

		// perhaps some broadcast here, 
		// sync slaves and master
	}

	return total;
}

double SBDPOmegaProcess::NonMPIMoveOmega(double tuning, int nrep)	{

	double total=0;

	for (int rep=0; rep<nrep; rep++)	{

		// yet to be implemented
		// total += OmegaIncrementalFiniteMove();
		/*
		MoveOccupiedCompAlloc(5);
		MoveAdjacentCompAlloc(5);
		*/

		UpdateOmegaSuffStat();

		total += MoveOmegaValuesAndHyperParameters(tuning,nrep);

		// Move kappa etc..

		ResampleOmegaWeights();
	}

	return total;
}

void SBDPOmegaProcess::SampleOmegaWeights()	{
}

void SBDPOmegaProcess::ResampleOmegaWeights()	{
}
// double MoveOccupiedCompAlloc(int nrep);
// double MoveAdjacentCompAlloc(int nrep);
// double MoveKappa


/*
double SBDPOmegaProcess::MixMoveOmega(int nmix, double tuning, int nsiterep, int nhyperrep)	{

	for (int mix=0; mix<nmix; mix++)	{

		// to be written:
		// reallocation move

		UpdateOmegaSuffStatLogProb();

	}

	UpdateOmega();
	return 1.0;
}

double SBDPOmegaProcess::GlobalMixMoveOmega(int nmix, double tuning, int nsiterep, int nhyperrep)	{

	MESSAGE signal = MIXMOVEOMEGA;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&nmix,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&nsiterep,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&tuning,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

	int* tmp = new int[GetNsite()];

	for (int mix=0; mix<nmix; mix++)	{

		// slave move site allocations

		// !!!
		// collect site allocations
		// in fact, could be done only after the mix/nmix loop
		MPI_Status stat;
		for(int i=1; i<GetNprocs(); i++) {
			MPI_Recv(tmp,GetNsite(),MPI_INT,i,TAG1,MPI_COMM_WORLD,&stat);
			for(int j=GetProcSiteMin(i); j<GetProcSiteMax(i); ++j) {
				if (ActiveSite(j))	{
					omegaalloc[j] = tmp[j];
				}
			}
		}

		GlobalUpdateOmegaSuffStat()

		// (1) pure MH
		MoveOmegas(tuning, nsiterep);

		// (2)
		// resample omegas (Gibbs)

		// prior
		// omega[k] ~ gamma(omegaalpha, omegabeta)
		// p(omega[k]) ~ omega[k]^(omegaalpha-1) * exp(-omegabeta*omega[k])

		// likelihood
		// p(D | omega[k]) ~ omega[k]^compomegasuffstatcount[k] * exp(- compomegasuffstatbeta[k]*omega[k])

		// posterior gamma
		// p(omega[k] | D) ~ omega[k]^(omegaalpha + compomegasuffstatcount[k] -1) * exp(- (omegabeta + compomegasuffstatbeta[k]) * omega[k])
		// omega[k] | D ~ gamma(omegaalpha + compsuffstatcount[k], omegabeta + compsuffstatbeta[k])

		// move hyperparam
		MoveOmegaHyper(tuning,nhyperrep);
		MoveOmegaHyper(0.1*tuning,nhyperrep);

		// (3)
		// ...

		// send omegas

		// send hyperparam
		// MPI_Bcast(&omegaalpha,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		// MPI_Bcast(&omegabeta,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	}
	// here, gather allocations
	// and broadcast them across slaves

	delete[] tmp;

	UpdateOmega();

	return 1.0;
}

void SBDPOmegaProcess::SlaveMixMoveOmega()	{

	double tuning = 1;
	int nrep = 1;
	int nmix = 1;
	MPI_Bcast(&nmix,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&nrep,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&tuning,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

	for (int mix=0; mix<nmix; mix++)	{

		// should move only those components that have been assigned as a job to that slave
		MoveOmegas(tuning, nrep);
		MPI_Send(omegaalloc,GetNsite(),MPI_INT,0,TAG1,MPI_COMM_WORLD);

		// get hyperparam
		MPI_Bcast(&omegaalpha,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(&omegabeta,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	}

	UpdateOmega();
}
*/


