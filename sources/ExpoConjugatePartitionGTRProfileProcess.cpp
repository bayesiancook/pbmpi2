
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "ExpoConjugatePartitionGTRProfileProcess.h"
#include "Random.h"
#include "Parallel.h"
	
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* ExpoConjugatePartitionGTRProfileProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------


void ExpoConjugatePartitionGTRProfileProcess::Create()	{
	if (! rrsuffstatcount)	{
		PartitionGTRProfileProcess::Create();
		allocrrsuffstatcount = new int[Nrr*Npart];
		allocrrsuffstatbeta = new double[Nrr*Npart];
		rrsuffstatcount = new int*[Npart];
		rrsuffstatbeta = new double*[Npart];
		for (int part=0; part<Npart; part++)	{
			rrsuffstatcount[part] = allocrrsuffstatcount + part*Nrr;
			rrsuffstatbeta[part] = allocrrsuffstatbeta + part*Nrr;
		}
	}
}

void ExpoConjugatePartitionGTRProfileProcess::Delete()	{
	if (rrsuffstatcount)	{
		delete[] rrsuffstatcount;
		delete[] rrsuffstatbeta;
		delete[] allocrrsuffstatcount;
		delete[] allocrrsuffstatbeta;
		rrsuffstatcount = 0;
		rrsuffstatbeta = 0;
		PartitionGTRProfileProcess::Delete();
	}
}

double ExpoConjugatePartitionGTRProfileProcess::MoveRR()	{
	GlobalUpdateRRSuffStat();
	for (int part=0; part<Npart; part++)	{
		for (int i=0; i<GetNrr(); i++)	{
			rr[part][i] = rnd::GetRandom().Gamma(1.0 + rrsuffstatcount[part][i], 1.0 + rrsuffstatbeta[part][i]);
		}
	}
	return 1;
}

void ExpoConjugatePartitionGTRProfileProcess::GlobalUpdateRRSuffStat()	{

	if (GetNprocs() > 1)	{
	MPI_Status stat;
	MESSAGE signal = UPDATE_RRATE;

	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	for (int part=0; part<Npart; part++)	{
		for(int i=0; i<GetNrr(); i++) {
			rrsuffstatcount[part][i] = 0;
			rrsuffstatbeta[part][i] = 0.0;
		}
	}

	int ivector[GetNrr()*Npart];
	double dvector[GetNrr()*Npart];
	for(int i=1; i<GetNprocs(); i++) {
		MPI_Recv(ivector,GetNrr()*Npart,MPI_INT,i,TAG1,MPI_COMM_WORLD,&stat);
		for(int j=0; j<GetNrr()*Npart; j++) {
			allocrrsuffstatcount[j] += ivector[j];
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	for(int i=1; i<GetNprocs(); i++) {
		MPI_Recv(dvector,GetNrr()*Npart,MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
		for(int j=0; j<GetNrr()*Npart; j++) {
			allocrrsuffstatbeta[j] += dvector[j];
		}
	}

	MPI_Bcast(allocrrsuffstatcount,GetNrr()*Npart,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(allocrrsuffstatbeta,GetNrr()*Npart,MPI_DOUBLE,0,MPI_COMM_WORLD);
	}
	else	{
		UpdateRRSuffStat();
	}
}

void ExpoConjugatePartitionGTRProfileProcess::SlaveUpdateRRSuffStat()	{

	UpdateRRSuffStat();
	MPI_Send(allocrrsuffstatcount,GetNrr()*Npart,MPI_INT,0,TAG1,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Send(allocrrsuffstatbeta,GetNrr()*Npart,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);

	MPI_Bcast(allocrrsuffstatcount,GetNrr()*Npart,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(allocrrsuffstatbeta,GetNrr()*Npart,MPI_DOUBLE,0,MPI_COMM_WORLD);
}
