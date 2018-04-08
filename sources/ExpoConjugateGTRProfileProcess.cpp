
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "ExpoConjugateGTRProfileProcess.h"
#include "Random.h"
#include "Parallel.h"
	
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* ExpoConjugateGTRProfileProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------


void ExpoConjugateGTRProfileProcess::Create()	{
	if (! rr)	{
		GTRProfileProcess::Create();
	}
	if (! rrsuffstatcount)	{
		rrsuffstatcount = new double[Nrr];
		rrsuffstatbeta = new double[Nrr];
	}
}

void ExpoConjugateGTRProfileProcess::Delete()	{
	if (rrsuffstatcount)	{
		delete[] rrsuffstatcount;
		delete[] rrsuffstatbeta;
		rrsuffstatcount = 0;
		rrsuffstatbeta = 0;
    }
    if (!rr)    {
		GTRProfileProcess::Delete();
	}
}

double ExpoConjugateGTRProfileProcess::MoveRR()	{
	GlobalUpdateRRSuffStat();
	for (int i=0; i<GetNrr(); i++)	{
		rr[i] = rnd::GetRandom().Gamma(1.0 + rrsuffstatcount[i], 1.0 + rrsuffstatbeta[i]);
	}
	return 1;
}

void ExpoConjugateGTRProfileProcess::GlobalUpdateRRSuffStat()	{

	if (GetNprocs() > 1)	{
	MPI_Status stat;
	MESSAGE signal = UPDATE_RRATE;

	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	for(int i=0; i<GetNrr(); i++) {
		rrsuffstatcount[i] = 0;
		rrsuffstatbeta[i] = 0.0;
	}

	double ivector[GetNrr()];
	double dvector[GetNrr()];
	for(int i=1; i<GetNprocs(); i++) {
		MPI_Recv(ivector,GetNrr(),MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
		for(int j=0; j<GetNrr(); j++) {
			rrsuffstatcount[j] += ivector[j];
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	for(int i=1; i<GetNprocs(); i++) {
		MPI_Recv(dvector,GetNrr(),MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
		for(int j=0; j<GetNrr(); j++) {
			rrsuffstatbeta[j] += dvector[j];
		}
	}

	MPI_Bcast(rrsuffstatcount,GetNrr(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(rrsuffstatbeta,GetNrr(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	}
	else	{
		UpdateRRSuffStat();
	}
}

void ExpoConjugateGTRProfileProcess::SlaveUpdateRRSuffStat()	{

	UpdateRRSuffStat();
	MPI_Send(rrsuffstatcount,GetNrr(),MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Send(rrsuffstatbeta,GetNrr(),MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);

	MPI_Bcast(rrsuffstatcount,GetNrr(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(rrsuffstatbeta,GetNrr(),MPI_DOUBLE,0,MPI_COMM_WORLD);
}
