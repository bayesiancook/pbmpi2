#include "AASubSelMixtureProfileProcess.h"
#include "Parallel.h"

void AASubSelMixtureProfileProcess::Create()    {

    GTRMixtureProfileProcess::Create();
    if (! rrsuffstatcount)  {
        rrsuffstatcount = new double[Nrr];
        rrsuffstatbeta = new double[Nrr];
    }

}

void AASubSelMixtureProfileProcess::Delete()    {

    if (rrsuffstatcount)   {
		delete[] rrsuffstatcount;
		delete[] rrsuffstatbeta;
		rrsuffstatcount = 0;
		rrsuffstatbeta = 0;
    }
    GTRMixtureProfileProcess::Delete();
}

double AASubSelMixtureProfileProcess::MoveRR()	{
	GlobalUpdateRRSuffStat();
	for (int i=0; i<GetNrr(); i++)	{
		rr[i] = rnd::GetRandom().Gamma(1.0 + rrsuffstatcount[i], 1.0 + rrsuffstatbeta[i]);
	}
	return 1;
}

void AASubSelMixtureProfileProcess::GlobalUpdateRRSuffStat()	{

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

void AASubSelMixtureProfileProcess::SlaveUpdateRRSuffStat()	{

	UpdateRRSuffStat();
	MPI_Send(rrsuffstatcount,GetNrr(),MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Send(rrsuffstatbeta,GetNrr(),MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);

	MPI_Bcast(rrsuffstatcount,GetNrr(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(rrsuffstatbeta,GetNrr(),MPI_DOUBLE,0,MPI_COMM_WORLD);
}
