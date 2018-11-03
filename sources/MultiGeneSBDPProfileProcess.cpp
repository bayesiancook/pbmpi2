
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include <cassert>

#include "Parallel.h"
#include "MultiGeneSBDPProfileProcess.h"

void MultiGeneSBDPProfileProcess::Create()	{

	if (! kappaarray)	{
		kappaarray = new double[Ngene];
		statalphaarray = new double[Ngene];
        nmodearray = new double[Ngene];
		MultiGeneProfileProcess::Create();
		SBDPProfileProcess::Create();
	}
}

void MultiGeneSBDPProfileProcess::Delete()	{

	if (kappaarray)	{
		SBDPProfileProcess::Delete();
		MultiGeneProfileProcess::Delete();
		delete[] kappaarray;
		kappaarray = 0;
		delete[] statalphaarray;
		statalphaarray = 0;
        delete[] nmodearray;
        nmodearray = 0;
	}
}

double MultiGeneSBDPProfileProcess::LogHyperPrior()	{

	if (kappaprior != 2)	{
		return 0;
	}
	// assumes that kappas have already been collected

	double beta = 1.0 / kappamean;
	double alpha = 1.0 / kapparelvar;
	double total = Ngene * (alpha*log(beta) - rnd::GetRandom().logGamma(alpha));
	for (int gene=0; gene<Ngene; gene++)	{
		total += (alpha-1)*log(kappaarray[gene]) - beta*kappaarray[gene];
	}

	// prior mean ? 
	total -= kappamean / 10;

	// prior relvar:
	total -= kapparelvar;
	return total;
}

void MultiGeneSBDPProfileProcess::SampleHyper()	{

	if (kappaprior == 2)	{
		kappamean = 10 * rnd::GetRandom().sExpo();
		kapparelvar = rnd::GetRandom().sExpo();
	}
}

void MultiGeneSBDPProfileProcess::PriorSampleHyper()	{

	if (kappaprior == 2)	{
		kappamean = 10 * rnd::GetRandom().sExpo();
		kapparelvar = rnd::GetRandom().sExpo();
	}
}

double MultiGeneSBDPProfileProcess::MoveHyper(double tuning, int nrep)	{

	if (kappaprior == 2)	{
		GlobalCollectKappas();
		for (int rep=0; rep<nrep; rep++)	{
			MoveKappaMean(tuning);
			MoveKappaRelVar(tuning);
		}
		GlobalUpdateParameters();
	}
	return 1.0;
}

double MultiGeneSBDPProfileProcess::MoveKappaMean(double tuning)	{

	double deltalogprob = - LogHyperPrior();
	double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
	double e = exp(m);
	kappamean *= e;
	deltalogprob += LogHyperPrior();
	deltalogprob += m;
	int accepted = (log(rnd::GetRandom().Uniform()) < deltalogprob);
	if (!accepted)	{
		kappamean /= e;
	}
	return ((double) accepted);
}

double MultiGeneSBDPProfileProcess::MoveKappaRelVar(double tuning)	{

	double deltalogprob = - LogHyperPrior();
	double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
	double e = exp(m);
	kapparelvar *= e;
	deltalogprob += LogHyperPrior();
	deltalogprob += m;
	int accepted = (log(rnd::GetRandom().Uniform()) < deltalogprob);
	if (!accepted)	{
		kapparelvar /= e;
	}
	return ((double) accepted);
}

void MultiGeneSBDPProfileProcess::GlobalCollectKappas()	{

	MESSAGE signal = COLLECTKAPPAS;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Status stat;
	double* tmpkappaarray = new double[Ngene];
	for(int i=1; i<GetNprocs(); i++) {
		MPI_Recv(tmpkappaarray,Ngene,MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
		for (int gene=0; gene<Ngene; gene++)	{
			if (genealloc[gene] == i)	{
				kappaarray[gene] = tmpkappaarray[gene];
			}
		}
	}

    MPI_Barrier(MPI_COMM_WORLD);

	for(int i=1; i<GetNprocs(); i++) {
		MPI_Recv(tmpkappaarray,Ngene,MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
		for (int gene=0; gene<Ngene; gene++)	{
			if (genealloc[gene] == i)	{
				statalphaarray[gene] = tmpkappaarray[gene];
			}
		}
	}

    MPI_Barrier(MPI_COMM_WORLD);

	for(int i=1; i<GetNprocs(); i++) {
		MPI_Recv(tmpkappaarray,Ngene,MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
		for (int gene=0; gene<Ngene; gene++)	{
			if (genealloc[gene] == i)	{
				nmodearray[gene] = tmpkappaarray[gene];
			}
		}
	}
	delete[] tmpkappaarray;
}

void MultiGeneSBDPProfileProcess::SlaveCollectKappas()	{

	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			kappaarray[gene] = GetSBDPProfileProcess(gene)->GetKappa();
			statalphaarray[gene] = GetSBDPProfileProcess(gene)->GetMeanDirWeight();
			nmodearray[gene] = GetSBDPProfileProcess(gene)->GetNDisplayedComponent();
		}
	}
	MPI_Send(kappaarray,Ngene,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

	MPI_Send(statalphaarray,Ngene,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

	MPI_Send(nmodearray,Ngene,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

double MultiGeneSBDPProfileProcess::GlobalGetMeanNcomponent()	{

	MESSAGE signal = MEANNCOMPONENT;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	int ncomp = 0;
	MPI_Status stat;
	for(int i=1; i<GetNprocs(); i++) {
		int tmp;
		MPI_Recv(&tmp,1,MPI_INT,i,TAG1,MPI_COMM_WORLD,&stat);
		ncomp += tmp;
	}
	return ((double) ncomp) / Ngene;
}

double MultiGeneSBDPProfileProcess::GlobalGetMeanStatEnt()	{

	MESSAGE signal = MEANSTATENT;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	double tot = 0;
	MPI_Status stat;
	for(int i=1; i<GetNprocs(); i++) {
		double tmp;
		MPI_Recv(&tmp,1,MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
		tot += tmp;
	}
	return tot / Ngene;
}

double MultiGeneSBDPProfileProcess::GlobalGetMeanStatAlpha()	{

	MESSAGE signal = MEANSTATALPHA;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	double tot = 0;
	MPI_Status stat;
	for(int i=1; i<GetNprocs(); i++) {
		double tmp;
		MPI_Recv(&tmp,1,MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
		tot += tmp;
	}
	return tot / Ngene;
}

double MultiGeneSBDPProfileProcess::GlobalGetMeanKappa()	{

	MESSAGE signal = MEANKAPPA;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	double tot = 0;
	MPI_Status stat;
	for(int i=1; i<GetNprocs(); i++) {
		double tmp;
		MPI_Recv(&tmp,1,MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
		tot += tmp;
	}
	return tot / Ngene;
}

void MultiGeneSBDPProfileProcess::SlaveGetMeanNcomponent() {
	int ncomp = 0;
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			ncomp += GetSBDPProfileProcess(gene)->GetNDisplayedComponent();
		}
	}
	MPI_Send(&ncomp,1,MPI_INT,0,TAG1,MPI_COMM_WORLD);
}

void MultiGeneSBDPProfileProcess::SlaveGetMeanStatEnt() {
	double tot = 0;
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			tot += GetSBDPProfileProcess(gene)->GetStatEnt();
		}
	}
	MPI_Send(&tot,1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

void MultiGeneSBDPProfileProcess::SlaveGetMeanStatAlpha() {
	double tot = 0;
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			tot += GetSBDPProfileProcess(gene)->GetMeanDirWeight();
		}
	}
	MPI_Send(&tot,1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

void MultiGeneSBDPProfileProcess::SlaveGetMeanKappa() {
	double tot = 0;
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			tot += GetSBDPProfileProcess(gene)->GetKappa();
		}
	}
	MPI_Send(&tot,1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

