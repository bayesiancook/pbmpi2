
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
#include <string.h>

#include "MultiGeneRateProcess.h"


void MultiGeneRateProcess::Create()	{

	DGamRateProcess::Create();
	if ((! genealpha) && (! GlobalAlpha()))	{
		genealpha = new double[Ngene];
		tmpgenealpha = new double[Ngene];
	}
}

void MultiGeneRateProcess::Delete()	{

	if (genealpha)	{
		delete[] genealpha;
		delete[] tmpgenealpha;
		genealpha = 0;
		tmpgenealpha = 0;
	}
	DGamRateProcess::Delete();
}

void MultiGeneRateProcess::SampleRate()	{
	if (GlobalAlpha())	{
		DGamRateProcess::SampleRate();
	}
	else	{
		meanalpha = 1;
		varalpha = 1;
		for (int gene=0; gene<Ngene; gene++)	{
			genealpha[gene] = rnd::GetRandom().Gamma(3,3);
		}
	}
}

void MultiGeneRateProcess::PriorSampleRate()	{
	if (GlobalAlpha())	{
		DGamRateProcess::PriorSampleRate();
	}
	else	{
		meanalpha = rnd::GetRandom().sExpo();
		varalpha = rnd::GetRandom().sExpo();
		double a = meanalpha * meanalpha / varalpha;
		double b = meanalpha / varalpha;
		for (int gene=0; gene<Ngene; gene++)	{
			genealpha[gene] = rnd::GetRandom().Gamma(a,b);
		}
	}
}

double MultiGeneRateProcess::GlobalGetMeanAlpha()	{

	MESSAGE signal = MEANALPHA;
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

void MultiGeneRateProcess::SlaveGetMeanAlpha() {
	double tot = 0;
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			tot += GetRateProcess(gene)->GetAlpha();
		}
	}
	MPI_Send(&tot,1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

void MultiGeneRateProcess::GlobalCollectGeneAlphas()	{

	MESSAGE signal = COLLECTALPHA;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Status stat;
	for(int i=1; i<GetNprocs(); i++) {
		MPI_Recv(tmpgenealpha,Ngene,MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
		for (int gene=0; gene<Ngene; gene++)	{
			if (genealloc[gene] == myid)	{
				genealpha[gene] = tmpgenealpha[gene];
			}
		}
	}
}

void MultiGeneRateProcess::SlaveCollectGeneAlphas() {
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			tmpgenealpha[gene] = GetRateProcess(gene)->GetAlpha();
		}
	}
	MPI_Send(tmpgenealpha,Ngene,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

double MultiGeneRateProcess::Move(double tuning, int nrep)	{

	double ret = 0;
	chronorate.Start();
	if (globalalpha)	{
		ret = MoveAlpha(tuning, nrep);
	}
	else	{
		GlobalCollectGeneAlphas();
		ret += MoveGeneRateHyperParams(tuning,nrep);
		ret += MoveGeneRateHyperParams(0.1*tuning,nrep);
	}
	chronorate.Stop();
}

double MultiGeneRateProcess::MoveGeneRateHyperParams(double tuning, int nrep)	{

	double nacc = 0;
	for (int rep=0; rep<nrep; rep++)	{

		double bkmeanalpha = meanalpha;
		double bkvaralpha = varalpha;

		double deltalogprob = - LogRatePrior();
		double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
		double e = exp(m);
		int choose = (int) (2 * rnd::GetRandom().Uniform());
		if (choose)	{
			meanalpha *= e;
		}
		else	{
			varalpha *= e;
		}
		deltalogprob += LogRatePrior();
		deltalogprob += m;
		int accepted = (log(rnd::GetRandom().Uniform()) < deltalogprob);

		if (accepted)	{
			nacc++;
		}
		else	{
			meanalpha = bkmeanalpha;
			varalpha = bkvaralpha;
		}
	}
	return ((double) nacc) / nrep;
}


double MultiGeneRateProcess::LogRatePrior()	{
	if (globalalpha)	{
		return DGamRateProcess::LogRatePrior();
	}
	double tot = -meanalpha - varalpha;
	double a = meanalpha * meanalpha / varalpha;
	double b = meanalpha / varalpha;
	for (int gene=0; gene<Ngene; gene++)	{
		tot += a*log(b) - rnd::GetRandom().logGamma(a) + (a-1)*log(genealpha[gene]) - b*genealpha[gene];
	}
	return tot;
}

void MultiGeneRateProcess::SlaveUpdateSiteRateSuffStat()	{

	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->UpdateSiteRateSuffStat();
		}
	}
}


void MultiGeneRateProcess::UpdateRateSuffStat() {

	for(int i=0; i<GetNcat(); i++) {
		ratesuffstatcount[i] = 0;
		ratesuffstatbeta[i] = 0.0;
	}

	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			GetRateProcess(gene)->UpdateRateSuffStat();
			const int* count = GetRateProcess(gene)->GetRateSuffStatCount();
			const double* beta = GetRateProcess(gene)->GetRateSuffStatBeta();
			for(int i=0; i<GetNcat(); i++) {
				ratesuffstatcount[i] += count[i];
				ratesuffstatbeta[i] += beta[i];
			}
		}
	}
}
