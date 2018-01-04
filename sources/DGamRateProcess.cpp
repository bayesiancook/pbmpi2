
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "DGamRateProcess.h"
#include "Random.h"
#include "IncompleteGamma.h"

#include "Parallel.h"
#include <string.h>

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* DGamRateProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------


void DGamRateProcess::Create()	{
	if (! rate)	{
		if (! GetNcat())	{
			cerr << "error: DGamRateProcess::Ncat has not been initialized\n";
			exit(1);
		}
		if (withpinv)	{
			cerr << "number of rate categories: " << Ncat << " + 1 for invariable sites\n";
			Ncat++;
		}
		RateProcess::Create();
		rate = new double[GetNcat()];
		ratesuffstatcount = new double[GetNcat()];
		ratesuffstatbeta = new double[GetNcat()];
	}
}

void DGamRateProcess::Delete() 	{
	delete[] rate;
	delete[] ratesuffstatcount;
	delete[] ratesuffstatbeta;
	rate = 0;
	ratesuffstatcount = 0;
	ratesuffstatbeta = 0;
	RateProcess::Delete();
}

void DGamRateProcess::ToStream(ostream& os)	{
	os << alpha << '\n';
	if (withpinv)	{
		os << pinv << '\n';
	}
}

void DGamRateProcess::FromStream(istream& is)	{
	double tmpalpha, tmppinv;
	is >> tmpalpha;
	tmppinv = 0.1;
	if (withpinv)	{
		is >> tmppinv;
	}
	SetRateParams(tmpalpha,tmppinv);
}

void DGamRateProcess::UpdateDiscreteCategories()	{

	if (withpinv)	{
		double* x = new double[GetNcat()-1];
		double* y = new double[GetNcat()-1];
		double lg = rnd::GetRandom().logGamma(alpha+1.0);
		for (int i=0; i<GetNcat()-1; i++)	{
			x[i] = PointGamma((i+1.0)/(GetNcat()-1),alpha,alpha);
		}
		for (int i=0; i<GetNcat()-2; i++)	{
			y[i] = IncompleteGamma(alpha*x[i],alpha+1,lg);
		}
		y[GetNcat()-2] = 1.0;

		rate[0] = 0;
	
		rate[1] = (GetNcat()-1) * y[0];
		for (int i=1; i<(GetNcat()-1); i++)	{
			rate[i+1] = (GetNcat()-1) * (y[i] - y[i-1]);
		}
		delete[] x;
		delete[] y;
	}
	else	{
		double* x = new double[GetNcat()];
		double* y = new double[GetNcat()];
		double lg = rnd::GetRandom().logGamma(alpha+1.0);
		for (int i=0; i<GetNcat(); i++)	{
			x[i] = PointGamma((i+1.0)/GetNcat(),alpha,alpha);
		}
		for (int i=0; i<GetNcat()-1; i++)	{
			y[i] = IncompleteGamma(alpha*x[i],alpha+1,lg);
		}
		y[GetNcat()-1] = 1.0;

		rate[0] = GetNcat() * y[0];
		for (int i=1; i<GetNcat(); i++)	{
			rate[i] = GetNcat() * (y[i] - y[i-1]);
		}
		delete[] x;
		delete[] y;
	}
}

void DGamRateProcess::SampleRate()	{
	if (! FixAlpha())	{
		// SetAlpha(1.0);
		alpha = 1.0;
	}
	if (withpinv  && (! FixPinv()))	{
		pinv = 0.01;
	}
	UpdateDiscreteCategories();
}

void DGamRateProcess::PriorSampleRate()	{
	double tmpalpha = alpha;
	if (! FixAlpha())	{
		double a = meanalpha * meanalpha / varalpha;
		double b = meanalpha / varalpha;
		int count = 0;
		while ((count < 1000) && (tmpalpha < alphamin))	{
			tmpalpha = rnd::GetRandom().Gamma(a,b);
			count++;
		}
		if (count == 1000)	{
			cerr << "error in DGamRateProcess::PriorSampleRate\n";
			exit(1);
		}
	}
	double tmppinv = pinv;
	if (withpinv && (! FixPinv()))	{
		double a = meanpinv * invconcpinv;
		double b = (1-meanpinv) * invconcpinv;
		double x = rnd::GetRandom().sGamma(a);
		double y = rnd::GetRandom().sGamma(b);
		tmppinv = x / (x+y);
	}
	SetRateParams(tmpalpha,tmppinv);
}

double DGamRateProcess::LogRatePrior()	{
	double a = meanalpha * meanalpha / varalpha;
	double b = meanalpha / varalpha;
	double ret = a*log(b) - rnd::GetRandom().logGamma(a) + (a-1)*log(alpha) - b*alpha;

	if (withpinv)	{
		double a = meanpinv * invconcpinv;
		double b = (1-meanpinv) * invconcpinv;
		ret += rnd::GetRandom().logGamma(a+b) - rnd::GetRandom().logGamma(a) - rnd::GetRandom().logGamma(b) + (a-1)*log(pinv) + (b-1)*log(1-pinv);
	}

	if (isnan(ret))	{
		cerr << "in DGamRateProcess::LogRatePrior: nan\n";
		exit(1);
	}
	if (isinf(ret))	{
		cerr << "in DGamRateProcess::LogRatePrior: inf\n";
		cerr << meanpinv << '\t' << invconcpinv << '\t' << pinv << '\n';
		exit(1);
	}
	return ret;
}

double DGamRateProcess::RateSuffStatLogProb()	{
	double total = 0;
	for (int k=0; k<GetNcat(); k++)	{
		if (ratesuffstatcount[k] > 0)	{
			total += ratesuffstatcount[k] * log(rate[k]);
		}
		total -= ratesuffstatbeta[k] * rate[k];
	}
	if (isnan(total))	{
		cerr << "in DGamRateProcess::RateSuffStatLogProb: nan log prob\n";
		for (int k=0; k<GetNcat(); k++)	{
			cerr << rate[k] << '\t' << log(rate[k]) << '\t' << ratesuffstatcount[k] << '\t' << ratesuffstatbeta[k] << '\t' << ratesuffstatcount[k] * log(rate[k]) << '\n';
		}
		exit(1);
	}
	if (isinf(total))	{
		cerr << "in DGamRateProcess::RateSuffStatLogProb: inf log prob\n";
		for (int k=0; k<GetNcat(); k++)	{
			cerr << rate[k] << '\t' << log(rate[k]) << '\t' << ratesuffstatcount[k] << '\t' << ratesuffstatbeta[k] << '\t' << ratesuffstatcount[k] * log(rate[k]) << '\n';
		}
		exit(1);
	}
	return total;
}

double DGamRateProcess::MoveRateParams(double tuning, int nrep)	{

	GlobalUpdateRateSuffStat();
	double naccepted = 0;
	for (int rep=0; rep<nrep; rep++)	{
		naccepted += MoveAlpha(tuning);
		if (withpinv)	{
			naccepted += MovePinv(tuning);
		}
	}
	return naccepted / nrep;
}

double DGamRateProcess::NonMPIMoveRateParams(double tuning, int nrep)	{

	UpdateRateSuffStat();
	double naccepted = 0;
	for (int rep=0; rep<nrep; rep++)	{
		naccepted += MoveAlpha(tuning);
		if (withpinv)	{
			naccepted += MovePinv(tuning);
		}
	}
	return naccepted / nrep;
}

double DGamRateProcess::MoveAlpha(double tuning)	{

	double bkalpha = alpha;
	double bkpinv = pinv;
	double deltalogprob = -LogRatePrior() - RateSuffStatLogProb();
	double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
	double e = exp(m);
	double newalpha = alpha * e;
	SetRateParams(newalpha,bkpinv);
	deltalogprob += m + LogRatePrior() + RateSuffStatLogProb();
	int accepted = (log(rnd::GetRandom().Uniform()) < deltalogprob);
	if (alpha < alphamin)	{
		accepted = 0;
	}
	if (!accepted)	{
		SetRateParams(bkalpha,bkpinv);
	}
	return ((double) accepted);
}

double DGamRateProcess::MovePinv(double tuning)	{

	double a = meanpinv * invconcpinv + Ninv;
	double b = (1-meanpinv) * invconcpinv + GetNsite() - Ninv;
	double x = rnd::GetRandom().sGamma(a);
	double y = rnd::GetRandom().sGamma(b);
	double bkalpha = alpha;
	double tmppinv = x / (x+y);
	SetRateParams(bkalpha,tmppinv);
	return 1.0;
}


void DGamRateProcess::GlobalUpdateRateSuffStat()	{

	if (GetNprocs() > 1)	{
	// MPI2
	// should ask the slaves to call their UpdateRateSuffStat
	// and then gather the statistics;
	MPI_Status stat;
	MESSAGE signal = UPDATE_RATE;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	for(int i=0; i<GetNcat(); i++) {
		ratesuffstatcount[i] = 0;
		ratesuffstatbeta[i] = 0.0;
	}
	double ivector[GetNcat()];
	double dvector[GetNcat()];
        for(int i=1; i<GetNprocs(); i++) {
                MPI_Recv(ivector,GetNcat(),MPI_DOUBLE,MPI_ANY_SOURCE,TAG1,MPI_COMM_WORLD,&stat);
                for(int j=0; j<GetNcat(); j++) {
                        ratesuffstatcount[j] += ivector[j];                      
                }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        for(int i=1; i<GetNprocs(); i++) {
                MPI_Recv(dvector,GetNcat(),MPI_DOUBLE,MPI_ANY_SOURCE,TAG1,MPI_COMM_WORLD,&stat);
                for(int j=0; j<GetNcat(); j++) {
                        ratesuffstatbeta[j] += dvector[j]; 
                }
        }
	if (withpinv)	{
		Ninv = 0;
		MPI_Barrier(MPI_COMM_WORLD);
		double tmp;
		for(int i=1; i<GetNprocs(); i++) {
			MPI_Recv(&tmp,1,MPI_DOUBLE,MPI_ANY_SOURCE,TAG1,MPI_COMM_WORLD,&stat);
			Ninv += tmp;
		}
	}
	}
	else	{
		UpdateRateSuffStat();
	}
	/*
	cerr << "update rate suff stat\n";
	for (int i=0; i<GetNcat(); i++)	{
		cerr << ratesuffstatcount[i] << '\t' << ratesuffstatbeta[i] << '\n';
	}
	*/
}

void DGamRateProcess::UpdateRateSuffStat()	{

	for (int i=0; i<GetNcat(); i++)	{
		ratesuffstatcount[i] = 0;
		ratesuffstatbeta[i] = 0.0;
	}
	Ninv = 0;
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i))	{
			ratesuffstatcount[ratealloc[i]] += GetSiteRateSuffStatCount(i);
			ratesuffstatbeta[ratealloc[i]] += GetSiteRateSuffStatBeta(i);
			if (withpinv)	{
				if (ratealloc[i] == 0)	{
					Ninv++;
				}
			}
		}
	}

}	

void DGamRateProcess::SlaveUpdateRateSuffStat()	{

	UpdateRateSuffStat();

	MPI_Send(ratesuffstatcount,GetNcat(),MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Send(ratesuffstatbeta,GetNcat(),MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
	if (withpinv)	{
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Send(&Ninv,1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
	}
}	
