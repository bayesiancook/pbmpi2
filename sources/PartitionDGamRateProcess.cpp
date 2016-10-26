
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "PartitionDGamRateProcess.h"
#include "Random.h"
#include "IncompleteGamma.h"

#include "Parallel.h"
#include <string.h>

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* PartitionDGamRateProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------


void PartitionDGamRateProcess::Create()	{
	if (! rate)	{
		if (! GetNcat())	{
			cerr << "error: PartitionDGamRateProcess::Ncat has not been initialized\n";
			exit(1);
		}
		if (withpinv)	{
			cerr << "number of rate categories: " << Ncat << " + 1 for invariable sites\n";
			Ncat++;
		}
		RateProcess::Create();
		rate = new double*[Npart];
		allocratesuffstatcount = new int[Npart*GetNcat()];
		allocratesuffstatbeta = new double[Npart*GetNcat()];
		ratesuffstatcount = new int*[Npart];
		ratesuffstatbeta = new double*[Npart];
		for (int part=0; part<Npart; part++)	{
			rate[part] = new double[GetNcat()];
			ratesuffstatcount[part] = allocratesuffstatcount + part*GetNcat();
			ratesuffstatbeta[part] = allocratesuffstatbeta + part*GetNcat();
		}
		Ninv = new int[Npart];
	}
}

void PartitionDGamRateProcess::Delete() 	{
	delete[] Ninv;
	for (int part=0; part<Npart; part++)	{
		delete[] rate[part];
	}
	delete[] rate;
	delete[] ratesuffstatcount;
	delete[] ratesuffstatbeta;
	delete[] allocratesuffstatcount;
	delete[] allocratesuffstatbeta;
	rate = 0;
	allocratesuffstatcount = 0;
	allocratesuffstatbeta = 0;
	ratesuffstatcount = 0;
	ratesuffstatbeta = 0;
	RateProcess::Delete();
}

void PartitionDGamRateProcess::ToStream(ostream& os)	{
	for (int part=0; part<Npart; part++)	{
		os << alpha[part] << '\t';
	}
	os << '\n';
	if (withpinv)	{
		for (int part=0; part<Npart; part++)	{
			os << pinv[part] << '\t';
		}
		os << '\n';
	}
}

void PartitionDGamRateProcess::FromStream(istream& is)	{
	double* tmpalpha = new double[Npart];
	double* tmppinv = new double[Npart];
	for (int part=0; part<Npart; part++)	{
		is >> tmpalpha[part];
	}
	for (int part=0; part<Npart; part++)	{
		tmppinv[part] = 0.1;
	}
	if (withpinv)	{
		for (int part=0; part<Npart; part++)	{
			is >> tmppinv[part];
		}
	}
	SetRateParams(tmpalpha,tmppinv);
	delete[] tmpalpha;
	delete[] tmppinv;
}

void PartitionDGamRateProcess::UpdateDiscreteCategories()	{

	for (int part=0; part<Npart; part++)	{
		UpdateDiscreteCategories(part);
	}
}

void PartitionDGamRateProcess::UpdateDiscreteCategories(int part)	{

	if (withpinv)	{
		double* x = new double[GetNcat()-1];
		double* y = new double[GetNcat()-1];
		double lg = rnd::GetRandom().logGamma(alpha[part]+1.0);
		for (int i=0; i<GetNcat()-1; i++)	{
			x[i] = PointGamma((i+1.0)/(GetNcat()-1),alpha[part],alpha[part]);
		}
		for (int i=0; i<GetNcat()-2; i++)	{
			y[i] = IncompleteGamma(alpha[part]*x[i],alpha[part]+1,lg);
		}
		y[GetNcat()-2] = 1.0;

		rate[part][0] = 0;
	
		rate[part][1] = (GetNcat()-1) * y[0];
		for (int i=1; i<(GetNcat()-1); i++)	{
			rate[part][i+1] = (GetNcat()-1) * (y[i] - y[i-1]);
		}
		delete[] x;
		delete[] y;
	}
	else	{
		double* x = new double[GetNcat()];
		double* y = new double[GetNcat()];
		double lg = rnd::GetRandom().logGamma(alpha[part]+1.0);
		for (int i=0; i<GetNcat(); i++)	{
			x[i] = PointGamma((i+1.0)/GetNcat(),alpha[part],alpha[part]);
		}
		for (int i=0; i<GetNcat()-1; i++)	{
			y[i] = IncompleteGamma(alpha[part]*x[i],alpha[part]+1,lg);
		}
		y[GetNcat()-1] = 1.0;

		rate[part][0] = GetNcat() * y[0];
		for (int i=1; i<GetNcat(); i++)	{
			rate[part][i] = GetNcat() * (y[i] - y[i-1]);
		}
		delete[] x;
		delete[] y;
	}
}

void PartitionDGamRateProcess::SampleRate()	{
	for (int part=0; part<Npart; part++)	{
		if (! FixAlpha())	{
			alpha[part] = 1.0;
		}
		if (withpinv  && (! FixPinv()))	{
			pinv[part] = 0.01;
		}
	}
	UpdateDiscreteCategories();
}

void PartitionDGamRateProcess::PriorSampleRate()	{
	
	double* tmpalpha = new double[Npart];
	double* tmppinv = new double[Npart];

	for (int part=0; part<Npart; part++)	{
		tmpalpha[part] = alpha[part];
		if (! FixAlpha())	{
			double a = meanalpha * meanalpha / varalpha;
			double b = meanalpha / varalpha;
			int count = 0;
			while ((count < 1000) && (tmpalpha[part] < alphamin))	{
				tmpalpha[part] = rnd::GetRandom().Gamma(a,b);
				count++;
			}
			if (count == 1000)	{
				cerr << "error in PartitionDGamRateProcess::PriorSampleRate\n";
				exit(1);
			}
		}
		tmppinv[part] = pinv[part];
		if (withpinv && (! FixPinv()))	{
			double a = meanpinv * invconcpinv;
			double b = (1-meanpinv) * invconcpinv;
			double x = rnd::GetRandom().sGamma(a);
			double y = rnd::GetRandom().sGamma(b);
			tmppinv[part] = x / (x+y);
		}
	}
	SetRateParams(tmpalpha,tmppinv);
	delete[] tmpalpha;
	delete[] tmppinv;
}

double PartitionDGamRateProcess::LogRatePrior()	{

	double tot = 0;
	for (int part=0; part<Npart; part++)	{
		tot += LogRatePrior(part);
	}
	return tot;
}

double PartitionDGamRateProcess::LogRatePrior(int part)	{
	double a = meanalpha * meanalpha / varalpha;
	double b = meanalpha / varalpha;
	double ret = a*log(b) - rnd::GetRandom().logGamma(a) + (a-1)*log(alpha[part]) - b*alpha[part];

	if (withpinv)	{
		double a = meanpinv * invconcpinv;
		double b = (1-meanpinv) * invconcpinv;
		ret += rnd::GetRandom().logGamma(a+b) - rnd::GetRandom().logGamma(a) - rnd::GetRandom().logGamma(b) + (a-1)*log(pinv[part]) + (b-1)*log(1-pinv[part]);
	}

	if (isnan(ret))	{
		cerr << "in PartitionDGamRateProcess::LogRatePrior: nan\n";
		exit(1);
	}
	if (isinf(ret))	{
		cerr << "in PartitionDGamRateProcess::LogRatePrior: inf\n";
		cerr << meanpinv << '\t' << invconcpinv << '\t' << pinv[part] << '\n';
		exit(1);
	}
	return ret;
}

double PartitionDGamRateProcess::RateSuffStatLogProb()	{

	double tot = 0;
	for (int part=0; part<Npart; part++)	{
		tot += RateSuffStatLogProb(part);
	}
	return tot;
}

double PartitionDGamRateProcess::RateSuffStatLogProb(int part)	{
	double total = 0;
	for (int k=0; k<GetNcat(); k++)	{
		if (ratesuffstatcount[part][k] > 0)	{
			total += ratesuffstatcount[part][k] * log(rate[part][k]);
		}
		total -= ratesuffstatbeta[part][k] * rate[part][k];
	}
	if (isnan(total))	{
		cerr << "in PartitionDGamRateProcess::RateSuffStatLogProb: nan log prob\n";
		for (int k=0; k<GetNcat(); k++)	{
			cerr << rate[part][k] << '\t' << log(rate[part][k]) << '\t' << ratesuffstatcount[part][k] << '\t' << ratesuffstatbeta[part][k] << '\t' << ratesuffstatcount[part][k] * log(rate[part][k]) << '\n';
		}
		exit(1);
	}
	if (isinf(total))	{
		cerr << "in PartitionDGamRateProcess::RateSuffStatLogProb: inf log prob\n";
		for (int k=0; k<GetNcat(); k++)	{
			cerr << rate[part][k] << '\t' << log(rate[part][k]) << '\t' << ratesuffstatcount[part][k] << '\t' << ratesuffstatbeta[part][k] << '\t' << ratesuffstatcount[part][k] * log(rate[part][k]) << '\n';
		}
		exit(1);
	}
	return total;
}

double PartitionDGamRateProcess::MoveRateParams(double tuning, int nrep)	{

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

double PartitionDGamRateProcess::NonMPIMoveRateParams(double tuning, int nrep)	{

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

double PartitionDGamRateProcess::MoveAlpha(double tuning)	{

	double nacc = 0;
	for (int part=0; part<Npart; part++)	{
		nacc += MoveAlpha(part,tuning);
	}
	return nacc / Npart;
}

double PartitionDGamRateProcess::MoveAlpha(int part, double tuning)	{

	double bkalpha = alpha[part];
	double bkpinv = pinv[part];
	double deltalogprob = -LogRatePrior(part) - RateSuffStatLogProb(part);
	double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
	double e = exp(m);
	double newalpha = alpha[part] * e;
	SetRateParams(part,newalpha,bkpinv);
	deltalogprob += m + LogRatePrior(part) + RateSuffStatLogProb(part);
	int accepted = (log(rnd::GetRandom().Uniform()) < deltalogprob);
	if (alpha[part] < alphamin)	{
		accepted = 0;
	}
	if (!accepted)	{
		SetRateParams(part,bkalpha,bkpinv);
	}
	return ((double) accepted);
}

double PartitionDGamRateProcess::MovePinv(double tuning)	{

	double nacc = 0;
	for (int part=0; part<Npart; part++)	{
		nacc += MovePinv(part,tuning);
	}
	return nacc / Npart;
}

double PartitionDGamRateProcess::MovePinv(int part, double tuning)	{

	double a = meanpinv * invconcpinv + Ninv[part];
	double b = (1-meanpinv) * invconcpinv + GetPartNsite(part) - Ninv[part];
	double x = rnd::GetRandom().sGamma(a);
	double y = rnd::GetRandom().sGamma(b);
	double bkalpha = alpha[part];
	double tmppinv = x / (x+y);
	SetRateParams(part,bkalpha,tmppinv);
	return 1.0;
}


void PartitionDGamRateProcess::GlobalUpdateRateSuffStat()	{

	if (GetNprocs() > 1)	{
	// MPI2
	// should ask the slaves to call their UpdateRateSuffStat
	// and then gather the statistics;
	MPI_Status stat;
	MESSAGE signal = UPDATE_RATE;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	for (int part=0; part<Npart; part++)	{
		for(int i=0; i<GetNcat(); i++) {
			ratesuffstatcount[part][i] = 0;
			ratesuffstatbeta[part][i] = 0.0;
		}
	}
	int ivector[Npart*GetNcat()];
	double dvector[Npart*GetNcat()];
        for(int i=1; i<GetNprocs(); i++) {
                MPI_Recv(ivector,Npart*GetNcat(),MPI_INT,MPI_ANY_SOURCE,TAG1,MPI_COMM_WORLD,&stat);
                for(int j=0; j<Npart*GetNcat(); j++) {
                        allocratesuffstatcount[j] += ivector[j];                      
                }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        for(int i=1; i<GetNprocs(); i++) {
                MPI_Recv(dvector,GetNcat(),MPI_DOUBLE,MPI_ANY_SOURCE,TAG1,MPI_COMM_WORLD,&stat);
                for(int j=0; j<Npart*GetNcat(); j++) {
                        allocratesuffstatbeta[j] += dvector[j]; 
                }
        }
	if (withpinv)	{
		for (int part=0; part<Npart; part++)	{
			Ninv[part] = 0;
		}
		MPI_Barrier(MPI_COMM_WORLD);
		int* tmp = new int[Npart];
		for(int i=1; i<GetNprocs(); i++) {
			MPI_Recv(tmp,Npart,MPI_INT,MPI_ANY_SOURCE,TAG1,MPI_COMM_WORLD,&stat);
			for (int part=0; part<Npart; part++)	{
				Ninv[part] += tmp[part];
			}
		}
	}
	}
	else	{
		UpdateRateSuffStat();
	}
}

void PartitionDGamRateProcess::UpdateRateSuffStat()	{

	for (int part=0; part<Npart; part++)	{
		for(int i=0; i<GetNcat(); i++) {
			ratesuffstatcount[part][i] = 0;
			ratesuffstatbeta[part][i] = 0.0;
		}
		Ninv[part] = 0;
	}
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i))	{
			ratesuffstatcount[partalloc[i]][ratealloc[i]] += GetSiteRateSuffStatCount(i);
			ratesuffstatbeta[partalloc[i]][ratealloc[i]] += GetSiteRateSuffStatBeta(i);
			if (withpinv)	{
				if (ratealloc[i] == 0)	{
					Ninv[partalloc[i]]++;
				}
			}
		}
	}

}	

void PartitionDGamRateProcess::SlaveUpdateRateSuffStat()	{

	UpdateRateSuffStat();

	MPI_Send(allocratesuffstatcount,Npart*GetNcat(),MPI_INT,0,TAG1,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Send(allocratesuffstatbeta,Npart*GetNcat(),MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
	if (withpinv)	{
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Send(Ninv,Npart,MPI_INT,0,TAG1,MPI_COMM_WORLD);
	}
}	

