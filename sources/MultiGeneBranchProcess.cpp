
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

#include "MultiGeneBranchProcess.h"

void MultiGeneBranchProcess::Create()	{

	GammaBranchProcess::Create();
	if ((! geneblarray) && (! GlobalBranchLengths()))	{
		allocgeneblarray = new double[Ngene*GetNbranch()];
		geneblarray = new double*[Ngene];
		alloctmpgeneblarray = new double[Ngene*GetNbranch()];
		tmpgeneblarray = new double*[Ngene];
		for (int gene=0; gene<Ngene; gene++)	{
			geneblarray[gene] = allocgeneblarray + gene*GetNbranch();
			tmpgeneblarray[gene] = alloctmpgeneblarray + gene*GetNbranch();
		}

		if (mappsuffstat)	{
		allocgeneblcount = new int[Ngene*GetNbranch()];
		geneblcount = new int*[Ngene];
		alloctmpgeneblcount = new int[Ngene*GetNbranch()];
		tmpgeneblcount  = new int*[Ngene];
		for (int gene=0; gene<Ngene; gene++)	{
			geneblcount[gene] = allocgeneblcount + gene*GetNbranch();
			tmpgeneblcount[gene] = alloctmpgeneblcount + gene*GetNbranch();
		}

		allocgeneblbeta = new double[Ngene*GetNbranch()];
		geneblbeta = new double*[Ngene];
		alloctmpgeneblbeta = new double[Ngene*GetNbranch()];
		tmpgeneblbeta = new double*[Ngene];
		for (int gene=0; gene<Ngene; gene++)	{
			geneblbeta[gene] = allocgeneblbeta + gene*GetNbranch();
			tmpgeneblbeta[gene] = alloctmpgeneblbeta + gene*GetNbranch();
		}
		}

		totlength = new double[GetNbranch()];
		totloglength = new double[GetNbranch()];
	}
}

void MultiGeneBranchProcess::Delete()	{

	if (geneblarray)	{
		delete[] allocgeneblarray;
		delete[] alloctmpgeneblarray;
		delete[] geneblarray;
		delete[] tmpgeneblarray;
		geneblarray = 0;
		delete[] allocgeneblcount;
		delete[] alloctmpgeneblcount;
		delete[] geneblcount;
		delete[] tmpgeneblcount;
		delete[] allocgeneblbeta;
		delete[] alloctmpgeneblbeta;
		delete[] geneblbeta;
		delete[] tmpgeneblbeta;

		delete[] totloglength;
		delete[] totlength;
	}
	GammaBranchProcess::Delete();
}

void MultiGeneBranchProcess::SampleLengthHyperParameters()	{

	if (GlobalBranchLengths())	{
		GammaBranchProcess::SampleLengthHyperParameters();
	}
	else	{
		meanbranchmean = 0.1;
		relvarbranchmean = 1;
		meanbranchrelvar = 1;
		relvarbranchrelvar = 1;
		double amean = 1.0 / relvarbranchmean;
		double bmean = 1.0 / meanbranchmean;
		double arelvar = 1.0 / relvarbranchrelvar;
		double brelvar = 1.0 / meanbranchrelvar;
		for (int j=1; j<GetNbranch(); j++)	{
			/*
			branchmean[j] = rnd::GetRandom().Gamma(amean,bmean);
			branchrelvar[j] = rnd::GetRandom().Gamma(arelvar,brelvar);
			*/
			branchmean[j] = 0.1;
			branchrelvar[j] = 1.0;
		}	
	}
}

void MultiGeneBranchProcess::PriorSampleLengthHyperParameters()	{

	if (GlobalBranchLengths())	{
		GammaBranchProcess::PriorSampleLengthHyperParameters();
	}
	else	{
		meanbranchmean = 0.1 * rnd::GetRandom().sExpo();
		relvarbranchmean = rnd::GetRandom().sExpo();
		meanbranchrelvar = rnd::GetRandom().sExpo();
		relvarbranchrelvar = rnd::GetRandom().sExpo();
		double amean = 1.0 / relvarbranchmean;
		double bmean = 1.0 / meanbranchmean;
		double arelvar = 1.0 / relvarbranchrelvar;
		double brelvar = 1.0 / meanbranchrelvar;
		for (int j=1; j<GetNbranch(); j++)	{
			branchmean[j] = rnd::GetRandom().Gamma(amean,bmean);
			branchrelvar[j] = rnd::GetRandom().Gamma(arelvar,brelvar);
		}	
	}
}

double MultiGeneBranchProcess::LogLengthPrior()	{

	if (! GlobalBranchLengths())	{
		return LogGeneLengthSuffStatPrior();
	}
	return GammaBranchProcess::LogLengthPrior();
}

double MultiGeneBranchProcess::LogLengthHyperPrior()	{

	if (! GlobalBranchLengths())	{
		double total = 0;
		for (int j=1; j<GetNbranch(); j++)	{
			total += LogLengthHyperPrior(j);
		}
		return total;
	}
	return GammaBranchProcess::LogLengthHyperPrior();
}

double MultiGeneBranchProcess::LogLengthHyperPrior(int j)	{

	double amean = 1.0 / relvarbranchmean;
	double bmean = 1.0 / meanbranchmean;
	double arelvar = 1.0 / relvarbranchrelvar;
	double brelvar = 1.0 / meanbranchrelvar;

	double total = 0;
	total += amean * log(bmean) - rnd::GetRandom().logGamma(amean) + (amean-1) * log(branchmean[j]) - bmean * branchmean[j];
	total += arelvar * log(brelvar) - rnd::GetRandom().logGamma(arelvar) + (arelvar-1) * log(branchrelvar[j]) - brelvar * branchrelvar[j];
	if (isnan(total))	{
		cerr << "multi gene log length hyper prior is nan\n";
		exit(1);
	}
	if (isinf(total))	{
		cerr << "multi gene log length hyper prior is inf\n";
		exit(1);
	}
	return total;
}

double MultiGeneBranchProcess::LogLengthHyperHyperPrior()	{

	double ret = - 10*meanbranchmean - relvarbranchmean - meanbranchrelvar - relvarbranchrelvar;
	if (isnan(ret))	{
		cerr << "mult gene log length hyper hyper prior is nan\n";
		exit(1);
	}
	if (isinf(ret))	{
		cerr << "mult gene log length hyper hyper prior is inf\n";
		exit(1);
	}
	return ret;
}

void MultiGeneBranchProcess::ComputeGeneLengthSuffStat()	{

	for (int j=1; j<GetNbranch(); j++)	{
		totlength[j] = 0;
		totloglength[j] = 0;
		for (int gene=0; gene<Ngene; gene++)	{
			double l = geneblarray[gene][j];
			if (l <= 0)	{
				cerr << "error in MultiGeneBranchProcess::ComputeGeneLengthSuffStat: invalid branch length : " << l << '\n';
				exit(1);
			}
			totlength[j] += geneblarray[gene][j];
			totloglength[j] += log(geneblarray[gene][j]);
		}
	}
}

double MultiGeneBranchProcess::LogGeneLengthSuffStatPrior()	{

	ComputeGeneLengthSuffStat();
	double total = 0;
	for (int j=1; j<GetNbranch(); j++)	{
		total += LogGeneLengthSuffStatPrior(j);
	}
	return total;
}

double MultiGeneBranchProcess::LogGeneLengthSuffStatPrior(int j)	{

	double alpha = 1.0 / branchrelvar[j];
	double beta = 1.0 / branchmean[j];
	double ret = Ngene*(alpha*log(beta) - rnd::GetRandom().logGamma(alpha)) + (alpha-1)*totloglength[j] - beta*totlength[j];
	if (isnan(ret))	{
		cerr << "multi gene log gene length suff stat prior is nan\n";
		exit(1);
	}
	if (isinf(ret))	{
		cerr << "multi gene log gene length suff stat prior is inf\n";
		cerr << Ngene << '\t' << alpha << '\t' << beta << '\t' << totloglength[j] << '\t' << totlength[j] << '\n';
		exit(1);
	}
	return ret;
}

double MultiGeneBranchProcess::LogGeneLengthMarginalSuffStatProb()	{

	double total = 0;
	for (int j=1; j<GetNbranch(); j++)	{
		total += LogGeneLengthMarginalSuffStatProb(j);
	}
	return total;
}

double MultiGeneBranchProcess::LogGeneLengthMarginalSuffStatProb(int j)	{

	double alpha = 1.0 / branchrelvar[j];
	double beta = 1.0 / branchmean[j];
	double total = Ngene*(alpha*log(beta) - rnd::GetRandom().logGamma(alpha));
	for (int gene=0; gene<Ngene; gene++)	{
		total -= (alpha + geneblcount[gene][j])*log(beta + geneblbeta[gene][j]) - rnd::GetRandom().logGamma(alpha + geneblcount[gene][j]);
	}
	if (isnan(total))	{
		cerr << "multi gene log gene length suff stat prob is nan\n";
		exit(1);
	}
	if (isinf(total))	{
		cerr << "multi gene log gene length suff stat prob is inf\n";
		exit(1);
	}
	return total;
}

double MultiGeneBranchProcess::Move(double tuning, int nrep)	{

	if (! GlobalBranchLengths())	{

		GlobalCollectGeneBranchLengths();
		ComputeGeneLengthSuffStat();
		if (mappsuffstat)	{
			GlobalCollectGeneLengthMappingSuffStat();
		}

		for (int rep=0; rep<nrep; rep++)	{
			MoveGeneLengthHyperParameters(tuning,10);
			MoveGeneLengthHyperParameters(0.1*tuning,10);
			MoveGeneLengthHyperHyperParameters(tuning,20);
			MoveGeneLengthHyperHyperParameters(0.1*tuning,20);
		}
		return 1.0;
	}
	return GammaBranchProcess::Move(tuning,nrep);
}

double MultiGeneBranchProcess::MoveGeneLengthHyperParameters(double tuning, int nrep)	{

	int Naccepted = 0;
	for (int j=1; j<GetNbranch(); j++)	{
		for (int rep=0; rep<nrep; rep++)	{

			double bkmean = branchmean[j];
			double bkrelvar = branchrelvar[j];

			double deltalogprob = - LogLengthHyperPrior(j);
			if (mappsuffstat)	{
				deltalogprob -= LogGeneLengthMarginalSuffStatProb(j);
			}
			else	{
				deltalogprob -= LogGeneLengthSuffStatPrior(j);
			}
			double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
			double e = exp(m);

			// branchmean[j] *= e;
			if (rnd::GetRandom().Uniform() < 0.5)	{
				branchmean[j] *= e;
			}
			else	{
				branchrelvar[j] *= e;
			}

			deltalogprob += LogLengthHyperPrior(j);
			if (mappsuffstat)	{
				deltalogprob += LogGeneLengthMarginalSuffStatProb(j);
			}
			else	{
				deltalogprob += LogGeneLengthSuffStatPrior(j);
			}
			deltalogprob += m;
			int accepted = (log(rnd::GetRandom().Uniform()) < deltalogprob);
			if (accepted)	{
				Naccepted ++;
			}
			else	{
				branchmean[j] = bkmean;
				branchrelvar[j] = bkrelvar;
			}
		}
	}
	return ((double) Naccepted) / nrep / (GetNbranch()-1);
}

double MultiGeneBranchProcess::MoveGeneLengthHyperHyperParameters(double tuning, int nrep)	{

	int Naccepted = 0;
	for (int rep=0; rep<nrep; rep++)	{

		double bkmeanbranchmean = meanbranchmean;
		double bkmeanbranchrelvar = meanbranchrelvar;
		double bkrelvarbranchmean = relvarbranchmean;
		double bkrelvarbranchrelvar = relvarbranchrelvar;

		double deltalogprob = - LogLengthHyperHyperPrior() - LogLengthHyperPrior();
		double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
		double e = exp(m);

		// int c = (int) (2 * rnd::GetRandom().Uniform());
		int c = (int) (4 * rnd::GetRandom().Uniform());
		if (c == 0)	{
			meanbranchmean *= e;
		}
		else if (c == 1)	{
			relvarbranchmean *= e;
		}
		else if (c == 2)	{
			meanbranchrelvar *= e;
		}
		else	{
			relvarbranchrelvar *= e;
		}

		deltalogprob += LogLengthHyperHyperPrior() + LogLengthHyperPrior();

		deltalogprob += m;
		int accepted = (log(rnd::GetRandom().Uniform()) < deltalogprob);
		if (accepted)	{
			Naccepted ++;
		}
		else	{
			meanbranchmean = bkmeanbranchmean;
			relvarbranchmean = bkrelvarbranchmean;
			meanbranchrelvar = bkmeanbranchrelvar;
			relvarbranchrelvar = bkrelvarbranchrelvar;
		}
	}
	return ((double) Naccepted) / nrep;
}

double MultiGeneBranchProcess::GlobalGetMeanTotalLength()	{

	MESSAGE signal = MEANTOTLENGTH;
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

void MultiGeneBranchProcess::SlaveGetMeanTotalLength() {
	double tot = 0;
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			double tmp = GetBranchProcess(gene)->GetTotalLength();
			tot += tmp;
		}
	}
	MPI_Send(&tot,1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

void MultiGeneBranchProcess::GlobalCollectGeneBranchLengths()	{

	MESSAGE signal = COLLECTLENGTHS;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Status stat;
	for(int i=1; i<GetNprocs(); i++) {
		MPI_Recv(alloctmpgeneblarray,Ngene*GetNbranch(),MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
		for (int gene=0; gene<Ngene; gene++)	{
			if (genealloc[gene] == i)	{
				for (int j=1; j<GetNbranch(); j++)	{
					geneblarray[gene][j] = tmpgeneblarray[gene][j];
					if (! geneblarray[gene][j])	{
						cerr << "error in global collect gene bl: null bl\n";
						exit(1);
					}
				}
			}
		}
	}
}

void MultiGeneBranchProcess::SlaveCollectGeneBranchLengths() {
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			const double* bl = GetBranchProcess(gene)->GetBranchLengths();
			for (int j=1; j<GetNbranch(); j++)	{
				if (isnan(bl[j]) || (! bl[j]))	{
					cerr << "error in slave collect; null bl\n";
					cerr << j << '\n';
					for (int j=0; j<GetNbranch(); j++)	{
						cerr << j << '\t' << bl[j] << '\n';
					}
					exit(1);
				}
				tmpgeneblarray[gene][j] = bl[j];
			}
		}
	}
	MPI_Send(alloctmpgeneblarray,Ngene*GetNbranch(),MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

void MultiGeneBranchProcess::GlobalCollectGeneLengthMappingSuffStat()	{

	MESSAGE signal = COLLECTLENGTHSUFFSTAT;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Status stat;
	for(int i=1; i<GetNprocs(); i++) {
		MPI_Recv(alloctmpgeneblcount,Ngene*GetNbranch(),MPI_INT,i,TAG1,MPI_COMM_WORLD,&stat);
		MPI_Recv(alloctmpgeneblbeta,Ngene*GetNbranch(),MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
		for (int gene=0; gene<Ngene; gene++)	{
			if (genealloc[gene] == i)	{
				for (int j=1; j<GetNbranch(); j++)	{
					geneblcount[gene][j] = tmpgeneblcount[gene][j];
					geneblbeta[gene][j] = tmpgeneblbeta[gene][j];
				}
			}
		}
	}
}

void MultiGeneBranchProcess::SlaveCollectGeneLengthMappingSuffStat() {
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			const double* beta = process[gene]->GetBranchLengthSuffStatBeta();
			const int* count = process[gene]->GetBranchLengthSuffStatCount();
			for (int j=1; j<GetNbranch(); j++)	{
				tmpgeneblcount[gene][j] = count[j];
				tmpgeneblbeta[gene][j] = beta[j];
			}
		}
	}
	MPI_Send(alloctmpgeneblcount,Ngene*GetNbranch(),MPI_INT,0,TAG1,MPI_COMM_WORLD);
	MPI_Send(alloctmpgeneblbeta,Ngene*GetNbranch(),MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

void MultiGeneBranchProcess::SlaveKnit()	{

	int arg;
	MPI_Bcast(&arg,1,MPI_INT,0,MPI_COMM_WORLD);
	LocalKnit(arg);
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			GetBranchProcess(gene)->LocalKnit(arg);
		}
	}
}


void MultiGeneBranchProcess::SlaveDetach(int n,int m) {

	LocalDetach(n,m);
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->LocalDetach(n,m);
		}
	}
}

void MultiGeneBranchProcess::SlaveAttach(int n,int m,int p,int q) {

	LocalAttach(n,m,p,q);
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->LocalAttach(n,m,p,q);
		}
	}
}

void MultiGeneBranchProcess::SlaveDetach1(int n,int m) {

	LocalDetach1(n,m);
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->LocalDetach1(n,m);
		}
	}
}

void MultiGeneBranchProcess::SlaveAttach1(int n,int m,int p,int q) {

	LocalAttach1(n,m,p,q);
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->LocalAttach1(n,m,p,q);
		}
	}
}

void MultiGeneBranchProcess::SlaveDetach2(int n,int m) {

	LocalDetach2(n,m);
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->LocalDetach2(n,m);
		}
	}
}

void MultiGeneBranchProcess::SlaveAttach2(int n,int m,int p,int q) {

	LocalAttach2(n,m,p,q);
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->LocalAttach2(n,m,p,q);
		}
	}
}

/*
void MultiGeneBranchProcess::SlaveSwapRoot()	{

	SwapRoot();
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->SwapRoot();
		}
	}
}
*/
