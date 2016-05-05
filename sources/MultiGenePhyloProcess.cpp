
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

#include "MultiGenePhyloProcess.h"

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
			tot += GetBranchProcess(gene)->GetTotalLength();
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
				if (! bl[j])	{
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


void MultiGenePhyloProcess::New()	{


	CreateMPI(0);
	// SetMPI(myid,nprocs);
	AllocateAlignments(datafile,treefile);
	SetProfileDim();
	SetTree(treefile);

	Create();

	if (! GetMyid())	{
		GlobalBroadcastTree();
		Sample();
		GlobalUpdateParameters();
		GlobalSample();
		GlobalUnfold();
	}
	else	{
		SlaveBroadcastTree();
	}

	if (BPP)	{
		BPP->RegisterWithTaxonSet(GetData()->GetTaxonSet());
	}
}


void MultiGenePhyloProcess::Create()	{

	PhyloProcess::Create();
	MultiGeneRateProcess::Create();
	MultiGeneBranchProcess::Create();
	MultiGeneMPIModule::Create();
}

void MultiGenePhyloProcess::Delete()	{

	MultiGeneMPIModule::Delete();
	MultiGeneBranchProcess::Delete();
	MultiGeneRateProcess::Delete();
	PhyloProcess::Delete();
}

void MultiGeneMPIModule::Create()	{

	if (! genelnL)	{
		genelnL = new double[Ngene];
		tmpgenelnL = new double[Ngene];
		process = new PhyloProcess*[Ngene];
		for (int gene=0; gene<Ngene; gene++)	{
			process[gene] = 0;
		}
	}
}

void MultiGeneMPIModule::Delete()	{

	if (genelnL)	{
		delete[] genelnL;
		delete[] tmpgenelnL;
		genelnL = 0;
		for (int gene=0; gene<Ngene; gene++)	{
			delete process[gene];
			process[gene] = 0;
		}
		delete[] process;
	}
}

void MultiGenePhyloProcess::AllocateAlignments(string datafile, string treefile)	{

	ifstream is(datafile.c_str());
	is >> Ngene;
	ifstream* tis = 0;
	genename = new string[Ngene];
	genesize = new int[Ngene];
	genealloc = new int[Ngene];
	int* geneweight = new int[Ngene];
	// genedata = new SequenceAlignment*[Ngene];
	for (int gene=0; gene<Ngene; gene++)	{
		is >> genename[gene];
		SequenceAlignment* localdata = new FileSequenceAlignment(genename[gene]);
		int nstate = localdata->GetNstate();
		if (! gene)	{
			data = localdata;
		}
		else	{
			if (nstate != data->GetNstate())	{
				cerr << "error: all data files do not have the same alphabet\n";
				cerr << nstate << '\t' << data->GetNstate() << '\n';
				exit(1);
			}
		}

		genesize[gene] = localdata->GetNsite();
		geneweight[gene] = localdata->GetNsite() * localdata->GetNtaxa();
		if (gene)	{
			delete localdata;
		}
	}
	delete tis;
	tis = 0;

	// sort alignments by decreasing size
	int permut[Ngene];
	for (int gene=0; gene<Ngene; gene++)	{
		permut[gene] = gene;
	}
	for (int i=0; i<Ngene; i++)	{
		for (int j=Ngene-1; j>i; j--)	{
			if (geneweight[permut[i]] < geneweight[permut[j]])	{
			// if (genesize[permut[i]] < genesize[permut[j]])	{
				int tmp = permut[i];
				permut[i] = permut[j];
				permut[j] = tmp;
			}
		}
	}

	int totsize[nprocs];
	for (int i=0; i<nprocs; i++)	{
		totsize[i] = 0;
	}

	for (int i=0; i<Ngene; i++)	{
		int gene = permut[i];
		int size = geneweight[gene];
		// int size = genesize[gene];

		int min = 0;
		int jmin = 0;
		for (int j=1; j<nprocs; j++)	{
			if ((j==1) || (min > totsize[j]))	{
				min = totsize[j];
				jmin = j;
			}
		}
		genealloc[gene] = jmin;
		totsize[jmin] += size;
	}

	if (totsize[0])	{
		cerr << "error in alloc\n";
		exit(1);
	}
	int total = 0;
	for (int i=1; i<nprocs; i++)	{
		int tot = 0;
		for (int gene=0; gene<Ngene; gene++)	{
			if (genealloc[gene] == i)	{
				tot += geneweight[gene];
				// tot += genesize[gene];
				total++;
			}
		}
		if (tot != totsize[i])	{
			cerr << "error in allocation\n";
			cerr << tot << '\t' << totsize[i] << '\n';
			exit(1);
		}
	}
	if (total != Ngene)	{
		cerr << "error in total allocation\n";
		exit(1);
	}

	globalnsite = new int[nprocs];
	for (int i=0; i<nprocs; i++)	{
		globalnsite[i] = 0;
	}
	for (int gene=0; gene<Ngene; gene++)	{
		if ((genealloc[gene] < 0) || (genealloc[gene] >= nprocs))	{
			cerr << "alloc : " << genealloc[gene] << '\t' << gene << '\n';
			exit(1);
		}
		globalnsite[0] += genesize[gene];
		globalnsite[genealloc[gene]] += genesize[gene];
	}
	if (! myid)	{
		cerr << '\n';
		cerr << "number of sites allocated to each slave:\n";
		for (int i=1; i<nprocs; i++)	{
			cerr << i << '\t' << globalnsite[i] << '\n';
		}
		cerr << '\n';
		// cerr << "total: " << GetGlobalNsite() << '\n';
	}
	
	// check total size
	if (! myid)	{
		int tot = 0;
		for (int i=1; i<nprocs; i++)	{
			tot += globalnsite[i];
		}
		if (tot != globalnsite[0])	{
			cerr << "error in total size during gene allocation\n";
			exit(1);
		}
		int tot2 = 0;
		for (int gene=0; gene<Ngene; gene++)	{
			tot2 += genesize[gene];
		}
		if (tot2 != tot)	{
			cerr << "error during alloc: total size does not match\n";
			exit(1);
		}
	}
	delete[] geneweight;
}

int MultiGenePhyloProcess::SpecialSlaveExecute(MESSAGE signal)	{

	switch(signal) {
	case SAMPLE:
		SlaveSample();
		return 1;
		break;
	case GENE_MOVE:
		SlaveGeneMove();
		return 1;
		break;
	case MEANALPHA:
		SlaveGetMeanAlpha();
		return 1;
		break;
	case COLLECTALPHA:
		SlaveCollectGeneAlphas();
		return 1;
		break;
	case MEANTOTLENGTH:
		SlaveGetMeanTotalLength();
		return 1;
		break;
	case COLLECTLENGTHS:
		SlaveCollectGeneBranchLengths();
		return 1;
		break;
	case COLLECTLENGTHSUFFSTAT:
		SlaveCollectGeneLengthMappingSuffStat();
		return 1;
		break;

	default:
		return 0;
		// PhyloProcess::SlaveExecute(signal);
	}
}

// should be defined at the level of the non specialized PhyloProcess class
void MultiGenePhyloProcess::GlobalSample()	{
	assert(myid == 0);
	MESSAGE signal = SAMPLE;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
}

void MultiGenePhyloProcess::SlaveSample()	{
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->Sample();
		}
	}
}

void MultiGenePhyloProcess::GlobalGeneMove()	{
	assert(myid == 0);
	MESSAGE signal = GENE_MOVE;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
}

void MultiGenePhyloProcess::SlaveGeneMove()	{
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->AugmentedMove();
		}
	}
}

void MultiGenePhyloProcess::SlaveBroadcastTree()	{

	int len;
	MPI_Bcast(&len,1,MPI_INT,0,MPI_COMM_WORLD);
	unsigned char* bvector = new unsigned char[len];
	MPI_Bcast(bvector,len,MPI_UNSIGNED_CHAR,0,MPI_COMM_WORLD);
	ostringstream os;
	for (int i=0; i<len; i++)	{
		os << bvector[i];
	}
	istringstream is(os.str());
	tree->ReadFromStream(is);
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			istringstream is(os.str());
			process[gene]->GetTree()->ReadFromStream(is);
			process[gene]->GetTree()->RegisterWith(GetData()->GetTaxonSet());
			process[gene]->CloneTree();
			process[gene]->GetTree2()->RegisterWith(GetData()->GetTaxonSet());
		}
	}
	delete[] bvector;
}

void MultiGenePhyloProcess::SlaveUnfold()	{
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->Unfold();
		}
	}
}

void MultiGenePhyloProcess::SlaveCollapse()	{
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->Collapse();
		}
	}
}

void MultiGenePhyloProcess::SlaveActivateSumOverRateAllocations()	{
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->ActivateSumOverRateAllocations();
			process[gene]->sumratealloc = 1;
		}
	}
	sumratealloc = 1;
}

void MultiGenePhyloProcess::SlaveInactivateSumOverRateAllocations()	{
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->InactivateSumOverRateAllocations(process[gene]->ratealloc);
			process[gene]->sumratealloc = 0;
		}
	}
	sumratealloc = 0;
}

void MultiGenePhyloProcess::SlaveActivateZip()	{
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->ActivateZip();
		}
	}
}

void MultiGenePhyloProcess::SlaveInactivateZip()	{
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->InactivateZip();
		}
	}
}

void MultiGenePhyloProcess::SlaveUpdateConditionalLikelihoods()	{
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->UpdateConditionalLikelihoods();
		}
	}
}

void MultiGenePhyloProcess::SlaveComputeNodeLikelihood(int fromindex,int auxindex) {
	double totlogl = 0;
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			genelnL[gene] = process[gene]->LocalComputeNodeLikelihood(fromindex,auxindex);
			totlogl += genelnL[gene];
		}
	}
	MPI_Send(&totlogl,1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

void MultiGenePhyloProcess::SlaveGetFullLogLikelihood()	{

	cerr << "in multigenephyloprocess::slavegetfullloglikelihood\n";
	exit(1);
	double totlogl = 0;
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			genelnL[gene] = process[gene]->GetFullLogLikelihood();
			totlogl += genelnL[gene];
		}
	}
	double sum[2];
	sum[0] = totlogl;
	// normally, should be the likelihood conditional on allocations
	sum[1] = totlogl;
	// sum[1] = logL;
	MPI_Send(&totlogl,2,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

void MultiGenePhyloProcess::SlaveReset(int n,bool v) {
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->SlaveReset(n,v);
		}
	}
}

void MultiGenePhyloProcess::SlaveMultiply(int n,int m,bool v) {
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->SlaveMultiply(n,m,v);
		}
	}
}

void MultiGenePhyloProcess::SlaveMultiplyByStationaries(int n,bool v) {
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->SlaveMultiplyByStationaries(n,v);
		}
	}
}

void MultiGenePhyloProcess::SlaveInitialize(int n,int m,bool v) {
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->SlaveInitialize(n,m,v);
		}
	}
}

void MultiGenePhyloProcess::SlavePropagate(int n,int m,bool v,double t) {
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->SlavePropagate(n,m,v,t);
		}
	}
}

void MultiGenePhyloProcess::SlaveProposeMove(int n,double x) {
	PhyloProcess::SlaveProposeMove(n,x);
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->SlaveProposeMove(n,x);
		}
	}
}

void MultiGenePhyloProcess::SlaveRestore(int n) {
	PhyloProcess::SlaveRestore(n);
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->SlaveRestore(n);
		}
	}
}

void MultiGenePhyloProcess::SlaveRoot(int n) {
	PhyloProcess::SlaveRoot(n);
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->SlaveRoot(n);
		}
	}
}

void MultiGenePhyloProcess::SlaveBackupTree()	{
	GetTree()->Backup();
	GetTree2()->Backup();
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->GetTree()->Backup();
			process[gene]->GetTree2()->Backup();
		}
	}
}

void MultiGenePhyloProcess::SlaveRestoreTree()	{
	GetTree()->Restore();
	GetTree2()->Restore();
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->GetTree()->Restore();
			process[gene]->GetTree2()->Restore();
		}
	}
}

void MultiGenePhyloProcess::SlaveSwapTree()	{
	GetTree()->Swap();
	GetTree2()->Swap();
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->GetTree()->Swap();
			process[gene]->GetTree2()->Swap();
		}
	}
}

void MultiGenePhyloProcess::SlaveGibbsSPRScan(int idown, int iup)	{
	for(int i=0; i<GetNbranch(); i++) {
		loglarray[i] = 0.0;
	}
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->LocalGibbsSPRScan(idown,iup);
			for(int i=0; i<GetNbranch(); i++) {
				loglarray[i] += process[gene]->loglarray[i];
			}
		}
	}
	MPI_Send(loglarray,GetNbranch(),MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

void MultiGenePhyloProcess::SlaveSetMinMax()	{

	double minmax[2];
	MPI_Bcast(minmax,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->SetMinMax(minmax[0],minmax[1]);
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

void MultiGeneBranchProcess::SlaveSwapRoot()	{

	SwapRoot();
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->SwapRoot();
		}
	}
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

void MultiGeneProfileProcess::GlobalUpdateSiteProfileSuffStat()	{

	if (GetNprocs() > 1)	{
		MPI_Status stat;
		MESSAGE signal = UPDATE_SPROFILE;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	}
	else	{
		UpdateSiteProfileSuffStat();
	}
}

void MultiGeneProfileProcess::SlaveUpdateSiteProfileSuffStat()	{

	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->UpdateSiteProfileSuffStat();
		}
	}
}

void MultiGenePhyloProcess::UpdateBranchLengthSuffStat() {

	for(int i=0; i<GetNbranch(); i++) {
		branchlengthsuffstatcount[i] = 0;
		branchlengthsuffstatbeta[i] = 0.0;
	}

	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->UpdateBranchLengthSuffStat();
			const int* count = process[gene]->GetBranchLengthSuffStatCount();
			const double* beta = process[gene]->GetBranchLengthSuffStatBeta();
			for(int i=0; i<GetNbranch(); i++) {
				branchlengthsuffstatcount[i] += count[i];
				branchlengthsuffstatbeta[i] += beta[i];
			}
		}
	}
}
