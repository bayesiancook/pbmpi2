
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "BiologicalSequences.h"
#include "FiniteProfileProcess.h"
#include "Random.h"
#include "Parallel.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* FiniteProfileProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------


void FiniteProfileProcess::Create()	{
	if (! weight)	{
		if (Ncomponent == -1)	{
			if (Npart)	{
				Ncomponent = Npart;
			}
			else	{
				Ncomponent = GetNsite();
			}
		}
		DirichletMixtureProfileProcess::Create();
		weight = new double[GetNmodeMax()];
	}
}

void FiniteProfileProcess::Delete()	{
	if (weight)	{
		if (statfix)	{
			for (int i=0; i<Nfixcomp; i++)	{
				delete[] statfix[i];
			}
			delete[] statfix;
			delete[] empweight;
		}
		delete[] weight;
		MixtureProfileProcess::Delete();
	}
}


double FiniteProfileProcess::Move(double tuning, int n, int nrep)	{

	double ret = 0;
	if (GetNprocs() > 1)	{
		ret = MPIMove(tuning,n,nrep);
	}
	else	{
		ret = NonMPIMove(tuning,n,nrep);
	}
	return ret;
}

double FiniteProfileProcess::MPIMove(double tuning, int n, int nrep)	{

	totchrono.Start();

	for (int rep=0; rep<nrep; rep++)	{

		GlobalParametersMove();

		if ((Ncomponent > 1) && (!Npart))	{
			// allocations
			GlobalUpdateParameters();
			GlobalUpdateSiteProfileSuffStat();

			incchrono.Start();
			GlobalIncrementalFiniteMove(1);
			incchrono.Stop();
		}

		if (empmix != 1)	{
			// profiles
			GlobalUpdateParameters();
			GlobalUpdateSiteProfileSuffStat();
			GlobalUpdateModeProfileSuffStat();

			profilechrono.Start();
			GlobalMoveProfile(1,1,100);
			GlobalMoveProfile(1,3,100);
			GlobalMoveProfile(0.1,3,100);
			profilechrono.Stop();

            GlobalUpdateParameters();
            GlobalUpdateSiteProfileSuffStat();
			MoveHyper(tuning,10);
		}

		if (! fixncomp)	{
			MoveNcomponent(10);
			MoveWeightAlpha(tuning,10);
		}
	}

	ResampleWeights();
	GlobalUpdateParameters();

	totchrono.Stop();
	return 1;
}

double FiniteProfileProcess::NonMPIMove(double tuning, int n, int nrep)	{

	totchrono.Start();

	for (int rep=0; rep<nrep; rep++)	{

		GlobalParametersMove();

		if ((Ncomponent > 1) && (!Npart))	{
			// allocations
			UpdateSiteProfileSuffStat();

			incchrono.Start();
			IncrementalFiniteMove(1);
			incchrono.Stop();
		}

		if (empmix != 1)	{
			// profiles
			UpdateSiteProfileSuffStat();
			UpdateModeProfileSuffStat();
			profilechrono.Start();
			MoveProfile(1,1,100);
			MoveProfile(1,3,100);
			MoveProfile(0.1,3,100);
			profilechrono.Stop();
			MoveHyper(tuning,10);
		}

		if (! fixncomp)	{
			MoveNcomponent(10);
			MoveWeightAlpha(tuning,10);
		}
	}

	ResampleWeights();

	totchrono.Stop();
	return 1;
}

double FiniteProfileProcess::GlobalIncrementalFiniteMove(int nrep)	{

	// send command and arguments
	MESSAGE signal = REALLOC_MOVE;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&nrep,1,MPI_INT,0,MPI_COMM_WORLD);

	int NAccepted = 0;
	for (int rep=0; rep<nrep; rep++)	{

		// resample component weights based on current site allocations
		// and send them to slaves
		UpdateOccupancyNumbers();
		ResampleWeights();
		MPI_Bcast(weight,Ncomponent,MPI_DOUBLE,0,MPI_COMM_WORLD);

		// receive new site allocations from slave
		MPI_Status stat;
		int tmpalloc[GetNsite()];
		for(int i=1; i<GetNprocs(); ++i) {
			MPI_Recv(tmpalloc,GetNsite(),MPI_INT,i,TAG1,MPI_COMM_WORLD,&stat);
			for(int j=GetProcSiteMin(i); j<GetProcSiteMax(i); ++j) {
				if (ActiveSite(j))	{
					alloc[j] = tmpalloc[j];
					if ((alloc[j] < 0) || (alloc[j] >= Ncomponent))	{
						cerr << "alloc overflow\n";
						exit(1);
					}
				}
			}
		}
	}
	
	// final cleanup
	UpdateOccupancyNumbers();
	// UpdateModeProfileSuffStat();

	// CHECK that: might be useful depending on the exact submodel
	// ResampleWeights();
	// UpdateMatrices();

	return 1;
}

double FiniteProfileProcess::SlaveIncrementalFiniteMove()	{

	// parse argument sent by master
	int nrep;
	MPI_Bcast(&nrep,1,MPI_INT,0,MPI_COMM_WORLD);

	int NAccepted = 0;
	int Ntried = 0;

	double* bigarray = new double[Ncomponent * GetNsite()];
	double* bigcumul = new double[Ncomponent * GetNsite()];

	for (int rep=0; rep<nrep; rep++)	{

		// receive weights sent by master
		MPI_Bcast(weight,Ncomponent,MPI_DOUBLE,0,MPI_COMM_WORLD);

		// do the incremental reallocation move on my site range
		for (int site=GetSiteMin(); site<GetSiteMax(); site++)	{

			if (ActiveSite(site))	{

				double* mLogSamplingArray = bigarray + site * Ncomponent;
				double* cumul = bigcumul + site * Ncomponent;

				if (alloc[site] == -1)	{
					cerr << "error in SlaveIncrementalFiniteMove: alloc == -1\n";
					exit(1);
				}
				int bk = alloc[site];

				double max = 0;
				for (int mode = 0; mode < Ncomponent; mode++)	{
					mLogSamplingArray[mode] =  LogStatProb(site,mode);
					if ((!mode) || (max < mLogSamplingArray[mode]))	{
						max = mLogSamplingArray[mode];
					}
				}

				double total = 0;

				for (int mode = 0; mode < Ncomponent; mode++)	{
					double p = weight[mode] * exp(mLogSamplingArray[mode] - max);
					total += p;
					cumul[mode] = total;
				}

				double q = total * rnd::GetRandom().Uniform();
				int mode = 0;
				while ( (mode<Ncomponent) && (q > cumul[mode])) mode++;
				if (mode == Ncomponent)	{
					cerr << "error in switch mode: gibbs overflow\n";
					exit(1);
				}

				int Accepted = (mode != bk);
				if (Accepted)	{
					NAccepted ++;
				}
				Ntried ++;
				alloc[site] = mode;
			}
		}

		// send new allocations to master
		MPI_Send(alloc,GetNsite(),MPI_INT,0,TAG1,MPI_COMM_WORLD);
	}
	
	delete[] bigarray;
	delete[] bigcumul;
	return ((double) NAccepted) / Ntried;
}

double FiniteProfileProcess::IncrementalFiniteMove(int nrep)	{

	if (GetMyid())	{
		cerr << "error: slave in FiniteProfileProcess::IncrementalFiniteMove\n";
		exit(1);
	}

	int NAccepted = 0;

	double* mLogSamplingArray = new double[Ncomponent];
	double* cumul = new double[Ncomponent];

	for (int rep=0; rep<nrep; rep++)	{

		ResampleWeights();

		for (int site=0; site<GetNsite(); site++)	{

			if (ActiveSite(site))	{

			if (alloc[site] == -1)	{
				cerr << "error in IncrementalFiniteMove: alloc == -1\n";
				exit(1);
			}

			int bk = alloc[site];
			RemoveSite(site,bk);

			double max = 0;
			for (int mode = 0; mode < Ncomponent; mode++)	{
				mLogSamplingArray[mode] =  LogStatProb(site,mode);
				if ((!mode) || (max < mLogSamplingArray[mode]))	{
					max = mLogSamplingArray[mode];
				}
			}

			double total = 0;

			for (int mode = 0; mode < Ncomponent; mode++)	{
				double p = weight[mode] * exp(mLogSamplingArray[mode] - max);
				total += p;
				cumul[mode] = total;
			}

			double q = total * rnd::GetRandom().Uniform();
			int mode = 0;
			while ( (mode<Ncomponent) && (q > cumul[mode])) mode++;
			if (mode == Ncomponent)	{
				cerr << "error in switch mode: gibbs overflow\n";
				exit(1);
			}

			int Accepted = (mode != bk);
			if (Accepted)	{
				NAccepted ++;
			}
			alloc[site] = mode;
			AddSite(site,mode);
			}

		}
	}
	
	delete[] cumul;
	delete[] mLogSamplingArray;
	// UpdateModeProfileSuffStat();
	return ((double) NAccepted) / GetNsite() / nrep;
}

void FiniteProfileProcess::ResampleWeights()	{
    if (GetNcomponent() > 1)    {
        UpdateOccupancyNumbers();
        double total = 0;
        for (int k=0; k<GetNcomponent(); k++)	{
            weight[k] = rnd::GetRandom().sGamma(weightalpha + occupancy[k]);
            if (weight[k] < 1e-10)	{
                weight[k] = 1e-10;
            }
            total += weight[k];
        }
        for (int k=0; k<GetNcomponent(); k++)	{
            weight[k] /= total;
        }
    }
}

void FiniteProfileProcess::SampleWeights()	{
	double total = 0;
	for (int k=0; k<GetNcomponent(); k++)	{
		weight[k] = rnd::GetRandom().sGamma(weightalpha);
		/*
		if (empmix == 2)	{
			weight[k] = rnd::GetRandom().sGamma(empweight[k] * weightalpha * GetNcomponent());
		}
		else	{
			weight[k] = rnd::GetRandom().sGamma(weightalpha);
		}
		*/
		total += weight[k];
	}
	for (int k=0; k<GetNcomponent(); k++)	{
		weight[k] /= total;
	}
}

void FiniteProfileProcess::SampleHyper()	{

	if (empmix != 2)	{
		DirichletProfileProcess::SampleHyper();
	}
	statfixalpha = 0.01;
	weightalpha = 1.0;
	SampleWeights();
}

void FiniteProfileProcess::PriorSampleHyper()	{

	if (empmix != 2)	{
		DirichletProfileProcess::SampleHyper();
	}
	statfixalpha = 0.01;
	weightalpha = 1.0;
	if (empmix == 2)	{
		statfixalpha = 10 * rnd::GetRandom().sExpo();
	}
	if (! fixncomp)	{
		weightalpha = rnd::GetRandom().sExpo();
	}
	SampleWeights();
}

void FiniteProfileProcess::SampleAlloc()	{
	if (!GetNcomponent())	{
		cerr << "error in sample alloc: " << GetNcomponent() << '\n';
		exit(1);
	}
	if (empmix)	{
		ReadStatFix(mixtype);
		SetStatFix();
		BroadcastStatFix();
	}
	SampleWeights();
	if (Npart)	{
		for (int i=0; i<GetNsite(); i++)	{
			AddSite(i,partalloc[i]);
		}
	}
	else if (GetNcomponent() == GetNsite())	{
		for (int i=0; i<GetNsite(); i++)	{
			AddSite(i,i);
		}	
	}
	else	{
		for (int i=0; i<GetNsite(); i++)	{
			if (ActiveSite(i))	{
				int choose = rnd::GetRandom().FiniteDiscrete(GetNcomponent(),weight);
				AddSite(i,choose);
			}
		}
	}
	ResampleWeights();
}

void FiniteProfileProcess::SampleStat()	{
	if (! empmix)	{
		MixtureProfileProcess::SampleStat();
	}
}

void FiniteProfileProcess::BroadcastStatFix()	{

	MESSAGE signal = STATFIX;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	int ncat = GetNcomponent();
	MPI_Bcast(&ncat,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(empweight,ncat,MPI_DOUBLE,0,MPI_COMM_WORLD);
	for (int i=0; i<ncat; i++)	{
		MPI_Bcast(statfix[i],GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	}
}

void FiniteProfileProcess::SlaveGetStatFix()	{

	MPI_Bcast(&Nfixcomp,1,MPI_INT,0,MPI_COMM_WORLD);
	empweight = new double[Nfixcomp];
	MPI_Bcast(empweight,Nfixcomp,MPI_DOUBLE,0,MPI_COMM_WORLD);
	statfix = new double*[Nfixcomp];
	for (int i=0; i<Nfixcomp; i++)	{
		statfix[i] = new double[GetDim()];
		MPI_Bcast(statfix[i],GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	}

	if (empmix == 3)	{
		SetStatCons();
	}
}

void FiniteProfileProcess::SampleStat(int cat)	{

	if (! empmix)	{
		MixtureProfileProcess::SampleStat(cat);
	}
	else if (empmix == 2)	{
		double statmin = stateps;
		double total = 0;
		int infreached = 0;
		for (int k=0; k<GetDim(); k++)	{
			profile[cat][k] = rnd::GetRandom().sGamma(statfix[cat][k] / statfixalpha * GetDim());
			if (profile[cat][k] < statmin)	{
				profile[cat][k] = statmin;
				infreached = 1;
			}
			total += profile[cat][k];
		}
		for (int k=0; k<GetDim(); k++)	{
			profile[cat][k] /= total;
		}
		if (infreached)	{
			statinfcount++;
		}
		totstatcount++;
	}
}


double FiniteProfileProcess::LogWeightPrior()	{

	UpdateOccupancyNumbers();
	double total = 0;
	double totnsite = 0;
	for (int k=0; k<GetNcomponent(); k++)	{
		total += rnd::GetRandom().logGamma(weightalpha + occupancy[k]);
		totnsite += occupancy[k];
	}
	total -= rnd::GetRandom().logGamma(GetNcomponent() * weightalpha + totnsite);
	total += rnd::GetRandom().logGamma(GetNcomponent()*weightalpha) - GetNcomponent() * rnd::GetRandom().logGamma(weightalpha);
	return total;
}

double FiniteProfileProcess::LogHyperPrior()	{

	double total = 0;
	double sum = 0;
	if (empmix == 0)	{
		total += DirichletProfileProcess::LogHyperPrior();
	}
	else if (empmix == 1)	{
	}
	else if (empmix == 2)	{
		total -= 10 * statfixalpha;
	}
	else if (empmix == 3)	{
	}
	total -= weightalpha;
	return total;
}

double FiniteProfileProcess::LogStatPrior(int cat)	{

	double total = 0;
	if (empmix == 0)	{
		total = MixtureProfileProcess::LogStatPrior(cat);
	}
	else if (empmix == 1)	{
		total = 0;
	}
	else if (empmix == 1)	{
		for (int k=0; k<GetDim(); k++)	{
			total += (1.0 / statfixalpha * GetDim() * statfix[cat][k] - 1) * log(profile[cat][k]) - rnd::GetRandom().logGamma(1.0 / statfixalpha * GetDim() * statfix[cat][k]);
		}
		total += rnd::GetRandom().logGamma(1.0 / statfixalpha * GetDim());
	}
	else if (empmix == 3)	{
		total = LogStatPriorConstrained(cat);
	}
	return total;
	
}

double FiniteProfileProcess::MoveHyper(double tuning, int nrep)	{

	double total = 0;
	if (empmix == 0)	{
		total += DirichletProfileProcess::MoveHyper(tuning,nrep);
	}
	else if (empmix == 2)	{
		total += MoveStatFixAlpha(tuning,nrep);
	}
	UpdateOccupancyNumbers();
	ResampleEmptyProfiles();
	return total;
}

double FiniteProfileProcess::MoveStatFixAlpha(double tuning, int nrep)	{
	UpdateOccupancyNumbers();
	double naccepted = 0;
	for (int rep=0; rep<nrep; rep++)	{
		double deltalogprob = - LogHyperPrior() - MixtureProfileProcess::LogStatPrior();
		double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
		double e = exp(m);
		statfixalpha *= e;
		deltalogprob += LogHyperPrior() + MixtureProfileProcess::LogStatPrior();
		deltalogprob += m;
		int accepted = (log(rnd::GetRandom().Uniform()) < deltalogprob);
		if (accepted)	{
			naccepted++;
		}
		else	{
			statfixalpha /= e;
		}
	}
	return naccepted / nrep;
}

double FiniteProfileProcess::MoveWeightAlpha(double tuning, int nrep)	{
	UpdateOccupancyNumbers();
	double naccepted = 0;
	for (int rep=0; rep<nrep; rep++)	{
		double deltalogprob = - LogHyperPrior() - LogWeightPrior();
		double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
		double e = exp(m);
		weightalpha *= e;
		deltalogprob += LogHyperPrior() + LogWeightPrior();
		deltalogprob += m;
		int accepted = (log(rnd::GetRandom().Uniform()) < deltalogprob);
		if (accepted)	{
			naccepted++;
		}
		else	{
			weightalpha /= e;
		}
	}
	return naccepted / nrep;
}

double FiniteProfileProcess::MoveNcomponent(int nrep)	{

	int Kmax = GetNmodeMax();

	int kmax = Ncomponent-1;
	int kmin = 0;
	do	{
		while (!occupancy[kmax])	{
			kmax--;
		}
		while (occupancy[kmin])	{
			kmin++;
		}
		if (kmin < kmax)	{
			SwapComponents(kmin,kmax);
		}
	} while (kmin < kmax);
	int Kmin = kmin;
	if (Kmin > Kmax)	{
		Kmin--;
	}
	int K = Ncomponent;
	int nacc = 0;

	for (int rep=0; rep<nrep; rep++)	{
		double ratio = 1;
		int BK = K;
		if (K == Kmin)	{
			ratio *= 0.5 * ((double) K+1);
			ratio *= (K - Kmin + 1);
			ratio /= (GetNsite() + K * weightalpha);
			ratio *= K * weightalpha;
			K++;
		}
		else	{
			if (rnd::GetRandom().Uniform() < 0.5)	{
				ratio *= ((double) (K+1)) / (K - Kmin + 1);
				ratio /= (GetNsite() + K * weightalpha);
				ratio *= K * weightalpha;
				K++;
			}
			else	{
				if (K == Kmin + 1)	{
					ratio *= 2.0 / K;
					ratio *= (GetNsite() + (K - 1) * weightalpha);
					ratio /= (K-1) * weightalpha;
					K--;
				}
				else	{
					ratio *= ((double) (K - Kmin)) / K;
					ratio *= (GetNsite() + (K - 1) * weightalpha);
					ratio /= (K-1) * weightalpha;
					K--;
				}
			}
		}
		bool accept = ((rnd::GetRandom().Uniform() < ratio) && (K<=Kmax));
		if (! accept)	{
			K = BK;
		}
		else	{
			nacc ++;
		}
	}
	if (K > Ncomponent)	{
		for (int k=Ncomponent; k<K; k++)	{
			CreateComponent(k);
		}
	}
	else if (K < Ncomponent)	{
		for (int k=K; k<Ncomponent; k++)	{
			DeleteComponent(k);
		}
	}
	Ncomponent = K;
	return ((double) nacc) / nrep;
}


void FiniteProfileProcess::ReadStatFix(string filename)	{
	mixtype = filename;
	int Nstate = GetDim();
	if ((filename == "WLSR5") || (filename == "wlsr5"))	{
		if (Nstate != 20)	{
			cerr << "error: WLRS5 is for aminoacids\n";
		}
		Nfixcomp = WLSR5N;
		statfix = new double*[Nfixcomp];
		empweight = new double[Nfixcomp];
		for (int i=0; i<Nfixcomp-1; i++)	{
			statfix[i] = new double[Nstate];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] = WLSR5StatFix[i][k];
				if (statfix[i][k]<stateps)	{
					statfix[i][k] = stateps;
				}
				total += statfix[i][k];
			}
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] /= total;
			}
		}
		statfix[Nfixcomp-1] = new double[Nstate];	
		double* tmp = GetEmpiricalFreq();
		for (int k=0; k<Nstate; k++)	{
			statfix[Nfixcomp-1][k] = tmp[k];
		}
		for (int i=0; i<Nfixcomp; i++)	{
			empweight[i] = 1.0 / WLSR5N;
		}
	}
	else if (filename == "empfreq")	{
		Nfixcomp = 1;
		statfix = new double*[Nfixcomp];
		empweight = new double[Nfixcomp];
		double* freq = GetEmpiricalFreq();
		statfix[0] = new double[Nstate];
		for (int k=0; k<Nstate; k++)	{
			statfix[0][k] = freq[k];
		}
		empweight[0] = 1.0;
	}
	else if ((filename == "CG6") || (filename == "cg6") || (filename == "c6") || (filename == "C6"))	{
		if (Nstate != 20)	{
			cerr << "error: CG6 is for aminoacids\n";
		}
		Nfixcomp = 6;
		statfix = new double*[Nfixcomp];
		empweight = new double[Nfixcomp];
		for (int i=0; i<Nfixcomp; i++)	{
			statfix[i] = new double[Nstate];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] = CG6StatFix[i][k];
				if (statfix[i][k]<stateps)	{
					statfix[i][k] = stateps;
				}
				total += statfix[i][k];
			}
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] /= total;
			}
			empweight[i] = 1.0 / Nfixcomp;
		}
	}
	else if ((filename == "CG10") || (filename == "cg10"))	{
		if (Nstate != 20)	{
			cerr << "error: CG10 is for aminoacids\n";
		}
		Nfixcomp = 10;
		statfix = new double*[Nfixcomp];
		empweight = new double[Nfixcomp];
		for (int i=0; i<Nfixcomp; i++)	{
			statfix[i] = new double[Nstate];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] = CG10StatFix[i][k];
				if (statfix[i][k]<stateps)	{
					statfix[i][k] = stateps;
				}
				total += statfix[i][k];
			}
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] /= total;
			}
			empweight[i] = 1.0 / Nfixcomp;
		}
	}
	else if ((filename == "CG20") || (filename == "cg20"))	{
		if (Nstate != 20)	{
			cerr << "error: CG20 is for aminoacids\n";
		}
		Nfixcomp = 20;
		statfix = new double*[Nfixcomp];
		empweight = new double[Nfixcomp];
		for (int i=0; i<Nfixcomp; i++)	{
			statfix[i] = new double[Nstate];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] = CG20StatFix[i][k];
				if (statfix[i][k]<stateps)	{
					statfix[i][k] = stateps;
				}
				total += statfix[i][k];
			}
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] /= total;
			}
			empweight[i] = 1.0 / Nfixcomp;
		}
	}
	else if ((filename == "CG30") || (filename == "cg30"))	{
		if (Nstate != 20)	{
			cerr << "error: CG30 is for aminoacids\n";
		}
		Nfixcomp = 30;
		statfix = new double*[Nfixcomp];
		empweight = new double[Nfixcomp];
		for (int i=0; i<Nfixcomp; i++)	{
			statfix[i] = new double[Nstate];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] = CG30StatFix[i][k];
				if (statfix[i][k]<stateps)	{
					statfix[i][k] = stateps;
				}
				total += statfix[i][k];
			}
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] /= total;
			}
			empweight[i] = 1.0 / Nfixcomp;
		}
	}
	else if ((filename == "CG40") || (filename == "cg40"))	{
		if (Nstate != 20)	{
			cerr << "error: CG40is for aminoacids\n";
		}
		Nfixcomp = 40;
		statfix = new double*[Nfixcomp];
		empweight = new double[Nfixcomp];
		for (int i=0; i<Nfixcomp; i++)	{
			statfix[i] = new double[Nstate];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] = CG40StatFix[i][k];
				if (statfix[i][k]<stateps)	{
					statfix[i][k] = stateps;
				}
				total += statfix[i][k];
			}
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] /= total;
			}
			empweight[i] = 1.0 / Nfixcomp;
		}
	}
	else if ((filename == "CG50") || (filename == "cg50"))	{
		if (Nstate != 20)	{
			cerr << "error: CG50 is for aminoacids\n";
		}
		Nfixcomp = 50;
		statfix = new double*[Nfixcomp];
		empweight = new double[Nfixcomp];
		for (int i=0; i<Nfixcomp; i++)	{
			statfix[i] = new double[Nstate];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] = CG50StatFix[i][k];
				if (statfix[i][k]<stateps)	{
					statfix[i][k] = stateps;
				}
				total += statfix[i][k];
			}
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] /= total;
			}
			empweight[i] = 1.0 / Nfixcomp;
		}
	}
	else if ((filename == "CG60") || (filename == "cg60"))	{
		if (Nstate != 20)	{
			cerr << "error: CG60 is for aminoacids\n";
		}
		Nfixcomp = 60;
		statfix = new double*[Nfixcomp];
		empweight = new double[Nfixcomp];
		for (int i=0; i<Nfixcomp; i++)	{
			statfix[i] = new double[Nstate];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] = CG60StatFix[i][k];
				if (statfix[i][k]<stateps)	{
					statfix[i][k] = stateps;
				}
				total += statfix[i][k];
			}
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] /= total;
			}
			empweight[i] = 1.0 / Nfixcomp;
		}
	}
	else if ((filename == "c10") || (filename == "C10"))	{
		if (Nstate != 20)	{
			cerr << "error: C10 is for aminoacids\n";
		}
		Nfixcomp = C10N;
		statfix = new double*[Nfixcomp];
		empweight = new double[Nfixcomp];
		for (int i=0; i<Nfixcomp; i++)	{
			statfix[i] = new double[Nstate];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] = C10StatFix[i][k];
				if (statfix[i][k]<stateps)	{
					statfix[i][k] = stateps;
				}
				total += statfix[i][k];
			}
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] /= total;
			}
			empweight[i] = C10StatWeight[i];
		}
	}
	else if ((filename == "c20") || (filename == "C20"))	{
		if (Nstate != 20)	{
			cerr << "error: C20 is for aminoacids\n";
		}
		Nfixcomp = C20N;
		statfix = new double*[Nfixcomp];
		empweight = new double[Nfixcomp];
		for (int i=0; i<Nfixcomp; i++)	{
			statfix[i] = new double[Nstate];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] = C20StatFix[i][k];
				if (statfix[i][k]<stateps)	{
					statfix[i][k] = stateps;
				}
				total += statfix[i][k];
			}
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] /= total;
			}
			empweight[i] = C20StatWeight[i];
		}
	}
	else if ((filename == "c30") || (filename == "C30"))	{
		if (Nstate != 20)	{
			cerr << "error: C30 is for aminoacids\n";
		}
		Nfixcomp = C30N;
		statfix = new double*[Nfixcomp];
		empweight = new double[Nfixcomp];
		for (int i=0; i<Nfixcomp; i++)	{
			statfix[i] = new double[Nstate];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] = C30StatFix[i][k];
				if (statfix[i][k]<stateps)	{
					statfix[i][k] = stateps;
				}
				total += statfix[i][k];
			}
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] /= total;
			}
			empweight[i] = C30StatWeight[i];
		}
	}
	else if ((filename == "c30") || (filename == "C30"))	{
		if (Nstate != 20)	{
			cerr << "error: C30 is for aminoacids\n";
		}
		Nfixcomp = C30N;
		statfix = new double*[Nfixcomp];
		empweight = new double[Nfixcomp];
		for (int i=0; i<Nfixcomp; i++)	{
			statfix[i] = new double[Nstate];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] = C30StatFix[i][k];
				if (statfix[i][k]<stateps)	{
					statfix[i][k] = stateps;
				}
				total += statfix[i][k];
			}
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] /= total;
			}
			empweight[i] = C30StatWeight[i];
		}
	}
	else if ((filename == "c40") || (filename == "C40"))	{
		if (Nstate != 20)	{
			cerr << "error: C40 is for aminoacids\n";
		}
		Nfixcomp = C40N;
		statfix = new double*[Nfixcomp];
		empweight = new double[Nfixcomp];
		for (int i=0; i<Nfixcomp; i++)	{
			statfix[i] = new double[Nstate];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] = C40StatFix[i][k];
				if (statfix[i][k]<stateps)	{
					statfix[i][k] = stateps;
				}
				total += statfix[i][k];
			}
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] /= total;
			}
			empweight[i] = C40StatWeight[i];
		}
	}
	else if ((filename == "c50") || (filename == "C50"))	{
		if (Nstate != 20)	{
			cerr << "error: C50 is for aminoacids\n";
		}
		Nfixcomp = C50N;
		statfix = new double*[Nfixcomp];
		empweight = new double[Nfixcomp];
		for (int i=0; i<Nfixcomp; i++)	{
			statfix[i] = new double[Nstate];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] = C50StatFix[i][k];
				if (statfix[i][k]<stateps)	{
					statfix[i][k] = stateps;
				}
				total += statfix[i][k];
			}
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] /= total;
			}
			empweight[i] = C50StatWeight[i];
		}
	}
	else if ((filename == "c60") || (filename == "C60"))	{
		if (Nstate != 20)	{
			cerr << "error: C60 is for aminoacids\n";
		}
		Nfixcomp = C60N;
		statfix = new double*[Nfixcomp];
		empweight = new double[Nfixcomp];
		for (int i=0; i<Nfixcomp; i++)	{
			statfix[i] = new double[Nstate];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] = C60StatFix[i][k];
				if (statfix[i][k]<stateps)	{
					statfix[i][k] = stateps;
				}
				total += statfix[i][k];
			}
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] /= total;
			}
			empweight[i] = C60StatWeight[i];
		}
	}
	else if ((filename == "lg") || (filename == "LG"))	{
		Nfixcomp = 1;
		statfix = new double*[Nfixcomp];
		empweight = new double[Nfixcomp];
		statfix[0] = new double[Nstate];
		double total = 0;
		for (int k=0; k<Nstate; k++)	{
			statfix[0][k] = LG_Stat[k];
			total += statfix[0][k];
		}
		for (int k=0; k<Nstate; k++)	{
			statfix[0][k] /= total;
		}
		empweight[0] = 1.0;
	}
	else if ((filename == "uniform") || (filename == "Uniform"))	{
		Nfixcomp = 1;
		statfix = new double*[Nfixcomp];
		empweight = new double[Nfixcomp];
		for (int i=0; i<Nfixcomp; i++)	{
			statfix[i] = new double[Nstate];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] = 1.0 / Nstate;
				if (statfix[i][k]<stateps)	{
					statfix[i][k] = stateps;
				}
				total += statfix[i][k];
			}
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] /= total;
			}
			empweight[i] = 1.0;
		}
	}
	else	{
		ifstream is(filename.c_str());
		if (!is)	{
			cerr << "error : unrecognized file for empirical mixture : " << filename << '\n';
			exit(1);
		}
		int tmp;
		is >> tmp;
		if (tmp != Nstate)	{
			cerr << "error when reading empirical mixture : bad number of states\n";
			exit(1);
		}
		// read alphabet
		int permut[Nstate];
		for (int k=0; k<Nstate; k++)	{
			char c;
			is >> c;
			int l=0;
			// while ((l<Nstate) && (c != GetStateSpace()->GetCharState(l))) l++;
			while ((l<Nstate) && (c != AminoAcids[l])) l++;
			if (l == Nstate)	{
				cerr << "error when reading empirical mixture in " << filename << ": does not recognise letter " << c << '\n';
				cerr << "file should be formatted as follows\n";
				cerr << "list of amino acids (e.g. A C D E ...)\n";
				cerr << "Ncat\n";
				cerr << "empweight_1 (1 number) profile_1 (20 numbers)\n";
				cerr << "empweight_2 (1 number) profile_2 (20 numbers)\n";
				cerr << "...\n";
				cerr << "empweight_Ncat (1 number) profile_Ncat (20 numbers)\n";
				for (int a=0; a<Nstate; a++)	{
					cerr << "::" << GetStateSpace()->GetState(a) << "::";
				}
				cerr << '\n';
				exit(1); 
			}
			permut[k] = l;
		}
		is >> Nfixcomp;
		statfix = new double*[Nfixcomp];
		empweight = new double[Nfixcomp];
		for (int i=0; i<Nfixcomp; i++)	{
			statfix[i] = new double[Nstate];
			is >> empweight[i];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				is >> statfix[i][permut[k]];
				if (statfix[i][permut[k]]<stateps)	{
					statfix[i][permut[k]] = stateps;
				}
				total += statfix[i][permut[k]];
			}
			for (int k=0; k<Nstate; k++)	{
				statfix[i][k] /= total;
			}
		}
		double total = 0;
		for (int i=0; i<Nfixcomp; i++)	{
			total += empweight[i];
		}
		for (int i=0; i<Nfixcomp; i++)	{
			empweight[i] /= total;
		}
	}
}

void FiniteProfileProcess::SetStatFix()	{

	if (! Nfixcomp)	{
		cerr << "error in set stat fix\n";
		exit(1);
	}

	for (int k=0; k<Ncomponent; k++)	{
		DeleteComponent(k);
	}
	Ncomponent = Nfixcomp;
	for (int k=0; k<Ncomponent; k++)	{
		CreateComponent(k);
	}
	
	for (int k=0; k<Ncomponent; k++)	{
		weight[k] = empweight[k];
		for (int i=0; i<GetDim(); i++)	{
			profile[k][i] = statfix[k][i];
		}
	}
	fixncomp = true;
	if (empmix == 3)	{
		SetStatCons();
	}
}


double FiniteProfileProcess::LogStatPriorConstrained(int cat)	{

	if (! Ncons)	{
		return MixtureProfileProcess::LogStatPrior(cat);
	}

	int ok = 1;
	for (int k=0; k<GetDim(); k++)	{
		ok &= (profile[cat][k] >= conscutoff[statcons[cat][k]]);
		ok &= (profile[cat][k] <= conscutoff[statcons[cat][k] + 1]);
	}

	if (ok)	{
		return 0;
	}
	return log(0);
}

void FiniteProfileProcess::SetStatCons()	{

	Ncons = 2;
	conscutoff = new double[Ncons+1];
	conscutoff[0] = 0;
	conscutoff[1] = 0.2;
	conscutoff[2] = 1.0;

	statcons = new int*[Nfixcomp];
	for (int i=0; i<Nfixcomp; i++)	{
		statcons[i] = new int[GetDim()];
		for (int k=0; k<GetDim(); k++)	{
			int j = 0;
			while ((j <Ncons) && ((statfix[i][k] < conscutoff[j]) || (statfix[i][k] > conscutoff[j+1]))) j++;
			if (j == Ncons)	{
				cerr << "error in set constraints\n";
				exit(1);
			}
			statcons[i][k] = j;
		}
	}
}

/*
void FiniteProfileProcess::ReadConstraints(string filename)	{

		ifstream is(filename.c_str());
		if (!is)	{
			cerr << "error : unrecognized file for empirical mixture : " << filename << '\n';
			exit(1);
		}
		int tmp;
		is >> tmp;
		if (tmp != Nstate)	{
			cerr << "error when reading empirical mixture : bad number of states\n";
			exit(1);
		}
		// read alphabet
		int permut[Nstate];
		for (int k=0; k<Nstate; k++)	{
			char c;
			is >> c;
			int l=0;
			// while ((l<Nstate) && (c != GetStateSpace()->GetCharState(l))) l++;
			while ((l<Nstate) && (c != AminoAcids[l])) l++;
			if (l == Nstate)	{
				cerr << "error when reading empirical mixture in " << filename << ": does not recognise letter " << c << '\n';
				cerr << "file should be formatted as follows\n";
				cerr << "list of amino acids (e.g. A C D E ...)\n";
				cerr << "Ncat\n";
				cerr << "empweight_1 (1 number) profile_1 (20 numbers)\n";
				cerr << "empweight_2 (1 number) profile_2 (20 numbers)\n";
				cerr << "...\n";
				cerr << "empweight_Ncat (1 number) profile_Ncat (20 numbers)\n";
				for (int a=0; a<Nstate; a++)	{
					cerr << "::" << GetStateSpace()->GetState(a) << "::";
				}
				cerr << '\n';
				exit(1); 
			}
			permut[k] = l;
		}
		is >> Ncons;
		conscutoff = new double[Ncons+1];
		for (int i=0; i<=Ncons; i++)	{
			is >> conscutoff[i];
		}
		
		is >> Ncat;
		statcons = new int*[Ncat];
		for (int i=0; i<Ncat; i++)	{
			statcons[i] = new int[Nstate];
			for (int k=0; k<Nstate; k++)	{
				is >> statcons[i][permut[k]];
			}
		}
	}
}
*/

double FiniteProfileProcess::SMCAddSites()	{

	double logw = 0;

	for (int site = GetSiteMin(); site<GetSiteMax(); site++)	{
	// for (int site = GetBKSiteMax(); site<GetSiteMax(); site++)	{

		if (NewlyActivated(site))	{

			if (Ncomponent > 1)	{

				double logl[Ncomponent];
				double max = 0;
				int inf = 0;
				for (int k=0; k<Ncomponent; k++)	{
					AddSite(site,k);
					logl[k] = SiteLogLikelihood(site);
					if (isinf(logl[k]))	{
						cerr << "inf\n";
					}
					if ((!k) || (max < logl[k]))	{
						max = logl[k];
					}
					RemoveSite(site,k);
				}

				if (isinf(max))	{
					cerr << "error in smc add sites\n";
					for (int k=0; k<Ncomponent; k++)	{
						cerr << logl[k] << '\n';
					}
					exit(1);
				}

				double p[Ncomponent];
				double cumul[Ncomponent];
				double tot = 0;
				for (int k=0; k<Ncomponent; k++)	{
					double tmp = 0;
					if (! isinf(logl[k]))	{
						tmp = weight[k] * exp(logl[k] - max);
					}
					p[k] = tmp;
					tot += tmp;
					cumul[k] = tot;
				}
				if (isinf(tot))	{
					cerr << "in add site: tot is inf\n";
					exit(1);
				}
				if (isnan(tot))	{
					cerr << "in add site: tot is nan\n";
					exit(1);
				}
				logw += log(tot) + max;
				double u = tot * rnd::GetRandom().Uniform();
				int k = 0;
				while ((k<Ncomponent) && (cumul[k] < u))	{
					k++;
				}
				if (k == Ncomponent)	{
					cerr << "error in FiniteProfileProcess::SMCAddSites:overflow\n";
					exit(1);
				}
				AddSite(site,k);
				SiteLogLikelihood(site);
				SampleSiteMapping(site);
			}
			else	{

				AddSite(site,0);
				logw += SiteLogLikelihood(site);
				SampleSiteMapping(site);
			}

		}
	}

	if (GetMyid())	{
		MPI_Send(&logw,1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
		if (Ncomponent > 1)	{
			MPI_Send(alloc,GetNsite(),MPI_INT,0,TAG1,MPI_COMM_WORLD);
		}
	}

	return logw;
}

double FiniteProfileProcess::GlobalSMCAddSites()	{

	double ret = ProfileProcess::GlobalSMCAddSites();

	// receive new site allocations from slave
	if (Ncomponent > 1)	{
		MPI_Status stat;
		int tmpalloc[GetNsite()];
		for(int i=1; i<GetNprocs(); ++i) {
			MPI_Recv(tmpalloc,GetNsite(),MPI_INT,i,TAG1,MPI_COMM_WORLD,&stat);
			for(int j=GetProcSiteMin(i); j<GetProcSiteMax(i); ++j) {
				if (ActiveSite(j))	{
					alloc[j] = tmpalloc[j];
					if ((tmpalloc[j] < 0) || (tmpalloc[j] >= Ncomponent))	{
						cerr << "in SMC add\n";
						cerr << "alloc overflow\n";
						cerr << tmpalloc[j] << '\n';
						exit(1);
					}
				}
			}
		}
		UpdateOccupancyNumbers();
		ResampleWeights();
	}
	return ret;
}
