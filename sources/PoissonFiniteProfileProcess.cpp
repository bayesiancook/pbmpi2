
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "PoissonFiniteProfileProcess.h"
#include "Random.h"
#include "Parallel.h"


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* PoissonFiniteProfileProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------


double PoissonFiniteProfileProcess::Move(double tuning, int n, int nrep)	{

	double ret = 0;
	if (GetNprocs() > 1)	{
		ret = MPIMove(tuning,n,nrep);
	}
	else	{
		ret = NonMPIMove(tuning,n,nrep);
	}
	return ret;
}

double PoissonFiniteProfileProcess::MPIMove(double tuning, int n, int nrep)	{

	for (int rep=0; rep<nrep; rep++)	{
		if (Ncomponent > 1)	{
			GlobalUpdateParameters();
			GlobalUpdateSiteProfileSuffStat();
			UpdateModeProfileSuffStat();
			GlobalIncrementalFiniteMove(1);
		}

		if (! empmix)	{
			GlobalUpdateParameters();
			GlobalUpdateSiteProfileSuffStat();
			UpdateModeProfileSuffStat();
			MoveProfile();
			GlobalUpdateParameters();
			GlobalUpdateSiteProfileSuffStat();
			MoveHyper(tuning,10);
		}
	}
	return 1;
}

double PoissonFiniteProfileProcess::NonMPIMove(double tuning, int n, int nrep)	{

	for (int rep=0; rep<nrep; rep++)	{

		if (Ncomponent > 1)	{
			UpdateSiteProfileSuffStat();
			UpdateModeProfileSuffStat();
			IncrementalFiniteMove(1);
		}

		if (! empmix)	{
			UpdateSiteProfileSuffStat();
			UpdateModeProfileSuffStat();
			MoveProfile();
			UpdateComponents();
			UpdateSiteProfileSuffStat();
			MoveHyper(tuning,10);
		}
	}
	return 1;
}


void PoissonFiniteProfileProcess::ToStream(ostream& os)	{

	os << Ncomponent << '\n';
	for (int j=0; j<GetDim(); j++)	{
		os << dirweight[j] << '\t';
	}
	os << '\n';
	os << '\n';
	for (int i=0; i<Ncomponent; i++)	{
		for (int j=0; j<GetDim(); j++)	{
			os << profile[i][j] << '\t';
		}
		os << '\n';
	}

	for (int i=0; i<GetNsite(); i++)	{
		if (ActiveSite(i))	{
			os << alloc[i] << '\t';
		}
		else	{
			os << -1 << '\t';
		}
	}
	os << '\n';
}

void PoissonFiniteProfileProcess::FromStream(istream& is)	{

	is >> Ncomponent;
	
	for (int i=0; i<GetDim(); i++)	{
		is >> dirweight[i];
	}

	for (int i=0; i<Ncomponent; i++)	{
		for (int j=0; j<GetDim(); j++)	{
			is >> profile[i][j];
		}
	}

	for (int i=0; i<GetNsite(); i++)	{
		is >> alloc[i];
	}

	ResampleWeights();
	// CHECK some update here ?
}


double PoissonFiniteProfileProcess::GlobalIncrementalFiniteMove(int nrep)	{

	// send command and arguments
	MESSAGE signal = REALLOC_MOVE;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&nrep,1,MPI_INT,0,MPI_COMM_WORLD);

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
			for(int j=GetProcSiteMin(i); j<GetProcSiteMax(i); j++) {
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
	UpdateModeProfileSuffStat();

	return 1;
}

double PoissonFiniteProfileProcess::SlaveIncrementalFiniteMove()	{

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
				Ntried++;
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

double PoissonFiniteProfileProcess::IncrementalFiniteMove(int nrep)	{

	if (GetMyid())	{
		cerr << "error: slave in PoissonFiniteProfileProcess::IncrementalFiniteMove\n";
		exit(1);
	}

	int NAccepted = 0;

	double* bigarray = new double[Ncomponent * GetNsite()];
	double* bigcumul = new double[Ncomponent * GetNsite()];

	for (int rep=0; rep<nrep; rep++)	{

		// do the incremental reallocation move on my site range
		for (int site=GetSiteMin(); site<GetSiteMax(); site++)	{

			if (ActiveSite(site))	{

				double* mLogSamplingArray = bigarray + site * Ncomponent;
				double* cumul = bigcumul + site * Ncomponent;

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
				alloc[site] = mode;

			}
		}
	}
	
	delete[] bigarray;
	delete[] bigcumul;
	return ((double) NAccepted) / GetNsite() / nrep;
}

