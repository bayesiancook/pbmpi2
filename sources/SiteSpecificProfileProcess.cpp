
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "SiteSpecificProfileProcess.h"
#include "Random.h"
#include "Parallel.h"


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* SiteSpecificProfileProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void SiteSpecificProfileProcess::Create()	{
	if (! profile)	{
		if (! GetDim())	{
			cerr << "error: ProfileProcess::dim has not been initialized\n";
			exit(1);
		}
		ProfileProcess::Create();
		allocprofile = new double[GetNsite() * GetDim()];
		profile = new double*[GetNsite()];
		for (int i=0; i<GetNsite(); i++)	{
			profile[i] = allocprofile + i*GetDim();
		}
		logstatprior = new double[GetNsite()];
		profilesuffstatlogprob = new double[GetNsite()];
	}
}

void SiteSpecificProfileProcess::Delete()	{
	if (profile)	{
		delete[] profilesuffstatlogprob;
		delete[] logstatprior;
		delete[] allocprofile;
		delete[] profile;
		profile = 0;
		ProfileProcess::Delete();
	}
}

double SiteSpecificProfileProcess::GetStatEnt()	{
	double total = 0;
	int totnsite = 0;
	for (int k=0; k<GetNsite(); k++)	{
		if (ActiveSite(k))	{
			total += GetStatEnt(k);
			totnsite ++;
		}
	}
	return total / totnsite;
}

double SiteSpecificProfileProcess::GetStatEnt(int k)	{
	double total = 0;
	for (int i=0; i<GetDim(); i++)	{
		if (profile[k][i] <= 0)	{
			cerr << "error: 0 entry in profile\n";
			cerr << profile[k][i] << '\n';
			exit(1);
		}
		total -= profile[k][i] * log(profile[k][i]);
	}
	if (isnan(total))	{
		cerr << "entropy is nan\n";
		exit(1);
	}
	return  total;
}

void SiteSpecificProfileProcess::RenormalizeProfiles()	{
	for (int i=0; i<GetNsite(); i++)	{
		double total = 0;
		for (int k=0; k<GetDim(); k++)	{
			total += profile[i][k];
		}
		for (int k=0; k<GetDim(); k++)	{
			profile[i][k] /= total;
		}
	}
}

void SiteSpecificProfileProcess::SampleProfile()	{
	SampleHyper();
	SampleStat();
}

void SiteSpecificProfileProcess::SampleStat()	{
	for (int i=0; i<GetNsite(); i++)	{
		SampleStat(i);
	}
}

void SiteSpecificProfileProcess::SampleStat(int i)	{
	SampleFrequencyStat(profile[i]);
}

double SiteSpecificProfileProcess::LogProfilePrior()	{
	double total = 0;
	total += LogHyperPrior();
	total += LogStatPrior();
	return total;
}

double SiteSpecificProfileProcess::LogStatPrior()	{

	for (int i=0; i<GetNsite(); i++)	{
		if (ActiveSite(i))	{
			LogStatPrior(i);
		}
	}
	double total = 0;
	for (int i=0; i<GetNsite(); i++)	{
		if (ActiveSite(i))	{
			total += logstatprior[i];
		}
	}
	return total;
}

double SiteSpecificProfileProcess::LogStatPrior(int site)	{
	double tmp = LogFrequencyStatPrior(profile[site]);
	logstatprior[site] = tmp;
	return tmp;
}

double SiteSpecificProfileProcess::ProfileSuffStatLogProb()	{

	// simply, sum over all components
	for (int i=0; i<GetNsite(); i++)	{
		if (ActiveSite(i))	{
			ProfileSuffStatLogProb(i);
		}
	}
	double total = 0;
	for (int i=0; i<GetNsite(); i++)	{
		if (ActiveSite(i))	{
			total += profilesuffstatlogprob[i];
		}
	}
	return total;
}

double SiteSpecificProfileProcess::Move(double tuning, int n, int nrep)	{

	double ret = 0;
	if (GetNprocs() > 1)	{
		ret = MPIMove(tuning,n,nrep);
	}
	else	{
		ret = NonMPIMove(tuning,n,nrep);
	}
	return ret;
}

double SiteSpecificProfileProcess::MPIMove(double tuning, int n, int nrep)	{

	totchrono.Start();

	for (int rep=0; rep<nrep; rep++)	{

		GlobalParametersMove();

		profilechrono.Start();

		GlobalUpdateParameters();
		GlobalUpdateSiteProfileSuffStat();
		GlobalMoveProfile(tuning,n,nrep);

		profilechrono.Stop();

		// hyperparameters
		GlobalUpdateParameters();
		GlobalUpdateSiteProfileSuffStat();
		MoveHyper(tuning,10);
	}
	totchrono.Stop();
	return 1;
}

double SiteSpecificProfileProcess::NonMPIMove(double tuning, int n, int nrep)	{

	totchrono.Start();

	for (int rep=0; rep<nrep; rep++)	{

		GlobalParametersMove();

		profilechrono.Start();

		UpdateSiteProfileSuffStat();
		MoveProfile(tuning,n,nrep);

		profilechrono.Stop();

		UpdateSiteProfileSuffStat();
		MoveHyper(tuning,10);
	}
	totchrono.Stop();
	return 1;
}

double SiteSpecificProfileProcess::GlobalMoveProfile(double tuning, int n, int nrep)	{

	// send PROFILE_MOVE Message with n and nrep and tuning
	MESSAGE signal = PROFILE_MOVE;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	int* itmp = new int[2];
	itmp[0] = n;
	itmp[1] = nrep;
	MPI_Bcast(itmp,2,MPI_INT,0,MPI_COMM_WORLD);
	delete[] itmp;

	MPI_Bcast(&tuning,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

	// collect final values of profiles (+ total acceptance rate) from slaves
	MPI_Status stat;
	double total = 0;
	for(int i=1; i<GetNprocs(); ++i) {
		MPI_Recv(allocprofile + GetProcSiteMin(i)*GetDim(),GetProcSiteNumber(i)*GetDim(),MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
		double tmp = 0;
		MPI_Recv(&tmp,1,MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
		total += tmp; // (sum all acceptance rates)
	}
	// return average acceptance rate
	return total / GetNsite() / nrep;
}

void SiteSpecificProfileProcess::SlaveMoveProfile()	{

	// parse arguments sent by master

	int* itmp = new int[2];
	MPI_Bcast(itmp,2,MPI_INT,0,MPI_COMM_WORLD);
	int n = itmp[0];
	int nrep = itmp[1];
	delete[] itmp;

	double tuning = 0;
	MPI_Bcast(&tuning,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

	// move components in the range just computed
	double total = 0;
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i))	{
			total += MoveProfile(i,tuning,n,nrep);
		}
	}

	MPI_Send(allocprofile + GetSiteMin()*GetDim(),GetSiteNumber()*GetDim(),MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
	MPI_Send(&total,1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

double SiteSpecificProfileProcess::MoveProfile(double tuning, int n, int nrep)	{
	double total = 0;
	int count = 0;
	for (int i=0; i<GetNsite(); i++)	{
		if (ActiveSite(i))	{
			total += MoveProfile(i,tuning,n,nrep);
			count++;
		}
	}
	return total / count;
}

double SiteSpecificProfileProcess::MoveProfile(int site, double tuning, int n, int nrep)	{

	int naccepted = 0;
	double* bk = new double[GetDim()];
	for (int k=0; k<GetDim(); k++)	{
		bk[k] = profile[site][k];
	} 
	for (int rep=0; rep<nrep; rep++)	{
		double deltalogprob = - LogStatPrior(site) - ProfileSuffStatLogProb(site);
		double loghastings = ProfileProposeMove(profile[site],tuning,n,0,site,0);
		UpdateSite(site);
		deltalogprob += LogStatPrior(site) + ProfileSuffStatLogProb(site);
		deltalogprob += loghastings;
		int accepted = (log(rnd::GetRandom().Uniform()) < deltalogprob);
		if (accepted)	{
			naccepted ++;
			for (int k=0; k<GetDim(); k++)	{
				bk[k] = profile[site][k];
			} 
		}
		else	{
			for (int k=0; k<GetDim(); k++)	{
				profile[site][k] = bk[k];
			} 
			UpdateSite(site);
		}
	}
	delete[] bk;
	return naccepted / nrep;
}
