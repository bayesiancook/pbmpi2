
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "ExpoConjugateGTRMixtureProfileProcess.h"
#include "Random.h"
#include "Parallel.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* ExpoConjugateGTRMixtureProfileProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------


void ExpoConjugateGTRMixtureProfileProcess::Create()	{
	if (! profilesuffstatcount)	{
		DirichletProfileProcess::Create();
		ExpoConjugateGTRProfileProcess::Create();
		GTRMixtureProfileProcess::Create();
		profilesuffstatcount = new int*[GetNmodeMax()];
		profilesuffstatbeta = new double*[GetNmodeMax()];
		allocprofilesuffstatcount = new int[GetNmodeMax() * GetDim()];
		allocprofilesuffstatbeta = new double[GetNmodeMax() * GetDim()];
		for (int i=0; i<GetNmodeMax(); i++)	{
			profilesuffstatcount[i] = allocprofilesuffstatcount + i*GetDim();
			profilesuffstatbeta[i] = allocprofilesuffstatbeta + i*GetDim();
		}
	}
}

void ExpoConjugateGTRMixtureProfileProcess::Delete() {
	if (profilesuffstatcount)	{
		delete[] profilesuffstatcount;
		delete[] profilesuffstatbeta;
		delete[] allocprofilesuffstatcount;
		delete[] allocprofilesuffstatbeta;
		profilesuffstatcount = 0;
		profilesuffstatbeta = 0;
		GTRMixtureProfileProcess::Delete();
		DirichletProfileProcess::Delete();
	}
}

void ExpoConjugateGTRMixtureProfileProcess::GlobalUpdateModeProfileSuffStat()	{

	UpdateModeProfileSuffStat();

	if (GetNprocs() > 1)	{
		MESSAGE signal = UPDATE_MPROFILE;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

		MPI_Bcast(allocprofilesuffstatcount,Ncomponent*GetDim(),MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(allocprofilesuffstatbeta,Ncomponent*GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	}

}

void ExpoConjugateGTRMixtureProfileProcess::SlaveUpdateModeProfileSuffStat()	{

	MPI_Bcast(allocprofilesuffstatcount,Ncomponent*GetDim(),MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(allocprofilesuffstatbeta,Ncomponent*GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);
}

void ExpoConjugateGTRMixtureProfileProcess::UpdateModeProfileSuffStat()	{

	if (GetMyid())	{
		cerr << "error: slave in ExpoConjugateGTRMixtureProfileProcess::UpdateModeProfileSuffStat()\n";
		exit(1);
	}
	for (int i=0; i<GetNcomponent(); i++)	{
		for (int k=0; k<GetDim(); k++)	{
			profilesuffstatcount[i][k] = 0;
			profilesuffstatbeta[i][k] = 0;
		}
	}
	for (int i=0; i<GetNsite(); i++)	{
		if (ActiveSite(i))	{
			const int* count = GetSiteProfileSuffStatCount(i);
			const double* beta = GetSiteProfileSuffStatBeta(i);
			int cat = alloc[i];
			for (int k=0; k<GetDim(); k++)	{
				profilesuffstatcount[cat][k] += count[k];
				profilesuffstatbeta[cat][k] += beta[k];
			}
		}
	}
}

double ExpoConjugateGTRMixtureProfileProcess::ProfileSuffStatLogProb(int cat)	{
	double total = 0;
	for (int k=0; k<GetDim(); k++)	{
		total += profilesuffstatcount[cat][k] * log(profile[cat][k]) - profilesuffstatbeta[cat][k] * profile[cat][k];
	}
	profilesuffstatlogprob[cat] = total;
	return total;
}

/*
double ExpoConjugateGTRMixtureProfileProcess::ProfileProposeMove(double* profile, double tuning, int n, int K, int cat, double statmin)	{

	if (! proposemode)	{
		return ProfileProcess::ProfileProposeMove(profile,tuning,n,K,cat,statmin);
	}

	if (! statmin)	{
		statmin = stateps;
	}

	double bkprofile[GetDim()];

	int totcount = 0;
	for (int k=0; k<GetDim(); k++)	{
		totcount += profilesuffstatcount[cat][k];
	}
	int nstep = 5 + 100 * totcount;
	// int nstep = 1000;

	double logHastings = 0;

	double logp1 = 0;
	for (int k=0; k<GetDim(); k++)	{
		logp1 -= profilesuffstatbeta[cat][k] * profile[k];
		logp1 += (profilesuffstatcount[cat][k] + dirweight[k] - 1) * log(profile[k]);
	}

	for (int step=nstep-1; step>=0; step--)	{
		double gamma = 1.0 / nstep * (step+0.5);
		// one MH step
		// using p_gamma as a target
		for (int k=0; k<GetDim(); k++)	{
			bkprofile[k] = profile[k];
		}
		double logratio = 0;
		for (int k=0; k<GetDim(); k++)	{
			logratio += gamma * profilesuffstatbeta[cat][k] * profile[k];
			logratio -= (profilesuffstatcount[cat][k] + dirweight[k] - 1) * log(profile[k]);
			// logratio -= gamma * (profilesuffstatcount[cat][k] + dirweight[k] - 1) * log(profile[k]);
		}
		double logh = ProfileProcess::ProfileProposeMove(profile,tuning,n,K,cat,statmin);
		for (int k=0; k<GetDim(); k++)	{
			logratio -= gamma * profilesuffstatbeta[cat][k] * profile[k];
			logratio += (profilesuffstatcount[cat][k] + dirweight[k] - 1) * log(profile[k]);
			// logratio += gamma * (profilesuffstatcount[cat][k] + dirweight[k] - 1) * log(profile[k]);
		}
		logratio += logh;
		int accept = (log(rnd::GetRandom().Uniform()) < logratio);
		if (accept)	{
			logHastings -= logratio - logh;
		}
		else	{
			for (int k=0; k<GetDim(); k++)	{
				profile[k] = bkprofile[k];
			}
		}
	}

	for (int k=0; k<GetDim(); k++)	{
		logHastings += (profilesuffstatcount[cat][k] + dirweight[k] - 1) * log(profile[k]);
	}
	// draw from prior
	double tot = 0;
	for (int k=0; k<GetDim(); k++)	{
		profile[k] = rnd::GetRandom().sGamma(profilesuffstatcount[cat][k] + dirweight[k]);
		if (profile[k] < statmin)	{
			profile[k] = statmin;
		}
		tot += profile[k];
	}
	for (int k=0; k<GetDim(); k++)	{
		profile[k] /= tot;
		logHastings -= (profilesuffstatcount[cat][k] + dirweight[k] - 1) * log(profile[k]);
	}

	for (int step=0; step<nstep; step++)	{
		double gamma = 1.0 / nstep * (step+0.5);
		// one MH step
		// using p_gamma as a target
		for (int k=0; k<GetDim(); k++)	{
			bkprofile[k] = profile[k];
		}
		double logratio = 0;
		for (int k=0; k<GetDim(); k++)	{
			logratio += gamma * profilesuffstatbeta[cat][k] * profile[k];
			logratio -= (profilesuffstatcount[cat][k] + dirweight[k] - 1) * log(profile[k]);
			// logratio -= gamma * (profilesuffstatcount[cat][k] + dirweight[k] - 1) * log(profile[k]);
		}
		double logh = ProfileProcess::ProfileProposeMove(profile,tuning,n,K,cat,statmin);
		for (int k=0; k<GetDim(); k++)	{
			logratio -= gamma * profilesuffstatbeta[cat][k] * profile[k];
			logratio += (profilesuffstatcount[cat][k] + dirweight[k] - 1) * log(profile[k]);
			// logratio += gamma * (profilesuffstatcount[cat][k] + dirweight[k] - 1) * log(profile[k]);
		}
		logratio += logh;
		int accept = (log(rnd::GetRandom().Uniform()) < logratio);
		if (accept)	{
			logHastings -= logratio - logh;
		}
		else	{
			for (int k=0; k<GetDim(); k++)	{
				profile[k] = bkprofile[k];
			}
		}
	}

	double logp2 = 0;
	for (int k=0; k<GetDim(); k++)	{
		logp2 -= profilesuffstatbeta[cat][k] * profile[k];
		logp2 += (profilesuffstatcount[cat][k] + dirweight[k] - 1) * log(profile[k]);
	}

	double logmarginal = logp2 - logp1 + logHastings;
	if (isnan(logmarginal))	{
		cerr << "error in propose profile: nan\n";
		cerr << logp1 << '\t' << logp2 << '\t' << logHastings << '\n';
		exit(1);
	}
	if (isinf(logmarginal))	{
		cerr << "error in propose profile: inf\n";
		cerr << logp1 << '\t' << logp2 << '\t' << logHastings << '\n';
		exit(1);
	}
	meanlogmarginal += logmarginal;
	varlogmarginal += logmarginal * logmarginal;
	countlogmarginal++;

	return logHastings;
}
*/

/*
double ExpoConjugateGTRMixtureProfileProcess::ProfileProposeMove(double* profile, double tuning, int n, int K, int cat, double statmin)	{

	if (! proposemode)	{
		return ProfileProcess::ProfileProposeMove(profile,tuning,n,K,cat,statmin);
	}

	if (! statmin)	{
		statmin = stateps;
	}

	double alpha[GetDim()];
	double totweight = 0;
	double totalpha = 0;
	for (int k=0; k<GetDim(); k++)	{
		double tmp =  (profilesuffstatcount[cat][k] + dirweight[k]);
		totweight += tmp;
		tmp /= profilesuffstatbeta[cat][k];
		alpha[k] = tmp;
		totalpha += tmp;
	}
	double targetweight = 0;
	for (int k=0; k<GetDim(); k++)	{
		// alpha[k] *= totweight / totalpha;
		alpha[k] *= tuning * totweight / totalpha;
		targetweight += alpha[k];
	}

	double logHastings = 0;
	for (int k=0; k<GetDim(); k++)	{
		logHastings += (alpha[k] - 1) * log(profile[k]);
		// logHastings += (alpha[k] - 1) * log(profile[k]) - rnd::GetRandom().logGamma(alpha[k]);
	}


	double total = 0;
	for (int k=0; k<GetDim(); k++)	{
		profile[k] = rnd::GetRandom().sGamma(alpha[k]);
		if (profile[k] < statmin)	{
			profile[k] = statmin;
		}
		total += profile[k];
	}

	for (int k=0; k<GetDim(); k++)	{
		profile[k] /= total;
		logHastings -= (alpha[k] - 1) * log(profile[k]);
		// logHastings -= (alpha[k] - 1) * log(profile[k]) - rnd::GetRandom().logGamma(alpha[k]);
	}

	return logHastings;

}
*/

void ExpoConjugateGTRMixtureProfileProcess::SwapComponents(int cat1, int cat2)	{

	MixtureProfileProcess::SwapComponents(cat1,cat2);
	for (int k=0; k<GetDim(); k++)	{

		int tmp = profilesuffstatcount[cat1][k];
		profilesuffstatcount[cat1][k] = profilesuffstatcount[cat2][k];
		profilesuffstatcount[cat2][k] = tmp;

		double temp = profilesuffstatbeta[cat1][k];
		profilesuffstatbeta[cat1][k] = profilesuffstatbeta[cat2][k];
		profilesuffstatbeta[cat2][k] = temp;
	}
}


//-------------------------------------------------------------------------
//	* udpate eq. frequency profiles based on sufficient statistics
//	(CPU level 3)
//-------------------------------------------------------------------------

double ExpoConjugateGTRMixtureProfileProcess::LogStatProb(int site, int cat)	{
	const int* count = GetSiteProfileSuffStatCount(site);
	const double* beta = GetSiteProfileSuffStatBeta(site);
	double* pi = profile[cat];
	double total = 0;
	for (int k=0; k<GetDim(); k++)	{
		total += count[k] * log(pi[k]) - beta[k] * pi[k];
	}
	return total;
}


