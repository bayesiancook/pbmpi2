
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "ExpoConjugatePartitionGTRMixtureProfileProcess.h"
#include "Random.h"
#include "Parallel.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* ExpoConjugatePartitionGTRMixtureProfileProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------


void ExpoConjugatePartitionGTRMixtureProfileProcess::Create()	{
	if (! profilesuffstatcount)	{
		DirichletProfileProcess::Create();
		ExpoConjugatePartitionGTRProfileProcess::Create();
		PartitionGTRMixtureProfileProcess::Create();
		profilesuffstatcount = new double*[GetNmodeMax()];
		profilesuffstatbeta = new double*[GetNmodeMax()];
		allocprofilesuffstatcount = new double[GetNmodeMax() * GetDim()];
		allocprofilesuffstatbeta = new double[GetNmodeMax() * GetDim()];
		tmpprofilesuffstatcount = new double*[GetNmodeMax()];
		tmpprofilesuffstatbeta = new double*[GetNmodeMax()];
		alloctmpprofilesuffstatcount = new double[GetNmodeMax() * GetDim()];
		alloctmpprofilesuffstatbeta = new double[GetNmodeMax() * GetDim()];
		for (int i=0; i<GetNmodeMax(); i++)	{
			profilesuffstatcount[i] = allocprofilesuffstatcount + i*GetDim();
			profilesuffstatbeta[i] = allocprofilesuffstatbeta + i*GetDim();
			tmpprofilesuffstatcount[i] = alloctmpprofilesuffstatcount + i*GetDim();
			tmpprofilesuffstatbeta[i] = alloctmpprofilesuffstatbeta + i*GetDim();
			// initialize to 0: this is not useful for profilesuffstat,
			// but very important for tmpprofilesuffstat (used in an MPI_Reduce call, and thus, should be initialized)
			for (int j=0; j<GetDim(); j++)	{
				profilesuffstatcount[i][j] = 0;
				profilesuffstatbeta[i][j] = 0;
				tmpprofilesuffstatcount[i][j] = 0;
				tmpprofilesuffstatbeta[i][j] = 0;
			}
		}
	}
}

void ExpoConjugatePartitionGTRMixtureProfileProcess::Delete() {
	if (profilesuffstatcount)	{
		delete[] profilesuffstatcount;
		delete[] profilesuffstatbeta;
		delete[] allocprofilesuffstatcount;
		delete[] allocprofilesuffstatbeta;
		delete[] tmpprofilesuffstatcount;
		delete[] tmpprofilesuffstatbeta;
		delete[] alloctmpprofilesuffstatcount;
		delete[] alloctmpprofilesuffstatbeta;
		profilesuffstatcount = 0;
		profilesuffstatbeta = 0;
		PartitionGTRMixtureProfileProcess::Delete();
		DirichletProfileProcess::Delete();
	}
}

void ExpoConjugatePartitionGTRMixtureProfileProcess::GlobalUpdateModeProfileSuffStat()	{

	if (GetNprocs() > 1)	{
		MESSAGE signal = UPDATE_MPROFILE;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);


		// Achtung! MPI_Reduce includes contribution of master alloctmp arrays in the total sum
		// for that reason, alloctmp arrays should be set to 0,
		// ... but this is done in the constructor (see above)

		MPI_Reduce(alloctmpprofilesuffstatcount,allocprofilesuffstatcount,Ncomponent*GetDim(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Reduce(alloctmpprofilesuffstatbeta,allocprofilesuffstatbeta,Ncomponent*GetDim(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

		MPI_Bcast(allocprofilesuffstatcount,Ncomponent*GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(allocprofilesuffstatbeta,Ncomponent*GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	}

	else	{
		UpdateModeProfileSuffStat();
	}
}

void ExpoConjugatePartitionGTRMixtureProfileProcess::SlaveUpdateModeProfileSuffStat()	{

	UpdateModeProfileSuffStat();
	for (int k=0; k<Ncomponent*GetDim(); k++)	{
		alloctmpprofilesuffstatcount[k] = allocprofilesuffstatcount[k];
		alloctmpprofilesuffstatbeta[k] = allocprofilesuffstatbeta[k];
	}

	MPI_Reduce(alloctmpprofilesuffstatcount,allocprofilesuffstatcount,Ncomponent*GetDim(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(alloctmpprofilesuffstatbeta,allocprofilesuffstatbeta,Ncomponent*GetDim(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

	MPI_Bcast(allocprofilesuffstatcount,Ncomponent*GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(allocprofilesuffstatbeta,Ncomponent*GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);

}

void ExpoConjugatePartitionGTRMixtureProfileProcess::UpdateModeProfileSuffStat()	{

	for (int i=0; i<GetNcomponent(); i++)	{
		for (int k=0; k<GetDim(); k++)	{
			profilesuffstatcount[i][k] = 0;
			profilesuffstatbeta[i][k] = 0;
		}
	}
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i))	{
			const double* count = GetSiteProfileSuffStatCount(i);
			const double* beta = GetSiteProfileSuffStatBeta(i);
			int cat = alloc[i];
			for (int k=0; k<GetDim(); k++)	{
				profilesuffstatcount[cat][k] += count[k];
				profilesuffstatbeta[cat][k] += beta[k];
			}
		}
	}
}

double ExpoConjugatePartitionGTRMixtureProfileProcess::ProfileSuffStatLogProb(int cat)	{
	double total = 0;
	for (int k=0; k<GetDim(); k++)	{
		total += profilesuffstatcount[cat][k] * log(profile[cat][k]) - profilesuffstatbeta[cat][k] * profile[cat][k];
	}
	profilesuffstatlogprob[cat] = total;
	return total;
}

void ExpoConjugatePartitionGTRMixtureProfileProcess::SwapComponents(int cat1, int cat2)	{

	MixtureProfileProcess::SwapComponents(cat1,cat2);
	for (int k=0; k<GetDim(); k++)	{

		double tmp = profilesuffstatcount[cat1][k];
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

double ExpoConjugatePartitionGTRMixtureProfileProcess::LogStatProb(int site, int cat)	{
	const double* count = GetSiteProfileSuffStatCount(site);
	const double* beta = GetSiteProfileSuffStatBeta(site);
	double* pi = profile[cat];
	double total = 0;
	for (int k=0; k<GetDim(); k++)	{
		total += count[k] * log(pi[k]) - beta[k] * pi[k];
	}
	return total;
}


