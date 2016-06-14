
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
		tmpprofilesuffstatcount = new int*[GetNmodeMax()];
		tmpprofilesuffstatbeta = new double*[GetNmodeMax()];
		alloctmpprofilesuffstatcount = new int[GetNmodeMax() * GetDim()];
		alloctmpprofilesuffstatbeta = new double[GetNmodeMax() * GetDim()];
		for (int i=0; i<GetNmodeMax(); i++)	{
			profilesuffstatcount[i] = allocprofilesuffstatcount + i*GetDim();
			profilesuffstatbeta[i] = allocprofilesuffstatbeta + i*GetDim();
			tmpprofilesuffstatcount[i] = alloctmpprofilesuffstatcount + i*GetDim();
			tmpprofilesuffstatbeta[i] = alloctmpprofilesuffstatbeta + i*GetDim();
		}
	}
}

void ExpoConjugateGTRMixtureProfileProcess::Delete() {
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
		GTRMixtureProfileProcess::Delete();
		DirichletProfileProcess::Delete();
	}
}

void ExpoConjugateGTRMixtureProfileProcess::GlobalUpdateModeProfileSuffStat()	{

	if (GetNprocs() > 1)	{
		MESSAGE signal = UPDATE_MPROFILE;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Reduce(alloctmpprofilesuffstatcount,allocprofilesuffstatcount,Ncomponent*GetDim(),MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Reduce(alloctmpprofilesuffstatbeta,allocprofilesuffstatbeta,Ncomponent*GetDim(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Bcast(allocprofilesuffstatcount,Ncomponent*GetDim(),MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(allocprofilesuffstatbeta,Ncomponent*GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	}

	else	{
		UpdateModeProfileSuffStat();
	}
}

void ExpoConjugateGTRMixtureProfileProcess::SlaveUpdateModeProfileSuffStat()	{

	UpdateModeProfileSuffStat();
	for (int k=0; k<Ncomponent*GetDim(); k++)	{
		alloctmpprofilesuffstatcount[k] = allocprofilesuffstatcount[k];
		alloctmpprofilesuffstatbeta[k] = allocprofilesuffstatbeta[k];
	}
	MPI_Reduce(alloctmpprofilesuffstatcount,allocprofilesuffstatcount,Ncomponent*GetDim(),MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(alloctmpprofilesuffstatbeta,allocprofilesuffstatbeta,Ncomponent*GetDim(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(allocprofilesuffstatcount,Ncomponent*GetDim(),MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(allocprofilesuffstatbeta,Ncomponent*GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);

}

void ExpoConjugateGTRMixtureProfileProcess::UpdateModeProfileSuffStat()	{

	for (int i=0; i<GetNcomponent(); i++)	{
		for (int k=0; k<GetDim(); k++)	{
			profilesuffstatcount[i][k] = 0;
			profilesuffstatbeta[i][k] = 0;
		}
	}
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
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


