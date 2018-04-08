
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "ExpoConjugatePartitionGTRSubstitutionProcess.h"


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* ExpoConjugate GTR Substitution Processes
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
//	* Gather sufficient statistics for resampling continuous parameters of the model
//	(branch lengths, relative rates, and site-specific rates and profiles)
//	just scanning through substitution mappings (singly linked lists of substitutions and waiting times)
//	counting them and summing over their exponential rates
//	(CPU level 1)
//-------------------------------------------------------------------------


void ExpoConjugatePartitionGTRSubstitutionProcess::AddRRSuffStat(double* allocrrsuffstatcount, double* allocrrsuffstatbeta, BranchSitePath** patharray, double branchlength, int* nonmissing)	{
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i) && (nonmissing[i] == 1))	{
		// if (ActiveSite(i))	{
			patharray[i]->AddRRSuffStat(allocrrsuffstatcount+partalloc[i]*GetNrr(),allocrrsuffstatbeta+partalloc[i]*GetNrr(),GetRate(i)*branchlength,GetMatrix(i), GetSiteRR(i),GetNstate(i));
			// patharray[i]->AddRRSuffStat(allocrrsuffstatcount+partalloc[i]*GetNrr(),allocrrsuffstatbeta+partalloc[i]*GetNrr(),GetRate(i)*branchlength,GetProfile(i),GetNstate(i));
		}
	}
}

void ExpoConjugatePartitionGTRSubstitutionProcess::AddSiteRateSuffStat(double* siteratesuffstatcount, double* siteratesuffstatbeta, BranchSitePath** patharray, double branchlength, int* nonmissing)	{
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i) && (nonmissing[i] == 1))	{
		// if (ActiveSite(i))	{
			patharray[i]->AddRateSuffStat(siteratesuffstatcount[i],siteratesuffstatbeta[i],branchlength,GetSiteRR(i),GetProfile(i),GetNstate(i));
		}
	}
}

void ExpoConjugatePartitionGTRSubstitutionProcess::AddBranchLengthSuffStat(double& count, double& beta, BranchSitePath** patharray, int* nonmissing)	{
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i) && (nonmissing[i] == 1))	{
		// if (ActiveSite(i))	{
			patharray[i]->AddRateSuffStat(count,beta,GetRate(i),GetSiteRR(i),GetProfile(i),GetNstate(i));
		}
	}
}

void ExpoConjugatePartitionGTRSubstitutionProcess::AddSiteProfileSuffStat(double** siteprofilesuffstatcount, double** siteprofilesuffstatbeta, BranchSitePath** patharray, double branchlength, int* nonmissing)	{

	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i))	{
			if (nonmissing[i] == 1)	{
				patharray[i]->AddProfileSuffStat(siteprofilesuffstatcount[i],siteprofilesuffstatbeta[i],GetRate(i)*branchlength,GetSiteRR(i),GetNstate(i));
			}
			else if (nonmissing[i] == 2)	{
				siteprofilesuffstatcount[i][patharray[i]->GetFinalState()]++;
			}
		}
	}
}

