
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "GeneralPathSuffStatMatrixPhyloProcess.h"
#include "Parallel.h"
#include <string.h>


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* General Path Suff Stat PhyloProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void GeneralPathSuffStatMatrixPhyloProcess::Unfold()	{

	DeleteMappings();
	ActivateSumOverRateAllocations();
	/*
	if (!sumratealloc)	{
		DrawAllocations(0);
		InactivateSumOverRateAllocations();
	}
	*/
	// this will in fact create only the matrices that did not already exist
	// CreateMatrices();
	// this one is important
	UpdateMatrices();
	// CreateConditionalLikelihoods();
	UpdateConditionalLikelihoods();
	activesuffstat = false;
}

void GeneralPathSuffStatMatrixPhyloProcess::Collapse()	{

	// if (sumratealloc)	{
	DrawAllocations(0);
	InactivateSumOverRateAllocations();
	// }
	SampleNodeStates();
	// DeleteConditionalLikelihoods();
	FillMissingMap();
	SampleSubstitutionMappings(GetRoot());
	// DeleteMatrices();
	activesuffstat = true;
}

void GeneralPathSuffStatMatrixPhyloProcess::GlobalUnfold()	{

	if (GetNprocs() > 1)	{
		GlobalUpdateParameters();

		// CreateMatrices();

		MESSAGE signal = UNFOLD;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

		GlobalUpdateConditionalLikelihoods();
	}
	else	{
		// CreateMatrices();
		Unfold();
		UpdateConditionalLikelihoods();
	}
}


void GeneralPathSuffStatMatrixPhyloProcess::CreateSuffStat()	{

	PhyloProcess::CreateSuffStat();
	if (sitepaircount)	{
		cerr << "error in PhyloProcess::CreateSuffStat\n";
		exit(1);
	}
	siterootstate = new int[GetNsite()];
	sitepaircount = new map<pair<int,int>, int>[GetNsite()];
	sitewaitingtime = new map<int,double>[GetNsite()];
}

void GeneralPathSuffStatMatrixPhyloProcess::DeleteSuffStat()	{

	delete[] siterootstate;
	delete[] sitepaircount;
	delete[] sitewaitingtime;
	siterootstate = 0;
	sitepaircount = 0;
	sitewaitingtime = 0;
	PhyloProcess::DeleteSuffStat();
}

void GeneralPathSuffStatMatrixPhyloProcess::UpdateSiteProfileSuffStat()	{

	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i))	{
			sitepaircount[i].clear();
			sitewaitingtime[i].clear();
		}
	}

	for (int j=0; j<GetNbranch(); j++)	{
		AddSiteProfileSuffStat(siterootstate,sitepaircount,sitewaitingtime,submap[j],blarray[j],missingmap[j]);
		// AddSiteProfileSuffStat(siterootstate,sitepaircount,sitewaitingtime,submap[j],blarray[j],(j == 0));
	}
}

void GeneralPathSuffStatMatrixPhyloProcess::UpdateSiteRateSuffStat()	{

	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i))	{
			siteratesuffstatcount[i] = 0;
			siteratesuffstatbeta[i] = 0;
		}
	}

	for (int j=1; j<GetNbranch(); j++)	{
		AddSiteRateSuffStat(siteratesuffstatcount,siteratesuffstatbeta,submap[j],blarray[j],missingmap[j]);
	}
}

void GeneralPathSuffStatMatrixPhyloProcess::UpdateBranchLengthSuffStat()	{

	branchlengthsuffstatcount[0] = 0;
	branchlengthsuffstatbeta[0] = 0;
	for (int j=1; j<GetNbranch(); j++)	{
		int& count = branchlengthsuffstatcount[j];
		double& beta = branchlengthsuffstatbeta[j];
		count = 0;
		beta = 0;
		AddBranchLengthSuffStat(count,beta,submap[j],missingmap[j]);
	}
}


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* MPI 
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void GeneralPathSuffStatMatrixPhyloProcess::GlobalUpdateSiteProfileSuffStat()	{

	if (GetNprocs() > 1)	{
		for (int i=0; i<GetNsite(); i++)	{
			sitepaircount[i].clear();
			sitewaitingtime[i].clear();
		}

		MPI_Status stat;
		MESSAGE signal = UPDATE_SPROFILE;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

		int inalloc = GetMaxSiteNumber() * (GetNstate()*GetNstate() + 1);
		int dnalloc = GetMaxSiteNumber() * GetNstate();
		int* ivector = new int[inalloc];
		double* dvector = new double[dnalloc];

		for(int i=1; i<GetNprocs(); i++) {
			MPI_Recv(ivector,GetProcSiteNumber(i) * (GetNstate()*GetNstate() + 1),MPI_INT,i,TAG1,MPI_COMM_WORLD,&stat);
			int m = 0;
			for(int j=GetProcSiteMin(i); j<GetProcSiteMax(i); j++) {
				if (ActiveSite(j))	{
					siterootstate[j] = ivector[m];
					m++;
					for(int k=0; k<GetNstate(); k++) {
						for(int l=0; l<GetNstate(); l++) {
							if (ivector[m])	{
								sitepaircount[j][pair<int,int>(k,l)] = ivector[m];
							}
							m++;
						}
					}
				}
			}
		}

		for(int i=1; i<GetNprocs(); i++) {
			MPI_Recv(dvector,GetProcSiteNumber(i)*GetNstate(),MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
			int m = 0;
			for(int j=GetProcSiteMin(i); j<GetProcSiteMax(i); j++) {
				if (ActiveSite(j))	{
					for(int k=0; k<GetNstate(); k++) {
						if (dvector[m])	{
							sitewaitingtime[j][k] = dvector[m];
						}
						m++;
					}
				}
			}
		}

		delete[] ivector;
		delete[] dvector;
	}
	else	{
		UpdateSiteProfileSuffStat();
	}
}

void GeneralPathSuffStatMatrixPhyloProcess::SlaveUpdateSiteProfileSuffStat()	{

	UpdateSiteProfileSuffStat();
	int nn = GetNstate()*GetNstate()+1;
	int iworkload = (GetSiteMax() - GetSiteMin())*(GetNstate()*GetNstate()+1);
	int* ivector = new int[iworkload];
	for (int i=0; i<iworkload; i++)	{
		ivector[i] = 0;
	}

	for(int j=GetSiteMin(); j<GetSiteMax(); j++) {
		if (ActiveSite(j))	{
			ivector[(j-GetSiteMin())*nn] = siterootstate[j];
			map<pair<int,int>, int>& paircount = GetSitePairCount(j);
			for (map<pair<int,int>, int>::iterator i = paircount.begin(); i!= paircount.end(); i++)	{
				ivector[(j-GetSiteMin())*nn + 1 + i->first.first*GetNstate() + i->first.second] = i->second;
			}
		}
	}
	MPI_Send(ivector,iworkload,MPI_INT,0,TAG1,MPI_COMM_WORLD);

	int dworkload = (GetSiteMax() - GetSiteMin())*GetNstate();
	double* dvector = new double[dworkload];
	for (int i=0; i<dworkload; i++)	{
		dvector[i] = 0;
	}
	for(int j=GetSiteMin(); j<GetSiteMax(); ++j) {
		if (ActiveSite(j))	{
			map<int,double>& waitingtime = GetSiteWaitingTime(j);
			for (map<int,double>::iterator i = waitingtime.begin(); i!= waitingtime.end(); i++)	{
				dvector[(j-GetSiteMin())*GetNstate() + i->first] = i->second;
			}
		}
	}
	MPI_Send(dvector,dworkload,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);

	delete[] ivector;
	delete[] dvector;
}


int GeneralPathSuffStatMatrixPhyloProcess::CountMapping(int i)	{
	cerr << "in count mapping\n";
	exit(1);
	int count = 0;
	for(int k=0; k<GetNstate(i); ++k) {
		for(int l=0; l<GetNstate(i); ++l) {
			count+=sitepaircount[i][pair<int,int>(k,l)];
		}
	}
	return count;
}

