
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

/*
void GeneralPathSuffStatMatrixPhyloProcess::Unfold()	{

	DeleteMappings();
	ActivateSumOverRateAllocations();
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
*/


void GeneralPathSuffStatMatrixPhyloProcess::CreateSuffStat()	{

	PhyloProcess::CreateSuffStat();
	if (sitepaircount)	{
		cerr << "error in PhyloProcess::CreateSuffStat\n";
		exit(1);
	}
	if (GetMyid() || (GetNprocs() == 1))	{
		siterootstate = new int[GetNsite()];
		sitepaircount = new map<pair<int,int>, int>[GetNsite()];
		sitewaitingtime = new map<int,double>[GetNsite()];
	}
}

void GeneralPathSuffStatMatrixPhyloProcess::DeleteSuffStat()	{

	if (sitepaircount)	{
		delete[] siterootstate;
		delete[] sitepaircount;
		delete[] sitewaitingtime;
	}
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
			siterootstate[i] = -1;
		}
	}

	for (int j=0; j<GetNbranch(); j++)	{
		AddSiteProfileSuffStat(siterootstate,sitepaircount,sitewaitingtime,submap[j],blarray[j],missingmap[j]);
	}

	// check
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i))	{
			if (siterootstate[i] == -1)	{
				if (! GetData()->AllMissingColumn(i))	{
					cerr << "error in GPSSMatrixPhyloProcess::UpdateSiteProfileSuffStat: site root state is -1\n";
					cerr << "site : " << i << '\n';
					GetData()->PrintColumn(cerr,i);
					exit(1);
				}
				else	{
					if (sitepaircount[i].size() || sitewaitingtime[i].size())	{
						cerr << "error in GPSSMatrixPhyloProcess:: all missing column has non empty suff stats\n";
						cerr << "site : " << i << '\n';
						GetData()->PrintColumn(cerr,i);
						exit(1);
					}
				}
			}
		}
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
		double& count = branchlengthsuffstatcount[j];
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
		/*
		for (int i=0; i<GetNsite(); i++)	{
			sitepaircount[i].clear();
			sitewaitingtime[i].clear();
		}
		*/

		MPI_Status stat;
		MESSAGE signal = UPDATE_SPROFILE;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	}
	else	{
		UpdateSiteProfileSuffStat();
	}
}

void GeneralPathSuffStatMatrixPhyloProcess::SlaveUpdateSiteProfileSuffStat()	{

	UpdateSiteProfileSuffStat();
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

