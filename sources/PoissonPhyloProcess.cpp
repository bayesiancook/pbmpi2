
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "Parallel.h"
#include <string>
#include "PoissonPhyloProcess.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* Poisson PhyloProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void PoissonPhyloProcess::Create()	{
	if (! zipdata)	{
		truedata = data;
		zipdata = new ZippedSequenceAlignment(truedata);
		data = zipdata;
		PhyloProcess::Create();
		truedata->GetEmpiricalFreq(empfreq);
		CreateZip();
	}
}

void PoissonPhyloProcess::Delete()	{
	if (zipdata)	{
		PhyloProcess::Delete();
		delete zipdata;
		zipdata = 0;
	}
}

void PoissonPhyloProcess::CreateSuffStat()	{

	PhyloProcess::CreateSuffStat();
	if (siteprofilesuffstatcount)	{
		cerr << "error in PoissonPhyloProcess::CreateSuffStat\n";
		exit(1);
	}
	if (GetMyid() || (GetNprocs() == 1))	{
		allocsiteprofilesuffstatcount = new int[(GetSiteMax() - GetSiteMin())*GetDim()];
		siteprofilesuffstatcount = new int*[GetNsite()];
		for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
			siteprofilesuffstatcount[i] = allocsiteprofilesuffstatcount + (i-GetSiteMin())*GetDim();
		}
	}
}

void PoissonPhyloProcess::DeleteSuffStat()	{

	if (siteprofilesuffstatcount)	{
		delete[] siteprofilesuffstatcount;
		siteprofilesuffstatcount = 0;
		delete[] allocsiteprofilesuffstatcount;
		allocsiteprofilesuffstatcount = 0;
	}
	PhyloProcess::DeleteSuffStat();
}

void PoissonPhyloProcess::Collapse()	{

	// if (sumratealloc)	{
	DrawAllocations(0);
	InactivateSumOverRateAllocations();
	// }
	if ((topobf == 1) || (topobf == 3))	{
		SetMinMax(bffrac,1);
		SampleNodeStates();
		FillMissingMap();
		SampleSubstitutionMappings(GetRoot());
		PoissonUpdateSiteProfileSuffStat();
		SlaveSwapTree();
		SetMinMax(0,bffrac);
		SampleNodeStates();
		FillMissingMap();
		SampleSubstitutionMappings(GetRoot());
		PoissonUpdateSiteProfileSuffStat();
		SlaveSwapTree();
		SetMinMax(0,1);
	}
	else	{
		SampleNodeStates();
		FillMissingMap();
		SampleSubstitutionMappings(GetRoot());
		PoissonUpdateSiteProfileSuffStat();
	}
	activesuffstat = true;
}

void PoissonPhyloProcess::UpdateBranchLengthSuffStat()	{

	branchlengthsuffstatbeta[0] = 0;
	branchlengthsuffstatcount[0] = 0;
	for (int j=1; j<GetNbranch(); j++)	{
		branchlengthsuffstatbeta[j] = 0;
		branchlengthsuffstatcount[j] = 0;
		AddBranchLengthSuffStat(branchlengthsuffstatcount[j],branchlengthsuffstatbeta[j],submap[j],missingmap[j]);
	}
}

void PoissonPhyloProcess::UpdateSiteRateSuffStat()	{
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i))	{
			siteratesuffstatcount[i] = 0;
			siteratesuffstatbeta[i] = 0;
		}
	}
	for (int j=1; j<GetNbranch(); j++)	{
		AddSiteRateSuffStat(siteratesuffstatcount,siteratesuffstatbeta,blarray[j],submap[j],missingmap[j]);
	}
}

int PoissonPhyloProcess::RecursiveUpdateSiteProfileSuffStat(const Link* from, int site)	{

	int state = -1;
	if (from->isLeaf())	{
		state = GetData(from)[site];
	}
	else	{
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			int tmp = RecursiveUpdateSiteProfileSuffStat(link->Out(),site);
			if (tmp != -1)	{
				if ((state != -1) && (state != tmp))	{
					cerr << "error in PoissonPhyloProcess::RecursiveUpdateSiteProfileSuffStat: state should be identical\n";
					cerr << "site : " << site << '\n';
					cerr << state << '\t' << tmp << '\t' << GetZipSize(site) << '\t' << GetOrbitSize(site) << '\n';
					cerr << "rate : " << GetRate(site) << '\n';
					truedata->PrintColumn(cerr,site);
					SampleTrueNodeStates(GetRoot());
					RecursiveUnzipBranchSitePath(GetRoot());
					WriteTreeMapping(cerr, GetRoot(), site);
					exit(1);
				}
				state = tmp;
			}
		}
	}
	if (state != -1)	{
		BranchSitePath* path = submap[GetBranchIndex(from->GetBranch())][site];
		if (from->isRoot() || path->GetNsub())	{
			if ((GetZipSize(site) != GetOrbitSize(site)) && (state == GetOrbitSize(site)))	{
				cerr << "error in PoissonPhyloProcess::RecursiveUpdateSiteProfileSuffStat: state should be observed at tips\n";
				exit(1);
			}
			int truestate = GetStateFromZip(site,state);
			siteprofilesuffstatcount[site][truestate]++;
			state = -1;
		}
	}
	return state;
}


void PoissonPhyloProcess::UpdateSiteProfileSuffStat()	{
}

void PoissonPhyloProcess::PoissonUpdateSiteProfileSuffStat()	{

	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i))	{
			for (int k=0; k<GetDim(); k++)	{
				siteprofilesuffstatcount[i][k] = 0;
			}
			RecursiveUpdateSiteProfileSuffStat(GetRoot(),i);
		}
	}
}

void PoissonPhyloProcess::GlobalUpdateSiteProfileSuffStat()	{

	if (GetNprocs() > 1)	{
		MPI_Status stat;
		MESSAGE signal = UPDATE_SPROFILE;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	}
	else	{
		UpdateSiteProfileSuffStat();
	}
}

void PoissonPhyloProcess::SlaveUpdateSiteProfileSuffStat()	{

	UpdateSiteProfileSuffStat();
}


void PoissonPhyloProcess::GlobalSetTestData()	{

	testnsite = testdata->GetNsite();
	int* tmp = new int[testnsite * GetNtaxa()];
	testdata->GetDataVector(tmp);

	if (GetNprocs() > 1)	{
		MESSAGE signal = SETTESTDATA;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(&testnsite,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(tmp,testnsite*GetNtaxa(),MPI_INT,0,MPI_COMM_WORLD);
	}
	else	{
		data->SetTestData(testnsite,0,0,testnsite,tmp);
	}
	delete[] tmp;

	// GlobalCollapse();
	DeleteZip();
	delete zipdata;

	zipdata = new ZippedSequenceAlignment(truedata);
	CreateZip();
}

void PoissonPhyloProcess::SlaveSetTestData()	{

	MPI_Bcast(&testnsite,1,MPI_INT,0,MPI_COMM_WORLD);
	int* tmp = new int[testnsite * GetNtaxa()];
	MPI_Bcast(tmp,testnsite*GetNtaxa(),MPI_INT,0,MPI_COMM_WORLD);
	
	SetTestSiteMinAndMax();
	truedata->SetTestData(testnsite,GetSiteMin(),testsitemin,testsitemax,tmp);

	delete[] tmp;

	Collapse();
	DeleteZip();
	delete zipdata;

	zipdata = new ZippedSequenceAlignment(truedata);
	CreateZip();
	Unfold();
}

void PoissonPhyloProcess::RecursiveUnzipBranchSitePath(const Link* from){
	for (const Link* link=from->Next(); link!=from; link=link->Next()){
		RecursiveUnzipBranchSitePath(link->Out());
		UnzipBranchSitePath(submap[GetBranchIndex(link->GetBranch())], GetStates(link->GetNode()), GetStates(link->Out()->GetNode()));
	}
}


void PoissonPhyloProcess::SlaveWriteMappings(){
	SampleTrueNodeStates(GetRoot());
	RecursiveUnzipBranchSitePath(GetRoot());
	PhyloProcess::SlaveWriteMappings();
}


void PoissonPhyloProcess::SampleTrueNodeStates(const Link* from)	{
	
	if (from->isRoot())	{
		ChooseTrueStates(submap[0],GetStates(from->Out()->GetNode()),GetStates(from->GetNode()),from->isRoot());
	}
	else{
		ChooseTrueStates(submap[GetBranchIndex(from->GetBranch())],GetStates(from->Out()->GetNode()),GetStates(from->GetNode()),from->isRoot());
	}

	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		SampleTrueNodeStates(link->Out());
	}
}

