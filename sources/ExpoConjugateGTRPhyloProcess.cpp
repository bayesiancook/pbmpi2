
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#include "ExpoConjugateGTRPhyloProcess.h"
#include "Parallel.h"
#include <string.h>

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* ExpoConjugateGTRPhyloProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void ExpoConjugateGTRPhyloProcess::CreateSuffStat()	{

	PhyloProcess::CreateSuffStat();
	if (siteprofilesuffstatcount)	{
		cerr << "error in ExpoConjugateGTRPhyloProcess::CreateSuffStat\n";
		exit(1);
	}
	if (GetMyid() || (GetNprocs() == 1))	{
		allocsiteprofilesuffstatcount = new int[(GetSiteMax() - GetSiteMin())*GetDim()];
		allocsiteprofilesuffstatbeta = new double[(GetSiteMax() - GetSiteMin())*GetDim()];
		siteprofilesuffstatcount = new int*[GetNsite()];
		siteprofilesuffstatbeta = new double*[GetNsite()];

		for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
			siteprofilesuffstatcount[i] = allocsiteprofilesuffstatcount + (i-GetSiteMin())*GetDim();
			siteprofilesuffstatbeta[i] = allocsiteprofilesuffstatbeta + (i-GetSiteMin())*GetDim();
		}
	}
}

void ExpoConjugateGTRPhyloProcess::DeleteSuffStat()	{

	if (siteprofilesuffstatcount)	{
		delete[] siteprofilesuffstatcount;
		delete[] siteprofilesuffstatbeta;
		siteprofilesuffstatcount = 0;
		siteprofilesuffstatbeta = 0;
		delete[] allocsiteprofilesuffstatcount;
		delete[] allocsiteprofilesuffstatbeta;
		allocsiteprofilesuffstatcount = 0;
		allocsiteprofilesuffstatbeta = 0;
	}
	PhyloProcess::DeleteSuffStat();
}

void ExpoConjugateGTRPhyloProcess::UpdateRRSuffStat()	{

	for (int k=0; k<GetNrr(); k++)	{
		rrsuffstatcount[k] = 0;
		rrsuffstatbeta[k] = 0;
	}
	for (int j=1; j<GetNbranch(); j++)	{
		AddRRSuffStat(rrsuffstatcount,rrsuffstatbeta,submap[j],blarray[j],missingmap[j]);
	}
}

void ExpoConjugateGTRPhyloProcess::UpdateSiteRateSuffStat()	{

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

void ExpoConjugateGTRPhyloProcess::UpdateBranchLengthSuffStat()	{

	
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

void ExpoConjugateGTRPhyloProcess::UpdateSiteProfileSuffStat()	{

	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i))	{
			for (int k=0; k<GetDim(); k++)	{
				siteprofilesuffstatcount[i][k] = 0;
				siteprofilesuffstatbeta[i][k] = 0;
			}
		}
	}
	for (int j=0; j<GetNbranch(); j++)	{
		AddSiteProfileSuffStat(siteprofilesuffstatcount,siteprofilesuffstatbeta,submap[j],blarray[j],missingmap[j]);
	}
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i))	{
			if (GetData()->AllMissingColumn(i))	{
				for (int k=0; k<GetDim(); k++)	{
					if (siteprofilesuffstatcount[i][k] != 0)	{
						cerr << "error in ExpoConjugateGTRPhyloProcess::UpdateSiteProfileSuffStat: all missing column has counts\n";
						cerr << "site : " << i << '\n';
						GetData()->PrintColumn(cerr,i);
						exit(1);
					}
					if (siteprofilesuffstatbeta[i][k] != 0)	{
						cerr << "error in ExpoConjugateGTRPhyloProcess::UpdateSiteProfileSuffStat: all missing column has positive beta\n";
						cerr << "site : " << i << '\n';
						GetData()->PrintColumn(cerr,i);
						exit(1);
					}
				}
			}
		}
	}
}

void ExpoConjugateGTRPhyloProcess::GlobalUpdateSiteProfileSuffStat()	{

	if (GetNprocs() > 1)	{
		MPI_Status stat;
		MESSAGE signal = UPDATE_SPROFILE;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	}
	else	{
		UpdateSiteProfileSuffStat();
	}
}

void ExpoConjugateGTRPhyloProcess::SlaveUpdateSiteProfileSuffStat()	{

	UpdateSiteProfileSuffStat();
}

int ExpoConjugateGTRPhyloProcess::GlobalCountMapping()	{

	GlobalUpdateSiteProfileSuffStat();
	return PhyloProcess::GlobalCountMapping();
}

int ExpoConjugateGTRPhyloProcess::CountMapping()	{

	int total = 0;	
	for(int i=GetSiteMin(); i<GetSiteMax(); i++){
		if (ActiveSite(i))	{
			total += CountMapping(i);
		}
	}
	return total;
}

int ExpoConjugateGTRPhyloProcess::CountMapping(int i)	{

	const int* tmp = GetSiteProfileSuffStatCount(i);
	int total = 0;
	for (int k=0; k<GetNstate(i); k++)	{
		total += tmp[k];
	}
	total--;
	return total;
}


void ExpoConjugateGTRPhyloProcess::ReadSiteProfiles(string name, int burnin, int every, int until)	{

	double** sitestat = new double*[GetNsite()];
	for (int i=0; i<GetNsite(); i++)	{
		sitestat[i] = new double[GetDim()];
		for (int k=0; k<GetDim(); k++)	{
			sitestat[i][k] = 0;
		}
	}
	ifstream is((name + ".chain").c_str());
	if (!is)	{
		cerr << "error: no .chain file found\n";
		exit(1);
	}

	cerr << "burnin : " << burnin << "\n";
	cerr << "until : " << until << '\n';
	int i=0;
	while ((i < until) && (i < burnin))	{
		FromStream(is);
		i++;
	}
	int samplesize = 0;

	while (i < until)	{
		cerr << ".";
		cerr.flush();
		samplesize++;
		FromStream(is);
		i++;

		for (int i=0; i<GetNsite(); i++)	{
			double* p = GetProfile(i);
			for (int k=0; k<GetDim(); k++)	{
				sitestat[i][k] += p[k];
			}
		}
		int nrep = 1;
		while ((i<until) && (nrep < every))	{
			FromStream(is);
			i++;
			nrep++;
		}
	}
	cerr << '\n';
	
	ofstream os((name + ".siteprofiles").c_str());
	for (int k=0; k<GetDim(); k++)	{
		os << GetStateSpace()->GetState(k) << ' ';
	}
	os << '\n';
	os << '\n';
	for (int i=0; i<GetNsite(); i++)	{
		os << i+1;
		for (int k=0; k<GetDim(); k++)	{
			sitestat[i][k] /= samplesize;
			os << '\t' << sitestat[i][k];
		}
		os << '\n';
	}
	cerr << "mean site-specific profiles in " << name << ".siteprofiles\n";
	cerr << '\n';
}
