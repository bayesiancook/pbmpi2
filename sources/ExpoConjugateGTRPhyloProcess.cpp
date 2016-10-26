
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


double ExpoConjugateGTRPhyloProcess::LengthRelRateMove(double tuning, int nrep)	{

	double naccept = 0;
	for (int rep=0; rep<nrep; rep++)	{
		double deltalogratio = - LogRRPrior() - LogLengthPrior();
		double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
		double e = exp(m);
		int nbranch = MoveAllBranches(e);
		for (int i=0; i<GetNrr(); i++)	{
			rr[i] /= e;
		}
		deltalogratio += LogRRPrior() + LogLengthPrior();
		deltalogratio += (nbranch-GetNrr()) * m;

		int accepted = (log(rnd::GetRandom().Uniform()) < deltalogratio);

		if (accepted)	{
			naccept++;
		}
		else	{
			MoveAllBranches(1.0/e);
			for (int i=0; i<GetNrr(); i++)	{
				rr[i] *= e;
			}
		}	
	}
	return naccept / nrep;
}
