
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "GeneralPathSuffStatMultipleMatrixMixtureProfileProcess.h"
#include "Random.h"
#include "Parallel.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* GeneralPathSuffStatMultipleMatrixMixtureProfileProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void GeneralPathSuffStatMultipleMatrixMixtureProfileProcess::Create()	{
	if (! profilepaircount)	{
		MultipleMatrixMixtureProfileProcess::Create();
		profilepaircount = new map< pair<int,int>, int>*[GetNmodeMax()];
		profilewaitingtime = new map<int,double>*[GetNmodeMax()];
		profilerootcount = new map<int,int>*[GetNmodeMax()];
		for (int k=0; k<GetNmodeMax(); k++)	{
			profilepaircount[k] = new map< pair<int,int>, int>[GetNSubAlloc()];
			profilewaitingtime[k] = new map<int,double>[GetNSubAlloc()];
			profilerootcount[k] = new map<int,int>[GetNSubAlloc()];
		}
	}
}

void GeneralPathSuffStatMultipleMatrixMixtureProfileProcess::Delete() {
	if (profilepaircount)	{
		for (int k=0; k<GetNmodeMax(); k++)	{
			delete[] profilepaircount[k];
			delete[] profilerootcount[k];
			delete[] profilewaitingtime[k];
		}
		delete[] profilepaircount;
		delete[] profilerootcount;
		delete[] profilewaitingtime;
		profilepaircount = 0;
		profilerootcount = 0;
		profilewaitingtime = 0;
		MultipleMatrixMixtureProfileProcess::Delete();
	}
}

void GeneralPathSuffStatMultipleMatrixMixtureProfileProcess::GlobalUpdateModeProfileSuffStat()	{

	UpdateModeProfileSuffStat();

	if (GetNprocs() > 1)	{
	MESSAGE signal = UPDATE_MPROFILE;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	int* pairvector = new int[GetNcomponent() * GetNSubAlloc() * GetNstate() * GetNstate()];
	for (int i=0; i<GetNcomponent() * GetNSubAlloc() * GetNstate() * GetNstate(); i++)	{
		pairvector[i] = 0;
	}

	int* rootvector = new int[GetNcomponent() * GetNSubAlloc() * GetNstate()];
	for (int i=0; i<GetNcomponent() * GetNSubAlloc() * GetNstate(); i++)	{
		rootvector[i] = 0;
	}

	double* waitvector = new double[GetNcomponent() * GetNSubAlloc() * GetNstate()];
	for (int i=0; i<GetNcomponent() * GetNSubAlloc() * GetNstate(); i++)	{
		waitvector[i] = 0;
	}

	for (int cat=0; cat<GetNcomponent(); cat++)	{
		for (int l=0; l<GetNSubAlloc(); l++)	{
			for (map<int,int>::iterator i = profilerootcount[cat][l].begin(); i!= profilerootcount[cat][l].end(); i++)	{
				rootvector[cat*GetNstate()*GetNSubAlloc() + l*GetNstate() + i->first] = i->second;
			}
			for (map<int,double>::iterator i = profilewaitingtime[cat][l].begin(); i!= profilewaitingtime[cat][l].end(); i++)	{
				waitvector[cat*GetNstate()*GetNSubAlloc() + l*GetNstate() + i->first] = i->second;
			}
			for (map<pair<int,int>, int>::iterator i = profilepaircount[cat][l].begin(); i!= profilepaircount[cat][l].end(); i++)	{
				pairvector[cat*GetNstate()*GetNstate()*GetNSubAlloc() + l*GetNstate() + i->first.first*GetNstate() + i->first.second] = i->second;
			}
		}
	}

	MPI_Bcast(pairvector,GetNcomponent()*GetNSubAlloc()*GetNstate()*GetNstate(),MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(rootvector,GetNcomponent()*GetNSubAlloc()*GetNstate(),MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(waitvector,GetNcomponent()*GetNSubAlloc()*GetNstate(),MPI_DOUBLE,0,MPI_COMM_WORLD);

	delete[] pairvector;
	delete[] rootvector;
	delete[] waitvector;
	}
}

void GeneralPathSuffStatMultipleMatrixMixtureProfileProcess::SlaveUpdateModeProfileSuffStat()	{

	int* pairvector = new int[GetNcomponent() * GetNSubAlloc() * GetNstate() * GetNstate()];
	for (int i=0; i<GetNcomponent() * GetNSubAlloc() * GetNstate() * GetNstate(); i++)	{
		pairvector[i] = 0;
	}

	int* rootvector = new int[GetNcomponent() * GetNSubAlloc() * GetNstate()];
	for (int i=0; i<GetNcomponent() * GetNSubAlloc() * GetNstate(); i++)	{
		rootvector[i] = 0;
	}

	double* waitvector = new double[GetNcomponent() * GetNSubAlloc() * GetNstate()];
	for (int i=0; i<GetNcomponent() * GetNSubAlloc() * GetNstate(); i++)	{
		waitvector[i] = 0;
	}

	MPI_Bcast(pairvector,GetNcomponent()*GetNSubAlloc()*GetNstate()*GetNstate(),MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(rootvector,GetNcomponent()*GetNSubAlloc()*GetNstate(),MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(waitvector,GetNcomponent()*GetNSubAlloc()*GetNstate(),MPI_DOUBLE,0,MPI_COMM_WORLD);

	int rm = 0;
	int pm = 0;
	int wm = 0;
	for (int cat=0; cat<GetNcomponent(); cat++)	{
		for (int l=0; l<GetNSubAlloc(); l++)	{
		profilepaircount[cat][l].clear();
		profilerootcount[cat][l].clear();
		profilewaitingtime[cat][l].clear();

		for (int i=0; i<GetNstate(); i++)	{
			if (rootvector[rm])	{
				profilerootcount[cat][l][i] = rootvector[rm];
			}
			rm++;
		}
		for (int i=0; i<GetNstate(); i++)	{
			if (waitvector[wm])	{
				profilewaitingtime[cat][l][i] = waitvector[wm];
			}
			wm++;
		}
		for (int i=0; i<GetNstate(); i++)	{
			for (int j=0; j<GetNstate(); j++)	{
				if (pairvector[pm])	{
					profilepaircount[cat][l][pair<int,int>(i,j)] = pairvector[pm];
				}
				pm++;
			}
		}
		}
	}

	delete[] pairvector;
	delete[] rootvector;
	delete[] waitvector;
}

void GeneralPathSuffStatMultipleMatrixMixtureProfileProcess::UpdateModeProfileSuffStat()	{

	if (GetMyid())	{
		cerr << "error: slave in GPSSMultipleMatrixMixtureProfileProcess::UpdateModeProfileSuffStat\n";
		exit(1);
	}

	for (int i=0; i<GetNcomponent(); i++)	{
		for (int l=0; l<GetNSubAlloc(); l++)	{
			profilepaircount[i][l].clear();
			profilerootcount[i][l].clear();
			profilewaitingtime[i][l].clear();
		}
	}
	for (int i=0; i<GetNsite(); i++)	{

		if (ActiveSite(i))	{
		
			map<pair<int,int>, int>& paircount = GetSitePairCount(i);
			map<int,double>& waitingtime = GetSiteWaitingTime(i);
			int rootstate = GetSiteRootState(i);
			int cat = alloc[i];
			int l = GetSubAlloc(i);
			if (rootstate != -1)	{
				profilerootcount[cat][l][rootstate]++;
				for (map<int,double>::iterator i = waitingtime.begin(); i!= waitingtime.end(); i++)	{
					profilewaitingtime[cat][l][i->first] += i->second;
				}
				for (map<pair<int,int>, int>::iterator i = paircount.begin(); i!= paircount.end(); i++)	{
					profilepaircount[cat][l][i->first] += i->second;
				}
			}
		}
	}
}

double GeneralPathSuffStatMultipleMatrixMixtureProfileProcess::ProfileSuffStatLogProb(int cat)	{

	double total = 0;
	for (int l=0; l<GetNSubAlloc(); l++)	{
		total += ProfileSuffStatLogProb(cat,l);
	}
	profilesuffstatlogprob[cat] = total;
	return total;
}

double GeneralPathSuffStatMultipleMatrixMixtureProfileProcess::ProfileSuffStatLogProb(int cat, int l)	{

	double total = 0;
	SubMatrix* mat = matrixarray[cat][l];
	if (! mat)	{
		cerr << "error : null matrix\n";
		cerr << cat << '\t' << Ncomponent << '\n';
		cerr << occupancy[cat] << '\n';
		exit(1);
	}
	const double* stat = mat->GetStationary();
	for (map<int,int>::iterator i = profilerootcount[cat][l].begin(); i!= profilerootcount[cat][l].end(); i++)	{
		total += i->second * log(stat[i->first]);
	}
	for (map<int,double>::iterator i = profilewaitingtime[cat][l].begin(); i!= profilewaitingtime[cat][l].end(); i++)	{
		total += i->second * (*mat)(i->first,i->first);
	}
	for (map<pair<int,int>, int>::iterator i = profilepaircount[cat][l].begin(); i!= profilepaircount[cat][l].end(); i++)	{
		total += i->second * log((*mat)(i->first.first, i->first.second));
	}
	// profilesuffstatlogprob[cat] = total;
	return total;
}

double GeneralPathSuffStatMultipleMatrixMixtureProfileProcess::LogStatProb(int site, int cat)	{

	double total = 0;
	int l = GetSubAlloc(site);
	SubMatrix* mat = matrixarray[cat][l];
	const double* stat = mat->GetStationary();
	int rootstate = GetSiteRootState(site);
	if (rootstate != -1)	{
		total += log(stat[GetSiteRootState(site)]);

		map<int,double>& waitingtime = GetSiteWaitingTime(site);
		for (map<int,double>::iterator i = waitingtime.begin(); i!= waitingtime.end(); i++)	{
			total += i->second * (*mat)(i->first,i->first);
		}

		map<pair<int,int>, int>& paircount = GetSitePairCount(site);
		for (map<pair<int,int>, int>::iterator i = paircount.begin(); i!= paircount.end(); i++)	{
			total += i->second * log((*mat)(i->first.first, i->first.second));
		}
	}
	return total;
}


void GeneralPathSuffStatMultipleMatrixMixtureProfileProcess::AddSite(int site, int cat)	{
	alloc[site] = cat;
	occupancy[cat] ++;
}

void GeneralPathSuffStatMultipleMatrixMixtureProfileProcess::RemoveSite(int site, int cat)	{
	if (cat != -1)	{
		occupancy[cat] --;
	}
}


