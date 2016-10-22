
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "GeneralPathSuffStatMatrixMixtureProfileProcess.h"
#include "Random.h"
#include "Parallel.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* GeneralPathSuffStatMatrixMixtureProfileProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void GeneralPathSuffStatMatrixMixtureProfileProcess::Create()	{
	if (! profilepaircount)	{
		MatrixMixtureProfileProcess::Create();
		profilepaircount = new map< pair<int,int>, int>[GetNmodeMax()];
		profilewaitingtime = new map<int,double>[GetNmodeMax()];
		profilerootcount = new map<int,int>[GetNmodeMax()];
	}
}

void GeneralPathSuffStatMatrixMixtureProfileProcess::Delete() {
	if (profilepaircount)	{
		delete[] profilepaircount;
		delete[] profilerootcount;
		delete[] profilewaitingtime;
		profilepaircount = 0;
		profilerootcount = 0;
		profilewaitingtime = 0;
		MatrixMixtureProfileProcess::Delete();
	}
}

void GeneralPathSuffStatMatrixMixtureProfileProcess::GlobalUpdateModeProfileSuffStat()	{

	if (GetNprocs() > 1)	{

		MESSAGE signal = UPDATE_MPROFILE;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

		int* pairvector = new int[GetNcomponent() * GetNstate() * GetNstate()];
		int* tmppairvector = new int[GetNcomponent() * GetNstate() * GetNstate()];
		for (int i=0; i<GetNcomponent() * GetNstate() * GetNstate(); i++)	{
			pairvector[i] = 0;
			tmppairvector[i] = 0;
		}

		int* rootvector = new int[GetNcomponent() * GetNstate()];
		int* tmprootvector = new int[GetNcomponent() * GetNstate()];
		for (int i=0; i<GetNcomponent() * GetNstate(); i++)	{
			rootvector[i] = 0;
			tmprootvector[i] = 0;
		}

		double* waitvector = new double[GetNcomponent() * GetNstate()];
		double* tmpwaitvector = new double[GetNcomponent() * GetNstate()];
		for (int i=0; i<GetNcomponent() * GetNstate(); i++)	{
			waitvector[i] = 0;
			tmpwaitvector[i] = 0;
		}

		for (int cat=0; cat<GetNcomponent(); cat++)	{
			for (map<int,int>::iterator i = profilerootcount[cat].begin(); i!= profilerootcount[cat].end(); i++)	{
				rootvector[cat*GetNstate() + i->first] = i->second;
			}
			for (map<int,double>::iterator i = profilewaitingtime[cat].begin(); i!= profilewaitingtime[cat].end(); i++)	{
				waitvector[cat*GetNstate() + i->first] = i->second;
			}
			for (map<pair<int,int>, int>::iterator i = profilepaircount[cat].begin(); i!= profilepaircount[cat].end(); i++)	{
				pairvector[cat*GetNstate()*GetNstate() + i->first.first*GetNstate() + i->first.second] = i->second;
			}
		}

		MPI_Reduce(tmppairvector,pairvector,GetNcomponent()*GetNstate()*GetNstate(),MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Reduce(tmprootvector,rootvector,GetNcomponent()*GetNstate(),MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Reduce(tmpwaitvector,waitvector,GetNcomponent()*GetNstate(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

		int rm = 0;
		int pm = 0;
		int wm = 0;
		for (int cat=0; cat<GetNcomponent(); cat++)	{
			profilepaircount[cat].clear();
			profilerootcount[cat].clear();
			profilewaitingtime[cat].clear();

			for (int i=0; i<GetNstate(); i++)	{
				if (rootvector[rm])	{
					profilerootcount[cat][i] = rootvector[rm];
				}
				rm++;
			}
			for (int i=0; i<GetNstate(); i++)	{
				if (waitvector[wm])	{
					profilewaitingtime[cat][i] = waitvector[wm];
				}
				wm++;
			}
			for (int i=0; i<GetNstate(); i++)	{
				for (int j=0; j<GetNstate(); j++)	{
					if (pairvector[pm])	{
						profilepaircount[cat][pair<int,int>(i,j)] = pairvector[pm];
					}
					pm++;
				}
			}
		}

		MPI_Bcast(pairvector,GetNcomponent()*GetNstate()*GetNstate(),MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(rootvector,GetNcomponent()*GetNstate(),MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(waitvector,GetNcomponent()*GetNstate(),MPI_DOUBLE,0,MPI_COMM_WORLD);

		delete[] pairvector;
		delete[] tmppairvector;
		delete[] rootvector;
		delete[] tmprootvector;
		delete[] waitvector;
		delete[] tmpwaitvector;
	}
	else	{
		UpdateModeProfileSuffStat();
	}
}


void GeneralPathSuffStatMatrixMixtureProfileProcess::SlaveUpdateModeProfileSuffStat()	{

	int* pairvector = new int[GetNcomponent() * GetNstate() * GetNstate()];
	for (int i=0; i<GetNcomponent() * GetNstate() * GetNstate(); i++)	{
		pairvector[i] = 0;
	}

	int* rootvector = new int[GetNcomponent() * GetNstate()];
	for (int i=0; i<GetNcomponent() * GetNstate(); i++)	{
		rootvector[i] = 0;
	}

	double* waitvector = new double[GetNcomponent() * GetNstate()];
	for (int i=0; i<GetNcomponent() * GetNstate(); i++)	{
		waitvector[i] = 0;
	}

	int* tmppairvector = new int[GetNcomponent() * GetNstate() * GetNstate()];
	for (int i=0; i<GetNcomponent() * GetNstate() * GetNstate(); i++)	{
		tmppairvector[i] = 0;
	}

	int* tmprootvector = new int[GetNcomponent() * GetNstate()];
	for (int i=0; i<GetNcomponent() * GetNstate(); i++)	{
		tmprootvector[i] = 0;
	}

	double* tmpwaitvector = new double[GetNcomponent() * GetNstate()];
	for (int i=0; i<GetNcomponent() * GetNstate(); i++)	{
		tmpwaitvector[i] = 0;
	}

	UpdateModeProfileSuffStat();

	for (int cat=0; cat<GetNcomponent(); cat++)	{
		for (map<int,int>::iterator i = profilerootcount[cat].begin(); i!= profilerootcount[cat].end(); i++)	{
			rootvector[cat*GetNstate() + i->first] = i->second;
			tmprootvector[cat*GetNstate() + i->first] = i->second;
		}
		for (map<int,double>::iterator i = profilewaitingtime[cat].begin(); i!= profilewaitingtime[cat].end(); i++)	{
			waitvector[cat*GetNstate() + i->first] = i->second;
			tmpwaitvector[cat*GetNstate() + i->first] = i->second;
		}
		for (map<pair<int,int>, int>::iterator i = profilepaircount[cat].begin(); i!= profilepaircount[cat].end(); i++)	{
			pairvector[cat*GetNstate()*GetNstate() + i->first.first*GetNstate() + i->first.second] = i->second;
			tmppairvector[cat*GetNstate()*GetNstate() + i->first.first*GetNstate() + i->first.second] = i->second;
		}
	}

	MPI_Reduce(tmppairvector,pairvector,GetNcomponent()*GetNstate()*GetNstate(),MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(tmprootvector,rootvector,GetNcomponent()*GetNstate(),MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(tmpwaitvector,waitvector,GetNcomponent()*GetNstate(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

	delete[] tmppairvector;
	delete[] tmprootvector;
	delete[] tmpwaitvector;

	MPI_Bcast(pairvector,GetNcomponent()*GetNstate()*GetNstate(),MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(rootvector,GetNcomponent()*GetNstate(),MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(waitvector,GetNcomponent()*GetNstate(),MPI_DOUBLE,0,MPI_COMM_WORLD);

	int rm = 0;
	int pm = 0;
	int wm = 0;
	for (int cat=0; cat<GetNcomponent(); cat++)	{
		profilepaircount[cat].clear();
		profilerootcount[cat].clear();
		profilewaitingtime[cat].clear();

		for (int i=0; i<GetNstate(); i++)	{
			if (rootvector[rm])	{
				profilerootcount[cat][i] = rootvector[rm];
			}
			rm++;
		}
		for (int i=0; i<GetNstate(); i++)	{
			if (waitvector[wm])	{
				profilewaitingtime[cat][i] = waitvector[wm];
			}
			wm++;
		}
		for (int i=0; i<GetNstate(); i++)	{
			for (int j=0; j<GetNstate(); j++)	{
				if (pairvector[pm])	{
					profilepaircount[cat][pair<int,int>(i,j)] = pairvector[pm];
				}
				pm++;
			}
		}
	}

	delete[] pairvector;
	delete[] rootvector;
	delete[] waitvector;
}

void GeneralPathSuffStatMatrixMixtureProfileProcess::UpdateModeProfileSuffStat()	{

	for (int i=0; i<GetNcomponent(); i++)	{
		profilepaircount[i].clear();
		profilerootcount[i].clear();
		profilewaitingtime[i].clear();
	}
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i))	{
			map<pair<int,int>, int>& paircount = GetSitePairCount(i);
			map<int,double>& waitingtime = GetSiteWaitingTime(i);
			int rootstate = GetSiteRootState(i);
			int cat = alloc[i];
			if (rootstate != -1)	{
				profilerootcount[cat][rootstate]++;
				for (map<int,double>::iterator i = waitingtime.begin(); i!= waitingtime.end(); i++)	{
					profilewaitingtime[cat][i->first] += i->second;
				}
				for (map<pair<int,int>, int>::iterator i = paircount.begin(); i!= paircount.end(); i++)	{
					profilepaircount[cat][i->first] += i->second;
				}
			}
		}
	}
}

/*
double GeneralPathSuffStatMatrixMixtureProfileProcess::ProfileSuffStatLogProb(int cat)	{

	double tot1 = ProfileSuffStatLogProb0(cat);
	SubMatrix* mat = matrixarray[cat];
	mat->CorruptMatrix();
	mat->UpdateMatrix();
	double tot2 = ProfileSuffStatLogProb0(cat);
	if (fabs(tot2-tot1) > 1e-6)	{
		cerr << "error in GPSS ProfileSuffStatLogProb\n";
		cerr << tot1 << '\t' << tot2 << '\n';
		exit(1);
	}
	return tot1;
}
*/

double GeneralPathSuffStatMatrixMixtureProfileProcess::ProfileSuffStatLogProb(int cat)	{
// double GeneralPathSuffStatMatrixMixtureProfileProcess::ProfileSuffStatLogProb0(int cat)	{

	double total = 0;
	SubMatrix* mat = matrixarray[cat];
	if (! mat)	{
		cerr << "error : null matrix\n";
		cerr << cat << '\t' << Ncomponent << '\n';
		cerr << occupancy[cat] << '\n';
		exit(1);
	}
	const double* stat = matrixarray[cat]->GetStationary();
	for (map<int,int>::iterator i = profilerootcount[cat].begin(); i!= profilerootcount[cat].end(); i++)	{
		total += i->second * log(stat[i->first]);
	}
	for (map<int,double>::iterator i = profilewaitingtime[cat].begin(); i!= profilewaitingtime[cat].end(); i++)	{
		total += i->second * (*mat)(i->first,i->first);
	}
	for (map<pair<int,int>, int>::iterator i = profilepaircount[cat].begin(); i!= profilepaircount[cat].end(); i++)	{
		total += i->second * log((*mat)(i->first.first, i->first.second));
	}
	profilesuffstatlogprob[cat] = total;
	return total;
}

double GeneralPathSuffStatMatrixMixtureProfileProcess::CountProfileSuffStatLogProb(int cat)	{

	double total = 0;
	SubMatrix* mat = matrixarray[cat];
	if (! mat)	{
		cerr << "error : null matrix\n";
		cerr << cat << '\t' << Ncomponent << '\n';
		cerr << occupancy[cat] << '\n';
		exit(1);
	}
	const double* stat = matrixarray[cat]->GetStationary();
	for (map<int,int>::iterator i = profilerootcount[cat].begin(); i!= profilerootcount[cat].end(); i++)	{
		total += i->second * log(stat[i->first]);
	}
	/*
	for (map<int,double>::iterator i = profilewaitingtime[cat].begin(); i!= profilewaitingtime[cat].end(); i++)	{
		total += i->second * (*mat)(i->first,i->first);
	}
	*/
	for (map<pair<int,int>, int>::iterator i = profilepaircount[cat].begin(); i!= profilepaircount[cat].end(); i++)	{
		total += i->second * log((*mat)(i->first.first, i->first.second));
	}
	profilesuffstatlogprob[cat] = total;
	return total;
}

double GeneralPathSuffStatMatrixMixtureProfileProcess::BetaProfileSuffStatLogProb(int cat)	{

	double total = 0;
	SubMatrix* mat = matrixarray[cat];
	if (! mat)	{
		cerr << "error : null matrix\n";
		cerr << cat << '\t' << Ncomponent << '\n';
		cerr << occupancy[cat] << '\n';
		exit(1);
	}
	for (map<int,double>::iterator i = profilewaitingtime[cat].begin(); i!= profilewaitingtime[cat].end(); i++)	{
		total += i->second * (*mat)(i->first,i->first);
	}
	return total;
}

double GeneralPathSuffStatMatrixMixtureProfileProcess::LogStatProb(int site, int cat)	{

	double total = 0;
	SubMatrix* mat = matrixarray[cat];
	const double* stat = matrixarray[cat]->GetStationary();
	int rootstate = GetSiteRootState(site);
	if (rootstate != -1)	{
		total += log(stat[rootstate]);

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


void GeneralPathSuffStatMatrixMixtureProfileProcess::AddSite(int site, int cat)	{
	alloc[site] = cat;
	occupancy[cat] ++;
}

void GeneralPathSuffStatMatrixMixtureProfileProcess::RemoveSite(int site, int cat)	{
	if (cat != -1)	{
		occupancy[cat] --;
	}
}


