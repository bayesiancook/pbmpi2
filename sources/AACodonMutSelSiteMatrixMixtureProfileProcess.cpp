
#include "Parallel.h"
#include "AACodonMutSelSiteMatrixMixtureProfileProcess.h"

void AACodonMutSelSiteMatrixMixtureProfileProcess::Create() {
	if (! omega0)	{
		omega0 = new double;
		*omega0 = 1;
		AACodonMutSelProfileProcess::Create();
		GeneralPathSuffStatSiteMatrixMixtureProfileProcess::Create();
		profilenonsynwaitingtime = new map<int,double>[GetNmodeMax()];
	}
}

void AACodonMutSelSiteMatrixMixtureProfileProcess::Delete() {
	if (omega0)	{
		delete omega0;
		omega0 = 0;
		delete[] profilenonsynwaitingtime;
		GeneralPathSuffStatSiteMatrixMixtureProfileProcess::Delete();
		AACodonMutSelProfileProcess::Delete();
	}
}

void AACodonMutSelSiteMatrixMixtureProfileProcess::GlobalUpdateModeProfileSuffStat()	{

	// UpdateMatrices();
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

		double* nonsynwaitvector = new double[GetNcomponent() * GetNstate()];
		double* tmpnonsynwaitvector = new double[GetNcomponent() * GetNstate()];
		for (int i=0; i<GetNcomponent() * GetNstate(); i++)	{
			nonsynwaitvector[i] = 0;
			tmpnonsynwaitvector[i] = 0;
		}

		for (int cat=0; cat<GetNcomponent(); cat++)	{
			for (map<int,int>::iterator i = profilerootcount[cat].begin(); i!= profilerootcount[cat].end(); i++)	{
				rootvector[cat*GetNstate() + i->first] = i->second;
			}
			for (map<int,double>::iterator i = profilewaitingtime[cat].begin(); i!= profilewaitingtime[cat].end(); i++)	{
				waitvector[cat*GetNstate() + i->first] = i->second;
			}
			for (map<int,double>::iterator i = profilenonsynwaitingtime[cat].begin(); i!= profilenonsynwaitingtime[cat].end(); i++)	{
				nonsynwaitvector[cat*GetNstate() + i->first] = i->second;
			}
			for (map<pair<int,int>, int>::iterator i = profilepaircount[cat].begin(); i!= profilepaircount[cat].end(); i++)	{
				pairvector[cat*GetNstate()*GetNstate() + i->first.first*GetNstate() + i->first.second] = i->second;
			}
		}

		MPI_Reduce(tmppairvector,pairvector,GetNcomponent()*GetNstate()*GetNstate(),MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
		// MPI_Barrier(MPI_COMM_WORLD);
		MPI_Reduce(tmprootvector,rootvector,GetNcomponent()*GetNstate(),MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
		// MPI_Barrier(MPI_COMM_WORLD);
		MPI_Reduce(tmpwaitvector,waitvector,GetNcomponent()*GetNstate(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		// MPI_Barrier(MPI_COMM_WORLD);
		MPI_Reduce(tmpnonsynwaitvector,nonsynwaitvector,GetNcomponent()*GetNstate(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		// MPI_Barrier(MPI_COMM_WORLD);

		int rm = 0;
		int pm = 0;
		int wm = 0;
		int nswm = 0;
		for (int cat=0; cat<GetNcomponent(); cat++)	{
			profilepaircount[cat].clear();
			profilerootcount[cat].clear();
			profilewaitingtime[cat].clear();
			profilenonsynwaitingtime[cat].clear();

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
				if (nonsynwaitvector[nswm])	{
					profilenonsynwaitingtime[cat][i] = nonsynwaitvector[nswm];
				}
				nswm++;
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
		// MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(rootvector,GetNcomponent()*GetNstate(),MPI_INT,0,MPI_COMM_WORLD);
		// MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(waitvector,GetNcomponent()*GetNstate(),MPI_DOUBLE,0,MPI_COMM_WORLD);
		// MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(nonsynwaitvector,GetNcomponent()*GetNstate(),MPI_DOUBLE,0,MPI_COMM_WORLD);
		// MPI_Barrier(MPI_COMM_WORLD);

		delete[] pairvector;
		delete[] tmppairvector;
		delete[] rootvector;
		delete[] tmprootvector;
		delete[] waitvector;
		delete[] tmpwaitvector;
		delete[] nonsynwaitvector;
		delete[] tmpnonsynwaitvector;
	}
	else	{
		UpdateModeProfileSuffStat();
	}
}


void AACodonMutSelSiteMatrixMixtureProfileProcess::SlaveUpdateModeProfileSuffStat()	{

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

	double* nonsynwaitvector = new double[GetNcomponent() * GetNstate()];
	for (int i=0; i<GetNcomponent() * GetNstate(); i++)	{
		nonsynwaitvector[i] = 0;
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

	double* tmpnonsynwaitvector = new double[GetNcomponent() * GetNstate()];
	for (int i=0; i<GetNcomponent() * GetNstate(); i++)	{
		tmpnonsynwaitvector[i] = 0;
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
		for (map<int,double>::iterator i = profilenonsynwaitingtime[cat].begin(); i!= profilenonsynwaitingtime[cat].end(); i++)	{
			nonsynwaitvector[cat*GetNstate() + i->first] = i->second;
			tmpnonsynwaitvector[cat*GetNstate() + i->first] = i->second;
		}
		for (map<pair<int,int>, int>::iterator i = profilepaircount[cat].begin(); i!= profilepaircount[cat].end(); i++)	{
			pairvector[cat*GetNstate()*GetNstate() + i->first.first*GetNstate() + i->first.second] = i->second;
			tmppairvector[cat*GetNstate()*GetNstate() + i->first.first*GetNstate() + i->first.second] = i->second;
		}
	}

	MPI_Reduce(tmppairvector,pairvector,GetNcomponent()*GetNstate()*GetNstate(),MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
	// MPI_Barrier(MPI_COMM_WORLD);
	MPI_Reduce(tmprootvector,rootvector,GetNcomponent()*GetNstate(),MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
	// MPI_Barrier(MPI_COMM_WORLD);
	MPI_Reduce(tmpwaitvector,waitvector,GetNcomponent()*GetNstate(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	// MPI_Barrier(MPI_COMM_WORLD);
	MPI_Reduce(tmpnonsynwaitvector,nonsynwaitvector,GetNcomponent()*GetNstate(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	// MPI_Barrier(MPI_COMM_WORLD);

	delete[] tmppairvector;
	delete[] tmprootvector;
	delete[] tmpwaitvector;
	delete[] tmpnonsynwaitvector;

	MPI_Bcast(pairvector,GetNcomponent()*GetNstate()*GetNstate(),MPI_INT,0,MPI_COMM_WORLD);
	// MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(rootvector,GetNcomponent()*GetNstate(),MPI_INT,0,MPI_COMM_WORLD);
	// MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(waitvector,GetNcomponent()*GetNstate(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	// MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(nonsynwaitvector,GetNcomponent()*GetNstate(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	// MPI_Barrier(MPI_COMM_WORLD);

	int rm = 0;
	int pm = 0;
	int wm = 0;
	int nswm = 0;
	for (int cat=0; cat<GetNcomponent(); cat++)	{
		profilepaircount[cat].clear();
		profilerootcount[cat].clear();
		profilewaitingtime[cat].clear();
		profilenonsynwaitingtime[cat].clear();

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
			if (nonsynwaitvector[nswm])	{
				profilenonsynwaitingtime[cat][i] = nonsynwaitvector[nswm];
			}
			nswm++;
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
	delete[] nonsynwaitvector;
}

void AACodonMutSelSiteMatrixMixtureProfileProcess::UpdateModeProfileSuffStat()	{

	// UpdateMatrices();
	for (int i=0; i<GetNcomponent(); i++)	{
		profilepaircount[i].clear();
		profilerootcount[i].clear();
		profilewaitingtime[i].clear();
		profilenonsynwaitingtime[i].clear();
	}
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i))	{
			int rootstate = GetSiteRootState(i);
			int cat = alloc[i];
			
			if (rootstate != -1)	{

				// double omega = GetSiteOmega(site);
				AACodonMutSelProfileSubMatrix* sitematrix = GetCodonMatrix(i);
				double omega = sitematrix->GetOmega();

				map<pair<int,int>, int>& paircount = GetSitePairCount(i);
				map<int,double>& waitingtime = GetSiteWaitingTime(i);

				profilerootcount[cat][rootstate]++;
				for (map<pair<int,int>, int>::iterator i = paircount.begin(); i!= paircount.end(); i++)	{
					profilepaircount[cat][i->first] += i->second;
				}
				for (map<int,double>::iterator i = waitingtime.begin(); i!= waitingtime.end(); i++)	{
					profilewaitingtime[cat][i->first] += i->second;
					profilenonsynwaitingtime[cat][i->first] += i->second * omega;
				}
			}
		}
	}
}

double AACodonMutSelSiteMatrixMixtureProfileProcess::LogStatProb(int site, int cat)	{

	double total = 0;
	int rootstate = GetSiteRootState(site);
	if (rootstate != -1)	{

		AACodonMutSelProfileSubMatrix* mat = GetCodonMatrix(site);
		mat->SetAAProfile(profile[cat]);
		const double* stat = mat->GetStationary();

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
	if (isnan(total))	{
		cerr << "in AACodonMutSelSiteMatrixMixtureProfileProcess::LogStatProb: nan\n";
		exit(1);
	}
	if (isinf(total))	{
		cerr << "in AACodonMutSelSiteMatrixMixtureProfileProcess::LogStatProb: inf\n";
		cerr << total << '\n';
		cerr << GetSiteOmega(site) << '\n';
		exit(1);
	}
	return total;
}

/*
double AACodonMutSelSiteMatrixMixtureProfileProcess::LogStatProb(int site, int cat)	{

	double total = 0;
	int rootstate = GetSiteRootState(site);
	if (rootstate != -1)	{

		AACodonMutSelProfileSubMatrix* mat = GetComponentCodonMatrix(cat);
		const double* stat = matrixarray[cat]->GetStationary();

		total += log(stat[rootstate]);

		// double omega = GetSiteOmega(site);
		AACodonMutSelProfileSubMatrix* sitematrix = GetCodonMatrix(site);
		double omega = sitematrix->GetOmega();

		map<int,double>& waitingtime = GetSiteWaitingTime(site);
		map<pair<int,int>, int>& paircount = GetSitePairCount(site);

		for (map<pair<int,int>, int>::iterator i = paircount.begin(); i!= paircount.end(); i++)	{
			total += i->second * log((*mat)(i->first.first, i->first.second));
		}

		for (map<int,double>::iterator i = waitingtime.begin(); i!= waitingtime.end(); i++)	{
			double tmp = - i->second * (mat->RateAwaySyn(i->first) + omega * mat->RateAwayNonsyn(i->first));
			total += tmp;
		}
	}
}
*/

double AACodonMutSelSiteMatrixMixtureProfileProcess::ProfileSuffStatLogProb(int cat)	{

	double total = 0;
	AACodonMutSelProfileSubMatrix* mat = GetComponentCodonMatrix(cat);
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
		total -= i->second * mat->RateAwaySyn(i->first);
	}
	for (map<int,double>::iterator i = profilenonsynwaitingtime[cat].begin(); i!= profilenonsynwaitingtime[cat].end(); i++)	{
		total -= i->second * mat->RateAwayNonsyn(i->first);
	}
	for (map<pair<int,int>, int>::iterator i = profilepaircount[cat].begin(); i!= profilepaircount[cat].end(); i++)	{
		total += i->second * log((*mat)(i->first.first, i->first.second));
	}
	profilesuffstatlogprob[cat] = total;
	return total;
}

double AACodonMutSelSiteMatrixMixtureProfileProcess::CountProfileSuffStatLogProb(int cat)	{

	double total = 0;
	AACodonMutSelProfileSubMatrix* mat = GetComponentCodonMatrix(cat);
	const double* stat = matrixarray[cat]->GetStationary();
	for (map<int,int>::iterator i = profilerootcount[cat].begin(); i!= profilerootcount[cat].end(); i++)	{
		total += i->second * log(stat[i->first]);
	}
	for (map<pair<int,int>, int>::iterator i = profilepaircount[cat].begin(); i!= profilepaircount[cat].end(); i++)	{
		total += i->second * log((*mat)(i->first.first, i->first.second));
	}
	return total;
}

double AACodonMutSelSiteMatrixMixtureProfileProcess::BetaProfileSuffStatLogProb(int cat)	{

	double total = 0;
	AACodonMutSelProfileSubMatrix* mat = GetComponentCodonMatrix(cat);
	for (map<int,double>::iterator i = profilewaitingtime[cat].begin(); i!= profilewaitingtime[cat].end(); i++)	{
		total -= i->second * mat->RateAwaySyn(i->first);
	}
	for (map<int,double>::iterator i = profilenonsynwaitingtime[cat].begin(); i!= profilenonsynwaitingtime[cat].end(); i++)	{
		total -= i->second * mat->RateAwayNonsyn(i->first);
	}
	return total;
}
