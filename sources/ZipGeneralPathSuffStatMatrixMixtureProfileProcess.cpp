
#include "ZipGeneralPathSuffStatMatrixMixtureProfileProcess.h"
#include "Parallel.h"

void ZipGeneralPathSuffStatMatrixMixtureProfileProcess::UpdateModeProfileSuffStat()	{

	if (GetMyid())	{
		cerr << "error: slave in GPSSMatrixMixtureProfileProcess::UpdateModeProfileSuffStat\n";
		exit(1);
	}

	for (int i=0; i<GetNcomponent(); i++)	{
		profilepaircount[i].clear();
		profilerootcount[i].clear();
		profilewaitingtime[i].clear();
	}
	for (int site=0; site<GetNsite(); site++)	{

		if (ActiveSite(site))	{
		
			map<pair<int,int>, int>& paircount = GetSitePairCount(site);
			map<int,double>& waitingtime = GetSiteWaitingTime(site);
			int rootstate = GetSiteRootState(site);
			int cat = alloc[site];
			int orbitsize = GetOrbitSize(site);

			if (rootstate != orbitsize)	{
				profilerootcount[cat][GetStateFromZip(site,rootstate)]++;
			}
			for (map<int,double>::iterator i = waitingtime.begin(); i!= waitingtime.end(); i++)	{
				if (i->first != orbitsize)	{
					profilewaitingtime[cat][GetStateFromZip(site,i->first)] += i->second;
				}
			}
			for (map<pair<int,int>, int>::iterator i = paircount.begin(); i!= paircount.end(); i++)	{
				if ((i->first.first != orbitsize) && (i->first.second != orbitsize))	{
					profilepaircount[cat][pair<int,int>(GetStateFromZip(site,i->first.first),GetStateFromZip(site,i->first.second))] += i->second;
				}
			}
		}
	}
}

double ZipGeneralPathSuffStatMatrixMixtureProfileProcess::LogStatProb(int site, int cat)	{

	// calculate based on zip matrix

	if (! GetMyid())	{
		cerr << "error: master in zip log stat prob\n";
		exit(1);
	}
	double total = 0;
	ZipSubMatrix* mat = zipmatrixarray[site];
	if (mat->GetFrom() != matrixarray[cat])	{
		mat->SetFrom(matrixarray[cat]);
		mat->CorruptMatrix();
	}

	const double* stat = mat->GetStationary();
	double** q = mat->GetQ();
	int nstate = mat->GetNstate();

	total += log(stat[GetSiteRootState(site)]);

	map<int,double>& waitingtime = GetSiteWaitingTime(site);
	for (map<int,double>::iterator i = waitingtime.begin(); i!= waitingtime.end(); i++)	{
		total += i->second * q[i->first][i->first];
	}

	map<pair<int,int>, int>& paircount = GetSitePairCount(site);
	for (map<pair<int,int>, int>::iterator i = paircount.begin(); i!= paircount.end(); i++)	{
		total += i->second * log(q[i->first.first][i->first.second]);
	}

	return total;
}

double ZipGeneralPathSuffStatMatrixMixtureProfileProcess::ProfileSuffStatLogProb()	{

	if (GetMyid())	{
		cerr << "error: slave in Zip ProfileSuffStatLogProb\n";
		exit(1);
	}
	GlobalProfileSuffStatLogProb();
	double total = 0;
	for (int k=0; k<Ncomponent; k++)	{
		total += profilesuffstatlogprob[k];
	}
	return total;
}

double ZipGeneralPathSuffStatMatrixMixtureProfileProcess::ProfileSuffStatLogProb(int cat)	{

	return profilesuffstatlogprob[cat];
}

void ZipGeneralPathSuffStatMatrixMixtureProfileProcess::GlobalProfileSuffStatLogProb()	{

	if (GetNprocs() > 1)	{
		sschrono.Start();
		MESSAGE signal = PROFILESSLOGPROB;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Status stat;

		// GlobalUpdateParameters();
		double* tmp = new double[Ncomponent];
		for (int k=0; k<Ncomponent; k++)	{
			profilesuffstatlogprob[k] = 0;
			tmp[k] = 0;
		}
		MPI_Reduce(tmp,profilesuffstatlogprob,Ncomponent,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		delete[] tmp;
		sschrono.Stop();
	}
	else	{
		for (int cat=0; cat<Ncomponent; cat++)	{
			double total = 0;
			for (int i=0; i<GetNsite(); i++)	{
				if (ActiveSite(i))	{
					if (alloc[i] == cat)	{
						total += LogStatProb(i,cat);
					}
				}
			}
			profilesuffstatlogprob[cat] = total;
		}
	}
}

void ZipGeneralPathSuffStatMatrixMixtureProfileProcess::SlaveProfileSuffStatLogProb()	{

	UpdateMatrices();
	double* tmp = new double[Ncomponent];
	for (int k=0; k<Ncomponent; k++)	{
		tmp[k] = 0;
	}
	for (int site=GetSiteMin(); site<GetSiteMax(); site++)	{
		if (ActiveSite(site))	{
			tmp[alloc[site]] += LogStatProb(site,alloc[site]);
		}
	}

	// send
	MPI_Reduce(tmp,profilesuffstatlogprob,Ncomponent,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	delete[] tmp;
}


double ZipGeneralPathSuffStatMatrixMixtureProfileProcess::GlobalMoveProfile(double tuning, int n, int nrep)	{

	UpdateOccupancyNumbers();

	// send PROFILE_MOVE Message with n and nrep and tuning
	MESSAGE signal = PROFILE_MOVE;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	int* itmp = new int[3+GetNsite()];
	itmp[0] = n;
	itmp[1] = nrep;
	itmp[2] = Ncomponent;
	for (int i=0; i<GetNsite(); i++)	{
		itmp[3+i] = alloc[i];
	}
	MPI_Bcast(itmp,3+GetNsite(),MPI_INT,0,MPI_COMM_WORLD);
	delete[] itmp;

	double* bkallocprofile = new double[Ncomponent * GetDim()];
	double** bkprofile = new double*[Ncomponent];
	for (int k=0; k<Ncomponent; k++)	{
		bkprofile[k] = bkallocprofile + k*GetDim();
	}
	double* profilelogratio = new double[Ncomponent];
	double* bkprofilesuffstatlogprob = new double[Ncomponent];
	
	int NAccepted = 0;

	// get current profile suffstat logprobs
	GlobalProfileSuffStatLogProb();

	for (int rep=0; rep<nrep; rep++)	{

		// backup profilesuffstatlogprobs
		for (int k=0; k<Ncomponent; k++)	{
			bkprofilesuffstatlogprob[k] = profilesuffstatlogprob[k];
		}

		// propose move to all profiles
		for (int k=0; k<Ncomponent; k++)	{
			if (occupancy[k])	{
				profilelogratio[k] = -profilesuffstatlogprob[k];
				for (int l=0; l<GetDim(); l++)	{
					bkprofile[k][l] = profile[k][l];
				}
				profilelogratio[k] += ProfileProposeMove(profile[k],tuning,n,GetDim(),k,0);
			}
		}

		// send new profiles
		MPI_Bcast(allocprofile,Ncomponent*GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);

		// get profile suff stat log probs
		GlobalProfileSuffStatLogProb();

		// decide acceptance on a per-mode basis 
		for (int k=0; k<Ncomponent; k++)	{
			if (occupancy[k])	{
				profilelogratio[k] += profilesuffstatlogprob[k];
				int accept = (log(rnd::GetRandom().Uniform()) < profilelogratio[k]);
				if (accept)	{
					NAccepted++;
				}
				else	{
					profilesuffstatlogprob[k] = bkprofilesuffstatlogprob[k];
					for (int l=0; l<GetDim(); l++)	{
						profile[k][l] = bkprofile[k][l];
					}
				}
			}
		}
	}

	// resample empty profiles
	ResampleEmptyProfiles();

	// resend all profiles
	MPI_Bcast(allocprofile,Ncomponent*GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);

	// check that profiles are normalized
	for (int k=0; k<Ncomponent; k++)	{
		double tot = 0;
		for (int i=0; i<GetDim(); i++)	{
			tot += profile[k][i];
		}
		if (fabs(tot - 1) > 1e-6)	{
			cerr << "normalization error : " << tot -1 << '\n';
			exit(1);
		}
	}

	UpdateMatrices();

	delete[] bkprofilesuffstatlogprob;
	delete[] profilelogratio;
	delete[] bkprofile;
	delete[] bkallocprofile;

	return ((double) NAccepted) / Ncomponent / nrep;
}

void ZipGeneralPathSuffStatMatrixMixtureProfileProcess::SlaveMoveProfile()	{

	// parse arguments sent by master

	int* itmp = new int[3+GetNsite()];
	MPI_Bcast(itmp,3+GetNsite(),MPI_INT,0,MPI_COMM_WORLD);
	int n = itmp[0];
	int nrep = itmp[1];
	Ncomponent = itmp[2];
	for (int i=0; i<GetNsite(); i++)	{
		alloc[i] = itmp[3+i];
	}
	delete[] itmp;

	// send profile suffstat logprobs
	MESSAGE signal;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	SlaveProfileSuffStatLogProb();

	for (int rep=0; rep<nrep; rep++)	{

		// receive new profiles
		MPI_Bcast(allocprofile,Ncomponent*GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);
		UpdateMatrices();

		// send profile suffstat logprobs
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
		SlaveProfileSuffStatLogProb();
	}

	// rereceive all profiles
	MPI_Bcast(allocprofile,Ncomponent*GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	UpdateMatrices();
}

