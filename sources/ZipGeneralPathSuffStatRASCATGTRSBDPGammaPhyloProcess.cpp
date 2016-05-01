
#include "ZipGeneralPathSuffStatRASCATGTRSBDPGammaPhyloProcess.h"
#include "Parallel.h"

void ZipGeneralPathSuffStatGTRSBDPProfileProcess::Create()	{
	if (! rrsuffstatcount)	{
		GeneralPathSuffStatGTRSBDPProfileProcess::Create();
		ZipGeneralPathSuffStatMatrixMixtureProfileProcess::Create();
		rrsuffstatcount = new int[Nrr];
		rrsuffstatbeta = new double[Nrr];
		profilesuffstatcount = new int*[GetNmodeMax()];
		profilesuffstatbeta = new double*[GetNmodeMax()];
		allocprofilesuffstatcount = new int[GetNmodeMax() * GetDim()];
		allocprofilesuffstatbeta = new double[GetNmodeMax() * GetDim()];
		for (int i=0; i<GetNmodeMax(); i++)	{
			profilesuffstatcount[i] = allocprofilesuffstatcount + i*GetDim();
			profilesuffstatbeta[i] = allocprofilesuffstatbeta + i*GetDim();
		}
	}
}

void ZipGeneralPathSuffStatGTRSBDPProfileProcess::Delete()	{
	if (rrsuffstatcount)	{
		ZipGeneralPathSuffStatMatrixMixtureProfileProcess::Delete();
		GeneralPathSuffStatGTRSBDPProfileProcess::Delete();
		delete[] profilesuffstatcount;
		delete[] profilesuffstatbeta;
		delete[] allocprofilesuffstatcount;
		delete[] allocprofilesuffstatbeta;
		delete[] rrsuffstatcount;
		delete[] rrsuffstatbeta;
		rrsuffstatcount = 0;
	}
}

void ZipGeneralPathSuffStatGTRSBDPProfileProcess::UpdateRRSuffStat()	{

	ZipGeneralPathSuffStatMatrixMixtureProfileProcess::UpdateModeProfileSuffStat();

	for (int k=0; k<GetNrr(); k++)	{
		rrsuffstatcount[k] = 0;
		rrsuffstatbeta[k] = 0;
	}
	for (int cat=0; cat<Ncomponent; cat++)	{

		for (map<int,double>::iterator i = profilewaitingtime[cat].begin(); i!= profilewaitingtime[cat].end(); i++)	{

			int initstate = i->first;
			for (int finalstate=0; finalstate<GetDim(); finalstate++)	{
				if (finalstate != initstate)	{
					int r = rrindex(initstate,finalstate,GetDim());
					rrsuffstatbeta[r] += i->second * profile[cat][finalstate];
				}
			}
		}

		for (map<pair<int,int>, int>::iterator i = profilepaircount[cat].begin(); i!= profilepaircount[cat].end(); i++)	{

			int initstate = i->first.first;
			int finalstate = i->first.second;
			int r = rrindex(initstate,finalstate,GetDim());
			rrsuffstatcount[r]++;
		}
	}
	/*
	for (int k=0; k<GetNrr(); k++)	{
		cerr << (1.0 + rrsuffstatcount[k]) / (1.0 + rrsuffstatbeta[k]) << '\t';
	}
	cerr << '\n';
	*/
}

void ZipGeneralPathSuffStatGTRSBDPProfileProcess::GlobalUpdateModeProfileSuffStat()	{

	ZipGeneralPathSuffStatMatrixMixtureProfileProcess::UpdateModeProfileSuffStat();
	ZipGeneralPathSuffStatGTRSBDPProfileProcess::UpdateModeProfileSuffStat();

	if (GetNprocs() > 1)	{
		MESSAGE signal = UPDATE_MPROFILE;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

		MPI_Bcast(allocprofilesuffstatcount,Ncomponent*GetDim(),MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(allocprofilesuffstatbeta,Ncomponent*GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	}
}

void ZipGeneralPathSuffStatGTRSBDPProfileProcess::SlaveUpdateModeProfileSuffStat()	{

	MPI_Bcast(allocprofilesuffstatcount,Ncomponent*GetDim(),MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(allocprofilesuffstatbeta,Ncomponent*GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);
}

void ZipGeneralPathSuffStatGTRSBDPProfileProcess::UpdateModeProfileSuffStat()	{

	if (GetMyid())	{
		cerr << "error: slave in ExpoConjugateGTRMixtureProfileProcess::UpdateModeProfileSuffStat()\n";
		exit(1);
	}
	for (int cat=0; cat<GetNcomponent(); cat++)	{
		for (int k=0; k<GetDim(); k++)	{
			profilesuffstatcount[cat][k] = 0;
			profilesuffstatbeta[cat][k] = 0;
		}

		for (map<int,double>::iterator i = profilewaitingtime[cat].begin(); i!= profilewaitingtime[cat].end(); i++)	{

			int initstate = i->first;
			for (int finalstate=0; finalstate<GetDim(); finalstate++)	{
				if (finalstate != initstate)	{
					int r = rrindex(initstate,finalstate,GetDim());
					profilesuffstatbeta[cat][finalstate] += i->second * rr[r];
				}
			}
		}

		for (map<pair<int,int>, int>::iterator i = profilepaircount[cat].begin(); i!= profilepaircount[cat].end(); i++)	{

			int finalstate = i->first.second;
			profilesuffstatcount[cat][finalstate]++;
		}

		for (map<int,int>::iterator i = profilerootcount[cat].begin(); i!= profilerootcount[cat].end(); i++)	{
			profilesuffstatcount[cat][i->first] += i->second;
		}
	}
}

double ZipGeneralPathSuffStatGTRSBDPProfileProcess::LogStatProb(int site, int cat)	{

	// return ZipMatrixMixtureProfileProcess::LogStatProb(site,cat);
	double total = 0;
	int zipsize = GetZipSize(site);
	int nstate = GetDim();
	int* zipindices = GetZipIndices(site);
	double zipstat[zipsize];
	for (int i=0; i<zipsize; i++)	{
		zipstat[i] = 0;
	}
	for (int i=0; i<nstate; i++)	{
		zipstat[zipindices[i]] += profile[cat][i];
	}

	total += log(zipstat[GetSiteRootState(site)]);

	map<pair<int,int>, int>& paircount = GetSitePairCount(site);
	for (map<pair<int,int>, int>::iterator i = paircount.begin(); i!= paircount.end(); i++)	{
		int i0 = i->first.first;
		int j0 = i->first.second;
		double tot = 0;
		for (int k=0; k<nstate; k++)	{
			int a = zipindices[k];
			if (a == i0)	{
				for (int l=0; l<nstate; l++)	{
					int b = zipindices[l];
					if (b == j0)	{
						tot += profile[cat][k] * profile[cat][l] * rr[rrindex(k,l,nstate)];
					}
				}
			}
		}
		tot /= zipstat[i0];
		total += i->second * log(tot);
	}

	map<int,double>& waitingtime = GetSiteWaitingTime(site);
	for (map<int,double>::iterator i = waitingtime.begin(); i!= waitingtime.end(); i++)	{
		int i0 = i->first;
		double tot = 0;
		for (int k=0; k<nstate; k++)	{
			if (zipindices[k] == i0)	{
				for (int l=0; l<nstate; l++)	{
					if (zipindices[l] != i0)	{
						tot += profile[cat][k] * profile[cat][l] * rr[rrindex(k,l,nstate)];
					}
				}
			}
		}
		tot /= zipstat[i0];
		total -= i->second * tot;
	}

	return total;
}
	
double ZipGeneralPathSuffStatGTRSBDPProfileProcess::ProfileSuffStatLogProb()	{

	if (proposemode)	{
		return MixtureProfileProcess::ProfileSuffStatLogProb();
	}
	return ZipGeneralPathSuffStatMatrixMixtureProfileProcess::ProfileSuffStatLogProb();
}

double ZipGeneralPathSuffStatGTRSBDPProfileProcess::ProfileSuffStatLogProb(int cat)	{

	if (! proposemode)	{
		cerr << "error : in ZipGTRSBSDP::ZipMatrixMixture::ProfileSuffStatLogProb(cat)\n";
	}
	return ApproxProfileSuffStatLogProb(cat);
}

double  ZipGeneralPathSuffStatGTRSBDPProfileProcess::ApproxSuffStatSampleRR()	{

	UpdateRRSuffStat();

	double logh = 0;
	for (int i=0; i<GetNrr(); i++)	{
		logh += rrsuffstatcount[i] * log(rr[i]) - rrsuffstatbeta[i] * rr[i];
		rr[i] = rnd::GetRandom().Gamma(1.0 + rrsuffstatcount[i], 1.0 + rrsuffstatbeta[i]);
		logh -= rrsuffstatcount[i] * log(rr[i]) - rrsuffstatbeta[i] * rr[i];
	}
	return logh;
}

double  ZipGeneralPathSuffStatGTRSBDPProfileProcess::MoveRR(double tuning, int n, int nrep)	{

	MESSAGE signal = RR_MOVE;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&nrep,1,MPI_INT,0,MPI_COMM_WORLD);

	int naccepted = 0;

	double bklogprob = ProfileSuffStatLogProb();
	double bk[Nrr]; 

	for (int i=0; i<nrep; i++)	{

		for (int j=0; j<Nrr; j++)	{
			bk[j] = rr[j];
		}
		
		double deltalogratio = - LogRRPrior() - bklogprob;

		double logh = 0;
		/*
		if (proposemode)	{
			logh = ApproxSuffStatSampleRR();
		}
		else	{
		*/
			int choose[n];
			rnd::GetRandom().DrawFromUrn(choose,n,Nrr);
			for (int j=0; j<n; j++)	{
				double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
				double e = exp(m);
				rr[choose[j]] *= e;
				logh += m;
			}
		// }

		deltalogratio += logh;

		MPI_Bcast(rr,Nrr,MPI_DOUBLE,0,MPI_COMM_WORLD);

		double logprob = ProfileSuffStatLogProb();

		deltalogratio += LogRRPrior() + logprob;
		// cerr << deltalogratio -logh << '\t' << logh << '\t' << deltalogratio << '\n';;
		int accepted = (log(rnd::GetRandom().Uniform()) < deltalogratio);
		if (accepted)	{
			naccepted++;
			bklogprob = logprob;
		}
		else	{
			for (int j=0; j<Nrr; j++)	{
				rr[j] = bk[j];
			}
		}
	}

	MPI_Bcast(rr,Nrr,MPI_DOUBLE,0,MPI_COMM_WORLD);
	UpdateMatrices();
	GlobalProfileSuffStatLogProb();

	return ((double) naccepted) / nrep;
}

void ZipGeneralPathSuffStatGTRSBDPProfileProcess::SlaveMoveRR()	{

	int nrep;
	MPI_Bcast(&nrep,1,MPI_INT,0,MPI_COMM_WORLD);

	// send profile suffstat logprobs
	MESSAGE signal;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	SlaveProfileSuffStatLogProb();

	for (int rep=0; rep<nrep; rep++)	{
		MPI_Bcast(rr,Nrr,MPI_DOUBLE,0,MPI_COMM_WORLD);

		// send profile suffstat logprobs
		MESSAGE signal;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
		SlaveProfileSuffStatLogProb();
	}
	MPI_Bcast(rr,Nrr,MPI_DOUBLE,0,MPI_COMM_WORLD);

	// send profile suffstat logprobs
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	SlaveProfileSuffStatLogProb();

	// UpdateMatrices();
}

double ZipGeneralPathSuffStatGTRSBDPProfileProcess::ApproxProfileSuffStatLogProb(int cat)	{

	double tot = 0;
	for (int k=0; k<GetDim(); k++)	{
		tot += profilesuffstatcount[cat][k] * log(profile[cat][k]) - profilesuffstatbeta[cat][k] * profile[cat][k];
	}
	return tot;
}

double ZipGeneralPathSuffStatGTRSBDPProfileProcess::ApproxSuffStatProfileProposeMove(int cat, int nstep)	{

	if (! proposemode)	{
		cerr << "error: approx suff stat propose move called when not in special mode\n";
		exit(1);
	}

	double logh = -ApproxProfileSuffStatLogProb(cat);

	double bkprofile[GetDim()];
	for (int k=0; k<GetDim(); k++)	{
		bkprofile[k] = profile[cat][k];
	}

	for (int rep=0; rep<nstep; rep++)	{

		double logratio = -ApproxProfileSuffStatLogProb(cat);

		double r = rnd::GetRandom().Uniform();
		if (r < 1.0 / 3)	{
			logratio += ProfileProposeMove(profile[cat],1.0,1,0,cat,0);
		}
		else if (r < 2.0 / 3)	{
			logratio += ProfileProposeMove(profile[cat],0.3,3,0,cat,0);
		}
		else	{
			logratio += ProfileProposeMove(profile[cat],0.1,5,0,cat,0);
		}

		logratio += ApproxProfileSuffStatLogProb(cat);

		if (! log(rnd::GetRandom().Uniform()) < logratio)	{
			for (int k=0; k<GetDim(); k++)	{
				profile[cat][k] = bkprofile[k];
			}
		}
	}

	logh += ApproxProfileSuffStatLogProb(cat);

	return logh;
}


double ZipGeneralPathSuffStatGTRSBDPProfileProcess::MultipleTryGlobalMixMove(int nrep, int nallocrep, double epsilon, int nprofilerep, int nmultryrep, double tuning, int nentry)	{

	if (Ncomponent != GetNmodeMax())	{
		cerr << "error in sbdp inc dp move: number of components\n";
		exit(1);
	}

	double NAccepted = 0;
	double NTried = 0;

	// define threshold between GIbbs and MH
	int K0 = GetNmodeMax();
	if (epsilon)	{
		double r = kappa / (1 + kappa);
		K0 = (int) (log(epsilon) / log(r));
		if (K0 >= GetNmodeMax())	{
			K0 = GetNmodeMax();
		}
	}
	// K0 = GetNmodeMax();

	// cerr << K0 << '\t' << GetNmodeMax() << '\t';

	// Chrono ch;
	// ch.Start();

	// send mixmove signal and tuning parameters
	MESSAGE signal = MTRYMIX_MOVE;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	int itmp[6];
	itmp[0] = nrep;
	itmp[1] = nallocrep;
	itmp[2] = K0;
	itmp[3] = nprofilerep;
	itmp[4] = nmultryrep;
	itmp[5] = nentry;

	MPI_Bcast(itmp,6,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&tuning,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

	int nocc = Ncomponent;
	int ntry = nmultryrep * (GetNprocs() - 1);

	double* allocbkprofile = new double[nocc * GetDim()];
	double** bkprofile = new double*[nocc];
	for (int k=0; k<nocc; k++)	{
		bkprofile[k] = allocbkprofile + k*GetDim();
	}

	double* allocnewprofile = new double[ntry * nocc * (GetDim() + 1)];
	double*** proposeprofile = new double**[ntry];
	double* tmp = allocnewprofile;
	for (int k=0; k<ntry; k++)	{
		proposeprofile[k] = new double*[nocc];
		for (int l=0; l<nocc; l++)	{
			proposeprofile[k][l] = tmp;
			tmp += GetDim() + 1;
		}
	}

	double* tmpallocsuffstatlogprob = new double[nocc*ntry];

	double* allocfwdsuffstatlogprob = new double[nocc*ntry];
	double* allocrevsuffstatlogprob = new double[nocc*ntry];
	double** fwdsuffstatlogprob = new double*[ntry];
	double** revsuffstatlogprob = new double*[ntry];
	for (int k=0; k<ntry; k++)	{
		fwdsuffstatlogprob[k] = allocfwdsuffstatlogprob + k*nocc;
		revsuffstatlogprob[k] = allocrevsuffstatlogprob + k*nocc;
	}
	
	for (int rep=0; rep<nrep; rep++)	{

		ResampleWeights();

		// mpi send weights
		MPI_Bcast(weight,Ncomponent,MPI_DOUBLE,0,MPI_COMM_WORLD);

		// here slaves do realloc moves

		// mpi receive new allocations
		MPI_Status stat;
		int tmpalloc[GetNsite()];
		for(int i=1; i<GetNprocs(); ++i) {
			MPI_Recv(tmpalloc,GetNsite(),MPI_INT,i,TAG1,MPI_COMM_WORLD,&stat);
			for(int j=GetProcSiteMin(i); j<GetProcSiteMax(i); j++) {
				if (ActiveSite(j))	{
					alloc[j] = tmpalloc[j];
					if ((alloc[j] < 0) || (alloc[j] >= Ncomponent))	{
						cerr << "alloc overflow\n";
						exit(1);
					}
				}
			}
		}

		MPI_Bcast(alloc,GetNsite(),MPI_INT,0,MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);

		UpdateOccupancyNumbers();

		/*
		int nocc = GetNmodeMax()-1;
		while (!occupancy[nocc])	{
			nocc--;
		}
		nocc++;

		MPI_Bcast(&nocc,1,MPI_INT,0,MPI_COMM_WORLD);
		*/

		// get current profile suffstat logprobs
		GlobalProfileSuffStatLogProb();

		for (int prep=0; prep<nprofilerep; prep++)	{

			// backup profiles
			for (int k=0; k<nocc; k++)	{
				if (occupancy[k])	{
					for (int l=0; l<GetDim(); l++)	{
						bkprofile[k][l] = profile[k][l];
					}
				}
			}

			// receive new profiles from slaves
			MPI_Status stat;
			double* tmpalloc = allocnewprofile;
			for(int i=1; i<GetNprocs(); ++i) {
				MPI_Recv(tmpalloc,nocc*nmultryrep*(GetDim() + 1),MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
				tmpalloc += nocc*nmultryrep*(GetDim() + 1);
			}

			// resend all profiles to slaves
			MPI_Bcast(allocnewprofile,nocc*ntry*(GetDim()+1),MPI_DOUBLE,0,MPI_COMM_WORLD);

			// get profilesuffstatlogprobs from slaves
			for (int k=0; k<nocc*ntry; k++)	{
				allocfwdsuffstatlogprob[k] = 0;
				tmpallocsuffstatlogprob[k] = 0;
			}
			MPI_Reduce(tmpallocsuffstatlogprob,allocfwdsuffstatlogprob,nocc*ntry,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

			// calculate probs and choose new profile
			double logweightfwd[nocc];
			int choose[nocc];

			for (int k=0; k<nocc; k++)	{
				if (occupancy[k])	{
					double max = 0;
					double logweight[ntry];
					for (int l=0; l<ntry; l++)	{

						if (isnan(fwdsuffstatlogprob[l][k]))	{
							cerr << "fwd ss log prob is nan\n";
							exit(1);
						}
						// log lambda + log proposal density + logprofilelikelihood
						logweight[l] = LogFrequencyStatPrior(proposeprofile[l][k]) + fwdsuffstatlogprob[l][k] - proposeprofile[l][k][GetDim()];
						if (isnan(logweight[l]))	{
							cerr << "fwd logweight is nan\n";
							cerr <<  LogFrequencyStatPrior(proposeprofile[l][k])  << '\t' <<  fwdsuffstatlogprob[l][k] << '\t' <<  proposeprofile[l][k][GetDim()] << '\n';
							double tot = 0;
							for (int j=0; j<GetDim(); j++)	{
								tot += proposeprofile[l][k][j];
								cerr << proposeprofile[l][k][j] << '\t';
							}
							cerr << '\n';
							cerr << tot << '\n';
							for (int j=0; j<GetDim(); j++)	{
								cerr << profile[k][j] << '\t';
							}
							cerr << '\n';
							exit(1);
						}
						if ((!l) || (max < logweight[l]))	{
							max = logweight[l];
						}
					}
					double cumul[ntry];
					double tot = 0;
					for (int l=0; l<ntry; l++)	{
						tot += exp(logweight[l] - max);
						cumul[l] = tot;
					}
					
					logweightfwd[k] = log(tot) + max;
					if (isnan(logweightfwd[k]))	{
						cerr << "logweight fwd is nan\n";
						exit(1);
					}

					double x = tot * rnd::GetRandom().Uniform();
					int m = 0;
					while ((m < ntry) && (x > cumul[m])) m++;
					if (m == ntry)	{
						cerr << "error in multiple try when choosing among the ntry profiles (forward)\n";
						exit(1);
					}
					choose[k] = m;

					for (int l=0; l<GetDim(); l++)	{
						profile[k][l] = proposeprofile[m][k][l];
					}
				}
			}

			// send array of chosen profiles to slaves
			MPI_Bcast(allocprofile,nocc*GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);

			// receive new profiles from slaves
			tmpalloc = allocnewprofile;
			for(int i=1; i<GetNprocs(); ++i) {
				MPI_Recv(tmpalloc,nocc*nmultryrep*(GetDim() + 1),MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
				tmpalloc += nocc*nmultryrep*(GetDim() + 1);
			}

			// resend all profiles to slaves
			MPI_Bcast(allocnewprofile,nocc*ntry*(GetDim()+1),MPI_DOUBLE,0,MPI_COMM_WORLD);

			// get profilesuffstatlogprobs from slaves
			for (int k=0; k<nocc*ntry; k++)	{
				allocrevsuffstatlogprob[k] = 0;
				tmpallocsuffstatlogprob[k] = 0;
			}
			MPI_Reduce(tmpallocsuffstatlogprob,allocrevsuffstatlogprob,nocc*ntry,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

			// calculate probs and choose new profile
			double logweightrev[nocc];

			// calculate reverse probs
			for (int k=0; k<nocc; k++)	{
				if (occupancy[k])	{
					double max = 0;
					double logweight[ntry];
					for (int l=0; l<ntry; l++)	{

						// log lambda + log proposal density + logprofilelikelihood
						logweight[l] = LogFrequencyStatPrior(proposeprofile[l][k]) + revsuffstatlogprob[l][k] - proposeprofile[l][k][GetDim()];
						if ((!l) || (max < logweight[l]))	{
							max = logweight[l];
						}
					}
					double tot = 0;
					for (int l=0; l<ntry; l++)	{
						tot += exp(logweight[l] - max);
					}
					
					logweightrev[k] = log(tot) + max;
					if (isnan(logweightrev[k]))	{
						cerr << "logweight is nan\n";
						exit(1);
					}
				}
			}

			// decide acceptance on a per-mode basis 
			for (int k=0; k<nocc; k++)	{
				if (occupancy[k])	{
					int accept = (log(rnd::GetRandom().Uniform()) < (logweightfwd[k] - logweightrev[k]));
					if (accept)	{
						NAccepted++;
						profilesuffstatlogprob[k] = fwdsuffstatlogprob[choose[k]][k];
					}
					else	{
						for (int l=0; l<GetDim(); l++)	{
							profile[k][l] = bkprofile[k][l];
						}
					}
					NTried++;
				}
			}

			// resample empty profiles
			ResampleEmptyProfiles();

			// resend all profiles to slaves
			MPI_Bcast(allocprofile,nocc*GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);

			// GlobalProfileSuffStatLogProb();

		}

		// resend all profiles
		// MPI_Bcast(allocprofile,Ncomponent*GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);

	}

	// delete arrays
	delete[] bkprofile;
	delete[] allocbkprofile;

	for (int k=0; k<ntry; k++)	{
		delete[] proposeprofile[k];
	}
	delete[] proposeprofile;
	delete[] allocnewprofile;

	delete[] fwdsuffstatlogprob;
	delete[] revsuffstatlogprob;

	delete[] tmpallocsuffstatlogprob;
	delete[] allocfwdsuffstatlogprob;
	delete[] allocrevsuffstatlogprob;

	// CreateMatrices();
	UpdateMatrices();
	GlobalProfileSuffStatLogProb();

	return NAccepted / NTried;
}

void ZipGeneralPathSuffStatGTRSBDPProfileProcess::MultipleTrySlaveMixMove()	{

	int itmp[6];
	MPI_Bcast(itmp,6,MPI_INT,0,MPI_COMM_WORLD);
	int nrep = itmp[0];
	int nallocrep = itmp[1];
	int K0 = itmp[2];
	int nprofilerep = itmp[3];
	int nmultryrep = itmp[4];
	int nentry = itmp[5];

	double tuning;
	MPI_Bcast(&tuning,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

	double* mLogSamplingArray = new double[Ncomponent];
	double* cumul = new double[Ncomponent];

	int nocc = Ncomponent;
	int ntry = nmultryrep * (GetNprocs() - 1);

	double* allocbkprofile = new double[nocc * GetDim()];
	double** bkprofile = new double*[nocc];
	for (int k=0; k<nocc; k++)	{
		bkprofile[k] = allocbkprofile + k*GetDim();
	}

	double* allocnewprofile = new double[ntry * nocc * (GetDim() + 1)];
	double*** proposeprofile = new double**[ntry];
	double* tmp = allocnewprofile;
	for (int k=0; k<ntry; k++)	{
		proposeprofile[k] = new double*[nocc];
		for (int l=0; l<nocc; l++)	{
			proposeprofile[k][l] = tmp;
			tmp += GetDim() + 1;
		}
	}

	double* tmpallocsuffstatlogprob = new double[nocc*ntry];

	double* allocfwdsuffstatlogprob = new double[nocc*ntry];
	double* allocrevsuffstatlogprob = new double[nocc*ntry];
	double** fwdsuffstatlogprob = new double*[ntry];
	double** revsuffstatlogprob = new double*[ntry];
	for (int k=0; k<ntry; k++)	{
		fwdsuffstatlogprob[k] = allocfwdsuffstatlogprob + k*nocc;
		revsuffstatlogprob[k] = allocrevsuffstatlogprob + k*nocc;
	}
	

	for (int rep=0; rep<nrep; rep++)	{

		// ResampleWeights();

		// receive weights
		MPI_Bcast(weight,Ncomponent,MPI_DOUBLE,0,MPI_COMM_WORLD);

		double totp = 0;
		for (int mode = 0; mode<K0; mode++)	{
			totp += weight[mode];
		}

		double totq = 0;
		for (int mode=K0; mode<GetNmodeMax(); mode++)	{
			totq += weight[mode];
		}

		// realloc move

		int NAccepted = 0;

		for (int allocrep=0; allocrep<nallocrep; allocrep++)	{

			for (int site=GetSiteMin(); site<GetSiteMax(); site++)	{

				if (ActiveSite(site))	{

					int bk = alloc[site];

					double max = 0;
					// double mean = 0;
					for (int mode = 0; mode<K0; mode++)	{
						mLogSamplingArray[mode] = LogStatProb(site,mode);
						if ((!mode) || (max < mLogSamplingArray[mode]))	{
							max = mLogSamplingArray[mode];
						}
						// mean += mLogSamplingArray[mode];
					}
					// mean /= K0;

					double total = 0;
					for (int mode = 0; mode<K0; mode++)	{
						double p = weight[mode] * exp(mLogSamplingArray[mode] - max);
						total += p;
						cumul[mode] = total;
					}
					if (isnan(total))	{
						cerr << "nan\n";
					}

					// double M = exp(mean- max);
					double M = 1;
					total += M * totq;
					double q = total * rnd::GetRandom().Uniform();
					int mode = 0;
					while ( (mode<K0) && (q > cumul[mode])) mode++;
					if (mode == K0)	{
						mode--;
						double r = (q - cumul[mode]) / M;
						while (r > 0)	{
							mode++;
							r -= weight[mode];
						}
					}

					// MH 
					double logratio = 0;
					if (mode >= K0)	{
						logratio += LogStatProb(site,mode) - max - log(M);
					}
					if (bk >= K0)	{
						logratio -= LogStatProb(site,bk) - max - log(M);
					}
					
					if (log(rnd::GetRandom().Uniform()) > logratio)	{
						mode = bk;
					}

					int Accepted = (mode != bk);
					if (Accepted)	{
						NAccepted ++;
					}
					alloc[site] = mode;
				}
			}
		}
		MPI_Send(alloc,GetNsite(),MPI_INT,0,TAG1,MPI_COMM_WORLD);
		MPI_Bcast(alloc,GetNsite(),MPI_INT,0,MPI_COMM_WORLD);

		UpdateOccupancyNumbers();

		MPI_Barrier(MPI_COMM_WORLD);


		// profile move

		/*
		int nocc;
		MPI_Bcast(&nocc,1,MPI_INT,0,MPI_COMM_WORLD);
		*/

		// send profile suffstat logprobs
		MESSAGE signal;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
		SlaveProfileSuffStatLogProb();

		for (int prep=0; prep<nprofilerep; prep++)	{

			// backup profiles
			for (int k=0; k<nocc; k++)	{
				if (occupancy[k])	{
					for (int l=0; l<GetDim(); l++)	{
						bkprofile[k][l] = profile[k][l];
					}
				}
			}

			for (int k=0; k<nocc; k++)	{
				for (int i=0; i<ntry; i++)	{
					for (int l=0; l<GetDim(); l++)	{
						proposeprofile[i][k][l] = 1.0 / GetDim();
					}
				}
			}
			for (int k=0; k<nocc; k++)	{
				if (occupancy[k])	{
					for (int i=0; i<nmultryrep; i++)	{
						int j = i + (GetMyid() - 1)*nmultryrep;
						SpecialProfileProposeMove(profile[k],proposeprofile[j][k],tuning,nentry);
						proposeprofile[j][k][GetDim()] = SpecialProfileProposeMoveLogProb(profile[k],proposeprofile[j][k],tuning,nentry);
					}
				}
			}

			// send profiles to master
			MPI_Send(allocnewprofile + (GetMyid()-1)*nocc*nmultryrep*(GetDim()+1), nocc*nmultryrep*(GetDim() + 1),MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
			// receive all profiles from master
			MPI_Bcast(allocnewprofile,nocc*ntry*(GetDim()+1),MPI_DOUBLE,0,MPI_COMM_WORLD);

			// calculate suff stat log probs
			for (int j=0; j<ntry; j++)	{
				for (int k=0; k<nocc; k++)	{
					if (occupancy[k])	{
						for (int l=0; l<GetDim(); l++)	{
							profile[k][l] = proposeprofile[j][k][l];
						}
					}
				}
				UpdateMatrices();

				for (int k=0; k<nocc; k++)	{
					fwdsuffstatlogprob[j][k] = 0;
				}
				for (int site=GetSiteMin(); site<GetSiteMax(); site++)	{
					if (ActiveSite(site))	{
						if (alloc[site] > nocc)	{
							cerr << "error: allocation out of bound in multiple try\n";
							exit(1);
						}
						fwdsuffstatlogprob[j][alloc[site]] += LogStatProb(site,alloc[site]);
					}
				}
			}

			// send suffstatlogprobs to master
			MPI_Reduce(tmpallocsuffstatlogprob,allocfwdsuffstatlogprob,nocc*ntry,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

			// get array of chosen profiles from master
			MPI_Bcast(allocprofile,nocc*GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);

			// propose new profiles
			// careful about exception for first series of first slave
			for (int k=0; k<nocc; k++)	{
				if (occupancy[k])	{
					for (int i=0; i<nmultryrep; i++)	{
						int j = i + (GetMyid() - 1)*nmultryrep;
						if ((GetMyid() == 1) && (!i))	{
							for (int l=0; l<GetDim(); l++)	{
								proposeprofile[j][k][l] = bkprofile[k][l];
							}
						}
						else	{
							SpecialProfileProposeMove(profile[k],proposeprofile[j][k],tuning,nentry);
						}
						proposeprofile[j][k][GetDim()] = SpecialProfileProposeMoveLogProb(profile[k],proposeprofile[j][k],tuning,nentry);
					}
				}
			}

			// send profiles to master
			MPI_Send(allocnewprofile + (GetMyid()-1)*nocc*nmultryrep*(GetDim()+1), nocc*nmultryrep*(GetDim()+1),MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
			// receive all profiles from master
			MPI_Bcast(allocnewprofile,nocc*ntry*(GetDim()+1),MPI_DOUBLE,0,MPI_COMM_WORLD);

			// calculate suff stat log probs
			for (int j=0; j<ntry; j++)	{
				for (int k=0; k<nocc; k++)	{
					for (int l=0; l<GetDim(); l++)	{
						profile[k][l] = proposeprofile[j][k][l];
					}
				}
				UpdateMatrices();

				for (int k=0; k<nocc; k++)	{
					revsuffstatlogprob[j][k] = 0;
				}
				for (int site=GetSiteMin(); site<GetSiteMax(); site++)	{
					if (ActiveSite(site))	{
						revsuffstatlogprob[j][alloc[site]] += LogStatProb(site,alloc[site]);
					}
				}
			}


			// send suffstatlogprobs to master
			MPI_Reduce(tmpallocsuffstatlogprob,allocrevsuffstatlogprob,nocc*ntry,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

			// get array of chosen profiles from master
			MPI_Bcast(allocprofile,nocc*GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);

			/*
			MESSAGE signal;
			MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
			SlaveProfileSuffStatLogProb();
			*/

		}

		// rereceive all profiles
		// MPI_Bcast(allocprofile,Ncomponent*GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);
		UpdateMatrices();

	}

	// delete arrays
	delete[] bkprofile;
	delete[] allocbkprofile;

	for (int k=0; k<ntry; k++)	{
		delete[] proposeprofile[k];
	}
	delete[] proposeprofile;
	delete[] allocnewprofile;

	delete[] fwdsuffstatlogprob;
	delete[] revsuffstatlogprob;

	delete[] tmpallocsuffstatlogprob;
	delete[] allocfwdsuffstatlogprob;
	delete[] allocrevsuffstatlogprob;

	// send profile suffstat logprobs
	MESSAGE signal;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	SlaveProfileSuffStatLogProb();

	delete[] cumul;
	delete[] mLogSamplingArray;
}

void ZipGeneralPathSuffStatGTRSBDPProfileProcess::SpecialProfileProposeMove(double* oldprofile, double* profile, double tuning, int n)	{

	double statmin = stateps;
	int K = GetDim();

	if (!n)	{ // dirichlet
		double total = 0;
		for (int i=0; i<K; i++)	{
			profile[i] = rnd::GetRandom().sGamma(tuning*oldprofile[i]);
			if (profile[i] < statmin)	{
				profile[i] = statmin;
			}
			total += profile[i];
		}
		for (int i=0; i<K; i++)	{
			profile[i] /= total;
		}
	}
	else	{
		for (int i=0; i<K; i++)	{
			profile[i] = oldprofile[i];
		}
		if (2*n > K)	{
			n = K / 2;
		}
		int* indices = new int[2*n];
		rnd::GetRandom().DrawFromUrn(indices,2*n,K);
		for (int i=0; i<n; i++)	{
			int i1 = indices[2*i];
			int i2 = indices[2*i+1];
			double tot = profile[i1] + profile[i2];
			double x = profile[i1];
			
			double h = tot * tuning * (rnd::GetRandom().Uniform() - 0.5);
			x += h;
			while ((x<0) || (x>tot))	{
				if (x<0)	{
					x = -x;
				}
				if (x>tot)	{
					x = 2*tot - x;
				}
			}
			profile[i1] = x;
			profile[i2] = tot - x;
		}
		delete[] indices;
	}
}

double ZipGeneralPathSuffStatGTRSBDPProfileProcess::SpecialProfileProposeMoveLogProb(double* oldprofile, double* profile, double tuning, int n)	{

	double ret = 0;
	if (!n)	{
		ret = rnd::GetRandom().logGamma(tuning);
		for (int i=0; i<GetDim(); i++)	{
			ret += (tuning * oldprofile[i] -1.0) * log(profile[i]) - rnd::GetRandom().logGamma(tuning*oldprofile[i]);
		}
	}
	return ret;
}

double ZipGeneralPathSuffStatGTRSBDPProfileProcess::ZipGlobalMixMove(int nrep, int nallocrep, double epsilon, int nprofilerep)	{

	if (Ncomponent != GetNmodeMax())	{
		cerr << "error in sbdp inc dp move: number of components\n";
		exit(1);
	}

	int NAccepted = 0;

	// define threshold between GIbbs and MH
	int K0 = GetNmodeMax();
	if (epsilon)	{
		double r = kappa / (1 + kappa);
		K0 = (int) (log(epsilon) / log(r));
		if (K0 >= GetNmodeMax())	{
			K0 = GetNmodeMax();
		}
	}
	// K0 = GetNmodeMax();

	// cerr << K0 << '\t' << GetNmodeMax() << '\t';

	// Chrono ch;
	// ch.Start();

	// send mixmove signal and tuning parameters
	MESSAGE signal = MIX_MOVE;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	int itmp[4];
	itmp[0] = nrep;
	itmp[1] = nallocrep;
	itmp[2] = K0;
	itmp[3] = nprofilerep;
	MPI_Bcast(itmp,4,MPI_INT,0,MPI_COMM_WORLD);


	double* bkallocprofile = new double[Ncomponent * GetDim()];
	double** bkprofile = new double*[Ncomponent];
	for (int k=0; k<Ncomponent; k++)	{
		bkprofile[k] = bkallocprofile + k*GetDim();
	}
	double* profilelogratio = new double[Ncomponent];
	double* bkprofilesuffstatlogprob = new double[Ncomponent];
	
	for (int rep=0; rep<nrep; rep++)	{

		ResampleWeights();

		// mpi send weights
		MPI_Bcast(weight,Ncomponent,MPI_DOUBLE,0,MPI_COMM_WORLD);

		// here slaves do realloc moves

		// mpi receive new allocations
		MPI_Status stat;
		int tmpalloc[GetNsite()];
		for(int i=1; i<GetNprocs(); ++i) {
			MPI_Recv(tmpalloc,GetNsite(),MPI_INT,i,TAG1,MPI_COMM_WORLD,&stat);
			for(int j=GetProcSiteMin(i); j<GetProcSiteMax(i); j++)	{
				if (ActiveSite(j))	{
					alloc[j] = tmpalloc[j];
					if ((alloc[j] < 0) || (alloc[j] >= Ncomponent))	{
						cerr << "alloc overflow\n";
						exit(1);
					}
				}
			}
		}

		UpdateOccupancyNumbers();

		MPI_Barrier(MPI_COMM_WORLD);

		// get current profile suffstat logprobs
		GlobalProfileSuffStatLogProb();

		int proposenstep = nprofilerep;
		int totalnstep = nprofilerep;
		if (proposemode)	{
			totalnstep = 1;
		}

		for (int prep=0; prep<totalnstep; prep++)	{

			// backup profilesuffstatlogprobs
			for (int k=0; k<Ncomponent; k++)	{
				bkprofilesuffstatlogprob[k] = profilesuffstatlogprob[k];
			}

			// propose move to all profiles
			for (int k=0; k<Ncomponent; k++)	{
				if (occupancy[k])	{
					profilelogratio[k] = - LogStatPrior(k) - profilesuffstatlogprob[k];
					for (int l=0; l<GetDim(); l++)	{
						bkprofile[k][l] = profile[k][l];
					}
					if (proposemode)	{
						ApproxSuffStatProfileProposeMove(k,proposenstep);
					}
					else	{
						double r = rnd::GetRandom().Uniform();
						if (r < 1.0 / 3)	{
							profilelogratio[k] += ProfileProposeMove(profile[k],1.0,1,0,k,0);
						}
						else if (r < 2.0 / 3)	{
							profilelogratio[k] += ProfileProposeMove(profile[k],0.3,3,0,k,0);
						}
						else	{
							profilelogratio[k] += ProfileProposeMove(profile[k],0.1,5,0,k,0);
						}
					}
				}
			}

			// send new profiles
			MPI_Bcast(allocprofile,Ncomponent*GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);

			// get profile suff stat log probs
			GlobalProfileSuffStatLogProb();

			// decide acceptance on a per-mode basis 
			int naccept = 0;
			for (int k=0; k<Ncomponent; k++)	{
				if (occupancy[k])	{
					profilelogratio[k] +=  LogStatPrior(k) + profilesuffstatlogprob[k];
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
	}

	UpdateMatrices();
	GlobalProfileSuffStatLogProb();

	delete[] bkprofilesuffstatlogprob;
	delete[] profilelogratio;
	delete[] bkprofile;
	delete[] bkallocprofile;

	// ch.Stop();
	// cerr << ch.GetTime() << '\t' << ch.GetTime() / K0 << '\n';
	
	return ((double) NAccepted) / GetNOccupiedComponent() / nrep;
}

void ZipGeneralPathSuffStatGTRSBDPProfileProcess::SlaveMixMove()	{
	if (proposemode)	{
		SBDPProfileProcess::SlaveMixMove();
	}
	else	{
		ZipSlaveMixMove();
	}
}

void ZipGeneralPathSuffStatGTRSBDPProfileProcess::ZipSlaveMixMove()	{

	int itmp[4];
	MPI_Bcast(itmp,4,MPI_INT,0,MPI_COMM_WORLD);
	int nrep = itmp[0];
	int nallocrep = itmp[1];
	int K0 = itmp[2];
	int nprofilerep = itmp[3];

	double* mLogSamplingArray = new double[Ncomponent];
	double* cumul = new double[Ncomponent];

	for (int rep=0; rep<nrep; rep++)	{

		// receive weights
		MPI_Bcast(weight,Ncomponent,MPI_DOUBLE,0,MPI_COMM_WORLD);

		double totp = 0;
		for (int mode = 0; mode<K0; mode++)	{
			totp += weight[mode];
		}

		double totq = 0;
		for (int mode=K0; mode<GetNmodeMax(); mode++)	{
			totq += weight[mode];
		}

		// realloc move

		int NAccepted = 0;

		for (int allocrep=0; allocrep<nallocrep; allocrep++)	{

			for (int site=GetSiteMin(); site<GetSiteMax(); site++)	{

				if (ActiveSite(site))	{

					int bk = alloc[site];

					double max = 0;
					// double mean = 0;
					for (int mode = 0; mode<K0; mode++)	{
						mLogSamplingArray[mode] = LogStatProb(site,mode);
						if ((!mode) || (max < mLogSamplingArray[mode]))	{
							max = mLogSamplingArray[mode];
						}
						// mean += mLogSamplingArray[mode];
					}
					// mean /= K0;

					double total = 0;
					for (int mode = 0; mode<K0; mode++)	{
						double p = weight[mode] * exp(mLogSamplingArray[mode] - max);
						total += p;
						cumul[mode] = total;
					}
					if (isnan(total))	{
						cerr << "nan in slave mix move\n";
						exit(1);
					}

					// double M = exp(mean- max);
					double M = 1;
					total += M * totq;
					double q = total * rnd::GetRandom().Uniform();
					int mode = 0;
					while ( (mode<K0) && (q > cumul[mode])) mode++;
					if (mode == K0)	{
						mode--;
						double r = (q - cumul[mode]) / M;
						while (r > 0)	{
							mode++;
							r -= weight[mode];
						}
					}

					// MH 
					double logratio = 0;
					if (mode >= K0)	{
						logratio += LogStatProb(site,mode) - max - log(M);
					}
					if (bk >= K0)	{
						logratio -= LogStatProb(site,bk) - max - log(M);
					}
					
					if (log(rnd::GetRandom().Uniform()) > logratio)	{
						mode = bk;
					}

					int Accepted = (mode != bk);
					if (Accepted)	{
						NAccepted ++;
					}
					alloc[site] = mode;
				}
			}
		}
		MPI_Send(alloc,GetNsite(),MPI_INT,0,TAG1,MPI_COMM_WORLD);

		MPI_Barrier(MPI_COMM_WORLD);

		// profile move

		// send profile suffstat logprobs
		MESSAGE signal;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
		SlaveProfileSuffStatLogProb();

		int proposenstep = nprofilerep;
		int totalnstep = nprofilerep;
		if (proposemode)	{
			totalnstep = 1;
		}

		for (int prep=0; prep<totalnstep; prep++)	{

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

	// send profile suffstat logprobs
	MESSAGE signal;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	SlaveProfileSuffStatLogProb();

	delete[] cumul;
	delete[] mLogSamplingArray;
}

void ZipGeneralPathSuffStatRASCATGTRSBDPGammaPhyloProcess::SlaveExecute(MESSAGE signal)	{

	switch(signal) {

	case PROFILESSLOGPROB:
		SlaveProfileSuffStatLogProb();
		break;

	case RR_MOVE:
		SlaveMoveRR();
		break;

	case MTRYMIX_MOVE:
		MultipleTrySlaveMixMove();
		break;

	default:
		GeneralPathSuffStatRASCATGTRSBDPGammaPhyloProcess::SlaveExecute(signal);
	}
}

