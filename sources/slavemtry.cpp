

void ZipGTRSBDPProfileProcess::MultipleTrySlaveMixMove()	{

	int itmp[5];
	MPI_Bcast(itmp,5,MPI_INT,0,MPI_COMM_WORLD);
	int nrep = itmp[0];
	int nallocrep = itmp[1];
	int K0 = itmp[2];
	int nprofilerep = itmp[3];
	int nmultryrep = itmp[4];

	double* mLogSamplingArray = new double[Ncomponent];
	double* cumul = new double[Ncomponent];

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
		MPI_Send(alloc,GetFinalNsite(),MPI_INT,0,TAG1,MPI_COMM_WORLD);

		// useful ? could just as well send the nocc
		// receive full set of allocs
		MPI_Bcast(alloc,GetFinalNsite(),MPI_INT,0,MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);

		// profile move

		UpdateOccupancyNumbers();

		int nocc = GetNmodeMax()-1;
		while (!occupancy[nocc])	{
			nocc--;
		}
		nocc++;

		int ntry = nmultryrep * (GetNprocs() - 1);

		double* allocnewprofile = new double[ntry * nocc * (GetDim() + 2)];
		double*** proposeprofile = new double**[ntry];
		double* tmp = allocnewprofile;
		for (int k=0; k<ntry; k++)	{
			proposeprofile[k] = new double*[nocc];
			for (int l=0; l<nocc; l++)	{
				proposeprofile[k][l] = tmp;
				tmp += GetDim() + 2;
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
		
		// send profile suffstat logprobs
		MESSAGE signal;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
		SlaveProfileSuffStatLogProb();

		for (int prep=0; prep<nprofilerep; prep++)	{

			// propose new profiles

			// send profiles to master
			MPI_Send(allocnewprofile, nocc*nmultryrep*(GetDim() + 2),MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);

			// receive all profiles from master
			MPI_Bcast(allocnewprofile,nocc*ntry*(GetDim()+2),MPI_DOUBLE,0,MPI_COMM_WORLD);

			// calculate suff stat log probs
			// some work here
			// UpdateMatrices();

			// send suffstatlogprobs to master
			MPI_Reduce(tmpallocsuffstatlogprob,allocfwdsuffstatlogprob,nocc*ntry,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

			// get array of chosen profiles from master
			MPI_Bcast(allocprofile,nocc*GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);

			// propose new profiles
			// careful about exception for first series of first slave

			// send profiles to master
			MPI_Send(allocnewprofile, nocc*nmultryrep*(GetDim() + 2),MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);

			// receive all profiles from master
			MPI_Bcast(allocnewprofile,nocc*ntry*(GetDim()+2),MPI_DOUBLE,0,MPI_COMM_WORLD);

			// calculate suff stat log probs
			// some work here
			// UpdateMatrices();

			// send suffstatlogprobs to master
			MPI_Reduce(tmpallocsuffstatlogprob,allocfwdsuffstatlogprob,nocc*ntry,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

			// get array of chosen profiles from master
			MPI_Bcast(allocprofile,nocc*GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);
		}

		// rereceive all profiles
		MPI_Bcast(allocprofile,Ncomponent*GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);
		UpdateMatrices();

		// delete arrays
		for (int k=0; k<ntry; k++)	{
			delete[] proposeprofile[k];
		}
		delete[] proposeprofile;
		delete[] allocnewprofile;

		for (int k=0; k<nocc; k++)	{
			delete[] fwdsuffstatlogprob[k];
			delete[] revsuffstatlogprob[k];
		}
		delete[] fwdsuffstatlogprob;
		delete[] revsuffstatlogprob;

		delete[] tmpallocsuffstatlogprob;
		delete[] allocfwdsuffstatlogprob;
		delete[] allocrevsuffstatlogprob;

	}

	// send profile suffstat logprobs
	MESSAGE signal;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	SlaveProfileSuffStatLogProb();

	delete[] cumul;
	delete[] mLogSamplingArray;
}
