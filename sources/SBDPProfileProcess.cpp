
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "SBDPProfileProcess.h"
#include "Random.h"
#include "Parallel.h"



//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* SBDPProfileProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------


void SBDPProfileProcess::Create()	{

	if (! V)	{
		DPProfileProcess::Create();
		V = new double[GetNmodeMax()];
		weight = new double[GetNmodeMax()];
	}
}

void SBDPProfileProcess::Delete()	{

	if (V)	{
		delete[] V;
		delete[] weight;
		DPProfileProcess::Delete();
	}
}

double SBDPProfileProcess::Move(double tuning, int nmix, int nrep, int nalloc)	{

	double ret = 0;
	if (GetNprocs() > 1)	{
		ret = MPIMove(tuning,nmix,nrep,nalloc);
	}
	else	{
		ret = NonMPIMove(tuning,nmix,nrep,nalloc);
	}
	return ret;
}

double SBDPProfileProcess::MPIMove(double tuning, int nmix, int nrep, int nalloc)	{

	totchrono.Start();

	for (int rep=0; rep<nrep; rep++)	{

		// relative rates
		GlobalParametersMove();

		incchrono.Start();
		GlobalUpdateParameters();
		GlobalUpdateSiteProfileSuffStat();
		GlobalUpdateModeProfileSuffStat();

		if ((!rep) && InitIncremental)	{
			cerr << "init incremental\n";
			InitIncremental--;
			IncrementalSampleAlloc();
			GlobalUpdateModeProfileSuffStat();
		}

		if (proposemode)	{
			if (allocmode)	{
				profacc += GlobalMixMove(nmix,nalloc,0.1,40);
			}
			else	{
				profacc += GlobalMixMove(nmix,nalloc,0.001,40);
			}
		}
		else	{
			if (allocmode)	{
				profacc += GlobalMixMove(nmix,nalloc,0.1,40);
			}
			else	{
				profacc += GlobalMixMove(nmix,nalloc,0.001,40);
			}
		}

		proftry ++;
		MoveOccupiedCompAlloc(5);
		MoveAdjacentCompAlloc(5);
		incchrono.Stop();

		// hyperparameters
		// globalupdates useless as long as not relying on sufficient statistics
		// themselves dependent on new parameter values
		GlobalUpdateParameters();
		GlobalUpdateSiteProfileSuffStat();
		MoveHyper(tuning,10);
	}
	totchrono.Stop();
	return 1;
}

double SBDPProfileProcess::NonMPIMove(double tuning, int nmix, int nrep, int nalloc)	{

	totchrono.Start();

	for (int rep=0; rep<nrep; rep++)	{

		GlobalParametersMove();

		incchrono.Start();
		UpdateSiteProfileSuffStat();
		UpdateModeProfileSuffStat();

		if ((!rep) && InitIncremental)	{
			cerr << "init incremental\n";
			InitIncremental--;
			IncrementalSampleAlloc();
			UpdateModeProfileSuffStat();
		}

		profacc += MixMove(nmix,nalloc,0.001,40);
		proftry ++;

		MoveOccupiedCompAlloc(5);
		MoveAdjacentCompAlloc(5);
		incchrono.Stop();

		UpdateSiteProfileSuffStat();
		MoveHyper(tuning,10);
	}
	totchrono.Stop();
	return 1;
}

void SBDPProfileProcess::SampleAlloc()	{

	for (int k=0; k<GetNmodeMax(); k++)	{
		CreateComponent(k);
	}
	Ncomponent = GetNmodeMax();
	SampleWeights();

	// SampleWeights();
	for (int i=0; i<GetNsite(); i++)	{
		if (ActiveSite(i))	{
			double U = rnd::GetRandom().Uniform();
			double total = weight[0];
			int k = 0;
			while ((k<GetNmodeMax()) && (total < U))	{
				k++;
				total += weight[k];
			}
			if (k == GetNmodeMax())	{
				cerr << "error in SBDPProfileProcess::SampleAlloc: overflow\n";
				exit(1);
			}
			AddSite(i,k);
		}
	}
}

void SBDPProfileProcess::IncrementalSampleAlloc()	{

	cerr << "in SBDPProfileProces::IncrementalSampleAlloc: check GetNsite()\n";
	exit(1);

	kappa = 0.1;

	for (int i=0; i<GetNsite(); i++)	{
		RemoveSite(i,alloc[i]);
	}

	AddSite(0,0);
	Ncomponent = 1;
	
	for (int i=0; i<GetNsite(); i++)	{

		int K = Ncomponent + 1;
		if (K > GetNmodeMax())	{
			K--;
		}
		double* p = new double[K];
		double total = 0;
		double max = 0;
		for (int k=0; k<K; k++)	{
			double w = occupancy[k];
			if (! w)	{
				w = kappa;
			}
			double tmp = log(w) * LogProxy(i,k);
			if ((!k) || (max < tmp))	{
				max = tmp;
			}
			p[k] = tmp;
		}
		for (int k=0; k<K; k++)	{
			double tmp = exp(p[k] - max);
			total += tmp;
			p[k] = total;
		}
		double q = total * rnd::GetRandom().Uniform();
		int k = 0;
		while ((k<K) && (q > p[k])) k++;
		if (k == K)	{
			cerr << "error in draw dp mode: overflow\n";
			exit(1);
		}
		if (k==Ncomponent)	{
			if (Ncomponent <= GetNmodeMax())	{
				Ncomponent++;
			}
		}
		AddSite(i,k);
		delete[] p;
	}

	Ncomponent = GetNmodeMax();
	ResampleWeights();
	cerr << "init incremental ok\n";
}

void SBDPProfileProcess::SampleWeights()	{

	double cumulProduct = 1.0;
	double totweight = 0;
	double v, x, y;
	for (int k=0; k<GetNcomponent(); k++)	{
		x = rnd::GetRandom().sGamma(1.0);
		y = rnd::GetRandom().sGamma(kappa);
		v = x / (x+y);
		V[k] = v;
		if (k == GetNcomponent() - 1)	{
			V[k] = 1;
			v = 1;
		}
		weight[k] = v * cumulProduct;
		cumulProduct *= (1 - v);	
		totweight += weight[k];
	}
}

void SBDPProfileProcess::ResampleWeights()	{

	UpdateOccupancyNumbers();

	int remainingOcc = GetNactiveSite();
	// int remainingOcc = GetNsite();
	double cumulProduct = 1.0;
	double totweight = 0;
	double v, x, y;
	for (int k=0; k<GetNcomponent(); k++)	{
		remainingOcc -= occupancy[k];
		x = rnd::GetRandom().sGamma(1 + occupancy[k]);
		y = rnd::GetRandom().sGamma(kappa + remainingOcc);
		v = x / (x+y);
		V[k] = v;
		if (k == GetNcomponent() - 1)	{
			double tmp = cumulProduct * (1 - v);
			if (maxweighterror < tmp)	{
				maxweighterror = tmp;
			}
			V[k] = 1;
			v = 1;
		}
		weight[k] = v * cumulProduct;
		cumulProduct *= (1 - v);
		totweight += weight[k];
	}
}

double SBDPProfileProcess::MoveOccupiedCompAlloc(int k0)	{

	int nrep = (int) (k0 * kappa);
	UpdateOccupancyNumbers();
	ResampleWeights();
	double total = 0.0;
	int Nocc = GetNOccupiedComponent();
	if (Nocc > 1)	{
		for (int i=0; i<nrep; i++)	{
			int* occupiedComponentIndices = new int[Nocc];
			int j=0;
			for (int k=0; k<GetNcomponent(); k++)	{
				if (occupancy[k] != 0)	{
					occupiedComponentIndices[j] = k;
					j++;
				}
			}
			if (j != Nocc)	{
				cerr << "error in MoveOccupiedCompAlloc.\n";
				exit(1);
			}
			int* indices = new int[2];
			rnd::GetRandom().DrawFromUrn(indices,2,Nocc);
			int cat1 = occupiedComponentIndices[indices[0]];
			int cat2 = occupiedComponentIndices[indices[1]];
			double logMetropolis = (occupancy[cat2] - occupancy[cat1]) * log(weight[cat1] / weight[cat2]);
			int accepted = (log(rnd::GetRandom().Uniform()) < logMetropolis);
			if (accepted)	{
				total += 1.0;
				// MixtureProfileProcess::SwapComponents(cat1,cat2);
				SwapComponents(cat1, cat2);
				double tempv = V[cat1];
				V[cat1] = V[cat2];
				V[cat2] = tempv;
				double tempw = weight[cat1];
				weight[cat1] = weight[cat2];
				weight[cat2] = tempw;
			
			}
			delete[] occupiedComponentIndices;
			delete[] indices; 
		}
		return total /= nrep;
	}
	return 0;
}

double SBDPProfileProcess::LogIntegratedAllocProb()	{
	int remainingOcc = GetNactiveSite();
	double total = 0;
	for (int k=0; k<GetNcomponent(); k++)	{
		if (remainingOcc)	{
			remainingOcc -= occupancy[k];
			total += log(kappa) + rnd::GetRandom().logGamma(1 + occupancy[k]) + rnd::GetRandom().logGamma(kappa + remainingOcc) - rnd::GetRandom().logGamma(1 + kappa + occupancy[k] + remainingOcc);
		}
	}
	if (remainingOcc)	{
		cerr << "error in allocation count\n";
		cerr << GetNOccupiedComponent() << '\n';
		cerr << GetNactiveSite() << '\n';
		exit(1);
	}
	return total;
}

double SBDPProfileProcess::MoveKappa(double tuning, int nrep)	{
	UpdateOccupancyNumbers();
	double naccepted = 0;
	for (int rep=0; rep<nrep; rep++)	{
		double deltalogprob = - LogHyperPrior() - LogIntegratedAllocProb();
		double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
		double e = exp(m);
		kappa *= e;
		deltalogprob += LogHyperPrior() + LogIntegratedAllocProb();
		deltalogprob += m;
		int accepted = (log(rnd::GetRandom().Uniform()) < deltalogprob);
		if (accepted)	{
			naccepted++;
		}
		else	{
			kappa /= e;
		}
	}
	return naccepted / nrep;
}

double SBDPProfileProcess::LogStatPrior()	{

	UpdateOccupancyNumbers();
	double total = 0;
	for (int i=0; i<GetNcomponent(); i++)	{
		if (occupancy[i])	{
			total += DPProfileProcess::LogStatPrior(i);
		}
	}
	return total;
}

void SBDPProfileProcess::SwapComponents(int cat1, int cat2)	{

	MixtureProfileProcess::SwapComponents(cat1,cat2);
	double tempv = V[cat1];
	V[cat1] = V[cat2];
	V[cat2] = tempv;
	double tempw = weight[cat1];
	weight[cat1] = weight[cat2];
	weight[cat2] = tempw;
}



double SBDPProfileProcess::MoveAdjacentCompAlloc(int k0)	{

	int nrep = (int) (k0 * kappa);
	ResampleWeights();
	
	double total = 0;

	for (int i=0; i<nrep; i++)	{
		int cat1 = (int)(rnd::GetRandom().Uniform() * (GetNcomponent()-2));  
		int cat2 = cat1 + 1;
		double logMetropolis = (occupancy[cat1] * log(1 - V[cat2])) - (occupancy[cat2] * log(1-V[cat1]));
		int accepted = (log(rnd::GetRandom().Uniform()) < logMetropolis);
		if (accepted)	{
			total += 1.0;
			SwapComponents(cat1,cat2);
		}
	}

	return total /= nrep;
}


double SBDPProfileProcess::GlobalMixMove(int nrep, int nallocrep, double epsilon, int nprofilerep)	{

	if (Ncomponent != GetNmodeMax())	{
		cerr << "error in sbdp inc dp move: number of components\n";
		exit(1);
	}

	// define threshold between GIbbs and MH
	int K0 = 0;
	if (allocmode)	{
		K0 = 1.0 / epsilon;
		if (K0 < 10)	{
			K0 = 10;
		}
	}
	else	{
		K0 = GetNmodeMax();
		if (epsilon)	{
			double r = kappa / (1 + kappa);
			K0 = (int) (log(epsilon) / log(r));
			if (K0 >= GetNmodeMax())	{
				K0 = GetNmodeMax();
			}
		}
	}

	// send mixmove signal and tuning parameters
	MESSAGE signal = MIX_MOVE;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	int itmp[4];
	itmp[0] = nrep;
	itmp[1] = nallocrep;
	itmp[2] = K0;
	itmp[3] = nprofilerep;
	MPI_Bcast(itmp,4,MPI_INT,0,MPI_COMM_WORLD);

	double* tmp = new double[Ncomponent * GetDim() + 1];

	double totalacc = 0;
	double totaltry = 0;

	for (int rep=0; rep<nrep; rep++)	{

		ResampleWeights();

		// mpi send message for realloc move
		// mpi send profiles and weights
		MPI_Bcast(weight,Ncomponent,MPI_DOUBLE,0,MPI_COMM_WORLD);

		// here slaves do realloc moves

		// mpi receive new allocations
		MPI_Status stat;
		int tmpalloc[GetNsite()+1];
		for(int i=1; i<GetNprocs(); ++i) {
			MPI_Recv(tmpalloc,GetNsite(),MPI_INT,i,TAG1,MPI_COMM_WORLD,&stat);
			for(int j=GetProcSiteMin(i); j<GetProcSiteMax(i); ++j) {
				if (ActiveSite(j))	{
					alloc[j] = tmpalloc[j];
					if ((alloc[j] < 0) || (alloc[j] >= Ncomponent))	{
						cerr << "in global mix move\n";
						cerr << "alloc overflow\n";
						exit(1);
					}
				}
			}
		}

		MPI_Barrier(MPI_COMM_WORLD);

		// broadcast new allocations
		MPI_Bcast(alloc,GetNsite(),MPI_INT,0,MPI_COMM_WORLD);

		// update mode profile suff stats and send them to slaves
		GlobalUpdateModeProfileSuffStat();

		// here slaves do profile moves

		// split Ncomponent items among GetNprocs() - 1 slaves
		UpdateOccupancyNumbers();
		int Nocc = GetNOccupiedComponent();
		int width = Nocc/(GetNprocs()-1);
		int cmin[GetNprocs()-1];
		int cmax[GetNprocs()-1];
		int dmin[GetNprocs()-1];
		int dmax[GetNprocs()-1];

		for(int i=0; i<GetNprocs()-1; ++i) {
			int ddmin = width * i;
			int ddmax = (i == GetNprocs() - 2) ? Nocc : width * (i+1);

			int k = -1;
			int ccmin = -1;
			while ((ccmin<Ncomponent) && (k<ddmin))	{
				ccmin++;
				if (ccmin == Ncomponent)	{
					cerr << "error in matmixslavemoveprofile: overflow\n";
					exit(1);
				}
				if (occupancy[ccmin])	{
					k++;
				}
			}
			int ccmax = ccmin;
			if (ddmax == Nocc)	{
				ccmax = Ncomponent;
			}
			else	{
				while ((ccmax<Ncomponent) && (k<ddmax))	{
					ccmax++;
					if (occupancy[ccmax])	{
						k++;
					}
				}
			}

			cmin[i] = ccmin;
			cmax[i] = ccmax;
			dmin[i] = ddmin;
			dmax[i] = ddmax;
			int nocc = 0;
			for (int j=ccmin; j<ccmax; j++)	{
				if (occupancy[j])	{
					nocc++;
				}
			}
			if (nocc != (dmax[i] - dmin[i]))	{
				cerr << "error: non matching numbers: " << i << '\t' << nocc << '\t' << ccmin << '\t' << ccmax << '\t' << ddmin << '\t' << ddmax << '\n';
				for (int j=ccmin; j<ccmax; j++)	{
					cerr << occupancy[j] << '\t';
				}
				cerr << '\n';
				exit(1);
			}
		}

		// collect final values of profiles (+ total acceptance rate) from slaves
		MPI_Status stat2;
		for(int i=1; i<GetNprocs(); ++i) {
			MPI_Recv(tmp,(dmax[i-1]-dmin[i-1])*GetDim()+1,MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat2);
			int l = 0;
			for(int j=cmin[i-1]; j<cmax[i-1]; ++j) {
				if (occupancy[j])	{
					double tot = 0;
					for (int k=0; k<GetDim(); k++)	{
						profile[j][k] = tmp[l];
						tot += profile[j][k];
						l++;
					}
					if (fabs(tot - 1) > 1e-6)	{
						cerr << "normalization error : " << tot -1 << '\n';
						cerr << "upon receiving\n";
						exit(1);
					}
				}
			}
			totalacc += tmp[l]; // (sum all acceptance rates)
			totaltry ++;
		}

		// resample empty profiles
		ResampleEmptyProfiles();

		MPI_Barrier(MPI_COMM_WORLD);
		// resend all profiles
		MPI_Bcast(allocprofile,Ncomponent*GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	}

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
	delete[] tmp;

	UpdateComponents();

	return totalacc / totaltry;
}

void SBDPProfileProcess::SlaveMixMove()	{

	int itmp[4];
	MPI_Bcast(itmp,4,MPI_INT,0,MPI_COMM_WORLD);
	int nrep = itmp[0];
	int nallocrep = itmp[1];
	int K0 = itmp[2];
	int nprofilerep = itmp[3];

	double* mLogSamplingArray = new double[Ncomponent];
	double* cumul = new double[Ncomponent];
	double* tmp = new double[Ncomponent * GetDim() + 1];

	for (int rep=0; rep<nrep; rep++)	{

		// realloc move

		MPI_Bcast(weight,Ncomponent,MPI_DOUBLE,0,MPI_COMM_WORLD);
		// MPI_Bcast(allocprofile,Ncomponent*GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);

		double totp = 0;
		double totq = 0;
		if (allocmode)	{
			cerr << "check alloc mode in SBDPProfileProcess::SlaveMixMove\n";
			exit(1);
		}
		if (! allocmode)	{
			for (int mode = 0; mode<K0; mode++)	{
				totp += weight[mode];
			}
			for (int mode=K0; mode<GetNmodeMax(); mode++)	{
				totq += weight[mode];
			}
		}

		int NAccepted = 0;

		for (int allocrep=0; allocrep<nallocrep; allocrep++)	{

			for (int site=GetSiteMin(); site<GetSiteMax(); site++)	{

				if (ActiveSite(site))	{

					if (alloc[site] == -1)	{
						cerr << "error in SBDPProfileProcess::SlaveMixMove: alloc == -1\n";
						exit(1);
					}

					if (allocmode)	{

						int bk = alloc[site];

						double selectweight[Ncomponent];
						double totweight = 0;
						for (int mode = 0; mode<Ncomponent; mode++)	{
							selectweight[mode] = weight[mode] * GetMinStat(profile[mode],site);
							totweight += selectweight[mode];
						}

						int select[Ncomponent];
						for (int mode=0; mode<Ncomponent; mode++)	{
							selectweight[mode] *= K0 / totweight;
							if (selectweight[mode] > 1)	{
								selectweight[mode] = 1;
							}
							if ((mode == bk) || (rnd::GetRandom().Uniform() < selectweight[mode]))	{
								select[mode] = 1;
							}
							else	{
								select[mode] = 0;
							}
						}
						double max = 0;
						// double mean = 0;
						int found = 0;
						for (int mode = 0; mode<Ncomponent; mode++)	{
							if (select[mode])	{
								mLogSamplingArray[mode] = LogStatProb(site,mode);
								if ((!found) || (max < mLogSamplingArray[mode]))	{
									max = mLogSamplingArray[mode];
									found = 1;
								}
							}
						}

						double total = 0;
						for (int mode = 0; mode<Ncomponent; mode++)	{
							double p = 0;
							if (select[mode])	{
								p = weight[mode] * exp(mLogSamplingArray[mode] - max) / selectweight[mode];
							}
							total += p;
							cumul[mode] = total;
						}
						if (isnan(total))	{
							cerr << "in alloc: nan\n";
							exit(1);
						}
						if (isinf(total))	{
							cerr << "in alloc: inf\n";
							exit(1);
						}

						double q = total * rnd::GetRandom().Uniform();
						int mode = 0;
						while ((mode < Ncomponent) && (q > cumul[mode]))	{
							mode++;
						}
						if (mode == Ncomponent)	{
							cerr << "error in alloc move: allocmode activated\n";
							exit(1);
						}

						int Accepted = (mode != bk);
						if (Accepted)	{
							NAccepted ++;
						}
						alloc[site] = mode;
					}
					else	{
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
							cerr << "in SBDPProfileProcess: nan\n";
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
		}
		MPI_Send(alloc,GetNsite(),MPI_INT,0,TAG1,MPI_COMM_WORLD);

		MPI_Barrier(MPI_COMM_WORLD);
		// profile move

		// receive new allocations
		MPI_Bcast(alloc,GetNsite(),MPI_INT,0,MPI_COMM_WORLD);

		// determine the range of components to move
		UpdateOccupancyNumbers();
		int Nocc = GetNOccupiedComponent();
		int width = Nocc/(GetNprocs()-1);
		int dmin = width * (GetMyid() - 1);
		int dmax = (GetMyid() == GetNprocs() - 1) ? Nocc : width * GetMyid();

		int k = -1;
		int cmin = -1;
		while ((cmin<Ncomponent) && (k<dmin))	{
			cmin++;
			if (cmin == Ncomponent)	{
				cerr << "error in matmixslavemoveprofile: overflow\n";
				exit(1);
			}
			if (occupancy[cmin])	{
				k++;
			}
		}
		int cmax = cmin;
		if (dmax == Nocc)	{
			cmax = Ncomponent;
		}
		else	{
			while ((cmax<Ncomponent) && (k<dmax))	{
				cmax++;
				if (occupancy[cmax])	{
					k++;
				}
			}
		}
		int nocc = 0;
		for (int j=cmin; j<cmax; j++)	{
			if (occupancy[j])	{
				nocc++;
			}
		}
		if (nocc != (dmax - dmin))	{
			cerr << "error : mismatch in nocc\n";
			exit(1);
		}

		// update sufficient statistics
		MESSAGE signal;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
		SlaveUpdateModeProfileSuffStat();
		// UpdateModeProfileSuffStat();

		// move components in the range just computed
		double accept = 0;
		double attempt = 0;
		for (int i=cmin; i<cmax; i++)	{
			if (occupancy[i])	{
				accept += MoveProfile(i,1,1,nprofilerep);
				accept += MoveProfile(i,1,3,nprofilerep);
				accept += MoveProfile(i,0.1,3,nprofilerep);
				attempt += 3;
			}
		}

		// send the new values of the profiles, plus the total success rate (total)
		int l = 0;
		for (int i=cmin; i<cmax; i++)	{
			if (occupancy[i])	{
				double tot = 0;
				for (int k=0; k<GetDim(); k++)	{
					tmp[l] = profile[i][k];
					tot += tmp[l];
					l++;
				}
				if (fabs(tot - 1) > 1e-6)	{
					cerr << "normalization error : " << tot -1 << '\n';
					cerr << "upon sending\n";
					cerr << cmin << '\t' << cmax << '\n';
					exit(1);
				}
			}
		}
		tmp[l] = accept / attempt;
		
		MPI_Send(tmp,(dmax-dmin)*GetDim()+1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);

		// rereceive all profiles
		MPI_Bcast(allocprofile,Ncomponent*GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);
		UpdateComponents();
	}

	delete[] cumul;
	delete[] mLogSamplingArray;
	delete[] tmp;
}


double SBDPProfileProcess::MixMove(int nrep, int nallocrep, double epsilon, int nprofilerep)	{

	if (GetMyid())	{
		cerr << "error: slave in SBDPProfileProcess::MixMove\n";
		exit(1);
	}

	if (Ncomponent != GetNmodeMax())	{
		cerr << "error in sbdp inc dp move: number of components\n";
		exit(1);
	}

	int NAccepted = 0;
	double profileaccept = 0;
	double profileattempt = 0;

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

	double* mLogSamplingArray = new double[Ncomponent];
	double* cumul = new double[Ncomponent];

	for (int rep=0; rep<nrep; rep++)	{

		ResampleWeights();

		// realloc move

		double totp = 0;
		for (int mode = 0; mode<K0; mode++)	{
			totp += weight[mode];
		}

		double totq = 0;
		for (int mode=K0; mode<GetNmodeMax(); mode++)	{
			totq += weight[mode];
		}

		int NAccepted = 0;

		for (int allocrep=0; allocrep<nallocrep; allocrep++)	{

			for (int site=0; site<GetNsite(); site++)	{

				if (ActiveSite(site))	{

					if (alloc[site] == -1)	{
						cerr << "error in SBDPProfileProcess::MixMove: alloc == -1\n";
						exit(1);
					}

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
						cerr << "in SBDPProfileProcess: nan\n";
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

		// update sufficient statistics
		UpdateModeProfileSuffStat();
		UpdateOccupancyNumbers();

		for (int i=0; i<GetNcomponent(); i++)	{
			if (occupancy[i])	{
				profileaccept += MoveProfile(i,1,1,nprofilerep);
				profileaccept += MoveProfile(i,1,3,nprofilerep);
				profileaccept += MoveProfile(i,0.1,3,nprofilerep);
				profileattempt += 3;
			}
		}

		UpdateComponents();
		ResampleEmptyProfiles();

	}

	delete[] cumul;
	delete[] mLogSamplingArray;

	return profileaccept / profileattempt;
	// return ((double) NAccepted) / GetNsite() / nrep;
}

double SBDPProfileProcess::IncrementalDPMove(int nrep, double epsilon)	{

	cerr << "in SBDPProfileProcess::IncrementalDPMove\n";
	cerr << "check GetNsite\n";
	exit(1);

	if (Ncomponent != GetNmodeMax())	{
		cerr << "error in sbdp inc dp move: number of components\n";
		exit(1);
	}

	ResampleEmptyProfiles();

	int NAccepted = 0;
	int K0 = GetNmodeMax();
	if (epsilon)	{
		double r = kappa / (1 + kappa);
		K0 = (int) (log(epsilon) / log(r));
		if (K0 >= GetNmodeMax())	{
			K0 = GetNmodeMax();
		}
	}

	double* mLogSamplingArray = new double[GetNmodeMax()];
	double* cumul = new double[GetNmodeMax()];


	for (int rep=0; rep<nrep; rep++)	{

		ResampleWeights();

		double totp = 0;
		for (int mode = 0; mode<K0; mode++)	{
			totp += weight[mode];
		}

		double totq = 0;
		for (int mode=K0; mode<GetNmodeMax(); mode++)	{
			totq += weight[mode];
		}

		for (int site=0; site<GetNsite(); site++)	{

			if (alloc[site] == -1)	{
				cerr << "error in SBDPProfileProcess::IncrementalDPMove: alloc == -1\n";
				exit(1);
			}

			int bk = alloc[site];

			double max = 0;
			double mean = 0;
			for (int mode = 0; mode<K0; mode++)	{
				mLogSamplingArray[mode] =  LogStatProb(site,mode);
				if ((!mode) || (max < mLogSamplingArray[mode]))	{
					max = mLogSamplingArray[mode];
				}
				mean += mLogSamplingArray[mode];
			}
			mean /= K0;

			double total = 0;
			for (int mode = 0; mode<K0; mode++)	{
				double p = weight[mode] * exp(mLogSamplingArray[mode] - max);
				total += p;
				cumul[mode] = total;
			}

			double M = exp(mean- max);
			M = 1;
			total += M * totq;
			double q = total * rnd::GetRandom().Uniform();
			int mode = 0;
			while ( (mode<K0) && (q > cumul[mode])) mode++;
			if (mode == K0)	{
				mode--;
				double r = (q - cumul[mode]) / M;
				// cerr << kappa << '\t' << 1 - totp << '\t' << r << '\t' << cumulProduct << '\n';
				// cerr << '\n';
				while (r > 0)	{
					mode++;
					/*
					if (mode >= Ncomponent)	{
						if (mode >= GetNmodeMax())	{
							cerr << "nmode max overflow\n";
							exit(1);
						}
						ResampleLastWeight();
						CreateComponent(Ncomponent);
						Ncomponent++;
					}
					*/
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
	
	delete[] cumul;
	delete[] mLogSamplingArray;

	// UpdateModeProfileSuffStat();
	return ((double) NAccepted) / GetNsite() / nrep;
}

double SBDPProfileProcess::SMCAddSites()	{

	// assumes weights have been resampled
	// SMC weights still to be correctly computed

	// int K0 = GetNmodeMax();
	int K0 = GetLastOccupiedComponent() + 50;
	if (K0 >= GetNmodeMax())	{
		K0 = GetNmodeMax();
	}

	double logw = 0;

	double* logl = new double[K0];
	double* p = new double[K0];
	double* cumul = new double[K0];

	// for (int site = GetBKSiteMax(); site<GetSiteMax(); site++)	{
	for (int site = GetSiteMin(); site<GetSiteMax(); site++)	{

		if (NewlyActivated(site))	{

			double max = 0;
			int inf = 0;
			for (int k=0; k<K0; k++)	{
				AddSite(site,k);
				logl[k] = SiteLogLikelihood(site);
				if (isinf(logl[k]))	{
					inf = 1;
				}
				if ((!k) || (max < logl[k]))	{
					max = logl[k];
				}
				RemoveSite(site,k);
			}
			if (isinf(max))	{
				cerr << "error in smc add sites\n";
				for (int k=0; k<K0; k++)	{
					cerr << logl[k] << '\n';
				}
				exit(1);
			}

			double tot = 0;
			for (int k=0; k<K0; k++)	{
				double tmp = 0;
				if (! isinf(logl[k]))	{
					tmp = weight[k] * exp(logl[k] - max);
				}
				p[k] = tmp;
				tot += tmp;
				cumul[k] = tot;
			}
			logw += log(tot) + max;
			if (isinf(logw))	{
				cerr << "in SBDPProfileProcess::SMCAddSites: lowg is inf\n";
				exit(1);
			}
			if (isnan(logw))	{
				cerr << "in SBDPProfileProcess::SMCAddSites: lowg is nan\n";
				exit(1);
			}
			double u = tot * rnd::GetRandom().Uniform();
			int k = 0;
			while ((k<K0) && (u > cumul[k]))	{
				k++;
			}
			if (k == K0)	{
				cerr << "error in FiniteProfileProcess::SMCAddSites:overflow\n";
				cerr << tot << '\t' << logw << '\n';
				cerr << max << '\n';
				cerr << u << '\n';
				exit(1);
			}
			AddSite(site,k);
			SiteLogLikelihood(site);
			/*
			if (inf)	{
				cerr << "in SBCPProfileProcess::SMCAddSites: inf\n";
				cerr << "chosen : " << k << '\t' << logl[k] << '\n';
				cerr << tot << '\t' << logw << '\n';
				cerr << max << '\n';
				cerr << '\n';
				for (int k=0; k<K0; k++)	{
					cerr << k << '\t' << weight[k] << '\t' << logl[k] << '\t' << p[k] << '\t' << cumul[k] << '\t' << u << '\n';
				}
		
				exit(1);
			}
			*/
			SampleSiteMapping(site);
		}
	}

	delete[] logl;
	delete[] p;
	delete[] cumul;

	// cerr << logw << '\n';
	// send back allocations to master
	if (GetMyid())	{
		MPI_Send(&logw,1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
		MPI_Send(alloc,GetNsite(),MPI_INT,0,TAG1,MPI_COMM_WORLD);
	}

	return logw;
}


double SBDPProfileProcess::GlobalSMCAddSites()	{

	double ret = ProfileProcess::GlobalSMCAddSites();

	// receive new site allocations from slave
	MPI_Status stat;
	int tmpalloc[GetNsite()];
	for(int i=1; i<GetNprocs(); ++i) {
		MPI_Recv(tmpalloc,GetNsite(),MPI_INT,i,TAG1,MPI_COMM_WORLD,&stat);
		for(int j=GetProcSiteMin(i); j<GetProcSiteMax(i); ++j) {
			if (ActiveSite(j))	{
				alloc[j] = tmpalloc[j];
				if ((tmpalloc[j] < 0) || (tmpalloc[j] >= Ncomponent))	{
					cerr << "in SMC add\n";
					cerr << "alloc overflow\n";
					cerr << tmpalloc[j] << '\n';
					exit(1);
				}
			}
		}
	}
	UpdateOccupancyNumbers();
	ResampleWeights();
	return ret;
}
