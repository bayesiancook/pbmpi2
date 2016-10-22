
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "MixtureProfileProcess.h"
#include "Random.h"
#include "Parallel.h"


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* MixtureProfileProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void MixtureProfileProcess::Create()	{
	if (! profile)	{
		if (! GetDim())	{
			cerr << "error: ProfileProcess::dim has not been initialized\n";
			exit(1);
		}
		ProfileProcess::Create();
		allocprofile = new double[GetNmodeMax() * GetDim()];
		profile = new double*[GetNmodeMax()];
		for (int i=0; i<GetNmodeMax(); i++)	{
			profile[i] = allocprofile + i*GetDim();
		}
		alloc = new int[GetNsite()];
		for (int i=0; i<GetNsite(); i++)	{
			if (ActiveSite(i))	{
				alloc[i] = 0;
			}
			else	{
				alloc[i] = -1;
			}
		}
		occupancy = new int[GetNmodeMax()];
		logstatprior = new double[GetNmodeMax()];
		profilesuffstatlogprob = new double[GetNmodeMax()];
	}
}

void MixtureProfileProcess::ActivateSumOverComponents()	{

	cerr <<  "allocate mtryalloc : " << mtryalloc << '\n';
	if (!mtryalloc)	{
		mtryalloc = new int*[GetNsite()];
		for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
			mtryalloc[i] = new int[sumovercomponents];
		}
		mtryweight = new double*[GetNsite()];
		for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
			mtryweight[i] = new double[sumovercomponents];
		}
	}
}

void MixtureProfileProcess::GlobalActivateSumOverComponents()	{

	if (GetNprocs() > 1)	{
		MESSAGE signal = ACTIVATEMTRY;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(&sumovercomponents,1,MPI_INT,0,MPI_COMM_WORLD);
	}
	else	{
		if (sumovercomponents > 0)	{
			ActivateSumOverComponents();
		}
	}
}

void MixtureProfileProcess::SlaveActivateSumOverComponents()	{

	MPI_Bcast(&sumovercomponents,1,MPI_INT,0,MPI_COMM_WORLD);
	if (sumovercomponents > 0)	{
		ActivateSumOverComponents();
	}
}

void MixtureProfileProcess::Delete()	{
	if (profile)	{
		if (mtryalloc)	{
			for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
				delete[] mtryalloc[i];
				delete[] mtryweight[i];
			}
			delete[] mtryalloc;
			delete[] mtryweight;
			mtryalloc = 0;
		}
		delete[] profilesuffstatlogprob;
		delete[] logstatprior;
		delete[] allocprofile;
		delete[] profile;
		delete[] alloc;
		delete[] occupancy;
		profile = 0;
		ProfileProcess::Delete();
	}
}

void MixtureProfileProcess::ChooseMultipleTryAlloc()	{

	int n = sumovercomponents;
	if (n > GetNcomponent())	{
		n = GetNcomponent();
	}

	double* probarray = new double[GetNcomponent()];
	double* cumul = new double[GetNcomponent()];
	double meandiv = 0;
	int nsite = 0;
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i))	{
			double tot = 0;
			for (int k=0; k<GetNcomponent(); k++)	{
				double tmp = GetWeight(k) * GetMinStat(profile[k],i);
				tot += tmp;
				probarray[k] = tmp;
				cumul[k] = tot;
			}
			for (int k=0; k<GetNcomponent(); k++)	{
				probarray[k] /= tot;
				cumul[k] /= tot;
			}
			// rnd::GetRandom().DrawFromUrn(mtryalloc[i],n,GetNcomponent(),probarray,mtryweight[i]);
			double div = 0;
			int jmin = 0;
			if (alloc[i] != -1)	{
				jmin = 1;
				mtryalloc[i][0] = alloc[i];
				mtryweight[i][0] = probarray[alloc[i]];
				div ++;
			}
			for (int j=jmin; j<n; j++)	{
				double u = rnd::GetRandom().Uniform();
				int k = 0;
				while ((k<GetNcomponent()) && (u > cumul[k]))	{
					k++;
				}
				if (k == GetNcomponent())	{
					cerr << "error in MixtureProfileProcess::ChooseMultipleTryAlloc: overflow\n";
					exit(1);
				}
				mtryalloc[i][j] = k;
				mtryweight[i][j] = probarray[k];

				int found = 0;
				for (int l=0; l<j; l++)	{
					if (mtryalloc[i][l] == k)	{
						found = 1;
					}
				}
				if (!found)	{
					div++;
				}
			}
			meandiv += div;
			nsite ++;
		}
	}

	meandiv /= nsite;

	// cerr << meandiv << '\n';

	delete[] probarray;
	delete[] cumul;
}

void MixtureProfileProcess::GlobalChooseMultipleTryAlloc()	{

	if (GetNprocs() > 1)	{
		MESSAGE signal = MTRYALLOC;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(&sumovercomponents,1,MPI_INT,0,MPI_COMM_WORLD);
	}
	else	{
		ChooseMultipleTryAlloc();
	}
}

void MixtureProfileProcess::SlaveChooseMultipleTryAlloc()	{

	MPI_Bcast(&sumovercomponents,1,MPI_INT,0,MPI_COMM_WORLD);
	ChooseMultipleTryAlloc();
}

double MixtureProfileProcess::GetStatEnt()	{
	double total = 0;
	UpdateOccupancyNumbers();
	int totnsite = 0;
	for (int k=0; k<GetNcomponent(); k++)	{
		total += occupancy[k] * GetStatEnt(k);
		totnsite += occupancy[k];
	}
	if (! totnsite)	{
		return log((double) GetDim());
	}
	return total / totnsite;
}

double MixtureProfileProcess::GetStatEnt(int k)	{
	double total = 0;
	for (int i=0; i<GetDim(); i++)	{
		if (profile[k][i] <= 0)	{
			cerr << "error: 0 entry in profile\n";
			cerr << profile[k][i] << '\n';
			exit(1);
		}
		total -= profile[k][i] * log(profile[k][i]);
	}
	if (isnan(total))	{
		cerr << "entropy is nan\n";
		exit(1);
	}
	return  total;
}

void MixtureProfileProcess::RenormalizeProfiles()	{
	for (int i=0; i<GetNcomponent(); i++)	{
		double total = 0;
		for (int k=0; k<GetDim(); k++)	{
			total += profile[i][k];
		}
		for (int k=0; k<GetDim(); k++)	{
			profile[i][k] /= total;
		}
	}
}

void MixtureProfileProcess::PriorSampleProfile()	{
	PriorSampleHyper();
	PriorSampleGlobalParameters();
	SampleAlloc();
	SampleStat();
}

void MixtureProfileProcess::SampleProfile()	{
	SampleHyper();
	SampleGlobalParameters();
	SampleAlloc();
	SampleStat();
}

void MixtureProfileProcess::SampleStat()	{
	for (int i=0; i<GetNcomponent(); i++)	{
		SampleStat(i);
	}
}

void MixtureProfileProcess::SampleStat(int i)	{
	SampleFrequencyStat(profile[i]);
}

double MixtureProfileProcess::ResampleEmptyProfiles()	{

	UpdateOccupancyNumbers();
	for (int i=0; i<GetNcomponent(); i++)	{
		if (! occupancy[i])	{
			SampleStat(i);
		}
	}
	return 1.0;
}

double MixtureProfileProcess::LogProfilePrior()	{
	double total = 0;
	total += LogHyperPrior();
	total += LogAllocPrior();
	total += LogStatPrior();
	return total;
}

void MixtureProfileProcess::UpdateOccupancyNumbers()	{
	for (int i=0; i<GetNcomponent(); i++)	{
		occupancy[i] = 0;
	}
	for (int i=0; i<GetNsite(); i++)	{
		if (ActiveSite(i))	{
			if (alloc[i] == -1)	{
				cerr << "error in MixtureProfileProcess::UpdateOccupancyNumber: alloc[i] == -1\n";
				exit(1);
			}
			occupancy[alloc[i]]++;
		}
	}
}

double MixtureProfileProcess::LogStatPrior()	{

	for (int i=0; i<GetNcomponent(); i++)	{
		LogStatPrior(i);
	}
	double total = 0;
	for (int i=0; i<GetNcomponent(); i++)	{
		total += logstatprior[i];
	}
	return total;
}

double MixtureProfileProcess::LogStatPrior(int cat)	{
	double tmp = LogFrequencyStatPrior(profile[cat]);
	logstatprior[cat] = tmp;
	return tmp;
}

double MixtureProfileProcess::ProfileSuffStatLogProb()	{

	// simply, sum over all components
	for (int i=0; i<GetNcomponent(); i++)	{
		ProfileSuffStatLogProb(i);
	}
	double total = 0;
	for (int i=0; i<GetNcomponent(); i++)	{
		total += profilesuffstatlogprob[i];
	}
	return total;
}

double MixtureProfileProcess::CountProfileSuffStatLogProb()	{

	// simply, sum over all components
	double total = 0;
	for (int i=0; i<GetNcomponent(); i++)	{
		total += CountProfileSuffStatLogProb(i);
	}
	return total;
}

double MixtureProfileProcess::BetaProfileSuffStatLogProb()	{

	// simply, sum over all components
	double total = 0;
	for (int i=0; i<GetNcomponent(); i++)	{
		total += BetaProfileSuffStatLogProb(i);
	}
	return total;
}

void MixtureProfileProcess::SwapComponents(int cat1, int cat2)	{

	int tmp = occupancy[cat1];
	occupancy[cat1] = occupancy[cat2];
	occupancy[cat2] = tmp;

	for (int k=0; k<GetDim(); k++)	{
		double tmp = profile[cat1][k];
		profile[cat1][k] = profile[cat2][k];
		profile[cat2][k] = tmp;
	}

	for (int i=0; i<GetNsite(); i++)	{
		if (ActiveSite(i))	{
			if (alloc[i] == -1)	{
				cerr << "error in MixtureProfileProcess::SwapComponents: alloc[i] == -1\n";
				exit(1);
			}
			if (alloc[i] == cat1)	{
				alloc[i] = cat2;
			}
			else if (alloc[i] == cat2)	{
				alloc[i] = cat1;
			}
		}
	}

	UpdateComponent(cat1);
	UpdateComponent(cat2);
}


double MixtureProfileProcess::GlobalMoveProfile(double tuning, int n, int nrep)	{

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

	int Nocc = GetNOccupiedComponent();

	double* dtmp = new double[1+Nocc*GetDim()];
	dtmp[0] = tuning;
	int k = 1;
	for (int i=0; i<Ncomponent; i++)	{
		if (occupancy[i])	{
			for (int j=0; j<GetDim(); j++)	{
				dtmp[k] = profile[i][j];
				k++;
			}
		}
	}
	MPI_Bcast(dtmp,1+Nocc*GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	delete[] dtmp;

	// split Ncomponent items among GetNprocs() - 1 slaves
	int width = Nocc/(GetNprocs()-1);
	int maxwidth = 0;
	int cmin[GetNprocs()-1];
	int cmax[GetNprocs()-1];
	int dmin[GetNprocs()-1];
	int dmax[GetNprocs()-1];

	for(int i=0; i<GetNprocs()-1; ++i) {

		int ddmin = width * i;
		int ddmax = (i == GetNprocs() - 2) ? Nocc : width * (i+1);
		if (maxwidth < (ddmax - ddmin))	{
			maxwidth = ddmax - ddmin;
		}

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
	MPI_Status stat;
	int bigdim = maxwidth * GetDim();
	double* tmp = new double[bigdim+1]; // (+1 for the acceptance rate)
	double total = 0;
	for(int i=1; i<GetNprocs(); ++i) {
		MPI_Recv(tmp,(dmax[i-1]-dmin[i-1])*GetDim()+1,MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
		int l = 0;
		for(int j=cmin[i-1]; j<cmax[i-1]; ++j) {
			if (occupancy[j])	{
				for (int k=0; k<GetDim(); k++)	{
					profile[j][k] = tmp[l];
					l++;
				}
			}
		}
		total += tmp[l]; // (sum all acceptance rates)
	}
	delete[] tmp;
	// return average acceptance rate
	return total / GetNOccupiedComponent();
	// return total / GetNOccupiedComponent();
}

void MixtureProfileProcess::SlaveMoveProfile()	{

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

	UpdateOccupancyNumbers();
	int Nocc = GetNOccupiedComponent();

	double* dtmp = new double[1 + Nocc*GetDim()];
	MPI_Bcast(dtmp,1+Nocc*GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	double tuning = dtmp[0];
	int k = 1;
	for (int i=0; i<Ncomponent; i++)	{
		if (occupancy[i])	{
			for (int j=0; j<GetDim(); j++)	{
				profile[i][j] = dtmp[k];
				k++;
			}
		}
	}
	delete[] dtmp;

	// determine the range of components to move
	int width = Nocc/(GetNprocs()-1);
	int dmin = width * (GetMyid() - 1);
	int dmax = (GetMyid() == GetNprocs() - 1) ? Nocc : width * GetMyid();

	k = -1;
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
	// UpdateModeProfileSuffStat();
	UpdateOccupancyNumbers();
	// move components in the range just computed
	double total = 0;
	for (int i=cmin; i<cmax; i++)	{
		if (occupancy[i])	{
			total += MoveProfile(i,tuning,n,nrep);
		}
	}

	// send the new values of the profiles, plus the total success rate (total)
	double* tmp = new double[(dmax - dmin) * GetDim() + 1];
	int l = 0;
	for (int i=cmin; i<cmax; i++)	{
		if (occupancy[i])	{
			for (int k=0; k<GetDim(); k++)	{
				tmp[l] = profile[i][k];
				l++;
			}
		}
	}
	tmp[l] = total;
	MPI_Send(tmp,(dmax-dmin)*GetDim()+1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
	delete[] tmp;
}

double MixtureProfileProcess::MoveProfile(double tuning, int n, int nrep)	{
	double total = 0;
	for (int i=0; i<GetNcomponent(); i++)	{
		total += MoveProfile(i,tuning,n,nrep);
	}
	return total / GetNcomponent();
}

double MixtureProfileProcess::MoveProfile(int cat, double tuning, int n, int nrep)	{

	double naccepted = 0;
	double* bk = new double[GetDim()];
	for (int k=0; k<GetDim(); k++)	{
		bk[k] = profile[cat][k];
	} 
	for (int rep=0; rep<nrep; rep++)	{
		double deltalogprob = - LogStatPrior(cat) - ProfileSuffStatLogProb(cat);
		double loghastings = ProfileProposeMove(profile[cat],tuning,n,0,cat,0);
		// double loghastings = ProfileProposeMove(profile[cat],tuning,n,0,g);
		UpdateComponent(cat);
		deltalogprob += LogStatPrior(cat) + ProfileSuffStatLogProb(cat);
		deltalogprob += loghastings;
		// cerr << "obtained : " << deltalogprob << '\t' << deltalogprob - loghastings << '\t' << loghastings << '\n';
		// int accepted = (g.RandU01() < exp(deltalogprob));
		int accepted = (log(rnd::GetRandom().Uniform()) < deltalogprob);
		if (accepted)	{
			naccepted ++;
			for (int k=0; k<GetDim(); k++)	{
				bk[k] = profile[cat][k];
			} 
		}
		else	{
			for (int k=0; k<GetDim(); k++)	{
				profile[cat][k] = bk[k];
			} 
			UpdateComponent(cat);
		}
	}
	delete[] bk;
	return naccepted / nrep;
}

double MixtureProfileProcess::GlobalSMCAddSites()	{

	double ret = ProfileProcess::GlobalSMCAddSites();

	// receive new site allocations from slave
	MPI_Status stat;
	int tmpalloc[GetNsite()];
	for(int i=1; i<GetNprocs(); ++i) {
		MPI_Recv(tmpalloc,GetNsite(),MPI_INT,i,TAG1,MPI_COMM_WORLD,&stat);
		for(int j=GetProcSiteMin(i); j<GetProcSiteMax(i); ++j) {
			if (ActiveSite(i))	{
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
	return ret;
}
