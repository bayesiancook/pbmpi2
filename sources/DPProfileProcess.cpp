
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "DPProfileProcess.h"
#include "Random.h"
#include "Parallel.h"


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* DPProfileProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

const double meandir[] = {0.499737,0.171262,0.183399,0.225593,0.197453,0.211819,0.173191,0.175454,0.3181,0.240008,0.187577,0.324778,0.205587,0.395097,0.162356,0.519427,0.526213,0.349177,0.0511527,0.130222};

void DPProfileProcess::SampleHyper()	{
	DirichletProfileProcess::SampleHyper();
	// kappa = 100;
	kappa = 1;
}

double DPProfileProcess::LogProxy(int site, int cat)	{
	return 0;
}

double DPProfileProcess::Move(double tuning, int n, int nrep)	{

	totchrono.Start();

	for (int rep=0; rep<nrep; rep++)	{

		// relative rates
		GlobalParametersMove();
		// for (int rep=0; rep<5; rep++)	{
		// allocations
		GlobalUpdateParameters();
		GlobalUpdateSiteProfileSuffStat();
		UpdateModeProfileSuffStat();
		// no parallel version of the allocation move for the moment
		incchrono.Start();
		IncrementalDPMove(5,1);
		incchrono.Stop();
		// profiles
		GlobalUpdateParameters();
		GlobalUpdateSiteProfileSuffStat();
		// UpdateModeProfileSuffStat();
		profilechrono.Start();
		GlobalMoveProfile(1,1,100);
		GlobalMoveProfile(1,3,100);
		GlobalMoveProfile(0.1,3,100);
		profilechrono.Stop();
		MoveHyper(tuning,10);
	}
	totchrono.Stop();
	return 1;
}

void DPProfileProcess::SampleAlloc()	{

	CreateComponent(0);
	int first = 0;
	while ((first < GetNsite()) && (! ActiveSite(first)))	{
		first++;
	}
	if (first == GetNsite())	{
		cerr << "error in DPProfileProcess::SampleAlloc: overflow\n";
		exit(1);
	}
	
	AddSite(first,0);
	Ncomponent = 1;
	
	for (int i=first+1; i<GetNsite(); i++)	{
		if (ActiveSite(i))	{
		double* p = new double[Ncomponent+1];
		double total = 0;
		double max = 0;
		for (int k=0; k<Ncomponent; k++)	{
			double tmp = log(occupancy[k]) * LogProxy(i,k);
			if ((!k) || (max < tmp))	{
				max = tmp;
			}
			p[k] = tmp;
		}
		p[Ncomponent] = log(kappa) + LogProxy(i,Ncomponent);
		if (max < p[Ncomponent])	{
			max = p[Ncomponent];
		}
		for (int k=0; k<=Ncomponent; k++)	{
			double tmp = exp(p[k] - max);
			total += tmp;
			p[k] = total;
		}
		double q = total * rnd::GetRandom().Uniform();
		int k = 0;
		while ((k<=Ncomponent) && (q > p[k])) k++;
		if (k == Ncomponent+1)	{
			cerr << "error in draw dp mode: overflow\n";
			exit(1);
		}
		if (k==Ncomponent)	{
			CreateComponent(k);
			Ncomponent++;
		}
		AddSite(i,k);
		delete[] p;
		}
	}
}

double DPProfileProcess::LogHyperPrior()	{
	double total = 0;
	if (kappaprior == 0)	{
		total = -kappa / 10.0;
	}
	else 	{
		total = -log(kappa);
		if ((kappa < 1e-4) || (kappa > 1e4))	{
			// total = InfProb;
			total -= 1.0 / 0;
		}
	}
	total += DirichletProfileProcess::LogHyperPrior();
	return total;
}

double DPProfileProcess::LogAllocPrior()	{
	double total = GetNOccupiedComponent() * log(kappa);
	for (int i=0; i<GetNsite(); i++)	{
		total -= log(kappa + i);
	}
	return total;
}

double DPProfileProcess::MoveHyper(double tuning, int nrep)	{
	double total = 0;
	total += MoveKappa(tuning,nrep);
	total += DirichletProfileProcess::MoveHyper(tuning,nrep);
	UpdateOccupancyNumbers();
	ResampleEmptyProfiles();
	return total;
}

double DPProfileProcess::MoveKappa(double tuning, int nrep)	{
	UpdateOccupancyNumbers();
	double naccepted = 0;
	for (int rep=0; rep<nrep; rep++)	{
		double deltalogprob = - LogHyperPrior() - LogAllocPrior();
		double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
		double e = exp(m);
		kappa *= e;
		deltalogprob += LogHyperPrior() + LogAllocPrior();
		deltalogprob += m;
		int accepted = (log(rnd::GetRandom().Uniform()) < deltalogprob);
		/*
		if (kappa < 1)	{
			accepted = 0;
		}
		*/
		if (accepted)	{
			naccepted++;
		}
		else	{
			kappa /= e;
		}
	}
	return naccepted / nrep;
}


double DPProfileProcess::IncrementalDPMove(int nrep, double epsilon)	{

	cerr << "in DPProfileProcess::IncrementalDPMove: should check active sites\n";
	exit(1);

	if (GetMyid())	{
		cerr << "error: slave in DPProfileProcess::IncrementalDPMove\n";
		exit(1);
	}
	// UpdateOccupancyNumbers();
	int NAccepted = 0;
	int Nrep = (GetNsite() * nrep )/ 10;

	for (int rep=0; rep<Nrep; rep++)	{

		int site = (int) (GetNsite() * rnd::GetRandom().Uniform());

	// for (int site=0; site<GetNsite(); site++)	{
		int bk = alloc[site];
		int k = occupancy[alloc[site]] > 1 ? Ncomponent : Ncomponent-1;
		int h = k + Nadd;

		double* mLogSamplingArray = new double[h];

		// draw a new matrix for Nmode <= i < h
		for (int i=Ncomponent; i<h ; i++)	{
			CreateComponent(i);
			// already called in CreateComponent
			// SampleStat(i);
		}

		RemoveSite(site,bk);

		// Gibbs

		double max = 0;
		for (int mode = 0; mode < h; mode++)	{
			mLogSamplingArray[mode] =  LogStatProb(site,mode);
			if ((!mode) || (max < mLogSamplingArray[mode]))	{
				max = mLogSamplingArray[mode];
			}
		}

		double* cumul = new double[h];
		double total = 0;

		for (int mode = 0; mode < h; mode++)	{

			double p = 0;			
			if (mode < Ncomponent)	{			
				if (occupancy[mode])	{
					p = (double) occupancy[mode];
				}
				else	{
					p = kappa / Nadd;
				}
				p *= exp(mLogSamplingArray[mode] - max);
			}
			else	{
				p = (kappa / Nadd) * exp(mLogSamplingArray[mode] - max);		
			}

			total += p;
			cumul[mode] = total;
		}

		double q = total * rnd::GetRandom().Uniform();
		int mode = 0;
		while ( (mode<h) && (q > cumul[mode])) mode++;
		if (mode == h)	{
			cerr << "error in switch mode: gibbs overflow\n";
			exit(1);
		}
		delete[] cumul;

		int Accepted = (mode != bk);
		if (Accepted)	{
			NAccepted ++;
		}
		AddSite(site,mode);

		if (mode >= Ncomponent)	{			// if it's a new one
			if (mode > Ncomponent)	{
				SwapComponents(mode, Ncomponent);
				mode = Ncomponent;
			}
			Ncomponent++;
		}
		if (! occupancy[bk])	{
			if (bk!= Ncomponent-1)	{
				SwapComponents(bk, Ncomponent-1);
			}
			Ncomponent--;
		}

		for (int k=Ncomponent; k<h; k++)	{
			DeleteComponent(k);
		}

		delete[] mLogSamplingArray;
	}
	
	return ((double) NAccepted) / GetNsite() / nrep * 10;
}

