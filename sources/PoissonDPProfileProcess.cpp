
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "PoissonDPProfileProcess.h"
#include "Random.h"


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* PoissonDPProfileProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void PoissonDPProfileProcess::ToStream(ostream& os)	{

	os << Ncomponent << '\n';
	os << kappa << '\n';
	for (int j=0; j<GetDim(); j++)	{
		os << dirweight[j] << '\t';
	}
	os << '\n';
	os << '\n';
	for (int i=0; i<Ncomponent; i++)	{
		for (int j=0; j<GetDim(); j++)	{
			os << profile[i][j] << '\t';
		}
		os << '\n';
	}

	for (int i=0; i<GetNsite(); i++)	{
		os << alloc[i] << '\t';
	}
	os << '\n';
}

void PoissonDPProfileProcess::FromStream(istream& is)	{

	is >> Ncomponent;
	is >> kappa;
	
	for (int i=0; i<GetDim(); i++)	{
		is >> dirweight[i];
	}

	for (int i=0; i<Ncomponent; i++)	{
		for (int j=0; j<GetDim(); j++)	{
			is >> profile[i][j];
		}
	}

	for (int i=0; i<GetNsite(); i++)	{
		is >> alloc[i];
	}

	// CHECK some update here ?
}

double PoissonDPProfileProcess::IncrementalDPMove(int nrep)	{

	cerr << "in PoissonDPProfileProcess::IncrementalDPMove: check active sites\n";
	exit(1);

	if (GetMyid())	{
		cerr << "error: slave in PoissonDPProfileProcess::IncrementalDPMove\n";
		exit(1);
	}

	UpdateOccupancyNumbers();
	int NAccepted = 0;
	int Nrep = (GetNsite() * nrep )/ 10;

	for (int rep=0; rep<Nrep; rep++)	{

		int site = (int) (GetNsite() * rnd::GetRandom().Uniform());

		int bk = alloc[site];
		int h = occupancy[bk] > 1 ? Ncomponent+1 : Ncomponent;

		// make a new mode Nmode <= i < h
		for (int i=Ncomponent; i<h ; i++)	{
			CreateComponent(i);
		}

		RemoveSite(site,bk);

		// Gibbs

		double total = 0;
		double* mModeGibbsGrid = new double[h];
		double* mLogSamplingArray = new double[h];

		double max = 0;
		for (int mode = 0; mode < h; mode++)	{
			// mLogSamplingArray[mode] =  LogStatProb(site,mode);
			mLogSamplingArray[mode] =  DiffLogSampling(mode,site);
			if ((!mode) || (max < mLogSamplingArray[mode]))	{
				max = mLogSamplingArray[mode];
			}
		}
		for (int mode = 0; mode < h; mode++)	{
			double temp=0;
			if (occupancy[mode])	{
				temp = occupancy[mode] * exp(mLogSamplingArray[mode] - max);
			}
			else	{
				temp = kappa * exp(mLogSamplingArray[mode] - max);
			}
			total += temp;
			mModeGibbsGrid[mode] = total;
		}

		double q = total * rnd::GetRandom().Uniform();
		int mode = 0;
		while ( (q > mModeGibbsGrid[mode]) && (mode < h)) mode++;
		if (mode == h)	{
			cerr << "error in switch mode integral\n";
			exit(1);
		}

		int Accepted = (mode != bk);

		if (Accepted)	{
			NAccepted++;
		}

		AddSite(site,mode);

		if (mode >= Ncomponent)	{
			if (mode > Ncomponent)	{
				SwapComponents(mode, Ncomponent);
			}
			Ncomponent++;
		}
		if (! occupancy[bk])	{
			if (bk != Ncomponent-1)	{
				SwapComponents(bk,Ncomponent-1);
			}
			Ncomponent--;
		}

		delete[] mModeGibbsGrid;
		delete[] mLogSamplingArray;
	}
	return ((double) NAccepted) / GetNsite() / nrep * 10;
}

/*

double PoissonDPProfileProcess::SMCAddSites()	{

	if (GetNprocs() > 1)	{
		cerr << "in DPProfileProcess::SMCAddSites: only in serial mode\n";
		exit(1);
	}

	double logw = 0;

	for (int site = GetBKSiteMax(); site<GetSiteMax(); site++)	{

		int h = Ncomponent+1;

		CreateComponent(Ncomponent);

		double logl[h];
		double cumul[h];
		double p[h];
		
		double max = 0;
		for (int k=0; k<h; k++)	{
			alloc[site] = k;
			// AddSite(site,k);
			logl[k] = SiteLogLikelihood(site);
			if ((!k) || (max < logl[k]))	{
				max = logl[k];
			}
			// RemoveSite(site,k);
		}

		double tot = 0;
		for (int k=0; k<K0; k++)	{
			double tmp = weight[k] * exp(logl[k] - max);
			p[k] = tmp;
			tot += tmp;
			cumul[k] = tmp;
		}
		logw += log(tot) + max;
		double u = tot * rnd::GetRandom().Uniform();
		int k = 0;
		while ((k<K0) && (cumul[k] > u))	{
			k++;
		}
		if (k == K0)	{
			cerr << "error in FiniteProfileProcess::SMCAddSites:overflow\n";
			exit(1);
		}
		AddSite(site,k);
		SiteLogLikelihood(site);
		SampleSiteMapping(site);
	}

	delete[] logl;
	delete[] p;
	delete[] cumul;

	// send back allocations to master
	if (GetMyid())	{
		MPI_Send(alloc,GetFinalNsite(),MPI_INT,0,TAG1,MPI_COMM_WORLD);
	}

	return logw;
}
*/

