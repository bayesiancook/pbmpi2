
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "RASCATSBDPGammaPhyloProcess.h"
#include "Parallel.h"
#include <string>

void RASCATSBDPGammaPhyloProcess::GlobalUpdateParameters()	{

	if (GetNprocs() > 1)	{
		RASCATGammaPhyloProcess::GlobalUpdateParameters();
		MPI_Bcast(V,GetNcomponent(),MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(weight,GetNcomponent(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	}
	else	{
		UpdateZip();
	}
}

void RASCATSBDPGammaPhyloProcess::SlaveUpdateParameters()	{

	RASCATGammaPhyloProcess::SlaveUpdateParameters();
	MPI_Bcast(V,GetNcomponent(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(weight,GetNcomponent(),MPI_DOUBLE,0,MPI_COMM_WORLD);
}

void RASCATSBDPGammaPhyloProcess::SlaveExecute(MESSAGE signal)	{

	switch(signal) {

	case MIX_MOVE:
		SlaveMixMove();
		break;
	case UPDATE_RATE:
		SlaveUpdateRateSuffStat();
		break;
	default:
		PhyloProcess::SlaveExecute(signal);
	}
}

void RASCATSBDPGammaPhyloProcess::SlaveComputeCVScore()	{

	int sitemin = GetSiteMin();
	int sitemax = GetSiteMin() + testsitemax - testsitemin;
	double** sitelogl = new double*[GetNsite()];
	for (int i=sitemin; i<sitemax; i++)	{
		sitelogl[i] = new double[GetNcomponent()];
	}
	
	// UpdateMatrices();

	for (int k=0; k<GetNcomponent(); k++)	{
		for (int i=sitemin; i<sitemax; i++)	{
			PoissonSBDPProfileProcess::alloc[i] = k;
			UpdateZip(i);
		}
		UpdateConditionalLikelihoods();
		for (int i=sitemin; i<sitemax; i++)	{
			sitelogl[i][k] = sitelogL[i];
		}
	}

	double total = 0;
	for (int i=sitemin; i<sitemax; i++)	{
		double max = 0;
		for (int k=0; k<GetNcomponent(); k++)	{
			if ((!k) || (max < sitelogl[i][k]))	{
				max = sitelogl[i][k];
			}
		}
		double tot = 0;
		double totweight = 0;
		for (int k=0; k<GetNcomponent(); k++)	{
			tot += weight[k] * exp(sitelogl[i][k] - max);
			totweight += weight[k];
		}
		total += log(tot) + max;
	}

	MPI_Send(&total,1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
	
	for (int i=sitemin; i<sitemax; i++)	{
		delete[] sitelogl[i];
	}
	delete[] sitelogl;
}

void RASCATSBDPGammaPhyloProcess::SlaveComputeSiteLogL()	{

	if (! SumOverRateAllocations())	{
		cerr << "rate error\n";
		exit(1);
	}

	double** sitelogl = new double*[GetNsite()];
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		sitelogl[i] = new double[GetNcomponent()];
	}
	
	// UpdateMatrices();

	for (int k=0; k<GetNcomponent(); k++)	{
		for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
			PoissonSBDPProfileProcess::alloc[i] = k;
		}
		UpdateConditionalLikelihoods();
		for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
			sitelogl[i][k] = sitelogL[i];
		}
	}

	double* meansitelogl = new double[GetNsite()];
	for (int i=0; i<GetNsite(); i++)	{
		meansitelogl[i] = 0;
	}
	double total = 0;
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		double max = 0;
		for (int k=0; k<GetNcomponent(); k++)	{
			if ((!k) || (max < sitelogl[i][k]))	{
				max = sitelogl[i][k];
			}
		}
		double tot = 0;
		double totweight = 0;
		for (int k=0; k<GetNcomponent(); k++)	{
			tot += weight[k] * exp(sitelogl[i][k] - max);
			totweight += weight[k];
		}
		meansitelogl[i] = log(tot) + max;
		total += meansitelogl[i] ;
	}

	MPI_Send(meansitelogl,GetNsite(),MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
	
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		delete[] sitelogl[i];
	}
	delete[] sitelogl;
	delete[] meansitelogl;

}


double RASCATSBDPGammaPhyloProcess::GetFullLogLikelihood()	{

	double** modesitelogL = new double*[GetNsite()];
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i))	{
			modesitelogL[i] = new double[GetNcomponent()];
		}
	}

	double totlogL = 0;

	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i))	{
			RemoveSite(i,SBDPProfileProcess::alloc[i]);
		}
	}

	for (int k=0; k<GetNcomponent(); k++)	{
		for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
			if (ActiveSite(i))	{
				AddSite(i,k);
				UpdateZip(i);
			}
		}
		UpdateConditionalLikelihoods();
		for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
			if (ActiveSite(i))	{
				modesitelogL[i][k] = sitelogL[i];
			}
		}
		for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
			if (ActiveSite(i))	{
				RemoveSite(i,k);
			}
		}
	}

	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i))	{
			double max = modesitelogL[i][0];
			for (int k=1; k<GetNcomponent(); k++)	{
				if (max < modesitelogL[i][k])	{
					max = modesitelogL[i][k];
				}
			}
			double total = 0;
			double cumul[GetNcomponent()];
			for (int k=0; k<GetNcomponent(); k++)	{
				double tmp = weight[k] * exp(modesitelogL[i][k] - max);
				total += tmp;
				cumul[k] = total;
			}

			double u = total * rnd::GetRandom().Uniform();
			int k = 0;
			while ((k<GetNcomponent()) && (u>cumul[k]))	{
				k++;
			}

			AddSite(i,k);
			UpdateZip(i);

			double sitetotlogL = log(total) + max;
			totlogL += sitetotlogL;
		}
	}

	// one last update so that cond likelihoods are in sync with new site allocations
	// this will also update logL
	UpdateConditionalLikelihoods();

	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i))	{
			delete[] modesitelogL[i];
		}
	}
	delete[] modesitelogL;

	return totlogL;
}

double RASCATSBDPGammaPhyloProcess::GlobalGetFullLogLikelihood()	{

	double totlogL = PhyloProcess::GlobalGetFullLogLikelihood();
	if (sumovercomponents && (Ncomponent > 1))	{
		// receive allocs from slaves
		MPI_Status stat;
		int tmpalloc[GetNsite()];
		for(int i=1; i<GetNprocs(); ++i) {
			MPI_Recv(tmpalloc,GetNsite(),MPI_INT,i,TAG1,MPI_COMM_WORLD,&stat);
			for(int j=GetProcSiteMin(i); j<GetProcSiteMax(i); ++j) {
				if (ActiveSite(j))	{
					SBDPProfileProcess::alloc[j] = tmpalloc[j];
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
		/*
		ResampleWeights();
		*/
	}
	GlobalUpdateParameters();
	return totlogL;
}

void RASCATSBDPGammaPhyloProcess::SlaveGetFullLogLikelihood()	{

	PhyloProcess::SlaveGetFullLogLikelihood();
	if (sumovercomponents && (Ncomponent > 1))	{
		MPI_Send(SBDPProfileProcess::alloc,GetNsite(),MPI_INT,0,TAG1,MPI_COMM_WORLD);
	}
}
