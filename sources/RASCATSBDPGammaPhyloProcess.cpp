
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

double RASCATSBDPGammaPhyloProcess::GlobalRestrictedTemperedMove()	{

	double tuning = 1.0;
	// important to start with that one
	// if marginal suff stat move is done before that in a multi gene context

	if (TemperedGene())	{
		PoissonSBDPProfileProcess::Move(1,1,1,1);
		GlobalUpdateParameters();
	}

	if (TemperedRate())	{
		DGamRateProcess::Move(tuning,10);
		DGamRateProcess::Move(0.3*tuning,10);
		DGamRateProcess::Move(0.03*tuning,10);
		DGamRateProcess::Move(0.3*tuning,10);
		DGamRateProcess::Move(tuning,10);
		GlobalUpdateParameters();
	}

	if (TemperedBL())	{
		GammaBranchProcess::Move(tuning,50);
		GlobalUpdateParameters();
	}

	if (TemperedRate())	{
		DGamRateProcess::Move(tuning,10);
		DGamRateProcess::Move(0.3*tuning,10);
		DGamRateProcess::Move(0.03*tuning,10);
		DGamRateProcess::Move(0.3*tuning,10);
		DGamRateProcess::Move(tuning,10);
		GlobalUpdateParameters();
	}

	if (TemperedGene())	{
		PoissonSBDPProfileProcess::Move(1,1,1,1);
		GlobalUpdateParameters();
	}
}

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
	case MTRYALLOC:
		SlaveChooseMultipleTryAlloc();
		break;
	case ACTIVATEMTRY:
		SlaveActivateSumOverComponents();
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

void RASCATSBDPGammaPhyloProcess::ReadStatMin(string name, int burnin, int every, int until)	{

	if (GetNprocs() > 1)	{
		cerr << "error: read stat min only in serial mode\n";
		exit(1);
	}

	ifstream is((name + ".chain").c_str());
	if (!is)	{
		cerr << "error: no .chain file found\n";
		exit(1);
	}

	cerr << "burnin : " << burnin << "\n";
	cerr << "until : " << until << '\n';
	int i=0;
	while ((i < until) && (i < burnin))	{
		FromStream(is);
		i++;
	}
	int samplesize = 0;

	double meanvar = 0;

	ofstream os((name + ".statmin").c_str());
	while (i < until)	{
		cerr << ".";
		cerr.flush();
		samplesize++;
		FromStream(is);
		i++;

		QuickUpdate();

		GlobalCollapse();

		GlobalUpdateSiteProfileSuffStat();
		// GlobalUpdateModeProfileSuffStat();

		/*
		for (int i=0; i<GetNsite(); i++)	{
			for (int k=0; k<GetNcomponent(); k++)	{
				os << log(weight[k]) << '\t' << log(GetMinStat(profile[k],i)) << '\t' << LogStatProb(i,k) << '\n';
			}
		}
		*/

		double* samplingprob = new double[GetNmodeMax()];
		double* fulllogp = new double[GetNmodeMax()];

		double logl = 0;
		double totrelvar = 0;
		for (int i=0; i<GetNsite(); i++)	{
			
			// what we want to estimate
			double max = 0;
			for (int k=0; k<GetNmodeMax(); k++)	{
				fulllogp[k] = LogStatProb(i,k);
				if ((! k) || (max < fulllogp[k]))	{
					max = fulllogp[k];
				}
			}
			double total = 0;
			for (int k=0; k<GetNmodeMax(); k++)	{
				total += weight[k] * exp(fulllogp[k] - max);
			}
			double m1 = log(total) + max;
			
			// compute sampling weights
			double tot = 0;
			for (int k=0; k<GetNmodeMax(); k++)	{
				double tmp = weight[k] * GetMinStat(profile[k],i);
				if (! tmp)	{
					cerr << "null min stat : " << weight[k] << '\t' << fulllogp[k] - max << '\t' << tmp << '\n';
					exit(1);
				}
				tot += tmp;
				samplingprob[k] = tmp;
			}
			for (int k=0; k<GetNmodeMax(); k++)	{
				samplingprob[k] /= tot;
			}

			double total2 = 0;
			for (int k=0; k<GetNmodeMax(); k++)	{
				total2 += weight[k] * weight[k] / samplingprob[k] * exp(2 * (fulllogp[k] - max));
			}
			total2 -= total*total;

			// double m2 = log(total2) + 2*max;
			// m2 -= m1*m1;

			// relative variance of estimator of L_i ~ variance of estimator of ln L_i
			double relvar = total2 / total / total;
			// cerr << total << '\t' << total2 << '\t' << relvar << '\n';

			logl += m1;
			totrelvar += relvar;
		}
			
		os << logl << '\t' << totrelvar << '\n';
		meanvar += totrelvar;

		GlobalUnfold();
		int nrep = 1;
		while ((i<until) && (nrep < every))	{
			FromStream(is);
			i++;
			nrep++;
		}
	}
	cerr << '\n';
	meanvar /= samplesize;
	cout << '\n';
	cout << "mean variance : " << meanvar << '\n';
	cout << '\n';
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
	int ncomp = GetNcomponent();
	if ((sumovercomponents > 0) && (sumovercomponents < GetNcomponent()))	{
		ncomp = sumovercomponents;
	}

	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i))	{
			modesitelogL[i] = new double[ncomp];
		}
	}

	double totlogL = 0;

	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i))	{
			RemoveSite(i,SBDPProfileProcess::alloc[i]);
		}
	}

	for (int k=0; k<ncomp; k++)	{
		for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
			if (ActiveSite(i))	{
				if (ncomp == GetNcomponent())	{
					AddSite(i,k);
				}
				else	{
					AddSite(i,mtryalloc[i][k]);
				}
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
				if (ncomp == GetNcomponent())	{
					RemoveSite(i,k);
				}
				else	{
					RemoveSite(i,mtryalloc[i][k]);
				}
			}
		}
	}

	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i))	{
			double max = modesitelogL[i][0];
			for (int k=1; k<ncomp; k++)	{
				if (max < modesitelogL[i][k])	{
					max = modesitelogL[i][k];
				}
			}
			double total = 0;
			double cumul[ncomp];
			for (int k=0; k<ncomp; k++)	{
				double w = 0;
				if (ncomp == GetNcomponent())	{
					w = weight[k];
				}
				else	{
					w = weight[mtryalloc[i][k]] / mtryweight[i][k];
				}
				double tmp = w * exp(modesitelogL[i][k] - max);
				total += tmp;
				cumul[k] = total;
			}

			double u = total * rnd::GetRandom().Uniform();
			int k = 0;
			while ((k<ncomp) && (u>cumul[k]))	{
				k++;
			}

			if (ncomp == GetNcomponent())	{
				AddSite(i,k);
			}
			else	{
				AddSite(i,mtryalloc[i][k]);
			}
			UpdateZip(i);

			if (ncomp < GetNcomponent())	{
				total /= ncomp;
			}
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

/*
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
*/

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
