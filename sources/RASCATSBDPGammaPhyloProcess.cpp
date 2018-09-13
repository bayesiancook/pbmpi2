
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
		PoissonSBDPProfileProcess::Move(1,5,1,1);
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
		PoissonSBDPProfileProcess::Move(1,5,1,1);
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

	int bksitemax = sitemax[myid];
	sitemax[myid] = GetSiteMin() + testsitemax - testsitemin;

	double** sitelogl = new double*[GetNsite()];
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		sitelogl[i] = new double[GetNcomponent()];
	}
	
	// double cutoff = 0;
	double cutoff = 1e-5;
	double total = 0;

	if (cutoff)	{
		int kmax = 0;
		double totw = weight[0];
		while ((kmax < GetNcomponent()) && (fabs(1-totw) > cutoff)) {
			kmax ++;
			totw += weight[kmax];
		}
		if (kmax == GetNcomponent())	{
			cerr << "error in SlaveComputeCVScore: overflow\n";
			exit(1);
		}
		kmax++;
		double remw = 1 - totw;
		int Ncomp = kmax;
		if ((kmax < GetNcomponent()) && (remw > 0))	{
			Ncomp++;
		}
		vector<double> w(Ncomp,0);
		for (int k=0; k<kmax; k++)	{
			for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
				PoissonSBDPProfileProcess::alloc[i] = k;
				UpdateZip(i);
			}
			UpdateConditionalLikelihoods();
			for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
				sitelogl[i][k] = sitelogL[i];
			}
			w[k] = weight[k];
		}

		if ((kmax < GetNcomponent()) && (remw > 0))	{
			int M = GetNcomponent() - kmax;
			vector<double> cumul(M,0);
			double tot = 0;
			for (int k=0; k<M; k++)	{
				tot += weight[kmax+k];
				cumul[k] = tot;
			}
			double q = tot * rnd::GetRandom().Uniform();
			int k = 0;
			while ((k<M) && (q>cumul[k]))	{
				k++;
			}
			if (k == M)	{
				cerr << "error in SlaveComputeCVScore: overflow when choosing low weight component\n";
				exit(1);
			}
			for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
				PoissonSBDPProfileProcess::alloc[i] = k;
				UpdateZip(i);
			}
			UpdateConditionalLikelihoods();
			for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
				sitelogl[i][kmax] = sitelogL[i];
			}
			w[kmax] = remw;
		}

		for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
			double max = 0;
			for (int k=0; k<Ncomp; k++)	{
				if ((!k) || (max < sitelogl[i][k]))	{
					max = sitelogl[i][k];
				}
			}
			double tot = 0;
			double totweight = 0;
			for (int k=0; k<Ncomp; k++)	{
				tot += w[k] * exp(sitelogl[i][k] - max);
				totweight += w[k];
			}
			total += log(tot) + max;
		}
	}
	else	{
		for (int k=0; k<GetNcomponent(); k++)	{
			for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
				PoissonSBDPProfileProcess::alloc[i] = k;
				UpdateZip(i);
			}
			UpdateConditionalLikelihoods();
			for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
				sitelogl[i][k] = sitelogL[i];
			}
		}

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
			total += log(tot) + max;
		}
	}

	MPI_Send(&total,1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
	
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		delete[] sitelogl[i];
	}
	delete[] sitelogl;
	sitemax[myid] = bksitemax;
}

/*
void RASCATSBDPGammaPhyloProcess::SlaveComputeCVScore()	{

	int sitemin = GetSiteMin();
	int sitemax = GetSiteMin() + testsitemax - testsitemin;

	int ncomp = GetNcomponent();
	if ((sumovercomponents > 0) && (sumovercomponents < GetNcomponent()))	{
		ncomp = sumovercomponents;
	}

	double* modesitelogL = new double[ncomp];

	double totlogL = 0;

	if (ncomp == GetNcomponent())	{

		for (int i=sitemin; i<sitemax; i++) {
			if (ActiveSite(i))	{

				// remove site
				RemoveSite(i,SBDPProfileProcess::alloc[i]);

				double max = 0;
				for (int k=0; k<GetNcomponent(); k++)	{
					AddSite(i,k);
					modesitelogL[k] = SiteLogLikelihood(i);
					if ((!k) || (max < modesitelogL[k]))	{
						max = modesitelogL[k];
					}
					RemoveSite(i,k);
				}

				double total = 0;
				double cumul[GetNcomponent()];
				for (int k=0; k<GetNcomponent(); k++)	{
					double tmp = weight[k] * exp(modesitelogL[k] - max);
					total += tmp;
					cumul[k] = total;
				}

				double u = total * rnd::GetRandom().Uniform();
				int k = 0;
				while ((k<GetNcomponent()) && (u>cumul[k]))	{
					k++;
				}

                AddSite(i,k);
                sitelogL[i] = modesitelogL[k];

				double sitetotlogL = log(total) + max;
				totlogL += sitetotlogL;
			}
		}
	}

	else	{

		for (int i=sitemin; i<sitemax; i++) {
			if (ActiveSite(i))	{

				// remove site
				RemoveSite(i,SBDPProfileProcess::alloc[i]);

				double max = 0;
				for (int k=0; k<ncomp; k++)	{
					if ((mtryalloc[i][k] < 0) || (mtryalloc[i][k] >= GetNcomponent()))	{
						cerr << "error in RASCATSBDPGammaPhyloProcess::GetFullLogLikelihood: wrong alloc : " << mtryalloc[i][k] << '\n';
						cerr << "site " << i << '\t' << "component: " << k << '\n';
						exit(1);
					}
					int found = 0;
					for (int l=0; l<k; l++)	{
						if (mtryalloc[i][k] == mtryalloc[i][l])	{
							found = 1;
							modesitelogL[k] = modesitelogL[l];
						}
					}
					if (! found)	{
						AddSite(i,mtryalloc[i][k]);
						modesitelogL[k] = SiteLogLikelihood(i);
						RemoveSite(i,mtryalloc[i][k]);
					}
					if ((!k) || (max < modesitelogL[k]))	{
						max = modesitelogL[k];
					}
				}

				double total = 0;
				double cumul[ncomp];
				for (int k=0; k<ncomp; k++)	{
					double tmp = exp(modesitelogL[k] - max) / mtryweight[i][k];
					// double tmp = weight[mtryalloc[i][k]] * exp(modesitelogL[k] - max);
					total += tmp;
					cumul[k] = total;
				}

				double u = total * rnd::GetRandom().Uniform();
				int k = 0;
				while ((k<ncomp) && (u>cumul[k]))	{
					k++;
				}

                AddSite(i,mtryalloc[i][k]);
                sitelogL[i] = modesitelogL[k];

				double sitetotlogL = log(total) + max;
				if (isinf(sitetotlogL))	{
					cerr << "sitetotlogL is inf: site " << i << '\n';
				}
				totlogL += sitetotlogL;
			}
		}
	}

	if (isinf(totlogL))	{
		cerr << "error in GetFullLogLikelihood: inf\n";
		exit(1);
	}

	MPI_Send(&totlogL,1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);

	delete[] modesitelogL;
}
*/

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

	int ncomp = GetNcomponent();
	if ((sumovercomponents > 0) && (sumovercomponents < GetNcomponent()))	{
		ncomp = sumovercomponents;
	}

	double* modesitelogL = new double[ncomp];

	double totlogL = 0;

	if (ncomp == GetNcomponent())	{

		for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
			if (ActiveSite(i))	{

				// remove site
				RemoveSite(i,SBDPProfileProcess::alloc[i]);

				double max = 0;
				for (int k=0; k<GetNcomponent(); k++)	{
					AddSite(i,k);
					modesitelogL[k] = SiteLogLikelihood(i);
					if ((!k) || (max < modesitelogL[k]))	{
						max = modesitelogL[k];
					}
					RemoveSite(i,k);
				}

				double total = 0;
				double cumul[GetNcomponent()];
				for (int k=0; k<GetNcomponent(); k++)	{
					double tmp = weight[k] * exp(modesitelogL[k] - max);
					total += tmp;
					cumul[k] = total;
				}

				double u = total * rnd::GetRandom().Uniform();
				int k = 0;
				while ((k<GetNcomponent()) && (u>cumul[k]))	{
					k++;
				}

				if (! reverseafterfull)	{
					AddSite(i,k);
					sitelogL[i] = modesitelogL[k];
					/*
					if (i == 33)	{
						cerr << "adding site " << i << "to " << k << '\n';
						cerr << modesitelogL[k] << '\n';
					}
					*/
				}

				double sitetotlogL = log(total) + max;
				totlogL += sitetotlogL;
			}
		}
	}

	else	{

		for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
			if (ActiveSite(i))	{

				// remove site
				RemoveSite(i,SBDPProfileProcess::alloc[i]);

				double max = 0;
				for (int k=0; k<ncomp; k++)	{
					if ((mtryalloc[i][k] < 0) || (mtryalloc[i][k] >= GetNcomponent()))	{
						cerr << "error in RASCATSBDPGammaPhyloProcess::GetFullLogLikelihood: wrong alloc : " << mtryalloc[i][k] << '\n';
						cerr << "site " << i << '\t' << "component: " << k << '\n';
						exit(1);
					}
					int found = 0;
					for (int l=0; l<k; l++)	{
						if (mtryalloc[i][k] == mtryalloc[i][l])	{
							found = 1;
							modesitelogL[k] = modesitelogL[l];
						}
					}
					if (! found)	{
						AddSite(i,mtryalloc[i][k]);
						modesitelogL[k] = SiteLogLikelihood(i);
						RemoveSite(i,mtryalloc[i][k]);
					}
					if ((!k) || (max < modesitelogL[k]))	{
						max = modesitelogL[k];
					}
				}

				double total = 0;
				double cumul[ncomp];
				for (int k=0; k<ncomp; k++)	{
					double tmp = exp(modesitelogL[k] - max) / mtryweight[i][k];
					// double tmp = weight[mtryalloc[i][k]] * exp(modesitelogL[k] - max);
					total += tmp;
					cumul[k] = total;
				}

				double u = total * rnd::GetRandom().Uniform();
				int k = 0;
				while ((k<ncomp) && (u>cumul[k]))	{
					k++;
				}

				if (! reverseafterfull)	{
					AddSite(i,mtryalloc[i][k]);
					sitelogL[i] = modesitelogL[k];
				}

				double sitetotlogL = log(total) + max;
				if (isinf(sitetotlogL))	{
					cerr << "sitetotlogL is inf: site " << i << '\n';
				}
				totlogL += sitetotlogL;
			}
		}
	}

	if (isinf(totlogL))	{
		cerr << "error in GetFullLogLikelihood: inf\n";
		exit(1);
	}

	delete[] modesitelogL;
	return totlogL;
}

double RASCATSBDPGammaPhyloProcess::GlobalGetFullLogLikelihood()	{

	double totlogL = PhyloProcess::GlobalGetFullLogLikelihood();
	if (! reverseafterfull)	{
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
		GlobalUpdateParameters();
	}
	return totlogL;
}

void RASCATSBDPGammaPhyloProcess::SlaveGetFullLogLikelihood()	{

	PhyloProcess::SlaveGetFullLogLikelihood();
	if (! reverseafterfull)	{
		MPI_Send(SBDPProfileProcess::alloc,GetNsite(),MPI_INT,0,TAG1,MPI_COMM_WORLD);
	}
}

void RASCATSBDPGammaPhyloProcess::ReadProfileDistribution(string name, int burnin, int every, int until, int ndisc, double cialpha, int nsample)	{

	double minhi = 100;
	double maxhi = -100;
	for (int a=0; a<Naa; a++)	{
		if (maxhi < HydrophobicityIndex_pH7[a])	{
			maxhi = HydrophobicityIndex_pH7[a];
		}
		if (minhi > HydrophobicityIndex_pH7[a])	{
			minhi = HydrophobicityIndex_pH7[a];
		}
	}

	double mines = 1;
	double maxes = 20;

	vector<double>* hicdf = new vector<double>[ndisc+1];
	double* higrid = new double[ndisc+1];
	for (int l=0; l<=ndisc; l++)	{
		higrid[l] = minhi + ((double) l) / ndisc * (maxhi - minhi);
	}
	double* tmphicdf = new double[ndisc+1];
	double* meanhicdf = new double[ndisc+1];
	for (int l=0; l<=ndisc; l++)	{
		meanhicdf[l] = 0;
	}

	vector<double>* escdf = new vector<double>[ndisc+1];
	double* esgrid = new double[ndisc+1];
	for (int l=0; l<=ndisc; l++)	{
		esgrid[l] = mines + ((double) l) / ndisc * (maxes - mines);
	}
	double* tmpescdf = new double[ndisc+1];
	double* meanescdf = new double[ndisc+1];
	for (int l=0; l<=ndisc; l++)	{
		meanescdf[l] = 0;
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

	ofstream mos((name + ".meandist").c_str());

	while (i < until)	{
		cerr << ".";
		cerr.flush();
		samplesize++;
		FromStream(is);
		ResampleWeights();
		i++;

		for (int l=0; l<=ndisc; l++)	{
			tmphicdf[l] = 0;
			tmpescdf[l] = 0;
		}

		// get spike distribution
		ostringstream s;
		s << name << "_" << samplesize << ".profilespikes";
		ofstream os(s.str().c_str());
		for (int i=0; i<GetNcomponent(); i++)	{
			os << weight[i] << '\t' << GetHI7(profile[i]) << '\t' << GetMolWeight(profile[i]) << '\t' << GetEffSize(profile[i]) << '\t' << GetTwoMaxRatio(profile[i]) << '\t' << GetMinMaxRatio(profile[i]);
			for (int k=0; k<GetDim(); k++)	{
				os << '\t' << profile[i][k];
			}
			os << '\n';

			double h = GetHI7(profile[i]);
			for (int l=0; l<=ndisc; l++)	{
				if (h < higrid[l])	{
					tmphicdf[l] += weight[i];
				}
			}

			double es = GetEffSize(profile[i]);
			for (int l=0; l<=ndisc; l++)	{
				if (h < esgrid[l])	{
					tmpescdf[l] += weight[i];
				}
			}
		}

		for (int l=0; l<=ndisc; l++)	{
			hicdf[l].push_back(tmphicdf[l]);
			meanhicdf[l] += tmphicdf[l];
		}
			
		for (int l=0; l<=ndisc; l++)	{
			escdf[l].push_back(tmpescdf[l]);
			meanescdf[l] += tmpescdf[l];
		}
			
		// this is the distribution of site-specific profiles, across sites, for current MCMC point
		ostringstream s2;
		s2 << name << "_" << samplesize << ".siteprofiledist";
		ofstream os2(s2.str().c_str());

		for (int i=0; i<GetNsite(); i++)	{
			double* p = GetProfile(i);
			os2 << i << '\t' << GetHI7(p) << '\t' << GetMolWeight(p) << '\t' << GetEffSize(p) << '\t' << GetTwoMaxRatio(p) << '\t' << GetMinMaxRatio(p);
			for (int k=0; k<GetDim(); k++)	{
				os2 << '\t' << p[k];
			}
			os2 << '\n';
		}

		// this is a random sample of nsample profiles from the current mixture
		// this random sample is saved in a separate file, for current MCMC point
		ostringstream s3;
		s3 << name << "_" << samplesize << ".profiledist";
		ofstream os3(s3.str().c_str());

		for (int i=0; i<nsample; i++)	{

			int k = rnd::GetRandom().DrawFromDiscreteDistribution(weight,GetNcomponent());
			double* p = profile[k];

			os3 << i << '\t' << GetHI7(p) << '\t' << GetMolWeight(p) << '\t' << GetEffSize(p) << '\t' << GetTwoMaxRatio(p) << '\t' << GetMinMaxRatio(p);
			for (int k=0; k<GetDim(); k++)	{
				os3 << '\t' << p[k];
			}
			os3 << '\n';
		}

		// another random sample, thinned, and pooled with all samples across the MCMC
		for (int i=0; i<nsample/10; i++)	{

			int k = rnd::GetRandom().DrawFromDiscreteDistribution(weight,GetNcomponent());
			double* p = profile[k];

			mos << i << '\t' << GetHI7(p) << '\t' << GetMolWeight(p) << '\t' << GetEffSize(p) << '\t' << GetTwoMaxRatio(p) << '\t' << GetMinMaxRatio(p);
			for (int k=0; k<GetDim(); k++)	{
				mos << '\t' << p[k];
			}
			mos << '\n';
		}

		int nrep = 1;
		while ((i<until) && (nrep < every))	{
			FromStream(is);
			i++;
			nrep++;
		}
	}
	cerr << '\n';

	for (int l=0; l<=ndisc; l++)	{
		meanhicdf[l] /= samplesize;
	}

	double* minhicdf = new double[ndisc+1];
	double* maxhicdf = new double[ndisc+1];

	int nmin = ((int) ( ((double) samplesize) * 0.5*cialpha));
	int nmax = samplesize - nmin;
	if (nmax == samplesize)	{
		nmax--;
	}

	for (int l=0; l<=ndisc; l++)	{
		sort(hicdf[l].begin(),hicdf[l].end());
		minhicdf[l] = hicdf[l][nmin];
		maxhicdf[l] = hicdf[l][nmax];
	}

	ofstream hos((name + ".hicdf").c_str());
	for (int l=0; l<=ndisc; l++)	{
		hos << higrid[l] << '\t' << meanhicdf[l] << '\t' << minhicdf[l] << '\t' << maxhicdf[l] << '\n';
	}

	double* minescdf = new double[ndisc+1];
	double* maxescdf = new double[ndisc+1];

	for (int l=0; l<=ndisc; l++)	{
		sort(escdf[l].begin(),escdf[l].end());
		minescdf[l] = escdf[l][nmin];
		maxescdf[l] = escdf[l][nmax];
	}

	ofstream esos((name + ".escdf").c_str());
	for (int l=0; l<=ndisc; l++)	{
		esos << esgrid[l] << '\t' << meanescdf[l] << '\t' << minescdf[l] << '\t' << maxescdf[l] << '\n';
	}

	delete[] meanhicdf;
	delete[] minhicdf;
	delete[] maxhicdf;
	delete[] tmphicdf;
	delete[] hicdf;
	delete[] higrid;

	delete[] meanescdf;
	delete[] minescdf;
	delete[] maxescdf;
	delete[] tmpescdf;
	delete[] escdf;
	delete[] esgrid;
}
