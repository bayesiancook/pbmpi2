
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "StringStreamUtils.h"
#include "RASCATFiniteGammaPhyloProcess.h"
#include "Parallel.h"
#include <string>

double RASCATFiniteGammaPhyloProcess::GlobalRestrictedTemperedMove()	{

	double tuning = 1.0;
	// important to start with that one
	// if marginal suff stat move is done before that in a multi gene context

	if (TemperedBL())	{
		// GammaBranchProcess::Move(tuning,0);
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
	// if (TemperedProfile())	{
		PoissonFiniteProfileProcess::Move(1,1,1);
		GlobalUpdateParameters();
	}
}


void RASCATFiniteGammaPhyloProcess::GlobalUpdateParameters()	{
	if (GetNprocs() > 1)	{
	// MPI2
	// should send the slaves the relevant information
	// about model parameters

	// for this model, should broadcast
	// double alpha
	// int Ncomponent
	// int* alloc
	// double* rr
	// double** profile
	// double* brancharray
	// (but should first call PutBranchLengthsIntoArray())
	// 
	// upon receiving this information
	// slave should 
	// store it in the local copies of the variables
	// and then call
	// SetBranchLengthsFromArray()
	// SetAlpha(inalpha)

	// ResampleWeights();
	RenormalizeProfiles();

	int i,j,nbranch = GetNbranch(),ni,nd,L1,L2;
	L1 = GetNmodeMax();
	L2 = GetDim();
	nd = 1 + nbranch + L1*L2 + GetDim() + 1;
	ni = 1 + GetNsite();
	int ivector[ni];
	double dvector[nd]; 
	MESSAGE signal = PARAMETER_DIFFUSION;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	// First we assemble the vector of doubles for distribution
	int index = 0;
	dvector[index] = GetAlpha();
	index++;
	for(i=0; i<nbranch; ++i) {
		dvector[index] = blarray[i];
		index++;
	}
	
	for(i=0; i<L1; ++i) {
		for(j=0; j<L2; ++j) {
			dvector[index] = profile[i][j];
			index++;
		}
	}
	for (int i=0; i<GetDim(); i++)	{
		dvector[index] = dirweight[i];
		index++;
	}

	// Now the vector of ints
	ivector[0] = GetNcomponent();
	for(i=0; i<GetNsite(); ++i) {
		ivector[1+i] = FiniteProfileProcess::alloc[i];
	}

	// Now send out the doubles and ints over the wire...
	MPI_Bcast(ivector,ni,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(dvector,nd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(weight,GetNcomponent(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	}
	else	{
		UpdateZip();
	}
}


double RASCATFiniteGammaPhyloProcess::GetFullLogLikelihood()	{

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
			RemoveSite(i,FiniteProfileProcess::alloc[i]);
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

			if (! reverseafterfull)	{
				if (ncomp == GetNcomponent())	{
					AddSite(i,k);
				}
				else	{
					AddSite(i,mtryalloc[i][k]);
				}
				UpdateZip(i);
			}

			if (ncomp < GetNcomponent())	{
				total /= ncomp;
			}
			double sitetotlogL = log(total) + max;
			totlogL += sitetotlogL;
		}
	}

	// one last update so that cond likelihoods are in sync with new site allocations
	// this will also update logL
	if (! reverseafterfull)	{
		UpdateConditionalLikelihoods();
	}

	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i))	{
			delete[] modesitelogL[i];
		}
	}
	delete[] modesitelogL;
	return totlogL;
}

double RASCATFiniteGammaPhyloProcess::GlobalGetFullLogLikelihood()	{

	double totlogL = PhyloProcess::GlobalGetFullLogLikelihood();
	// receive allocs from slaves
	if (! reverseafterfull)	{
		MPI_Status stat;
		int tmpalloc[GetNsite()];
		for(int i=1; i<GetNprocs(); ++i) {
			MPI_Recv(tmpalloc,GetNsite(),MPI_INT,i,TAG1,MPI_COMM_WORLD,&stat);
			for(int j=GetProcSiteMin(i); j<GetProcSiteMax(i); ++j) {
				if (ActiveSite(j))	{
					FiniteProfileProcess::alloc[j] = tmpalloc[j];
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

void RASCATFiniteGammaPhyloProcess::SlaveGetFullLogLikelihood()	{

	PhyloProcess::SlaveGetFullLogLikelihood();
	if (! reverseafterfull)	{
		MPI_Send(FiniteProfileProcess::alloc,GetNsite(),MPI_INT,0,TAG1,MPI_COMM_WORLD);
	}
}

void RASCATFiniteGammaPhyloProcess::SlaveExecute(MESSAGE signal)	{

	switch(signal) {

	case UPDATE_RATE:
		SlaveUpdateRateSuffStat();
		break;
	case REALLOC_MOVE:
		SlaveIncrementalFiniteMove();
		break;
	case STATFIX:
		SlaveGetStatFix();
		break;
	case ACTIVATEMTRY:
		SlaveActivateSumOverComponents();
		break;
	case MTRYALLOC:
		SlaveChooseMultipleTryAlloc();
		break;
	default:
		PhyloProcess::SlaveExecute(signal);
	}
}


void RASCATFiniteGammaPhyloProcess::SlaveUpdateParameters()	{
	int i,j,L1,L2,ni,nd,nbranch = GetNbranch();
	L1 = GetNmodeMax();
	L2 = GetDim();
	nd = 1 + nbranch + L1*L2 + GetDim() + 1;
	ni = 1 + GetNsite();
	int* ivector = new int[ni];
	double* dvector = new double[nd];
	MPI_Bcast(ivector,ni,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(dvector,nd,MPI_DOUBLE,0,MPI_COMM_WORLD);

	int index = 0;
	SetAlpha(dvector[index]);
	index++;
	for(i=0; i<nbranch; ++i) {
		blarray[i] = dvector[index];
		index++;
	}

	for(i=0; i<L1; ++i) {
		for(j=0; j<L2; ++j) {
			profile[i][j] = dvector[index];
			index++;
		}
	}
	for (int i=0; i<GetDim(); i++)	{
		dirweight[i] = dvector[index];
		index++;
	}

	Ncomponent = ivector[0];
	for(i=0; i<GetNsite(); ++i) {
		FiniteProfileProcess::alloc[i] = ivector[1+i];
	}

	MPI_Bcast(weight,GetNcomponent(),MPI_DOUBLE,0,MPI_COMM_WORLD);

	UpdateZip();

	delete[] dvector;
	delete[] ivector;

	// some upate here ?
}

void RASCATFiniteGammaPhyloProcess::ReadPB(int argc, char* argv[])	{

	string name = "";

	int burnin = -1;
	int every = 1;
	int until = -1;
	int ppred = 0;
	// 1 : plain ppred (outputs simulated data)
	// 2 : diversity statistic
	// 3 : compositional statistic
	int cv = 0;
	int sitelogl = 0;
	int map = 0;
	int rates = 0;
	string testdatafile = "";

	int rateprior = 0;
	int profileprior = 0;
	int rootprior = 1;

	string taxon1 = "None";
	string taxon2 = "None";
	string taxon3 = "None";
	string taxon4 = "None";
	int toponfrac = 100;
	int toponstep = 10;
	int bf = 0;
	double bfprop = 0.5;

	int sis = 0;
	double sisprop = 0.5;

	int temperedbl = 1;
	int temperedgene = 0;
	int temperedrate = 0;

	int sumcomp = 0;

	int bfl = 0;
	roottax1 = "None";
	roottax2 = "None";

	try	{

		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];
			if (s == "-div")	{
				ppred = 2;
			}
			else if (s == "-comp")	{
				ppred = 3;
			}
			else if (s == "-ppred")	{
				ppred = 1;
			}
			else if (s == "-ppredrate")	{
				i++;
				string tmp = argv[i];
				if (tmp == "prior")	{
					rateprior = 1;
				}
				else if ((tmp == "posterior") || (tmp == "post"))	{
					rateprior = 0;
				}
				else	{
					cerr << "error after ppredrate: should be prior or posterior\n";
					throw(0);
				}
			}
			else if (s == "-ppredprofile")	{
				i++;
				string tmp = argv[i];
				if (tmp == "prior")	{
					profileprior = 1;
				}
				else if ((tmp == "posterior") || (tmp == "post"))	{
					profileprior = 0;
				}
				else	{
					cerr << "error after ppredprofile: should be prior or posterior\n";
					throw(0);
				}
			}
			else if (s == "-ppredroot")	{
				i++;
				string tmp = argv[i];
				if (tmp == "prior")	{
					rootprior = 1;
				}
				else if ((tmp == "posterior") || (tmp == "post"))	{
					rootprior = 0;
				}
				else	{
					cerr << "error after ppredroot: should be prior or posterior\n";
					throw(0);
				}
			}
			else if (s == "-bf")	{
				bf = 1;
				i++;
				bfprop = atof(argv[i]);
			}
			else if (s == "-bl")	{
				bf = 2;
				i++;
				bfprop = atof(argv[i]);
			}
			else if (s == "-ibl")	{
				bf = 3;
				i++;
				bfprop = atof(argv[i]);
			}
			else if (s == "-sis")	{
				sis = 1;
				i++;
				sisprop = atof(argv[i]);
			}
			else if (s == "-bfl")	{
				bfl = 1;
				i++;
				taxon1 = argv[i];
				i++;
				taxon2 = argv[i];
				i++;
				taxon3 = argv[i];
				i++;
				taxon4 = argv[i];
				i++;
				toponfrac= atoi(argv[i]);
				i++;
				toponstep = atoi(argv[i]);
				i++;
				blfactor = atof(argv[i]);
			}
			else if (s == "-blfile")	{
				i++;
				blfile = argv[i];
			}
			else if (s == "-blfactor")	{
				i++;
				blfactor = atof(argv[i]);
			}
			/*
			else if (s == "-bf")	{
				bf = 1;
				i++;
				taxon1 = argv[i];
				i++;
				taxon2 = argv[i];
				i++;
				taxon3 = argv[i];
				i++;
				taxon4 = argv[i];
				i++;
				toponfrac= atoi(argv[i]);
				i++;
				toponstep = atoi(argv[i]);
			}
			*/
			else if (s == "-sumcomp")	{
				i++;
				sumcomp = atoi(argv[i]);
			}
			else if (s == "-fullsumcomp")	{
				sumcomp = -1;
			}
			else if (s == "+tmpbl")	{
				temperedbl = 1;
			}
			else if (s == "-tmpbl")	{
				temperedbl = 0;
			}
			else if (s == "+tmprate")	{
				temperedrate = 1;
			}
			else if (s == "-tmprate")	{
				temperedrate = 0;
			}
			else if (s == "+tmpprofile")	{
				temperedgene = 1;
			}
			else if (s == "-tmpprofile")	{
				temperedgene = 0;
			}
			else if (s == "-sitelogl")	{
				sitelogl = 1;
			}
			else if (s == "-r")	{
				rates = 1;
			}
			else if (s == "-map")	{
				map = 1;
			}
			else if (s == "-cv")	{
				cv = 1;
				i++;
				testdatafile = argv[i];
			}
			else if ( (s == "-x") || (s == "-extract") )	{
				i++;
				if (i == argc) throw(0);
				burnin = atoi(argv[i]);
				i++;
				if (i == argc) throw(0);
				s = argv[i];
				if (IsInt(s))	{
					every = atoi(argv[i]);
					i++;
					if (i == argc) throw(0);
					string tmp = argv[i];
					if (IsInt(tmp))	{
						until = atoi(argv[i]);
					}
					else	{
						i--;
					}
				}
				else {
					i--;
				}
			}
			else	{
				if (i != (argc -1))	{
					throw(0);
				}
				name = argv[i];
			}
			i++;
		}
	}
	catch(...)	{
		cerr << "error in command\n";
		cerr << '\n';
		exit(1);
	}

	if (until == -1)	{
		until = GetSize();
	}
	if (burnin == -1)	{
		burnin = GetSize() / 5;
	}

	if ((GetNprocs() == 1) && (ppred || cv || sitelogl))	{
		cerr << "error : should run readpb_mpi in mpi mode, with at least 2 processes\n";
		MPI_Finalize();
		exit(1);
	}

	if (cv)	{
		ReadCV(testdatafile,name,burnin,every,until);
	}
	else if (sis)	{
		FastReadSIS(name,burnin,every,until,sisprop);
	}
	else if (bf == 1)	{
		sumovercomponents = sumcomp;
		if (sumcomp)	{
			GlobalActivateSumOverComponents();
			ReadTopoBF(name,burnin,every,until,bfprop);
		}
		else	{
			FastReadTopoBF(name,burnin,every,until,bfprop);
		}
	}
	else if (bf == 2)	{
		FastReadTopoBL(name,burnin,every,until,bfprop);
	}
	else if (bf == 3)	{
		ReadTopoBL(name,burnin,every,until,bfprop);
	}
	else if (bfl)	{
		ReadTopoBL(name,burnin,every,until,taxon1,taxon2,taxon3,taxon4,toponfrac,toponstep);
	}
	/*
	else if (bf)	{
		SetTemperedBL(temperedbl);
		SetTemperedGene(temperedgene);
		SetTemperedRate(temperedrate);
		sumovercomponents = sumcomp;
		if (sumcomp > 0)	{
			GlobalActivateSumOverComponents();
		}
		ReadTopoBF(name,burnin,every,until,taxon1,taxon2,taxon3,taxon4,toponfrac,toponstep);
	}
	*/
	else if (sitelogl)	{
		ReadSiteLogL(name,burnin,every,until);
	}
	else if (rates)	{
		ReadSiteRates(name,burnin,every,until);
	}
	else if (ppred)	{
		PostPred(ppred,name,burnin,every,until,rateprior,profileprior,rootprior);
	}
	else if (map)	{
		ReadMap(name,burnin,every,until);
	}
	else	{
		Read(name,burnin,every,until);
	}
}

void RASCATFiniteGammaPhyloProcess::SlaveComputeCVScore()	{

	int sitemin = GetSiteMin();
	int sitemax = GetSiteMin() + testsitemax - testsitemin;
	double** sitelogl = new double*[GetNsite()];
	for (int i=sitemin; i<sitemax; i++)	{
		sitelogl[i] = new double[GetNcomponent()];
	}
	
	// UpdateMatrices();

	for (int k=0; k<GetNcomponent(); k++)	{
		for (int i=sitemin; i<sitemax; i++)	{
			PoissonFiniteProfileProcess::alloc[i] = k;
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

void RASCATFiniteGammaPhyloProcess::SlaveComputeSiteLogL()	{

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
			PoissonFiniteProfileProcess::alloc[i] = k;
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
	}

	MPI_Send(meansitelogl,GetNsite(),MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
	
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		delete[] sitelogl[i];
	}
	delete[] sitelogl;
	delete[] meansitelogl;
}

