
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

#include "RASCATGTRSBDPGammaPhyloProcess.h"
#include "Parallel.h"
#include <string>
#include <list>
#include "TexTab.h"

void RASCATGTRSBDPGammaPhyloProcess::GlobalUpdateParameters()	{

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

	int i,j,nrr,nbranch = GetNbranch(),ni,nd,L1,L2;
	nrr = GetNrr();
	L1 = GetNmodeMax();
	L2 = GetDim();
	nd = 1 + nbranch + nrr + L1*L2 + GetDim() + 1;
	ni = 1 + GetNsite();
	int ivector[ni];
	double dvector[nd]; 
	MESSAGE signal = PARAMETER_DIFFUSION;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	// GlobalBroadcastTree();

	// First we assemble the vector of doubles for distribution
	int index = 0;
	dvector[index] = GetAlpha();
	index++;
	
	for(i=0; i<nbranch; ++i) {
		dvector[index] = blarray[i];
		index++;
	}
	
	for(i=0; i<nrr ; ++i) {
		dvector[index] = rr[i];
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
	dvector[index] = kappa;
	index++;


	// Now the vector of ints
	ivector[0] = GetNcomponent();
	for(i=0; i<GetNsite(); ++i) {
		ivector[1+i] = SBDPProfileProcess::alloc[i];
	}

	// Now send out the doubles and ints over the wire...
	MPI_Bcast(ivector,ni,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(dvector,nd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(V,GetNcomponent(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(weight,GetNcomponent(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	}
	else	{
		UpdateMatrices();
	}
}


void RASCATGTRSBDPGammaPhyloProcess::SlaveExecute(MESSAGE signal)	{

	switch(signal) {

	case UPDATE_RATE:
		SlaveUpdateRateSuffStat();
		break;
	case UPDATE_RRATE:
		SlaveUpdateRRSuffStat();
		break;
	case PROFILE_MOVE:
		SlaveMoveProfile();
		break;
	case MIX_MOVE:
		SlaveMixMove();
		break;
	default:
		PhyloProcess::SlaveExecute(signal);
	}
}

void RASCATGTRSBDPGammaPhyloProcess::SlaveUpdateParameters()	{

	// SlaveBroadcastTree();

	int i,j,L1,L2,ni,nd,nbranch = GetNbranch(),nrr = GetNrr();
	L1 = GetNmodeMax();
	L2 = GetDim();
	nd = 1 + nbranch + nrr + L1*L2 + GetDim() + 1;
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

	for(i=0; i<nrr; ++i) {
		rr[i] = dvector[index];
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
	kappa = dvector[index];
	index++;

	Ncomponent = ivector[0];
	for(i=0; i<GetNsite(); ++i) {
		SBDPProfileProcess::alloc[i] = ivector[1+i];
	}
	delete[] dvector;
	delete[] ivector;

	MPI_Bcast(V,GetNcomponent(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(weight,GetNcomponent(),MPI_DOUBLE,0,MPI_COMM_WORLD);

	UpdateMatrices();
}


void RASCATGTRSBDPGammaPhyloProcess::ReadPB(int argc, char* argv[])	{

	string name = "";

	int burnin = -1;
	int every = 1;
	int until = -1;
	int ppred = 0;
	// 1 : plain ppred (outputs simulated data)
	// 2 : diversity statistic
	// 3 : compositional statistic

	int ss = 0;
	int nocc = 0;
	int cv = 0;
	int sitelogl = 0;
	int rr = 0;
	int rates = 0;
	int map = 0;
	string testdatafile = "";
	int testprofile = 0;
	double tuning = 1;

	int rateprior = 0;
	int profileprior = 0;
	int rootprior = 1;

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
			else if (s == "-sitelogl")	{
				sitelogl = 1;
			}
			else if (s == "-cv")	{
				cv = 1;
				i++;
				testdatafile = argv[i];
			}
			else if (s == "-testprof")	{
				i++;
				testprofile = atoi(argv[i]);
				i++;
				tuning = atof(argv[i]);
			}
			else if (s == "-ss")	{
				ss = 1;
			}
			else if (s == "-rr")	{
				rr = 1;
			}
			else if (s == "-r")	{
				rates = 1;
			}
			else if (s == "-map")	{
				map = 1;
			}
			else if (s == "-m")	{
				nocc = 1;
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

	if (nocc)	{
		ReadNocc(name,burnin,every,until);
	}
	else if (testprofile)	{
		ReadTestProfile(name,tuning,testprofile,burnin,every,until);
	}
	else if (cv)	{
		ReadCV(testdatafile,name,burnin,every,until);
	}
	else if (rates)	{
		ReadSiteRates(name,burnin,every,until);
	}
	else if (sitelogl)	{
		ReadSiteLogL(name,burnin,every,until);
	}
	else if (ss)	{
		ReadSiteProfiles(name,burnin,every,until);
	}
	else if (rr)	{
		ReadRelRates(name,burnin,every,until);
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

void RASCATGTRSBDPGammaPhyloProcess::ReadRelRates(string name, int burnin, int every, int until)	{

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

	double* meanrr = new double[GetNrr()];
	for (int k=0; k<GetNrr(); k++)	{
		meanrr[k] = 0;
	}

	while (i < until)	{
		cerr << ".";
		cerr.flush();
		samplesize++;
		FromStream(is);
		i++;

		for (int k=0; k<GetNrr(); k++)	{
			meanrr[k] += rr[k];
		}

		int nrep = 1;
		while ((i<until) && (nrep < every))	{
			FromStream(is);
			i++;
			nrep++;
		}
	}
	cerr << '\n';
	for (int k=0; k<GetNrr(); k++)	{
		meanrr[k] /= samplesize;
	}
	ofstream os((name + ".meanrr").c_str());
	for (int k=0; k<GetDim(); k++)	{
		os << GetStateSpace()->GetState(k) << ' ';
	}
	os << '\n';
	os << '\n';
	for (int i=0; i<GetDim(); i++)	{
		for (int j=i+1; j<GetDim(); j++)	{
			os << GetStateSpace()->GetState(i) << '\t' << GetStateSpace()->GetState(j) << '\t' << meanrr[rrindex(i,j,GetDim())] << '\n';
		}
	}
	cerr << "mean relative exchangeabilities in " << name << ".meanrr\n";

}

void RASCATGTRSBDPGammaPhyloProcess::ReadNocc(string name, int burnin, int every, int until)	{

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

	list<double> nocclist;

	while (i < until)	{
		cerr << ".";
		cerr.flush();
		samplesize++;
		FromStream(is);
		i++;

		UpdateOccupancyNumbers();
		double nocc = GetNOccupiedComponent();
		nocclist.push_back(nocc);

		int nrep = 1;
		while ((i<until) && (nrep < every))	{
			FromStream(is);
			i++;
			nrep++;
		}
	}
	cerr << '\n';
	cerr << "number of occupied component (credibility interval) : ";
	printCI(nocclist,cerr);
	cerr << '\n';
}

void RASCATGTRSBDPGammaPhyloProcess::SlaveComputeCVScore()	{

	if (! SumOverRateAllocations())	{
		cerr << "rate error\n";
		exit(1);
	}

	int sitemin = GetSiteMin();
	int sitemax = GetSiteMin() + testsitemax - testsitemin;
	double** sitelogl = new double*[GetNsite()];
	for (int i=sitemin; i<sitemax; i++)	{
		sitelogl[i] = new double[GetNcomponent()];
	}
	
	// UpdateMatrices();

	for (int k=0; k<GetNcomponent(); k++)	{
		for (int i=sitemin; i<sitemax; i++)	{
			ExpoConjugateGTRSBDPProfileProcess::alloc[i] = k;
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

void RASCATGTRSBDPGammaPhyloProcess::SlaveComputeSiteLogL()	{

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
			ExpoConjugateGTRSBDPProfileProcess::alloc[i] = k;
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

void RASCATGTRSBDPGammaPhyloProcess::ReadTestProfile(string name, int nrep, double tuning, int burnin, int every, int until)	{

	proposemode = 1;
	
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

	ofstream os((name + ".testprofile").c_str());

	while (i < until)	{
		cerr << ".";
		cerr.flush();
		samplesize++;
		FromStream(is);
		i++;

		// quick update and mapping on the fly
		MPI_Status stat;
		MESSAGE signal = BCAST_TREE;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
		GlobalBroadcastTree();
		GlobalUpdateConditionalLikelihoods();
		GlobalCollapse();

		// collect suff stats
		GlobalUpdateSiteProfileSuffStat();
		UpdateModeProfileSuffStat();
		UpdateOccupancyNumbers();

		// for each component
		// propose try new profiles
		// calculate the acceptance probability
		// as if in a gibbs
		// average over nrep replicates
		for (int k=0; k<Ncomponent; k++)	{
			double a = 0;
			for (int rep=0; rep<nrep; rep++)	{
				double logratio = - ProfileSuffStatLogProb(k) - LogStatPrior(k);
				double logh = ProfileProposeMove(profile[k],tuning,0,0,k,0);
				logratio += ProfileSuffStatLogProb(k) + LogStatPrior(k);
				logratio += logh;
				if (logratio < 0)	{
					a += exp(logratio);
				}
			}
			a /= nrep;
			os << occupancy[k] << '\t' << a << '\n';
		}

		int nrep = 1;
		while ((i<until) && (nrep < every))	{
			FromStream(is);
			i++;
			nrep++;
		}
	}
	cerr << '\n';

	// delete arrays

}
