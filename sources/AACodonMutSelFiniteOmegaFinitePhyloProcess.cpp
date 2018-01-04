
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
#include "AACodonMutSelFiniteOmegaFinitePhyloProcess.h"
#include "Parallel.h"
#include <string.h>


// MPI: these two functions are responsible for broadcasting/receiving the current state of the parameter vector
// are model dependent
// should be implemented in .cpp file
void AACodonMutSelFiniteOmegaFinitePhyloProcess::SlaveUpdateParameters()	{

	int i,j,L1,L2,ni,nd,nbranch = GetNbranch(),nnucrr = GetNnucrr(),nnucstat = 4, nomega = GetNomega();
	L1 = GetNmodeMax();
	L2 = GetDim();
	int nstate = GetData()->GetNstate();
	//nd = nbranch + nnucrr + nnucstat + L2 + L1*(L2+1) + nstate + 1; // check if these last terms are correct in this context...
	nd = nbranch + nnucrr + nnucstat + L2 + L1*(L2+1) + nstate + nomega*2; // check if these last terms are correct in this context...
	//ni = 1 + ProfileProcess::GetNsite();
	ni = 1 + ProfileProcess::GetNsite() * 2;
	int* ivector = new int[ni];
	double* dvector = new double[nd];
	MPI_Bcast(ivector,ni,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(dvector,nd,MPI_DOUBLE,0,MPI_COMM_WORLD);

	int index = 0;
	for(i=0; i<nbranch; ++i) {
		blarray[i] = dvector[index];
		index++;
	}
	for(i=0; i<nnucrr; ++i) {
		nucrr[i] = dvector[index];
		index++;
	}
	for(i=0; i<nnucstat; ++i) {
		nucstat[i] = dvector[index];
		index++;
	}
	for (i=0; i<nstate; i++)	{
		codonprofile[i] = dvector[index];
		index++;
	}
	for(i=0; i<L1; ++i) {
		for(j=0; j<L2; ++j) {
			profile[i][j] = dvector[index];
			index++;
		}
		weight[i] = dvector[index];
		index++;
	}
	for (int i=0; i<GetDim(); i++)	{
		dirweight[i] = dvector[index];
		index++;
	}
	for (int i=0; i<nomega; i++)	{
		omega[i] = dvector[index];
		index++;
		omegaweight[i] = dvector[index];
		index++;
	}
	int iindex = 0;
	Ncomponent = ivector[iindex];
	iindex++;
	for(i=0; i<ProfileProcess::GetNsite(); ++i) {
		FiniteProfileProcess::alloc[i] = ivector[iindex];
		iindex++;
		FiniteOmegaProcess::omegaalloc[i] = ivector[iindex];
		iindex++; 
	}
	//GetBranchLengthsFromArray();
	delete[] dvector;
	delete[] ivector;
	// this one is really important
	// in those cases where new components have appeared, or some old ones have disappeared
	// during allocation move on the master node.
	// 
	// note that CreateMatrices() in fact creates only those that are not yet allocated
	// CreateMatrices();
	UpdateMatrices();
}


void AACodonMutSelFiniteOmegaFinitePhyloProcess::SlaveExecute(MESSAGE signal)	{

	switch(signal) {

		case UPDATE_SITEOMEGA:
			SlaveUpdateSiteOmegaSuffStat();
			break;
		case UPDATE_OMEGA:
			SlaveUpdateOmegaSuffStat();
			break;
		case REALLOCOMEGA_MOVE:
			SlaveOmegaIncrementalFiniteMove();
			break;
		case STATFIX:
			SlaveGetStatFix();
			break;
		case REALLOC_MOVE:
			SlaveIncrementalFiniteMove();
			break;
		case PROFILE_MOVE:
			SlaveMoveProfile();
			break;
		default:
			PhyloProcess::SlaveExecute(signal);
	}
}

void AACodonMutSelFiniteOmegaFinitePhyloProcess::GlobalUpdateParameters() {

	if (GetNprocs() > 1)	{
	// MPI2
	// should send the slaves the relevant information
	// about model parameters
	// for this model, should broadcast
	// (but should first call PutBranchLengthsIntoArray())
	// 
	// upon receiving this information
	// slave should 
	// store it in the local copies of the variables
	// and then call
	// SetBranchLengthsFromArray()
	int i,j,nnucrr,nnucstat,nbranch = GetNbranch(),ni,nd,L1,L2,nomega;
	nnucrr = GetNnucrr();
	nnucstat = 4;
	L1 = GetNmodeMax();
	L2 = GetDim();
	int nstate = GetData()->GetNstate();
	nomega = GetNomega();	
	//nd = nbranch + nnucrr + nnucstat + L2 + L1*(L2+1) + nstate + 1;  // check if these last terms are correct in this context...
	nd = nbranch + nnucrr + nnucstat + L2 + L1*(L2+1) + nstate + nomega*2;  // check if these last terms are correct in this context...
	ni = 1 + ProfileProcess::GetNsite() * 2; // 1 for the number of componenets, and the rest for allocations
	int ivector[ni];
	double dvector[nd]; 
	MESSAGE signal = PARAMETER_DIFFUSION;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	
	int index = 0;
	// First we assemble the vector of doubles for distribution
	for(i=0; i<nbranch; ++i) {
		dvector[index] = blarray[i];
		index++;
	}

	for(i=0; i<nnucrr; ++i) {
		dvector[index] = nucrr[i];
		index++;
	}
	for(i=0; i<nnucstat; ++i) {
		dvector[index] = nucstat[i];
		index++;
	}

	for (i=0; i<nstate; i++)	{
		dvector[index] = codonprofile[i];
		index++;
	}

	for(i=0; i<L1; ++i) {
		for(j=0; j<L2; ++j) {
			dvector[index] = profile[i][j];
			index++;
		}
		dvector[index] = weight[i];
		index++;
	}
	for (int i=0; i<GetDim(); i++)	{
		dvector[index] = dirweight[i];
		index++;
	}
	for (int i=0; i<nomega; i++)	{
		dvector[index] = omega[i];
		index++;
		dvector[index] = omegaweight[i];
		index++;
	}

	// Now the vector of ints
	int iindex = 0;
	ivector[iindex] = GetNcomponent();
	iindex++;
	for(i=0; i<ProfileProcess::GetNsite(); ++i) {
		ivector[iindex] = FiniteProfileProcess::alloc[i];
		iindex++;
		ivector[iindex] = FiniteOmegaProcess::omegaalloc[i];
		iindex++;
	}

	// Now send out the doubles and ints over the wire...
	MPI_Bcast(ivector,ni,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(dvector,nd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	}
	// else	{
		UpdateMatrices();
	// }
	
}


void AACodonMutSelFiniteOmegaFinitePhyloProcess::SlaveComputeCVScore()	{

	int sitemin = GetSiteMin();
	int sitemax = GetSiteMin() + testsitemax - testsitemin;
	double** sitelogl = new double*[ProfileProcess::GetNsite()];
	for (int i=sitemin; i<sitemax; i++)	{
		sitelogl[i] = new double[GetNcomponent()];
	}
	
	for (int k=0; k<GetNcomponent(); k++)	{
		for (int i=sitemin; i<sitemax; i++)	{
			AACodonMutSelFiniteOmegaFiniteProfileProcess::alloc[i] = k;
		}
		UpdateComponent(k);
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

void AACodonMutSelFiniteOmegaFinitePhyloProcess::ReadPB(int argc, char* argv[])	{


	// Needs updating!

	string name = "";

	int burnin = 0;
	int every = 1;
	int until = -1;
	int ppred = 0;
	// 1 : plain ppred (outputs simulated data)
	// 2 : diversity statistic
	// 3 : compositional statistic

	int cv = 0;
	int sel = 0;
	int map = 0;
	string testdatafile = "";

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
			if (s == "-sel")	{
				sel = 1;
			}
			else if (s == "-map")	{
				map = 1;
			}
			else if (s == "-cv")	{
				cv = 1;
				i++;
				testdatafile = argv[i];
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

	if (map)	{
		ReadMap(name,burnin,every,until);
	}
	else if (cv)	{
		ReadCV(testdatafile,name,burnin,every,until,1,codetype);
	}
	//if (sel)	{
	//	ReadSDistributions(name,burnin,every,until);
	//}
	//else if (ppred)	{
	//	PostPred(ppred,name,burnin,every,until);
	//}
	else	{
		Read(name,burnin,every,until);
	}
}

void AACodonMutSelFiniteOmegaFinitePhyloProcess::Read(string name, int burnin, int every, int until)	{

	ifstream is((name + ".chain").c_str());
	if (! is)	{
		cerr << "error: did not find " << name << ".chain\n";
		exit(1);
	}

	int* siteomegagreaterthanonecount = new int[ProfileProcess::GetNsite()];
	for (int site = 0; site < ProfileProcess::GetNsite(); site++)        {
		siteomegagreaterthanonecount[site] = 0;
	}

	cerr << "burnin : " << burnin << "\n";
	cerr << "until : " << until << '\n';
	int i=0;
	while ((i < until) && (i < burnin))	{
		cerr << ".";
		cerr.flush();
		FromStream(is);
		i++;
	}
	cerr << "\nburnin complete\n";
	cerr.flush();
	int samplesize = 0;
	while (i < until)	{
		cerr << ".";
		cerr.flush();
		samplesize++;
		FromStream(is);
		i++;
		QuickUpdate();
		//UpdateMatrices();
		for (int site = 0; site < ProfileProcess::GetNsite(); site++)        {
			siteomegagreaterthanonecount[site] = 0;
			if (MixtureOmegaProcess::GetSiteOmega(site) > 1.0)	{
				siteomegagreaterthanonecount[site]++;
			}
		}
	}

	ofstream pomegagtone_os( (name + ".ppomegagtone").c_str(), std::ios::out);
	for (int site = 0; site < ProfileProcess::GetNsite(); site++)	{
		pomegagtone_os << site+1 << "\t" << (double) (siteomegagreaterthanonecount[site])/samplesize << "\n";
	}
	delete siteomegagreaterthanonecount;
}