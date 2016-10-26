
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
#include "AACodonMutSelSBDPOmegaSBDPPhyloProcess.h"
#include "Parallel.h"
#include <string.h>


// MPI: these two functions are responsible for broadcasting/receiving the current state of the parameter vector
// are model dependent
// should be implemented in .cpp file
void AACodonMutSelSBDPOmegaSBDPPhyloProcess::SlaveUpdateParameters()	{

	// SlaveBroadcastTree();

	int i,j,L1,L2,ni,nd,nbranch = GetNbranch(),nnucrr = GetNnucrr(),nnucstat = 4;
	L1 = GetNmodeMax();
	L2 = GetDim();
	int nstate = GetData()->GetNstate();

	// 2 for branchalpha branchbeta
	// blarray
	// nucrr
	// nucstat
	// profiles
	// dirweight
	// 1 for kappa
	// codonprofile
	// 1 for omega
	nd = nbranch + nnucrr + nnucstat + L1*L2 + GetDim() + 1 + nstate + ProfileProcess::GetNsite() + 2;
	ni = 1 + ProfileProcess::GetNsite();
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
	}
	for (int i=0; i<GetDim(); i++)	{
		dirweight[i] = dvector[index];
		index++;
	}
	kappa = dvector[index];
	index++;
	for(int i=0; i<ProfileProcess::GetNsite(); ++i) {
		omegaarray[i] = dvector[index];
		index++;
	}
	omegaalpha = dvector[index];
	index++;
	omegabeta = dvector[index];
	index++;
	
	Ncomponent = ivector[0];
	for(i=0; i<ProfileProcess::GetNsite(); ++i) {
		SBDPProfileProcess::alloc[i] = ivector[1+i];
	}
	//GetBranchLengthsFromArray();
	delete[] dvector;
	delete[] ivector;
	// this one is really important
	// in those cases where new components have appeared, or some old ones have disappeared
	// during allocation move on the master node.
	// 
	// note that CreateMatrices() in fact creates only those that are not yet allocated
	// and also deletes those that are now obsolete
	// CreateMatrices();
	UpdateMatrices();
}


void AACodonMutSelSBDPOmegaSBDPPhyloProcess::SlaveExecute(MESSAGE signal)	{

	switch(signal) {

		case PROFILE_MOVE:
			SlaveMoveProfile();
			break;
		case MIX_MOVE:
			SlaveMixMove();
			break;
		case NONSYNMAPPING:
			SlaveNonSynMapping();
			break;
		case UPDATE_SITEOMEGA:
			SlaveUpdateSiteOmegaSuffStat();
			break;
		case UPDATE_OMEGA:
			SlaveUpdateOmegaSuffStat();
			break;
		case MIXMOVEOMEGA:
			SlaveMixMoveOmega();
			break;
		//case COLLECTOMEGASUFFSTAT:
		//	SlaveCollectSiteOmegaSuffStats();
		//	break;
		default:
			PhyloProcess::SlaveExecute(signal);
	}
}

void AACodonMutSelSBDPOmegaSBDPPhyloProcess::GlobalUpdateParameters() {

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
	int i,j,nnucrr,nnucstat,nbranch = GetNbranch(),ni,nd,L1,L2;
	nnucrr = GetNnucrr();
	nnucstat = 4;	
	L1 = GetNmodeMax();
	L2 = GetDim();
	int nstate = GetData()->GetNstate();
	//nd = 2 + nbranch + nnucrr + + nnucstat + L1*L2 + GetDim() + 1;
	nd = nbranch + nnucrr + + nnucstat + L1*L2 + GetDim() + 1 + nstate + ProfileProcess::GetNsite() + 2;
	ni = 1 + ProfileProcess::GetNsite(); // 1 for the number of componenets, and the rest for allocations
	int ivector[ni];
	double dvector[nd]; 
	MESSAGE signal = PARAMETER_DIFFUSION;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	
	// GlobalBroadcastTree();
	// First we assemble the vector of doubles for distribution
	int index = 0;
	
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
	}
	for (int i=0; i<GetDim(); i++)	{
		dvector[index] = dirweight[i];
		index++;
	}
	dvector[index] = kappa;
	index++;

	for(int i=0; i<ProfileProcess::GetNsite(); ++i) {
		dvector[index] = omegaarray[i];
		index++;
	}
	dvector[index] = omegaalpha;
	index++;
	dvector[index] = omegabeta;
	index++;

	// Now the vector of ints
	ivector[0] = GetNcomponent();
	for(i=0; i<ProfileProcess::GetNsite(); ++i) {
		ivector[1+i] = SBDPProfileProcess::alloc[i];
	}

	// Now send out the doubles and ints over the wire...
	MPI_Bcast(ivector,ni,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(dvector,nd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	}
}

void AACodonMutSelSBDPOmegaSBDPPhyloProcess::SlaveComputeCVScore()	{

	int sitemin = GetSiteMin();
	int sitemax = GetSiteMin() + testsitemax - testsitemin;
	double** sitelogl = new double*[ProfileProcess::GetNsite()];
	for (int i=sitemin; i<sitemax; i++)	{
		sitelogl[i] = new double[GetNcomponent()];
	}
	
	for (int k=0; k<GetNcomponent(); k++)	{
		for (int i=sitemin; i<sitemax; i++)	{
			AACodonMutSelSBDPOmegaSBDPProfileProcess::alloc[i] = k;
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

void AACodonMutSelSBDPOmegaSBDPPhyloProcess::ReadPB(int argc, char* argv[])	{

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
	int mapstats = 0;
	string testdatafile = "";

	int rateprior = 0;
	int profileprior = 0;
	int rootprior = 0;

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
			else if (s == "-mapstats")	{
				mapstats = 1;
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
			else if (s == "-div")	{
				ppred = 2;
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
	else if (mapstats)	{
		ReadMapStats(name,burnin,every,until);
	}
	else if (ppred)	{
		PostPred(ppred,name,burnin,every,until,rateprior,profileprior,rootprior);
	}
	else	{
		Read(name,burnin,every,until);
	}
}

void AACodonMutSelSBDPOmegaSBDPPhyloProcess::ReadMapStats(string name, int burnin, int every, int until){
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

	
	ofstream ospost((name + ".nonsynpost").c_str());
	ofstream ospred((name + ".nonsynpred").c_str());
	ofstream ospvalue((name + ".nonsynpvalue").c_str());

	int samplesize = 0;
	int pvalue=0;
	int obs, pred;
	while (i < until)	{
		cerr << ".";
		cerr.flush();
		samplesize++;
		FromStream(is);
		i++;

		MPI_Status stat;
		MESSAGE signal = BCAST_TREE;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
		GlobalBroadcastTree();
		GlobalUpdateConditionalLikelihoods();
		GlobalCollapse();

		GlobalUpdateSiteProfileSuffStat();

		// write posterior
		obs = GlobalNonSynMapping();
		ospost << (double) (obs) / AACodonMutSelProfileProcess::GetNsite() << "\n";
		cerr << (double) (obs) / AACodonMutSelProfileProcess::GetNsite() << "\t";

		GlobalUnfold();

		//Posterior Prededictive Mappings
		GlobalUnclamp();
		GlobalCollapse();
		GlobalUpdateSiteProfileSuffStat();
		GlobalSetDataFromLeaves();

		// write posterior predictive
		pred = GlobalNonSynMapping();
		ospred << (double) (pred) / AACodonMutSelProfileProcess::GetNsite() << "\n";
		cerr << (double) (pred) / AACodonMutSelProfileProcess::GetNsite() << "\n";
	
		if (pred > obs) pvalue++;

		GlobalRestoreData();
		GlobalUnfold();

		int nrep = 1;
		while ((i<until) && (nrep < every))	{
			FromStream(is);
			i++;
			nrep++;
		}
	}

	ospvalue << (double) (pvalue) / samplesize << "\n";
	cerr << '\n';
}


void AACodonMutSelSBDPOmegaSBDPPhyloProcess::Read(string name, int burnin, int every, int until)	{
}

int AACodonMutSelSBDPOmegaSBDPPhyloProcess::CountNonSynMapping()	{

	int total = 0;	
	for(int i =GetSiteMin(); i <GetSiteMax(); i++){
		total += CountNonSynMapping(i);
	}
	return total;
}

int AACodonMutSelSBDPOmegaSBDPPhyloProcess::CountNonSynMapping(int i)	{
	int count = 0;
	for(int k=0; k<GetNstate(); ++k) {
		for(int l=0; l<GetNstate(); ++l) {
			count+=sitepaircount[i][pair<int,int>(k,l)];
		}
	}
	return count;
}

int AACodonMutSelSBDPOmegaSBDPPhyloProcess::GlobalNonSynMapping()	{

	MESSAGE signal = NONSYNMAPPING;
	MPI_Status stat;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	int totalcount=0;
	for (int i=1; i<GetNprocs(); i++)	{
		int count;
		MPI_Recv(&count,1,MPI_INT,MPI_ANY_SOURCE,TAG1,MPI_COMM_WORLD, &stat);
		totalcount += count;
	}
	return totalcount;

}

void AACodonMutSelSBDPOmegaSBDPPhyloProcess::SlaveNonSynMapping()	{

	int nonsyn = CountNonSynMapping();
	MPI_Send(&nonsyn,1,MPI_INT,0,TAG1,MPI_COMM_WORLD);

}
