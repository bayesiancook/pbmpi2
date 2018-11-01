
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

#include "RASCATGTRFiniteGammaPhyloProcess.h"
#include "Parallel.h"
#include <string>


void RASCATGTRFiniteGammaPhyloProcess::GlobalUpdateParameters()	{

	if (GetNprocs() > 1)	{

        // ResampleWeights();
        RenormalizeProfiles();

        int nd = 2 + GetNbranch() + GetNrr() + GetNmodeMax()*(GetDim()+1) + Nstatcomp*(GetDim()+1) + 1;
        int ni = 1 + GetNsite();
        int* ivector = new int[ni];
        double* dvector = new double[nd]; 

        MESSAGE signal = PARAMETER_DIFFUSION;
        MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

        // GlobalBroadcastTree();

        // First we assemble the vector of doubles for distribution
        int index = 0;
        dvector[index] = GetAlpha();
        index ++;
        dvector[index] = GetPinv();
        index ++;

        for (int i=0; i<GetNbranch(); i++)  {
            dvector[index] = blarray[i];
            index++;
        }
        
        for (int i=0; i<GetNrr(); i++)  {
            dvector[index] = rr[i];
            index++;
        }

        for (int i=0; i<GetNmodeMax(); i++) {
            for (int j=0; j<GetDim(); j++)  {
                dvector[index] = profile[i][j];
                index++;
            }
            dvector[index] = weight[i];
            index++;
        }

        for (int k=0; k<Nstatcomp; k++) {
            dvector[index] = statweight[k];
            index++;
            for (int i=0; i<GetDim(); i++)	{
                dvector[index] = dirweight[k][i];
                index++;
            }
        }

        dvector[index] = statfixalpha;
        index++;

        if (index != nd)    {
            cerr << "error in globalupdateparams: non matching double vector size\n";
            exit(1);
        }

        // Now the vector of ints
        ivector[0] = GetNcomponent();
        for(int i=0; i<GetNsite(); i++) {
            ivector[1+i] = FiniteProfileProcess::alloc[i];
        }

        // Now send out the doubles and ints over the wire...
        MPI_Bcast(ivector,ni,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(dvector,nd,MPI_DOUBLE,0,MPI_COMM_WORLD);

        delete[] dvector;
        delete[] ivector;
    }
    else    {
        UpdateMatrices();
    }
}

void RASCATGTRFiniteGammaPhyloProcess::SlaveUpdateParameters()	{

	// SlaveBroadcastTree();

    int nd = 2 + GetNbranch() + GetNrr() + GetNmodeMax()*(GetDim()+1) + Nstatcomp*(GetDim()+1) + 1;
    int ni = 1 + GetNsite();
    int* ivector = new int[ni];
    double* dvector = new double[nd]; 

	MPI_Bcast(ivector,ni,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(dvector,nd,MPI_DOUBLE,0,MPI_COMM_WORLD);

	int index = 0;
	SetRateParams(dvector[index],dvector[index+1]);
	index +=2;

    for (int i=0; i<GetNbranch(); i++)  {
		blarray[i] = dvector[index];
		index++;
	}

    for (int i=0; i<GetNrr(); i++)  {
		rr[i] = dvector[index];
		index++;
	}

    for (int i=0; i<GetNmodeMax(); i++) {
        for (int j=0; j<GetDim(); j++)  {
			profile[i][j] = dvector[index];
			index++;
		}
		weight[i] = dvector[index];
		index++;
	}

    for (int k=0; k<Nstatcomp; k++) {
        statweight[k] = dvector[index];
        index++;
        for (int i=0; i<GetDim(); i++)	{
            dirweight[k][i] = dvector[index];
            index++;
        }
    }

    statfixalpha = dvector[index];
    index++;

    if (index != nd)    {
        cerr << "error in slave update params: non matching dim\n";
        exit(1);
    }

	Ncomponent = ivector[0];
	for(int i=0; i<GetNsite(); i++) {
		FiniteProfileProcess::alloc[i] = ivector[1+i];
	}

	delete[] dvector;
	delete[] ivector;

	UpdateMatrices();
}


void RASCATGTRFiniteGammaPhyloProcess::SlaveExecute(MESSAGE signal)	{

	switch(signal) {

	case UPDATE_RATE:
		SlaveUpdateRateSuffStat();
		break;
	case UPDATE_RRATE:
		SlaveUpdateRRSuffStat();
		break;
	case REALLOC_MOVE:
		SlaveIncrementalFiniteMove();
		break;
	case PROFILE_MOVE:
		SlaveMoveProfile();
		break;
	case STATFIX:
		SlaveGetStatFix();
		break;
	default:
		PhyloProcess::SlaveExecute(signal);
	}
}

void RASCATGTRFiniteGammaPhyloProcess::ReadPB(int argc, char* argv[])	{

	string name = "";

	int burnin = -1;
	int every = 1;
	int until = -1;
	int ppred = 0;
	// 1 : plain ppred (outputs simulated data)
	// 2 : diversity statistic
	// 3 : compositional statistic

	int ss = 0;
	double cialpha = 0;
	string trueprofiles = "None";
	int ms = 0;
	int cv = 0;
	int sitelogl = 0;
	int map = 0;
	int rates = 0;
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
			else if (s == "-r")	{
				rates = 1;
			}
			else if (s == "-ss")	{
				ss = 1;
			}
			else if (s == "-ci")	{
				i++;
				cialpha = atof(argv[i]);
			}
			else if (s == "-true")	{
				i++;
				trueprofiles = argv[i];
			}
			else if (s == "-ms")	{
				ms = 1;
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
	else if (sitelogl)	{
		ReadSiteLogL(name,burnin,every,until);
	}
	else if (rates)	{
		ReadSiteRates(name,burnin,every,until);
	}
	else if (ss)	{
		ReadSiteProfiles(name,burnin,every,until,cialpha,trueprofiles);
	}
	else if (ms)	{
		ReadModeProfiles(name,burnin,every,until);
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

void RASCATGTRFiniteGammaPhyloProcess::ReadModeProfiles(string name, int burnin, int every, int until)	{

	double* modeweight = new double[GetNcomponent()];
	for (int i=0; i<GetNcomponent(); i++)	{
		modeweight[i] = 0;
	}
	double** modestat = new double*[GetNcomponent()];
	for (int i=0; i<GetNcomponent(); i++)	{
		modestat[i] = new double[GetDim()];
		for (int k=0; k<GetDim(); k++)	{
			modestat[i][k] = 0;
		}
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

	while (i < until)	{
		cerr << ".";
		cerr.flush();
		samplesize++;
		FromStream(is);
		i++;

		for (int i=0; i<GetNcomponent(); i++)	{
			double* p = profile[i];
			for (int k=0; k<GetDim(); k++)	{
				modestat[i][k] += p[k];
			}
			modeweight[i] += weight[i];
		}
		int nrep = 1;
		while ((i<until) && (nrep < every))	{
			FromStream(is);
			i++;
			nrep++;
		}
	}
	cerr << '\n';
	
	ofstream os((name + ".modeprofiles").c_str());
	/*
	for (int k=0; k<GetDim(); k++)	{
		os << GetStateSpace()->GetState(k) << ' ';
	}
	os << '\n';
	os << '\n';
	*/
	os << GetNcomponent() << '\t' << GetDim() << '\n';
	for (int i=0; i<GetNcomponent(); i++)	{
		modeweight[i] /= samplesize;
		os << modeweight[i];
		// os << i+1;
		for (int k=0; k<GetDim(); k++)	{
			modestat[i][k] /= samplesize;
			os << '\t' << modestat[i][k];
		}
		os << '\n';
	}
	cerr << "mean mixture profiles in " << name << ".modeprofiles\n";
	cerr << '\n';
}

void RASCATGTRFiniteGammaPhyloProcess::SlaveComputeCVScore()	{

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
			ExpoConjugateGTRFiniteProfileProcess::alloc[i] = k;
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

void RASCATGTRFiniteGammaPhyloProcess::SlaveComputeSiteLogL()	{

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
			ExpoConjugateGTRFiniteProfileProcess::alloc[i] = k;
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

