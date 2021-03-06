
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
#include "RASCATGammaPhyloProcess.h"
#include "Parallel.h"
#include <string>

void RASCATGammaPhyloProcess::GlobalUpdateParameters()	{

	if (GetNprocs() > 1)	{

        // ResampleWeights();
        RenormalizeProfiles();

        int nd = 2 + GetNbranch() + GetNmodeMax()*GetDim() + Nstatcomp*(GetDim()+1) + 1;
        int ni = 1 + GetNsite();
        double* dvector = new double[nd];
        int* ivector = new int[ni];

        MESSAGE signal = PARAMETER_DIFFUSION;
        MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

        // First we assemble the vector of doubles for distribution
        int index = 0;
        dvector[index] = GetAlpha();
        index++;
        dvector[index] = GetPinv();
        index++;

        for (int i=0; i<GetNbranch(); i++)  {
            dvector[index] = blarray[i];
            index++;
        }
        
        for (int i=0; i<GetNmodeMax(); i++) {
            for (int j=0; j<GetDim(); j++)  {
                dvector[index] = profile[i][j];
                index++;
            }
        }

        for (int k=0; k<Nstatcomp; k++) {
            dvector[index] = statweight[k];
            index++;
            for (int i=0; i<GetDim(); i++)	{
                dvector[index] = dirweight[k][i];
                index++;
            }
        }

        dvector[index] = kappa;
        index++;

        if (index != nd)    {
            cerr << "error in global update params: non matching dim\n";
            exit(1);
        }

        // Now the vector of ints
        ivector[0] = GetNcomponent();
        for(int i=0; i<GetNsite(); i++) {
            ivector[1+i] = DPProfileProcess::alloc[i];
        }

        // Now send out the doubles and ints over the wire...
        MPI_Bcast(ivector,ni,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(dvector,nd,MPI_DOUBLE,0,MPI_COMM_WORLD);

        delete[] dvector;
        delete[] ivector;
    }
    else	{
        UpdateZip();
    }
}

void RASCATGammaPhyloProcess::SlaveUpdateParameters()	{

    int nd = 2 + GetNbranch() + GetNmodeMax()*GetDim() + Nstatcomp*(GetDim()+1) + 1;
    int ni = 1 + GetNsite();
    double* dvector = new double[nd];
    int* ivector = new int[ni];

	MPI_Bcast(ivector,ni,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(dvector,nd,MPI_DOUBLE,0,MPI_COMM_WORLD);

	int index = 0;
	SetRateParams(dvector[index],dvector[index+1]);
	index += 2;

    for (int i=0; i<GetNbranch(); i++)  {
		blarray[i] = dvector[index];
		index++;
	}

    for (int i=0; i<GetNmodeMax(); i++) {
        for (int j=0; j<GetDim(); j++)  {
			profile[i][j] = dvector[index];
			index++;
		}
	}

    for (int k=0; k<Nstatcomp; k++) {
        statweight[k] = dvector[index];
        index++;
        for (int i=0; i<GetDim(); i++)	{
            dirweight[k][i] = dvector[index];
            index++;
        }
    }

	kappa = dvector[index];
	index++;

    if (index != nd)    {
        cerr << "error in slave update params: non matching dim\n";
        exit(1);
    }

	Ncomponent = ivector[0];
	for(int i=0; i<GetNsite(); i++) {
		DPProfileProcess::alloc[i] = ivector[1+i];
	}

	UpdateZip();

	delete[] dvector;
	delete[] ivector;
}

void RASCATGammaPhyloProcess::SlaveExecute(MESSAGE signal)	{

	switch(signal) {

	case UPDATE_RATE:
		SlaveUpdateRateSuffStat();
		break;
	default:
		PhyloProcess::SlaveExecute(signal);
	}
}


void RASCATGammaPhyloProcess::ReadPB(int argc, char* argv[])	{

	string name = "";

	int burnin = -1;
	int every = 1;
	int until = -1;
	int ppred = 0;
	int ss = 0;
	double cialpha = 0;
	string trueprofiles = "None";
	int map = 0;
	// 1 : plain ppred (outputs simulated data)
	// 2 : diversity statistic
	// 3 : compositional statistic
	int cv = 0;
	int sitelogl = 0;
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

	int statmin = 0;

	int bfl = 0;
	roottax1 = "None";
	roottax2 = "None";

	int ssdist = 0;
	// int meandirweight = 0;
	int ndisc = 100;
	int nsample = 1000;

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
			else if (s == "-statmin")	{
				statmin = 1;
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
			else if (s == "-reroot")	{
				i++;
				roottax1 = argv[i];
				i++;
				roottax2 = argv[i];
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
			else if (s == "-cv")	{
				cv = 1;
				i++;
				testdatafile = argv[i];
			}
            /*
			else if (s == "-meandirweight")	{
				meandirweight = 1;
			}
            */
			else if (s == "-ssdist")	{
				ssdist = 1;
				i++;
				ndisc = atoi(argv[i]);
				i++;
				cialpha = atof(argv[i]);
				i++;
				nsample = atoi(argv[i]);
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
			else if (s == "-map")	{
				map = 1;
			}
			else if (s == "-o")	{
				i++;
				outputname = argv[i];
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

	if (outputname == "")	{
		outputname = name;
	}

	if (ss)	{
		ReadSiteProfiles(name,burnin,every,until,cialpha,trueprofiles);
	}
	else if (ssdist)	{
		ReadProfileDistribution(name,burnin,every,until,ndisc,cialpha,nsample);
	}
    /*
	else if (meandirweight)	{
		ReadMeanDirWeight(name,burnin,every,until);
	}
    */
	else if (statmin)	{
		ReadStatMin(name,burnin,every,until);
	}
	else if (map)	{
		ReadMap(name,burnin,every,until);
	}
	else if (cv)	{
		ReadCV(testdatafile,name,burnin,every,until);
	}
	else if (sis)	{
		FastReadSIS(name,burnin,every,until,sisprop);
	}
	else if (bf == 1)	{
		FastReadTopoBF(name,burnin,every,until,bfprop);
		/*
		sumovercomponents = sumcomp;
		if (sumcomp)	{
			GlobalActivateSumOverComponents();
			ReadTopoBF(name,burnin,every,until,bfprop);
		}
		else	{
			FastReadTopoBF(name,burnin,every,until,bfprop);
		}
		*/
	}
	else if (bf == 2)	{
		FastReadTopoBL(name,burnin,every,until,bfprop);
	}
	/*
	else if (bf == 3)	{
		ReadTopoBL(name,burnin,every,until,bfprop);
	}
	*/
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
	else	{
		Read(name,burnin,every,until);
	}
}

/*
void RASCATGammaPhyloProcess::ReadMeanDirWeight(string name, int burnin, int every, int until)	{

	double* meandirweight = new double[GetDim()];
	for (int k=0; k<GetDim(); k++)	{
		meandirweight[k] = 0;
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

		for (int k=0; k<GetDim(); k++)	{
			meandirweight[k] += dirweight[k];
		}
		int nrep = 1;
		while ((i<until) && (nrep < every))	{
			FromStream(is);
			i++;
			nrep++;
		}
	}
	cerr << '\n';
	
	ofstream dos((name + ".dirweight").c_str());
	double total = 0;
	for (int k=0; k<GetDim(); k++)	{
		meandirweight[k] /= samplesize;
		total += meandirweight[k];
	}
	dos << total << '\n';
	for (int k=0; k<GetDim(); k++)	{
		meandirweight[k] /= total;
		dos << meandirweight[k] << '\t';
	}
	dos << '\n';
}
*/

void RASCATGammaPhyloProcess::ReadProfileDistribution(string name, int burnin, int every, int until, int ndisc, double cialpha, int nsample)	{

	vector<double>* hidist = 0;
	vector<double>** dist = 0;

	hidist = new vector<double>[GetNsite()];
	dist = new vector<double>*[GetNsite()];
	for (int i=0; i<GetNsite(); i++)	{
		dist[i] = new vector<double>[GetDim()];
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

		for (int i=0; i<GetNsite(); i++)	{
			double* p = GetProfile(i);
			double tmphi = GetHI7(p);
			for (int k=0; k<GetDim(); k++)	{
				dist[i][k].push_back(p[k]);
			}
			hidist[i].push_back(tmphi);
		}
		int nrep = 1;
		while ((i<until) && (nrep < every))	{
			FromStream(is);
			i++;
			nrep++;
		}
	}
	cerr << '\n';
	
	ofstream hos((name + ".samplehidist").c_str());
	for (int i=0; i<GetNsite(); i++)	{
		for (int j=0; j<samplesize; j++)	{
			hos << hidist[i][j] << '\t';
		}
		hos << '\n';
	}
}

void RASCATGammaPhyloProcess::ReadStatMin(string name, int burnin, int every, int until)	{

	cerr << "in RASCATGammaPhyloProcess::ReadStatMin\n";
	exit(1);
}

