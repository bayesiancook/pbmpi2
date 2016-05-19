
#include "PhyloProcess.h"
#include "Parallel.h"
#include "StringStreamUtils.h"
#include "TexTab.h"

void PhyloProcess::GlobalSetRatePrior(int inrateprior)	{

	rateprior = inrateprior;
	MESSAGE signal = SETRATEPRIOR;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&rateprior,1,MPI_INT,0,MPI_COMM_WORLD);
}

void PhyloProcess::SlaveSetRatePrior()	{

	MPI_Bcast(&rateprior,1,MPI_INT,0,MPI_COMM_WORLD);
}

void PhyloProcess::GlobalSetProfilePrior(int inprofileprior)	{

	profileprior = inprofileprior;
	MESSAGE signal = SETPROFILEPRIOR;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&profileprior,1,MPI_INT,0,MPI_COMM_WORLD);
}

void PhyloProcess::SlaveSetProfilePrior()	{

	MPI_Bcast(&profileprior,1,MPI_INT,0,MPI_COMM_WORLD);
}

void PhyloProcess::GlobalSetRootPrior(int inrootprior)	{

	rootprior = inrootprior;
	MESSAGE signal = SETROOTPRIOR;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&rootprior,1,MPI_INT,0,MPI_COMM_WORLD);
}

void PhyloProcess::SlaveSetRootPrior()	{

	MPI_Bcast(&rootprior,1,MPI_INT,0,MPI_COMM_WORLD);
}

void PhyloProcess::GlobalGetMeanSiteRate()	{

	if (GetNprocs() > 1)	{
		if (! meansiterate)	{
			meansiterate = new double[GetNsite()];
		}

		MPI_Status stat;
		MESSAGE signal = SITERATE;

		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

		for(int i=1; i<GetNprocs(); i++) {
			MPI_Recv(meansiterate+GetProcSiteMin(i),GetProcSiteMax(i) - GetProcSiteMin(i),MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
		}
	}
}

void PhyloProcess::SlaveSendMeanSiteRate()	{

	MPI_Send(meansiterate+GetSiteMin(),GetSiteMax()-GetSiteMin(),MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

void PhyloProcess::ReadPB(int argc, char* argv[])	{

	string name = "";

	int burnin = -1;
	int every = 1;
	int until = -1;
	int ppred = 0;
	int cv = 0;
	int sitelogl = 0;
	string testdatafile = "";
	int rateprior = 0;
	int profileprior = 0;
	int rootprior = 0;

	// 1 : plain ppred (outputs simulated data)
	// 2 : diversity statistic
	// 3 : compositional statistic

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
			else if (s == "-nsub")	{
				ppred = 4;
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
			else if (s == "-cv")	{
				cv = 1;
				i++;
				testdatafile = argv[i];
			}
			else if (s == "-sitelogl")	{
				sitelogl = 1;
			}
			else if ( (s == "-x") || (s == "-extract") )	{
				i++;
				if (i == argc) throw(0);
				burnin = atoi(argv[i]);
				i++;
				if (i == argc) throw(0);
				every = atoi(argv[i]);
				i++;
				if (i == argc) throw(0);
				string tmp = argv[i];
				if (IsFloat(tmp))	{
					until = atoi(argv[i]);
				}
				else	{
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
		MESSAGE signal = KILL;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Finalize();
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

	if (ppred)	{
		PostPred(ppred,name,burnin,every,until,rateprior,profileprior,rootprior);
	}
	else if (cv)	{
		ReadCV(testdatafile,name,burnin,every,until);
	}
	else if (sitelogl)	{
		ReadSiteLogL(name,burnin,every,until);
	}
	else	{
		Read(name,burnin,every,until);
	}
}

void PhyloProcess::Read(string name, int burnin, int every, int until)	{

	ifstream is((name + ".chain").c_str());
	if (!is)	{
		cerr << "error: no .chain file found\n";
		exit(1);
	}

	cerr << '\n';
	cerr << "burnin : " << burnin << "\n";
	cerr << "every  : " << every << '\n'; 
	cerr << "until  : " << until << '\n';
	cerr << '\n';

	cerr << "burnin\n";

	int i=0;
	while ((i < until) && (i < burnin))	{
		cerr << '.';
		FromStream(is);
		i++;
	}
	cerr << '\n';
	int samplesize = 0;

	list<double> lengthlist;
	list<double> alphalist;

	cerr << "sample\n";
	while (i < until)	{
		cerr << ".";
		cerr.flush();
		samplesize++;
		FromStream(is);
		i++;
		double alpha = GetAlpha();
		alphalist.push_back(alpha);
		double length = GetRenormTotalLength();
		lengthlist.push_back(length);
		int nrep = 1;
		while ((i<until) && (nrep < every))	{
			FromStream(is);
			i++;
			nrep++;
		}
	}
	cerr << '\n';
	cerr << '\n';
	cerr << "		   post. mean (95 % CI)\n";
	cerr << "tree length      : ";
	printCI(lengthlist, cerr);
	cerr << '\n';
	cerr << "alpha paarameter : ";
	printCI(alphalist, cerr);
	cerr << '\n';
	cerr << '\n';
}

void PhyloProcess::ReadSiteRates(string name, int burnin, int every, int until)	{

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

	double* meanrate = new double[GetNsite()];
	for (int i=0; i<GetNsite(); i++)	{
		meanrate[i] = 0;
	}

	while (i < until)	{
		cerr << ".";
		cerr.flush();
		samplesize++;
		FromStream(is);
		i++;

		QuickUpdate();

		GlobalGetMeanSiteRate();

		// double length = GetRenormTotalLength();
		for (int i=0; i<GetNsite(); i++)	{
			// meansiterate[i] *= length;
			meanrate[i] += meansiterate[i];
		}

		int nrep = 1;
		while ((i<until) && (nrep < every))	{
			FromStream(is);
			i++;
			nrep++;
		}
	}
	cerr << '\n';
	ofstream os((name + ".meansiterates").c_str());
	for (int i=0; i<GetNsite(); i++)	{
		meanrate[i] /= samplesize;
		os << i << '\t' << meanrate[i] << '\n';
	}
	cerr << "posterior mean relative site rates in " << name << ".meansiterates\n";

	delete[] meanrate;

}

void PhyloProcess::ReadTopoBF(string name, int burnin, int every, int until, string intaxon1, string intaxon2, string intaxon3, string intaxon4, int nstep)	{

	SetSpecialSPR(intaxon1,intaxon2,intaxon3,intaxon4);
	fixroot = 1;
	
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

	vector<double> deltalogp;
	vector<double> logbf;
	double meandeltalogp = 0;
	double vardeltalogp = 0;
	double meanlogbf = 0;
	double varlogbf = 0;

	int b = 100;

	ofstream logos((name + ".logbflog").c_str());

	while (i < until)	{
		cerr << ".";
		cerr.flush();
		samplesize++;
		FromStream(is);
		i++;

		QuickUpdate();

		logos << i << '\n';
		Trace(logos);
		for (int k=0; k<b; k++)	{
			Move(1.0);
			Trace(logos);
		}
		logos << '\n';
		logos.flush();

		double tmpdeltalogp = 0;
		double tmplogbf = 0;
		TemperedGibbsSPR(0,0,nstep,1,2,tmpdeltalogp,tmplogbf);
		deltalogp.push_back(tmpdeltalogp);
		logbf.push_back(tmplogbf);
		meandeltalogp += tmpdeltalogp;
		vardeltalogp += tmpdeltalogp * tmpdeltalogp;
		meanlogbf += tmplogbf;
		varlogbf += tmplogbf * tmplogbf;

		logos << tmplogbf << '\n' << '\n';
		logos.flush();

		int nrep = 1;
		while ((i<until) && (nrep < every))	{
			FromStream(is);
			i++;
			nrep++;
		}
	}
	cerr << '\n';

	meandeltalogp /= samplesize;
	vardeltalogp /= samplesize;
	vardeltalogp -= meandeltalogp * meandeltalogp;
	meanlogbf /= samplesize;
	varlogbf /= samplesize;
	varlogbf -= meanlogbf * meanlogbf;

	if (logbf.size() != samplesize)	{
		cerr << "error in read bf: non matching size\n";
		exit(1);
	}

	double max = logbf[0];
	for (int i=1; i<samplesize; i++)	{
		if (max < logbf[i])	{
			max = logbf[i];
		}
	}
	double mean = 0;
	double var = 0;
	for (int i=0; i<samplesize; i++)	{
		double tmp = exp(logbf[i] - max);
		mean += tmp;
		var += tmp*tmp;
	}
	mean /= samplesize;
	var /= samplesize;
	var -= mean*mean;
	double effsize = mean*mean / var;
	cout << taxon1 << '\t' << taxon2 << '\t' << taxon3 << '\t' << taxon4 << '\n';
	cout << "logbf1: " << log(mean) + max << '\t' << effsize << '\n';
	cout << "logbf2: " << meanlogbf << '\t' << varlogbf << '\n';
	cout << "dlogp: " << meandeltalogp << '\t' << vardeltalogp << '\n';

	ofstream os((name + ".logbf").c_str());
	os << taxon1 << '\t' << taxon2 << '\t' << taxon3 << '\t' << taxon4 << '\n';
	os << "logbf1: " << log(mean) + max << '\t' << effsize << '\n';
	os << "logbf2: " << meanlogbf << '\t' << varlogbf << '\n';
	os << "dlogp: " << meandeltalogp << '\t' << vardeltalogp << '\n';
	os << '\n';
	for (int i=0; i<samplesize; i++)	{
		os << logbf[i] << '\t' << deltalogp[i] << '\n';
	}
}

void PhyloProcess::PostPred(int ppredtype, string name, int burnin, int every, int until, int inrateprior, int inprofileprior, int inrootprior)	{

	GlobalSetRatePrior(inrateprior);
	GlobalSetProfilePrior(inprofileprior);
	GlobalSetRootPrior(inrootprior);

	ifstream is((name + ".chain").c_str());
	if (!is)	{
		cerr << "error: no .chain file found\n";
		exit(1);
	}

	double* obstaxstat = new double[GetNtaxa()];
	SequenceAlignment* datacopy  = new SequenceAlignment(GetData());
	double obs = 0;
	double obs2 = 0;
	if (ppredtype == 2)	{
		obs = data->GetMeanDiversity();
		// obs = GlobalGetMeanDiversity();
	}
	else if (ppredtype == 3)	{
		obs = GetObservedCompositionalHeterogeneity(obstaxstat,obs2);
	}

	cerr << "burnin: " << burnin << '\n';
	cerr << "every " << every << " points until " << until << '\n';
	// cerr << "number of points : " << (until - burnin)/every << '\n';
	int i=0;
	while ((i < until) && (i < burnin))	{
		FromStream(is);
		i++;
	}
	int samplesize = 0;
	double meanstat = 0;
	double varstat = 0;
	double ppstat = 0;
	double meanstat2 = 0;
	double varstat2 = 0;
	double ppstat2 = 0;
	double* meantaxstat = new double[GetNtaxa()];
	double* vartaxstat = new double[GetNtaxa()];
	double* pptaxstat = new double[GetNtaxa()];
	double* taxstat = new double[GetNtaxa()];
	for (int j=0; j<GetNtaxa(); j++)	{
		meantaxstat[j] = 0;
		vartaxstat[j] = 0;
		pptaxstat[j] = 0;
	}
	while (i < until)	{
		cerr << ".";
		samplesize++;
		FromStream(is);
		i++;

		// output tree
		ostringstream s;
		s << name << "_ppred" << samplesize << ".tree";
		ofstream os(s.str().c_str());
		SetNamesFromLengths();
		RenormalizeBranchLengths();
		GetTree()->ToStream(os);
		DenormalizeBranchLengths();
		os.close();

		MPI_Status stat;
		MESSAGE signal = BCAST_TREE;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
		GlobalBroadcastTree();
		GlobalUpdateConditionalLikelihoods();
		GlobalUnclamp();
		GlobalCollapse();
		GlobalSetDataFromLeaves();

		if (ppredtype > 1)	{
			double stat = 0;
			double stat2 = 0;
			if (ppredtype == 2)	{
				stat = data->GetMeanDiversity();
			}
			else if (ppredtype == 3)	{
				stat = GetCompositionalHeterogeneity(taxstat,stat2);
				for (int j=0; j<GetNtaxa(); j++)	{
					meantaxstat[j] += taxstat[j];
					vartaxstat[j] += taxstat[j] * taxstat[j];
					if (taxstat[j] > obstaxstat[j])	{
						pptaxstat[j] ++;
					}
				}
			}
			meanstat += stat;
			varstat += stat * stat;
			if (stat < obs)	{
				ppstat++;
			}
			meanstat2 += stat2;
			varstat2 += stat2 * stat2;
			if (stat2 < obs2)	{
				ppstat2++;
			}
		}
		else	{
			// write datafile
			ostringstream s;
			s << name << "_ppred" << samplesize << ".ali";
			ofstream os(s.str().c_str());
			data->ToStream(os);
			os.close();
		}

		GlobalRestoreData();
		GlobalUnfold();

		int nrep = 1;
		while ((i<until) && (nrep < every))	{
			FromStream(is);
			i++;
			nrep++;
		}
	}
	cerr << '\n';
	cerr << '\n';
	if (ppredtype > 1)	{
		meanstat /= samplesize;
		varstat /= samplesize;
		varstat -= meanstat * meanstat;
		ppstat /= samplesize;

		meanstat2 /= samplesize;
		varstat2 /= samplesize;
		varstat2 -= meanstat2 * meanstat2;
		ppstat2 /= samplesize;
	}

	if (ppredtype == 1)	{
		cerr << "datasets in " << name << "_ppred<rep>.ali\n";
	}
	if (ppredtype == 2)	{
		ofstream os((name + ".div").c_str());
		os << "diversity test\n";
		os << "obs div : " << obs << '\n';
		os << "mean div: " << meanstat << " +/- " << sqrt(varstat) << '\n';
		os << "z-score : " << (meanstat - obs) / sqrt(varstat) << '\n';
		os << "pp      : " << ppstat << '\n';
		cerr << "result of diversity test in " << name << ".div\n";
	}
	else if (ppredtype == 3)	{
		ofstream os((name + ".comp").c_str());
		os << "compositional homogeneity test\n";

		os << '\n';
		os << "max heterogeneity across taxa\n";
		os << "obs comp : " << obs << '\n';
		os << "mean comp: " << meanstat << " +/- " << sqrt(varstat) << '\n';
		os << "z-score : " << (obs - meanstat) / sqrt(varstat) << '\n';
		os << "pp      : " << (1 - ppstat) << '\n';

		os << '\n';
		os << "mean squared heterogeneity across taxa\n";
		os << "obs comp : " << obs2 << '\n';
		os << "mean comp: " << meanstat2 << " +/- " << sqrt(varstat2) << '\n';
		os << "z-score : " << (obs2 - meanstat2) / sqrt(varstat2) << '\n';
		os << "pp      : " << (1 - ppstat2) << '\n';

		os << '\n';
		os << "taxonname\tobs\tmean pred\tz-score\tpp\n";
		for (int j=0; j<GetNtaxa(); j++)	{
			meantaxstat[j] /= samplesize;
			vartaxstat[j] /= samplesize;
			pptaxstat[j] /= samplesize;
			vartaxstat[j] -= meantaxstat[j] * meantaxstat[j];
			os << GetTaxonSet()->GetTaxon(j) << '\t' << obstaxstat[j] << '\t' << meantaxstat[j] << '\t' << (obstaxstat[j] - meantaxstat[j])/sqrt(vartaxstat[j]) << '\t' << pptaxstat[j] << '\n';
		}
		cerr << "result of compositional homogeneity test in " << name << ".comp\n";
	}
	cerr << '\n';
}

void PhyloProcess::GlobalSetTestData()	{

	testnsite = testdata->GetNsite();
	int* tmp = new int[testnsite * GetNtaxa()];
	testdata->GetDataVector(tmp);

	if (GetNprocs() > 1)	{
		MESSAGE signal = SETTESTDATA;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(&testnsite,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(tmp,testnsite*GetNtaxa(),MPI_INT,0,MPI_COMM_WORLD);
	}
	else	{
		GetData()->SetTestData(testnsite,0,0,testnsite,tmp);
	}

	delete[] tmp;
}

void PhyloProcess::SlaveSetTestData()	{

	MPI_Bcast(&testnsite,1,MPI_INT,0,MPI_COMM_WORLD);
	int* tmp = new int[testnsite * GetNtaxa()];
	MPI_Bcast(tmp,testnsite*GetNtaxa(),MPI_INT,0,MPI_COMM_WORLD);
	
	SetTestSiteMinAndMax();
	GetData()->SetTestData(testnsite,GetSiteMin(),testsitemin,testsitemax,tmp);

	delete[] tmp;
}

void PhyloProcess::ReadCV(string testdatafile, string name, int burnin, int every, int until, int iscodon, GeneticCodeType codetype)	{
	
	ifstream is((name + ".chain").c_str());
	if (!is)	{
		cerr << "error: no .chain file found\n";
		exit(1);
	}

	if (iscodon)	{
		SequenceAlignment* tempdata = new FileSequenceAlignment(testdatafile);
		testdata = new CodonSequenceAlignment(tempdata,true,codetype);
	}
	else	{
		testdata = new FileSequenceAlignment(testdatafile);
	}
	GlobalSetTestData();

	cerr << "burnin: " << burnin << '\n';
	cerr << "every " << every << " points until " << until << '\n';
	int i=0;
	while ((i < until) && (i < burnin))	{
		cout << "before FromStream...\n";
		cout.flush();
		FromStream(is);
		cout << "after FromStream...\n";
		cout.flush();
		i++;
	}
	int samplesize = 0;
	vector<double> scorelist;


	while (i < until)	{
		cerr << ".";
		samplesize++;
		FromStream(is);
		i++;
		QuickUpdate();
		// Trace(cerr);
		MPI_Status stat;
		MESSAGE signal = CVSCORE;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

		double tmp = 0;
		double score = 0;
		for(int i=1; i<GetNprocs(); ++i) {
			MPI_Recv(&tmp,1,MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
			score += tmp;
		}
		scorelist.push_back(score);
		
		int nrep = 1;
		while ((i<until) && (nrep < every))	{
			FromStream(is);
			i++;
			nrep++;
		}
	}

	cerr << '\n';

	double max = 0;
	for (int j=0; j<samplesize; j++)	{
		if ((!j) || (max < scorelist[j]))	{
			max = scorelist[j];
		}
	}

	double tot = 0;
	for (int j=0; j<samplesize; j++)	{
		tot += exp(scorelist[j] - max);
	}
	tot /= samplesize;
	
	double meanscore = log(tot) + max;
	
	ofstream os((name + ".cv").c_str());
	os << meanscore << '\n';
	cerr << meanscore << '\n';
}

void PhyloProcess::ReadSiteLogL(string name, int burnin, int every, int until)	{

	ifstream is((name + ".chain").c_str());
	if (!is)	{
		cerr << "error: no .chain file found\n";
		exit(1);
	}

	cerr << "burnin: " << burnin << '\n';
	cerr << "every " << every << " points until " << until << '\n';
	int i=0;
	while ((i < until) && (i < burnin))	{
		FromStream(is);
		i++;
	}
	int samplesize = 0;

	double* tmp = new double[GetNsite()];
	double* mean = new double[GetNsite()];
	vector<double>* logl = new vector<double>[GetNsite()];

	for (int i=0; i<GetNsite(); i++)	{
		mean[i] = 0;
	}

	while (i < until)	{
		cerr << ".";
		samplesize++;
		FromStream(is);
		i++;
		QuickUpdate();
		MPI_Status stat;
		MESSAGE signal = SITELOGL;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

		double total = 0;
		for(int i=1; i<GetNprocs(); ++i) {
			MPI_Recv(tmp,GetNsite(),MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
			for (int j=GetProcSiteMin(i); j<GetProcSiteMax(i); j++)	{
				logl[j].push_back(tmp[j]);
				mean[j] += tmp[j];
				total += tmp[j];
			}
		}
		
		int nrep = 1;
		while ((i<until) && (nrep < every))	{
			FromStream(is);
			i++;
			nrep++;
		}
	}

	double meancpo = 0;
	double varcpo = 0;
	double* cpo = new double[GetNsite()];
	for (int i=0; i<GetNsite(); i++)	{
		double min = 0;
		for (vector<double>::iterator j=logl[i].begin(); j != logl[i].end(); j++)	{
			if (min > (*j))	{
				min = *j;
			}
		}
		double hmean = 0;
		int count = 0;
		for (vector<double>::iterator j=logl[i].begin(); j != logl[i].end(); j++)	{
			hmean += exp(min - (*j));
			count++;
		}
		hmean /= count;
		cpo[i] = min - log(hmean);
		meancpo += cpo[i];
		varcpo += cpo[i] * cpo[i];
	}
	meancpo /= GetNsite();
	varcpo /= GetNsite();
	varcpo -= meancpo * meancpo;

	ofstream os((name + ".sitelogl").c_str());
	double total = 0;
	for (int i=0; i<GetNsite(); i++)	{
		mean[i] /= samplesize;
		total += mean[i];
		os << i+1 << '\t' << mean[i] << '\t' << cpo[i] << '\n';
	}

	ofstream cos((name + ".cpo").c_str());
	cos << "posterior mean ln L : " << total << '\n';
	cos << "CPO : " << GetNsite() * meancpo << '\t' << meancpo << '\t' << sqrt(varcpo) << '\n';
	
	cerr << '\n';
	cerr << "posterior mean ln L : " << total << '\n';
	cerr << "site-specific posterior mean ln L in " << name << ".sitelogl\n";
	cerr << "CPO: " << GetNsite() * meancpo << '\t' << meancpo << '\t' << sqrt(varcpo) << '\n';
	cerr << '\n';

}


void PhyloProcess::ReadMap(string name, int burnin, int every, int until){
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
	double meandiff = 0;
	double vardiff = 0;
	double meanobs = 0;
	for(int i = 0; i < GetNsite(); i++){
		stringstream osfmap;
		osfmap << name << '_' << i << ".map";
		ofstream osmap((osfmap.str()).c_str());
		osmap.close();
	}
	while (i < until)	{
		cerr << ".";
		// cerr << i << '\t' << rnd::GetRandom().Uniform() << '\n';

		cerr.flush();
		samplesize++;
		FromStream(is);
		i++;

		// prepare file for ancestral node states
		ostringstream s;
		s << name << "_" << samplesize << ".nodestates";
		ofstream sos(s.str().c_str());

		// quick update and mapping on the fly
		MPI_Status stat;
		MESSAGE signal = BCAST_TREE;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
		GlobalBroadcastTree();
		GlobalUpdateConditionalLikelihoods();
		GlobalCollapse();

		// write posterior mappings
		GlobalWriteMappings(name);

		// write posterior ancestral node states
		GlobalSetNodeStates();
		WriteNodeStates(sos,GetRoot());
		sos << '\n';

		double obs = GlobalCountMapping();

		//Posterior Predictive Mappings
		GlobalUnfold();
		GlobalUnclamp();
		GlobalCollapse();

		GlobalSetDataFromLeaves();

		// write posterior predictive mappings
		GlobalWriteMappings(name);

		// write posterior predictive ancestral node states
		GlobalSetNodeStates();
		WriteNodeStates(sos,GetRoot());

		double pred = GlobalCountMapping();

		obs /= GetNsite();
		pred /= GetNsite();

		meandiff += obs - pred;
		vardiff += (obs-pred)*(obs-pred);
		meanobs += obs;

		GlobalRestoreData();
		GlobalUnfold();

		for(int i = 0; i < GetNsite(); i++){
			stringstream osfmap;
			osfmap << name << '_' << i << ".map";
			ofstream osmap((osfmap.str()).c_str(), ios_base::app);
			osmap << '\n';
			osmap.close();
		}

		int nrep = 1;
		while ((i<until) && (nrep < every))	{
			FromStream(is);
			i++;
			nrep++;
		}
	}
	cerr << '\n';
	meandiff /= samplesize;
	vardiff /= samplesize;
	vardiff -= meandiff*meandiff;
	meanobs /= samplesize;
	cerr << "mean obs : " << meanobs << '\n';
	cerr << meandiff << '\t' << sqrt(vardiff) << '\n';
}

void PhyloProcess::GlobalWriteMappings(string name){

	if (GetNprocs() > 1)	{
		MPI_Status stat;
		MESSAGE signal = WRITE_MAPPING;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

		 //send the chain name
		ostringstream os;
		os << name;
		string s = os.str();
		unsigned int len = s.length();
		unsigned char* bvector = new unsigned char[len];
		for (unsigned int i=0; i<len; i++)	{
			bvector[i] = s[i];
		}
		MPI_Bcast(&len,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(bvector,len,MPI_UNSIGNED_CHAR,0,MPI_COMM_WORLD);
		delete[] bvector;
	}
	else	{
		for(int i=0; i<GetNsite(); i++)	{
			stringstream osfmap;
			osfmap << name << '_' << i << ".map";
			ofstream osmap((osfmap.str()).c_str(), ios_base::app);
			WriteTreeMapping(osmap, GetRoot(), i);
			osmap.close();
		}
	}
}

void PhyloProcess::SlaveWriteMappings(){

	int len;
	MPI_Bcast(&len,1,MPI_INT,0,MPI_COMM_WORLD);
	unsigned char* bvector = new unsigned char[len];
	MPI_Bcast(bvector,len,MPI_UNSIGNED_CHAR,0,MPI_COMM_WORLD);
	ostringstream os;
	for (int i=0; i<len; i++)	{
		os << bvector[i];
	}
	string name = os.str();
	delete[] bvector;

	for(int i = GetSiteMin(); i < GetSiteMax(); i++){
		stringstream osfmap;
		osfmap << name << '_' << i << ".map";
		ofstream osmap((osfmap.str()).c_str(), ios_base::app);
		WriteTreeMapping(osmap, GetRoot(), i);
		osmap.close();
	}
}


void PhyloProcess::WriteTreeMapping(ostream& os, const Link* from, int i){
	if(from->isLeaf()){
		os << from->GetNode()->GetName();
	}
	else{
		os << '(';
		for (const Link* link=from->Next(); link!=from; link=link->Next()){
			WriteTreeMapping(os, link->Out(), i);
			if (link->Next() != from)       {
				os << ',';
			}
		}
		os << ')';
	}
	if(from->isRoot()){
		BranchSitePath* mybsp = submap[GetBranchIndex(from->Next()->GetBranch())][i];
		os << '_' << GetStateSpace()->GetState(mybsp->Init()->GetState()) << ";\n";     
	}
	else{
		BranchSitePath* mybsp = submap[GetBranchIndex(from->GetBranch())][i];
		double l = GetLength(from->GetBranch());
		os << '_' << GetStateSpace()->GetState(mybsp->Last()->GetState());
		for(Plink* plink = mybsp->Last(); plink ; plink = plink->Prev()){
			os << ':' << plink->GetRelativeTime() * l << ':' << GetStateSpace()->GetState(plink->GetState());
		}
	}
}

void PhyloProcess::WriteNodeStates(ostream& os, const Link* from)	{
	os << GetLeftMost(from) << '\t' << GetRightMost(from) << '\t';
	int nodelabel = GetNodeIndex(from->GetNode());
	for (int i=0; i<GetNsite(); i++)	{
		os << GetStateSpace()->GetState(nodestate[nodelabel][i]);
	}
	os << '\n';
		
	for (const Link* link=from->Next(); link!=from; link=link->Next()){
		WriteNodeStates(os,link->Out());
	}
}

int PhyloProcess::CountMapping()	{

	int total = 0;	
	for(int i=GetSiteMin(); i<GetSiteMax(); i++){
		total += CountMapping(i);
	}
	return total;
}

int PhyloProcess::CountMapping(int i)	{
	return 0;
}

int PhyloProcess::GlobalCountMapping()	{

	int totalcount=0;

	if (GetNprocs() > 1)	{
		MESSAGE signal = COUNTMAPPING;
		MPI_Status stat;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

		for (int i=1; i<GetNprocs(); i++)	{
			int count;
			MPI_Recv(&count,1,MPI_INT,MPI_ANY_SOURCE,TAG1,MPI_COMM_WORLD, &stat);
			totalcount += count;
		}
	}
	else	{
		totalcount = CountMapping();
	}

	return totalcount;
}

void PhyloProcess::SlaveCountMapping()	{

	int count = CountMapping();
	MPI_Send(&count,1,MPI_INT,0,TAG1,MPI_COMM_WORLD);

}

void PhyloProcess::GlobalUnclamp()	{

	if (GetNprocs() > 1)	{
		MESSAGE signal = UNCLAMP;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	}
	else	{
		GetData()->Unclamp();
	}
}

void PhyloProcess::GlobalRestoreData()	{

	if (GetNprocs() > 1)	{
		MESSAGE signal = RESTOREDATA;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	}
	else	{
		GetData()->Restore();
	}
}

void PhyloProcess::GlobalSetDataFromLeaves()	{

	if (GetNprocs() > 1)	{
		MESSAGE signal = SETDATA;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

		MPI_Status stat;

		int* tmp = new int[GetMaxSiteNumber() * GetNtaxa()];

		for(int i=1; i<GetNprocs(); i++) {
			MPI_Recv(tmp,GetProcSiteNumber(i)*GetNtaxa(),MPI_INT,i,TAG1,MPI_COMM_WORLD,&stat);
			int k = 0;
			for (int l=0; l<GetNtaxa(); l++)	{
				for (int j=GetProcSiteMin(i); j<GetProcSiteMax(i); j++)	{
					GetData()->SetState(l,j,tmp[k]);
					k++;
				}
			}
		}

		delete[] tmp;
	}
	else	{
		SetDataFromLeaves();
	}
}

void PhyloProcess::GlobalSetNodeStates()	{

	if (GetNprocs() > 1)	{
		MESSAGE signal = SETNODESTATES;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

		MPI_Status stat;

		int* tmp = new int[GetMaxSiteNumber() * GetNnode()];

		for(int i=1; i<GetNprocs(); i++) {
			MPI_Recv(tmp,GetProcSiteNumber(i)*GetNnode(),MPI_INT,i,TAG1,MPI_COMM_WORLD,&stat);
			int k = 0;
			for (int l=0; l<GetNnode(); l++)	{
				for (int j=GetProcSiteMin(i); j<GetProcSiteMax(i); j++)	{
					nodestate[l][j] = tmp[k];
					k++;
				}
			}
			delete[] tmp;
		}
	}
}


double PhyloProcess::GlobalGetMeanDiversity()	{

	double total = 0;
	if (GetNprocs() > 1)	{
		MESSAGE signal = GETDIV;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

		MPI_Status stat;

		for(int i=1; i<GetNprocs(); ++i) {
			double tmp;
			MPI_Recv(&tmp,1,MPI_DOUBLE,MPI_ANY_SOURCE,TAG1,MPI_COMM_WORLD,&stat);
			total += tmp;
		}
	}
	else	{
		total = GetData()->GetTotalDiversity(GetSiteMin(),GetSiteMax());
	}
	return total / GetNsite();
}

void PhyloProcess::SlaveRestoreData()	{
	GetData()->Restore();
}

void PhyloProcess::SlaveUnclamp()	{
	GetData()->Unclamp();
}

void PhyloProcess::SlaveSetDataFromLeaves()	{
	SetDataFromLeaves();

	// mpi send the array
	int* tmp = new int[(GetSiteMax() - GetSiteMin()) * GetNtaxa()];
	int k = 0;
	for (int i=0; i<GetNtaxa(); i++)	{
		for (int j=GetSiteMin(); j<GetSiteMax(); j++)	{
			tmp[k] = GetData()->GetState(i,j);
			k++;
		}
	}
	MPI_Send(tmp,(GetSiteMax()-GetSiteMin())*GetNtaxa(),MPI_INT,0,TAG1,MPI_COMM_WORLD);

	delete[] tmp;

}

void PhyloProcess::SlaveSetNodeStates()	{

	// mpi send the array
	int* tmp = new int[(GetSiteMax() - GetSiteMin()) * GetNnode()];
	int k = 0;
	for (int i=0; i<GetNnode(); i++)	{
		for (int j=GetSiteMin(); j<GetSiteMax(); j++)	{
			tmp[k] = nodestate[i][j];
			k++;
		}
	}
	MPI_Send(tmp,(GetSiteMax() - GetSiteMin())*GetNnode(),MPI_INT,0,TAG1,MPI_COMM_WORLD);

	delete[] tmp;

}

void PhyloProcess::SlaveGetMeanDiversity()	{

	double div = GetData()->GetTotalDiversity(GetSiteMin(),GetSiteMax());
	MPI_Send(&div,1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}
