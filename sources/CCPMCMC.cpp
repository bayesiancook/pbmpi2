
#include "CCP.h"

#include "Parallel.h"
MPI_Datatype Propagate_arg;

UnrootedBPP* ccp;
Tree* tree;

void RecursiveGibbsSPRScan(Link* from, Link* fromup, Link* down, Link* up, map<pair<Link*,Link*>,double>& loglmap)	{

	if (! from->isRoot())	{
		tree->Attach(down,up,from,fromup);
		loglmap[pair<Link*,Link*>(from,fromup)] = ccp->GetLogProb(tree);
		tree->Detach(down,up);
	}
	Link* trailer = from;
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		RecursiveGibbsSPRScan(link->Out(),trailer,down,up,loglmap);
		trailer = trailer->Next();
	}
}

int GibbsSPR()	{

	Link* down = 0;
	Link* up = 0;
	tree->DrawSubTree(down,up);
	Link* fromdown = tree->Detach(down,up);

	map<pair<Link*,Link*>,double> loglmap;
	RecursiveGibbsSPRScan(tree->GetRoot(),tree->GetRoot(),down,up,loglmap);

	// check for numerical problems
	for (map<pair<Link*,Link*>,double>::iterator i=loglmap.begin(); i!=loglmap.end(); i++)	{
		if (isnan(i->second))	{
			cerr << "nan log prob in gibbs\n";
			exit(1);
		}
		if (isinf(i->second))	{
			cerr << "inf log prob in gibbs\n";
			exit(1);
		}
	}

	// choose among all target positions including current position
	double max = 0;
	int notfound = 1;
	for (map<pair<Link*,Link*>,double>::iterator i=loglmap.begin(); i!=loglmap.end(); i++)	{
		if ((notfound) || (max < i->second))	{
			max = i->second;
			notfound = 0;
		}
	}

	// compute total prob mass
	double total = 0;
	int nchoice = 0;
	for (map<pair<Link*,Link*>,double>::iterator i=loglmap.begin(); i!=loglmap.end(); i++)	{
		nchoice++;
		double tmp = exp(i->second - max);
		total += tmp;
		if (isinf(total))	{
			cerr << "error in gibbs: inf\n";
			cerr << "total1\n";
			exit(1);
		}
		if (isnan(total))	{
			cerr << "error in gibbs: nan\n";
		}
	}
	if (! nchoice)	{
		cerr << "error in gibbs: no possible choice\n";
		exit(1);
	}

	// randomly choose final position
	double u = total * rnd::GetRandom().Uniform();
	map<pair<Link*,Link*>, double>::iterator i = loglmap.begin();
	double cumul = exp(i->second-max);
	while ((i!=loglmap.end()) && (cumul < u))	{
		i++;
		if (i == loglmap.end())	{
			cerr << "error in gibbs spr: overflow\n";
			exit(1);
		}
		cumul += exp(i->second-max);
	}
	
	// store values of target position
	Link* todown = i->first.first;
	Link* toup = i->first.second;

	int accepted = (todown != fromdown);

	tree->Attach(down,up,todown,toup);
	return accepted;
}

int MHGibbsSPR()	{

	double logp0 = ccp->GetLogProb(tree);
	Link* down = 0;
	Link* up = 0;
	tree->DrawSubTree(down,up);
	Link* fromdown = tree->Detach(down,up);

	map<pair<Link*,Link*>,double> loglmap;
	RecursiveGibbsSPRScan(tree->GetRoot(),tree->GetRoot(),down,up,loglmap);

	// check for numerical problems
	for (map<pair<Link*,Link*>,double>::iterator i=loglmap.begin(); i!=loglmap.end(); i++)	{
		if (isnan(i->second))	{
			cerr << "nan log prob in gibbs\n";
			exit(1);
		}
		if (isinf(i->second))	{
			cerr << "inf log prob in gibbs\n";
			exit(1);
		}
	}

	// choose among all target positions different from current position
	double max1 = 0;
	int notfound = 1;
	Link* fromup = 0;
	double logp1 = 0;
	for (map<pair<Link*,Link*>,double>::iterator i=loglmap.begin(); i!=loglmap.end(); i++)	{
		// grep values of current position
		if (i->first.first == fromdown)	{
			fromup = i->first.second;
			logp1 = i->second;
		}
		// find max, excluding current position
		if ((i->first.first != fromdown) && ((notfound) || (max1 < i->second)))	{
			max1 = i->second;
			notfound = 0;
		}
	}

	/*
	if (fromup->Next()->Out() == fromdown)	{
		cerr << "fromup next out == fromdown\n";
	}
	else	{
		cerr << "fromup next out != fromdown\n";
	}
	*/

	if (fabs(logp1 - logp0) > 1e-6)	{
		cerr << "error in mh: non matching probs\n";
		cerr << logp1 << '\t' << logp0 << '\n';
		exit(1);
	}

	// compute total prob mass
	double total1 = 0;
	int nchoice = 0;
	for (map<pair<Link*,Link*>,double>::iterator i=loglmap.begin(); i!=loglmap.end(); i++)	{
		if (i->first.first != fromdown)	{
			nchoice++;
			double tmp = exp(i->second - max1);
			total1 += tmp;
		}
		if (isinf(total1))	{
			cerr << "error in gibbs: inf\n";
			cerr << "total1\n";
			exit(1);
		}
		if (isnan(total1))	{
			cerr << "error in gibbs: nan\n";
		}
	}
	if (! nchoice)	{
		cerr << "error in gibbs: no possible choice\n";
		exit(1);
	}

	// randomly choose final position, duly excluding current position
	double u = total1 * rnd::GetRandom().Uniform();
	map<pair<Link*,Link*>, double>::iterator i = loglmap.begin();
	double cumul = 0;
	if (i->first.first != fromdown)	{
		cumul += exp(i->second-max1);
	}
	while ((i!=loglmap.end()) && (cumul < u))	{
		i++;
		if (i == loglmap.end())	{
			cerr << "error in gibbs spr: overflow\n";
			exit(1);
		}
		if (i->first.first != fromdown)	{
			cumul += exp(i->second-max1);
		}
	}
	
	// store values of target position
	Link* todown = i->first.first;
	Link* toup = i->first.second;
	double logp2 = i->second;

	// calculate probability of reverse move
	double max2 = 0;
	notfound = 1;
	for (map<pair<Link*,Link*>,double>::iterator i=loglmap.begin(); i!=loglmap.end(); i++)	{
		if ((i->first.first != todown) && ((notfound) || (max2 < i->second)))	{
			max2 = i->second;
			notfound = 0;
		}
	}

	double total2 = 0;
	for (map<pair<Link*,Link*>,double>::iterator i=loglmap.begin(); i!=loglmap.end(); i++)	{
		if (i->first.first != todown)	{
			double tmp = exp(i->second - max2);
			total2 += tmp;
		}
		if (isinf(total2))	{
			cerr << "error in gibbs: inf\n";
		}
		if (isnan(total2))	{
			cerr << "error in gibbs: nan\n";
		}
	}

	// hastings log ratio
	double logpfwd = logp2 - max1 - log(total1);
	double logprev = logp1 - max2 - log(total2);
	double logh = logprev - logpfwd;

	// MH log ratio
	double logratio = logh + logp2 - logp1;

	int accepted = (log(rnd::GetRandom().Uniform()) < logratio);
	if (accepted)	{
		tree->Attach(down,up,todown,toup);
		double logp3 = ccp->GetLogProb(tree);
		if (fabs(logp3 - logp2) > 1e-6)	{
			cerr << "error in mh accepted: non matching probs\n";
			cerr << logp2 << '\t' << logp3 << '\n';
			exit(1);
		}
	}
	else	{
		tree->Attach(down,up,fromdown,fromup);
		double logp3 = ccp->GetLogProb(tree);
		if (fabs(logp3 - logp1) > 1e-6)	{
			cerr << "error in mh accepted: non matching probs\n";
			cerr << logp1 << '\t' << logp3 << '\n';
			exit(1);
		}
	}
	return accepted;
}

int main(int argc, char* argv[])	{

	// get ccp from file
	string ccpfile = argv[1];
	double cutoff = atof(argv[2]);
	double beta = atof(argv[3]);
	int nrep = atoi(argv[4]);
	string name = argv[5];
	string type = argv[6];
	string ccptype = argv[7];

	int subsample = 100;

	if (ccptype == "bpp")	{
		ccp = new UnrootedBPP(ccpfile,cutoff,beta);
	}
	else if (ccptype == "ccp")	{
		ccp = new UnrootedCCP(ccpfile,cutoff,beta);
	}
	else	{
		cerr << "unrecognised ccp type\n";
		exit(1);
	}
	
	// get TaxonSet
	const TaxonSet* taxset = ccp->GetTaxonSet();

	// make random tree based on taxon set
	tree = new Tree(taxset);
	// ifstream is("ref.tre");
	// tree->ReadFromStream(is);
	tree->MakeRandomTree();
	tree->RegisterWith(taxset);
	tree->EraseInternalNodeName();

	BpList* mcmcbp = new BpList(taxset);
	// TpList* mcmctp = new TpList(taxset);

	double logprob = ccp->GetLogProb(tree);
	cerr << "initial log prob : " << logprob << '\n';
	// make a MCMC
	// gibbs based on uniform random draws of subtree to be pruned

	double naccept = 0;
	double ntry = 0;
	ofstream os((name + ".trace").c_str());
	ofstream tos((name + ".treelist").c_str());

	int size = 0;

	ccp->InactivateTP();

	while (1)	{

		for (int k=0; k<subsample; k++)	{
		for (int rep=0; rep<nrep; rep++)	{

			tree->RootAtRandom();

			// int accept = ccp->FullGibbsSPR(tree);
			int accept = 0;
			if (type == "mh")	{
				accept = MHGibbsSPR();
			}
			else if (type == "gibbs")	{
				accept = GibbsSPR();
			}
			else if (type == "fullgibbs")	{
				accept = ccp->FullGibbsSPR(tree);
			}
			
			naccept += accept;
			ntry++;
		}

		size++;
		if (size >= 100)	{
			ccp->ActivateTP();
		}
		mcmcbp->Add(tree,1);

		os << naccept/ntry << '\t' << ccp->GetLogProb(tree) << '\n';
		os.flush();

		tree->ToStream(tos);
		tos.flush();
		}

		ofstream bos((name + ".bpdiff").c_str());
		mcmcbp->PrintDiff(bos,ccp->bplist,cutoff);
		bos.flush();
	}
}

