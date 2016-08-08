
#include "CCP.h"

#include "Parallel.h"
MPI_Datatype Propagate_arg;

UnrootedBPP** ccp;
int Ngene;
Tree* tree;
double pi;

double GetLogProb(Tree* intree)	{

	double total = 0;
	for (int gene=0; gene<Ngene; gene++)	{
		// make a copy of the tree
		Tree* copytree = new Tree(intree);

		// prune the tree 
		copytree->RegisterWith(ccp[gene]->GetTaxonSet());
		int s1 = copytree->GetSize();
		int s2 = copytree->GetFullSize(copytree->GetRoot());
		if (s2 != 2*s1 - 2)	{
			cerr << "error when pruning tree\n";
			exit(1);
		}
		double tmp = ccp[gene]->GetLogProb(copytree);
		delete copytree;
		total += log(pi + (1-pi) * exp(tmp));
	}
	return total;
}


void RecursiveGibbsSPRScan(Link* from, Link* fromup, Link* down, Link* up, map<pair<Link*,Link*>,double>& loglmap)	{

	if (! from->isRoot())	{
		tree->Attach(down,up,from,fromup);
		loglmap[pair<Link*,Link*>(from,fromup)] = GetLogProb(tree);
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

int MovePi(double tuning)	{

	double bk = pi;
	double deltalogratio = - GetLogProb(tree);

	double u = tuning * (rnd::GetRandom().Uniform() - 0.5);
	pi += u;
	while ((pi < 0) || (pi > 1))	{
		if (pi < 0)	{
			pi = -pi;
		}
		if (pi > 1)	{
			pi = 2 - pi;
		}
	}
	deltalogratio += GetLogProb(tree);
	int accept = (log(rnd::GetRandom().Uniform()) < deltalogratio);
	if (! accept)	{
		pi = bk;
	}
	return accept;
}


int main(int argc, char* argv[])	{

	// get ccp from file
	string ccplist = argv[1];
	double cutoff = atof(argv[2]);
	double beta = atof(argv[3]);
	int nrep = atoi(argv[4]);
	string name = argv[5];
	string ccptype = argv[6];
	string treefile = argv[7];

	ifstream is(ccplist.c_str());
	is >> Ngene;
	ccp = new UnrootedBPP*[Ngene];
	for (int gene=0; gene<Ngene; gene++)	{
		string ccpfile;
		is >> ccpfile;
		if (ccptype == "bpp")	{
			ccp[gene] = new UnrootedBPP(ccpfile,cutoff,beta);
		}
		else if (ccptype == "ccp")	{
			ccp[gene] = new UnrootedCCP(ccpfile,cutoff,beta);
		}
		else	{
			cerr << "unrecognised ccp type\n";
			exit(1);
		}
	}
	
	tree = new Tree(treefile);
	tree->MakeTaxonSet();
	const TaxonSet* taxset = tree->GetTaxonSet();
	tree->RegisterWith(taxset);
	tree->EraseInternalNodeName();

	pi = 0;

	double naccept = 0;
	double ntry = 0;
	ofstream os((name + ".trace").c_str());
	ofstream tos((name + ".treelist").c_str());

	int size = 0;
	for (int gene=0; gene<Ngene; gene++)	{
		ccp[gene]->InactivateTP();
	}

	double logprob = GetLogProb(tree);
	cerr << "initial log prob : " << logprob << '\n';
	// make a MCMC
	// gibbs based on uniform random draws of subtree to be pruned

	while (1)	{

		for (int rep=0; rep<nrep; rep++)	{

			tree->RootAtRandom();

			int accept = GibbsSPR();
			naccept += accept;
			ntry++;
		}

		size++;
		if (size > 100)	{
			for (int gene=0; gene<Ngene; gene++)	{
				ccp[gene]->InactivateTP();
			}

			// MovePi(1);

		}

		os << naccept/ntry << '\t' << GetLogProb(tree) << '\t' << pi << '\n';
		os.flush();

		tree->ToStream(tos);
		tos.flush();
	}
}

