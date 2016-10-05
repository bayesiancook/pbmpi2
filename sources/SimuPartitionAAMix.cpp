
#include "Tree.h"
#include "SequenceAlignment.h"
#include "Random.h"
#include <map>
#include <iostream>
#include <sstream>

#include "StringStreamUtils.h"
#include "GTRSubMatrix.h"

#include "Parallel.h"
MPI_Datatype Propagate_arg;

class Simulator : public NewickTree {

	public:

	Simulator(string datafile, string treefile, string partitionfile, string paramfile, int inmask, string inbasename)	{

		mask = inmask;
		// take tree

		// take a protein datafile
		// calculate eq freqs with pseudocount

		// take GTR nucleotide mutation rate matrix
		// relative echange rates: fixed
		// equilibrium gc: either fixed or random
		// normalize GTR matrix

		// calculate an AA replacement rate normalizer (so that tree is non-syn length)

		// mu: homothetic factor for substitution rate (around 1)
		// Ne: relative effective population size (around 1)

		// apply Halpern and Bruno formalism (as in Holder)

		basename = inbasename; 

		tree = new Tree(treefile);

		protdata = new FileSequenceAlignment(datafile);
		if (protdata->GetNstate() != Naa)	{
			cerr << "error: should be protein datafile\n";
			exit(1);
		}
		statespace = protdata->GetStateSpace();

		Nsite = protdata->GetNsite();
		currentseq = new int[Nsite];

		ifstream pis(partitionfile.c_str());
		pis >> Ngene;
		genesize = new int[Ngene];
		genefirst = new int[Ngene];
		genedata = new SequenceAlignment*[Ngene];
		genealloc = new int[Nsite];
		int count = 0;
		for (int gene=0; gene<Ngene; gene++)	{
			pis >> genesize[gene];
			genefirst[gene] = count;
			count += genesize[gene];
			genedata[gene] = new SequenceAlignment(protdata,genefirst[gene],genesize[gene]);
			for (int i=0; i<genesize[gene]; i++)	{
				genealloc[genefirst[gene] + i] = gene;
			}
		}
		if (count != Nsite)	{
			cerr << "error: non matching total size\n";
			exit(1);
		}

		// get parameters from file
		ifstream prmis(paramfile.c_str());

		cerr << "Param file : " << paramfile << '\n';
		// read file
		string tmp;
		prmis >> tmp;

		if (tmp != "Alpha")	{
			cerr << "error: missing Alpha keyword\n";
			cerr << tmp << '\n';
			exit(1);
		}
		prmis >> alpha;

		prmis >> tmp;
		if (tmp != "mu")	{
			cerr << "error: missing mu keyword\n";
			cerr << tmp << '\n';
			exit(1);
		}
		prmis >> mu;


		int withempfreq = 1;
		prmis >> tmp;
		if (tmp != "Mixture")	{
			cerr << "error: missing Mixture keyword\n";
			cerr << tmp << '\n';
			exit(1);
		}
		prmis >> Ncat; 
		if (Ncat == -1)	{
			Ncat = Nsite;

			double pseudocount;
			int focus;
			prmis >> pseudocount >> focus;

			stat = new double*[Ncat];
			for (int k=0; k<Ncat; k++)	{
				stat[k] = new double[Naa];
			}

			for (int gene=0; gene<Ngene; gene++)	{
				genedata[gene]->GetSiteEmpiricalFreq(stat + genefirst[gene],pseudocount,focus);
			}

			alloc = new int[Nsite];
			for (int i=0; i<Nsite; i++)	{
				alloc[i] = i;
			}
		}
		else	{
			prmis >> tmp;
			if (tmp == "+F")	{
				withempfreq = 1;
			}
			else if (tmp == "-F")	{
				withempfreq = 0;
			}
			else	{
				cerr << "error: does not recognize mixture type\n";
				cerr << tmp << '\n';
				exit(1);
			}

			if (withempfreq)	{
				Ncat++;
			}
			mixweight = new double[Ncat];
			prmis >> tmp;
			if (tmp != "Weights")	{
				cerr << "error: missing Weights keyword\n";
				cerr << tmp << '\n';
				exit(1);
			}
			double totweight = 0;
			for (int k=0; k<Ncat; k++)	{
				prmis >> mixweight[k];
				totweight += mixweight[k];
			}
			for (int k=0; k<Ncat; k++)	{
				mixweight[k] /= totweight;
			}
			
			prmis >> tmp;
			if (tmp != "Stats")	{
				cerr << "error: missing Weights keyword\n";
				cerr << tmp << '\n';
				exit(1);
			}
			stat = new double*[Ncat];
			for (int k=0; k<Ncat; k++)	{
				stat[k] = new double[Naa];
			}
			for (int k=withempfreq; k<Ncat; k++)	{
				double tot = 0;
				for (int i=0; i<Naa; i++)	{
					prmis >> stat[k][i];
					tot += stat[k][i];
				}
				for (int i=0; i<Naa; i++)	{
					stat[k][i] /= tot;
				}
			}
			if (withempfreq)	{
				protdata->GetEmpiricalFreq(stat[0]);
			}
			alloc = new int[Nsite];
			for (int i=0; i<Nsite; i++)	{
				alloc[i] = rnd::GetRandom().FiniteDiscrete(Ncat,mixweight);
			}

		}

		// relative exchangeabilites
		Nrrcat = 5;
		Nrr = Naa * (Naa-1) / 2;
		rr = new double*[Nrrcat];
		for (int k=0; k<Nrrcat; k++)	{
			rr[k] = new double[Nrr];
			for (int i=0; i<Nrr; i++)	{
				if (k == 0)	{
					rr[k][i] = WAG_RR[i];
				}
				else if (k == 1)	{
					rr[k][i] = BLOSUM62_RR[i];
				}
				else if (k == 2)	{
					rr[k][i] = JTT_RR[i];
					if (JTT_RR[i] != JTT2_RR[i])	{
						cerr << "error: check on JTT\n";
						exit(1);
					}
				}
				else if (k == 3)	{
					rr[k][i] = DCmut_RR[i];
				}
				else if (k == 4)	{
					rr[k][i] = mtREV_RR[i];
				}
				else	{
					cerr << "error: only 5 categories\n";
					exit(1);
				}
			}
		}

		generralloc = new int[Ngene];
		for (int gene=0; gene<Ngene; gene++)	{
			generralloc[gene] = (int) (Nrrcat * rnd::GetRandom().Uniform());
		}

		ofstream pos((basename + ".param").c_str());

		pos << "rr categories\n";
		pos << "0 : WAG\n";
		pos << "1: BLOSUM62\n";
		pos << "2: JTT\n";
		pos << "3: DCmut\n";
		pos << "4: mtREV\n";
		pos << '\n';

		pos << "gene rr categories\n";
		for (int gene=0; gene<Ngene; gene++)	{
			pos << generralloc[gene];
			pos << '\t';
		}
		pos << '\n';

		if (Ncat != Nsite)	{
			pos << "site-specific mixture allocations\n";
			for (int i=0; i<Nsite; i++)	{
				pos << alloc[i] << '\t';
			}
			pos << '\n';
		}
		else	{
			pos << "site-specific eq freq\n";
			for (int i=0; i<Nsite; i++)	{
				for (int k=0; k<Naa; k++)	{
					pos << stat[i][k] << '\t';
				}
				pos << '\n';
			}
		}
		
		rate = new double[Nsite];
		pos << "site-specific rates\n";
		for (int i=0; i<Nsite; i++)	{
			rate[i] = rnd::GetRandom().Gamma(alpha,alpha);
			pos << rate[i] << '\t';
		}
		pos << '\n';
		
		taxonset = protdata->GetTaxonSet();

		// check whether tree and data fit together
		tree->RegisterWith(taxonset);
		Ntaxa = taxonset->GetNtaxa();

		cerr << "Nsite: " << Nsite << '\n';
		cerr << "Ntaxa: " << Ntaxa << '\n';

		CreateMatrices();
		UpdateMatrices();

		totlength = 0;
		RecursiveSetBranchLengths(GetRoot());

		cerr << "length  : " << totlength << '\n';
	}

	void CreateMatrices()	{
		Q = new GTRSubMatrix*[Nsite];
		for (int i=0; i<Nsite; i++)	{
			int gene = genealloc[i];
			Q[i] = new GTRSubMatrix(Naa,rr[generralloc[gene]],stat[alloc[i]],true); 
		}
	}

	void UpdateMatrices()	{
		for (int i=0; i<Nsite; i++)	{
			Q[i]->CorruptMatrix();
			Q[i]->UpdateMatrix();
		}
	}

	const Tree* GetTree()	const {
		return tree;
	}

	Tree* GetTree()	{
		return tree;
	}

	Link* GetRoot()	{
		return GetTree()->GetRoot();
	}

	const Link* GetRoot()	const {
		return GetTree()->GetRoot();
	}

	string GetNodeName(const Link* link)	const {
		return GetTree()->GetNodeName(link);
	}

	string GetBranchName(const Link* link)	const {
		if (link->isRoot())	{
			return "";
		}
		return GetTree()->GetBranchName(link);
	}

	void RecursiveSetBranchLengths(const Link* from)	{
		if (! from->isRoot())	{
			double l = atof(from->GetBranch()->GetName().c_str());
			if (l<=0)	{
				cerr << "null branch length : " << from->GetBranch()->GetName() << '\n';
				exit(1);
			}
			bl[from->GetBranch()] = l;
			totlength += l;
		}
		else	{
			bl[from->GetBranch()] = 0;
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveSetBranchLengths(link->Out());
		}
	}

	void Simulate()	{
		count = 0;
		RecursiveSimulate(GetRoot());
		cerr << '\n';
	}

	void RecursiveSimulate(const Link* from)	{

		cerr << '.';
		if (from->isRoot())	{
			RootSimulate();
		}
		else	{
			BranchSimulate(from);
		}
		StoreCurrentValues(from);

		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveSimulate(link->Out());
		}
	}

	void StoreCurrentValues(const Link* from)	{

		// store new values
		int* newseq = new int[Nsite];
		for (int i=0; i<Nsite; i++)	{
			newseq[i] = currentseq[i];
		}
		nodeseq[from->GetNode()] = newseq;
	}

	void RootSimulate()	{
		for (int i=0; i<Nsite; i++)	{
			currentseq[i] = Q[i]->DrawFromStationary();
		}
	}

	void BranchSimulate(const Link* from)	{


		int* ancseq = nodeseq[from->Out()->GetNode()];
		for (int i=0; i<Nsite; i++)	{
			currentseq[i] = ancseq[i];
		}

		for (int i=0; i<Nsite; i++)	{
			double l = bl[from->GetBranch()] * mu * rate[i];
			int state = currentseq[i];
			double t = 0;

			while (t < l)	{

				double dt = Q[i]->DrawWaitingTime(state);
				t += dt;
				if (t < l)	{
					int newstate = Q[i]->DrawOneStep(state);
					count++;
					state = newstate;
				}
			}
			currentseq[i] = state;
		}
	}

	void WriteSimu()	{

		// dataset
		int** data = new int*[Ntaxa];
		string* names = new string[Ntaxa];
		int n = 0;
		RecursiveMakeData(GetRoot(),data,names,n);
		if (n != Ntaxa)	{
			cerr << "error when making sequence alignment: wrong number of taxa : " << n << '\t' << " instead of " << Ntaxa << '\n';
			exit(1);
		}

		cerr << "make new ali\n";
		SequenceAlignment* protali = new SequenceAlignment(data,names,Nsite,statespace,taxonset);
		if (mask)	{
			protali->Mask(protdata);
		}

		cerr << "write ali\n";
		ofstream callos((basename + ".ali").c_str());
		protali->ToStream(callos);
	
		ofstream sos((basename + ".summary").c_str());

		sos << "tot number of subs per site : " << ((double) count) / Nsite << '\n';
		sos << '\n';
		sos << "fraction inv colums: " << ((double) protali->GetNumberConstantColumns()) / Nsite << '\n';
		sos << "ref frac inv colums: " << ((double) protdata->GetNumberConstantColumns()) / Nsite << '\n';
		sos << '\n';
		sos << "mean diversity     : " << protali->GetMeanDiversity() << '\n';
		sos << "ref  diversity     : " << protdata->GetMeanDiversity() << '\n';
		sos << '\n';
		sos << "mean pairwise diff : " << protali->GetMeanPairwiseDiff() << '\n';
		sos << "ref  pairwise diff : " << protdata->GetMeanPairwiseDiff() << '\n';
		sos << '\n';

		delete[] data;
		delete[] names;
	}

	void RecursiveMakeData(const Link* from, int** data, string* names, int& n)	{
		if (from->isLeaf())	{
			names[n] = from->GetNode()->GetName();
			data[n] = nodeseq[from->GetNode()];
			n++;
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveMakeData(link->Out(), data, names, n);
		}
	}

	int count;
	int mask;

	int Ntaxa;
	int Nsite;
	int* currentseq;
	int Nrr;

	int Ngene;
	int* genealloc;
	int* genesize;
	int* genefirst;

	int* alloc;
	double* rate;

	int withempfreq;
	double* mixweight;
	// double** geneweight;
	int* generralloc;

	GTRSubMatrix** Q;

	int Nrrcat;
	int Ncat;
	
	double alpha;

	double** rr;
	double** stat;

	double mu;

	Tree* tree;

	// with branch lengths, measured in time
	map<const Branch*, double> bl;
	double totlength;

	map<const Node*, int*> nodeseq;

	SequenceAlignment* protdata;
	SequenceAlignment** genedata;
	const TaxonSet* taxonset;
	StateSpace* statespace;

	string basename;

};

int main(int argc, char* argv[])	{

	string datafile = "";
	string treefile = "";
	string partitionfile = "";
	string paramfile = "";
	string basename = "";
	int mask = 0;

	try	{

		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];
			if ((s == "-d") || (s == "-data"))	{
				i++;
				datafile = argv[i];
			}
			else if (s == "-m")	{
				mask = 1;
			}
			else if ((s == "-t") || (s == "-tree"))	{
				i++;
				treefile = argv[i];
			}
			else if (s == "-part")	{
				i++;
				partitionfile = argv[i];
			}
			else if ((s == "-p") || (s == "-param"))	{
				i++;
				paramfile = argv[i];
			}
			else	{
				if (i != (argc -1))	{
					throw(0);
				}
				basename = argv[i];
			}
			i++;
		}
		if ((datafile == "") || (treefile == "") || (paramfile == "") || (basename == ""))	{
			throw(0);
		}
	}
	catch(...)	{
		cerr << '\n';
		cerr << "simucodon -p <paramfile> -t <treefile> -d <datafile> [-m] <basename>\n";
		cerr << '\n';
		exit(1);
	}

	cerr << "new sim\n";
	Simulator* sim = new Simulator(datafile,treefile,partitionfile,paramfile,mask,basename);

	cerr << "simulate\n";
	sim->Simulate();

	cerr << "write simu\n";
	sim->WriteSimu();

	cerr << "simu written in " << basename << '\n';
	cerr << '\n';

}

