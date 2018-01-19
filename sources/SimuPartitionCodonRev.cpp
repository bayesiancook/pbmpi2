
#include "Tree.h"
#include "SequenceAlignment.h"
#include "CodonSequenceAlignment.h"
#include "ProteinSequenceAlignment.h"
#include "Random.h"
#include <map>
#include <iostream>
#include <sstream>

#include "StringStreamUtils.h"
#include "CodonSubMatrix.h"

#include "Parallel.h"
MPI_Datatype Propagate_arg;

class Simulator : public NewickTree {

	public:

	Simulator(string datafile, string treefile, string partitionfile, string paramfile, string profilefile, int inmask, string inbasename)	{
	Simulator(string datafile, string treefile, string partitionfile, string rrfile, string paramfile, string profilefile, int inmask, int inrandomprofiles, int innfold, int ngene, int targetnsite, string inbasename)	{

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

		mask = inmask;
		randomprofiles = inrandomprofiles;
		nfold = innfold;

		basename = inbasename; 

		tree = new Tree(treefile);

		protdata = new FileSequenceAlignment(datafile);
		if (protdata->GetNstate() != Naa)	{
			cerr << "error: should be protein datafile\n";
			exit(1);
		}

		taxonset = protdata->GetTaxonSet();
		codonstatespace = new CodonStateSpace(Universal);

		// check whether tree and data fit together
		tree->RegisterWith(taxonset);
		Ntaxa = taxonset->GetNtaxa();

		if (partitionfile != "None")	{

            double meangenesize = 0;
			ifstream pis(partitionfile.c_str());
			int fromngene;
			pis >> fromngene;
			int* fromgenesize = new int[fromngene];
			int* fromgenefirst = new int[fromngene];
			SequenceAlignment** fromgenedata = new SequenceAlignment*[fromngene];
			int count = 0;
			for (int gene=0; gene<fromngene; gene++)	{
				pis >> fromgenesize[gene];
                meangenesize += fromgenesize[gene];
				fromgenefirst[gene] = count;
				count += fromgenesize[gene];
				fromgenedata[gene] = new SequenceAlignment(protdata,fromgenefirst[gene],fromgenesize[gene]);
			}
            meangenesize /= fromngene;

            if (targetnsite != -1)  {

                cerr << "random assignment until reaching target number of sites: " << targetnsite << '\n';
                Nsite = 0;
                Ngene = 0;
                while (Nsite < targetnsite) {
                    int alloc = (int) (rnd::GetRandom().Uniform() * fromngene);
                    genesize.push_back(fromgenesize[alloc]);
                    genefirst.push_back(Nsite);
                    genedata.push_back(new SequenceAlignment(fromgenedata[alloc]));
                    Nsite += genesize[Ngene];
                    Ngene++;
                }
            }


            else if (ngene != -1) {
                cerr << "random assignment until reaching target number of genes: " << ngene << '\n';
                Ngene = ngene;
                genesize.assign(Ngene,0);
                genefirst.assign(Ngene,0);
                genedata.assign(Ngene,(SequenceAlignment*)0);
                Nsite = 0;
                for (int gene=0; gene<Ngene; gene++)	{
                    int alloc = (int) (rnd::GetRandom().Uniform() * fromngene);
                    genesize[gene] = fromgenesize[alloc];
                    genefirst[gene] = Nsite;
                    Nsite += genesize[gene];
                    genedata[gene] = new SequenceAlignment(fromgenedata[alloc]);
                }
            }

            else    {
                Nsite = protdata->GetNsite() * nfold;
                Ngene = fromngene*nfold;
                if (nfold > 1)  {
                    cerr << "rolling circle: " << nfold << " x " << fromngene << " = " << Ngene << " genes in total\n";
                }
                else    {
                    cerr << "simulation based on template alignment and partition: same size, same number of genes\n";
                }
                genesize.assign(Ngene,0);
                genefirst.assign(Ngene,0);
                genedata.assign(Ngene,(SequenceAlignment*)0);
                int count = 0;
                int genecount = 0;
                for (int gene=0; gene<Ngene; gene++)	{
                    genesize[gene] = fromgenesize[genecount];
                    genefirst[gene] = count;
                    count += genesize[gene];
                    genedata[gene] = new SequenceAlignment(fromgenedata[genecount]);
                    genecount++;
                    if (genecount == fromngene)	{
                        genecount = 0;
                    }
                }
                if (count != Nsite)	{
                    cerr << "error: non matching total size\n";
                    exit(1);
                }

            }

            genealloc = new int[Nsite];
            for (int gene=0; gene<Ngene; gene++)	{
                for (int i=0; i<genesize[gene]; i++)	{
                    genealloc[genefirst[gene] + i] = gene;
                }
            }

            delete[] fromgenesize;
            delete[] fromgenefirst;
            for (int gene=0; gene<fromngene; gene++)	{
                delete fromgenedata[gene];
            }
            delete[] fromgenedata;
		}
		else	{
			Ngene = 1;
            Nsite = protdata->GetNsite() * nfold;
            genesize.assign(Ngene,0);
            genefirst.assign(Ngene,0);
            genedata.assign(Ngene,(SequenceAlignment*)0);
			genealloc = new int[Nsite];
			genesize[0] = Nsite;
			genefirst[0] = 0;
			genedata[0] = protdata;
			for (int i=0; i<Nsite; i++)	{
				genealloc[i] = 0;
			}
		}

        cerr << "total number of genes: " << Ngene << '\n';
        cerr << "total number of sites: " << Nsite << '\n';
		cerr << "total number of taxa : " << Ntaxa << '\n';

		currentseq = new int[Nsite];

		// get parameters from file
		ifstream prmis(paramfile.c_str());

		cerr << "Param file : " << paramfile << '\n';
		// read file
		string tmp;
		prmis >> tmp;

		if (tmp != "PseudoCount")	{
			cerr << "error: missing PseudoCount keyword\n";
			cerr << tmp << '\n';
			exit(1);
		}
		prmis >> pseudocount;
		prmis >> focus;

		aafreq = new double*[Nsite];
		for (int i=0; i<Nsite; i++)	{
			aafreq[i] = new double[Naa];
		}
		if (profilefile != "None")	{
            ifstream is(profilefile.c_str());
            int nprofile;
            is >> nprofile;
            if (randomprofiles) {
                profile = new double*[Nsite];
                for (int i=0; i<Nsite; i++)	{
                    aafreq[i] = new double[Naa];
                }
                for (int k=0; k<Nsite; k++)	{
                    int tmp;
                    is >> tmp;
                    for (int i=0; i<Naa; i++)	{
                        is >> aafreq[k][i];
                    }
                }
            }
            else    {
                ifstream is(profilefile.c_str());
                for (int k=0; k<Nsite; k++)	{
                    int tmp;
                    is >> tmp;
                    if (tmp < Nsite)    {
                        cerr << "error: not enough profiles in file\n";
                        exit(1);
                    }
                    for (int i=0; i<Naa; i++)	{
                        is >> aafreq[k][i];
                    }
                }
            }
		}
		else	{
			// make gene-specific empirical freqs
			for (int gene=0; gene<Ngene; gene++)	{
				genedata[gene]->GetSiteEmpiricalFreq(aafreq + genefirst[gene],pseudocount,focus);
			}
		}

		// nucleotide matrix
		Nrr = Nnuc * (Nnuc-1) / 2;
		nucrr = new double[Nrr];
		nucstat = new double[Nnuc];

		prmis >> tmp;
		if (tmp != "RelRates")	{
			cerr << "error: missing RelRates keyword\n";
			cerr << tmp << '\n';
			exit(1);
		}
		for (int i=0; i<Nrr; i++)	{
			prmis >> nucrr[i];
		}
		prmis >> rralpha;

		prmis >> tmp;
		if (tmp != "EqGC")	{
			cerr << "error: missing EqGC keyword\n";
			cerr << tmp << '\n';
			exit(1);
		}
		prmis >> gc;
		prmis >> gcalpha;
		if ((gc - gcalpha < 0) || (gc + gcalpha > 1))	{
			cerr << "error: gc overflow\n";
			exit(1);
		}
		nucstat[0] = 0.5 * (1-gc);
		nucstat[1] = 0.5 * gc;
		nucstat[2] = 0.5 * gc;
		nucstat[3] = 0.5 * (1-gc);

		double* meannucrr = new double[Nrr];
		double* varnucrr = new double[Nrr];
		double* minnucrr = new double[Nrr];
		double* maxnucrr = new double[Nrr];
		for (int i=0; i<Nrr; i++)	{
			meannucrr[i] = 0;
			varnucrr[i] = 0;
			minnucrr[i] = 100;
			maxnucrr[i] = 0;
		}
		double meangc = 0;
		double mingc = 1;
		double maxgc = 0;
		double vargc = 0;

		genenucrr = new double*[Ngene];
		genenucstat = new double*[Ngene];
		for (int gene=0; gene<Ngene; gene++)	{
			genenucrr[gene] = new double[Nrr];
			genenucstat[gene] = new double[Nnuc];
			for (int i=0; i<Nrr; i++)	{
				if (rralpha)	{
					double tmp = rnd::GetRandom().sGamma(rralpha * nucrr[i]);
					genenucrr[gene][i] = tmp;
					if (minnucrr[i] > tmp)	{
						minnucrr[i] = tmp;
					}
					if (maxnucrr[i] < tmp)	{
						maxnucrr[i] = tmp;
					}
					meannucrr[i] += tmp;
					varnucrr[i] += tmp * tmp;
				}
				else	{
					genenucrr[gene][i] = nucrr[i];
				}
			}
			double tot = 0;
			double tmpgc = gc + 2 * gcalpha * (rnd::GetRandom().Uniform() - 0.5);
			genenucstat[gene][1] = genenucstat[gene][2] = 0.5 * tmpgc;
			genenucstat[gene][0] = genenucstat[gene][3] = 0.5 * (1 - tmpgc);
			/*
			for (int i=0; i<Nnuc; i++)	{
				genenucstat[gene][i] = rnd::GetRandom().sGamma(gcalpha * nucstat[i]);
				tot += genenucstat[gene][i];
			}
			for (int i=0; i<Nnuc; i++)	{
				genenucstat[gene][i] /= tot;
			}
			double tmpgc = genenucstat[gene][1] + genenucstat[gene][2];
			*/
			meangc += tmpgc;
			vargc += tmpgc * tmpgc;
			if (mingc > tmpgc)	{
				mingc = tmpgc;
			}
			if (maxgc < tmpgc)	{
				maxgc = tmpgc;
			}
		}

		meangc /= Ngene;
		vargc /= Ngene;
		vargc -= meangc * meangc;
		cerr << "GC across genes: " << meangc << " +/- " << sqrt(vargc) << '\t' << mingc << '\t' << maxgc << '\n';
		for (int i=0; i<Nrr; i++)	{
			meannucrr[i] /= Ngene;
			varnucrr[i] /= Ngene;
			varnucrr[i] -= meannucrr[i] * meannucrr[i];
		}
		if (rralpha)	{
			cerr << "RR across genes:\n";
			for (int i=0; i<Nrr; i++)	{
				cerr << meannucrr[i] << " +/- " << sqrt(varnucrr[i]) << '\t' << minnucrr[i] << '\t' << maxnucrr[i] << '\n';
			}
			cerr << '\n';
		}

		omega = 1;

		prmis >> tmp;
		if (tmp != "Ne")	{
			cerr << "error: missing Ne keyword\n";
			cerr << tmp << '\n';
			exit(1);
		}
		prmis >> Ne;

		CreateMatrices();
		UpdateMatrices();
		nonsynrate = ComputeAverageNonSynRate();
		cerr << "average dN:" << nonsynrate << '\n';

		prmis >> tmp;
		if (tmp != "mu")	{
			cerr << "error: missing Scale keyword\n";
			cerr << tmp << '\n';
			exit(1);
		}
		prmis >> mu;

		totlength = 0;
		RecursiveSetBranchLengths(GetRoot());

		cerr << "length  : " << totlength << '\n';
	}

	void CreateMatrices()	{
		Q = new HBAACodonMutSelProfileSubMatrix*[Nsite];
		for (int i=0; i<Nsite; i++)	{
			Q[i] = new HBAACodonMutSelProfileSubMatrix(codonstatespace,genenucrr[genealloc[i]],genenucstat[genealloc[i]],aafreq[i],&omega,&Ne,false);
		}
	}

	void UpdateMatrices()	{
		for (int i=0; i<Nsite; i++)	{
			Q[i]->CorruptMatrix();
			Q[i]->UpdateMatrix();
		}
	}

	double ComputeAverageNonSynRate()	{
		double total = 0;
		for (int i=0; i<Nsite; i++)	{
			total += Q[i]->NonSynRate();
		}
		return total / Nsite;
	}

	static int rrindex(int i, int j, int nstate)	{
		return (i<j) ? (2 * nstate - i - 1) * i / 2 + j - i - 1 : (2 * nstate - j - 1) * j / 2 + i - j - 1 ;
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
		dscount = 0;
		dncount = 0;
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

		double l = bl[from->GetBranch()] * mu / nonsynrate;

		int* ancseq = nodeseq[from->Out()->GetNode()];
		for (int i=0; i<Nsite; i++)	{
			currentseq[i] = ancseq[i];
		}

		for (int i=0; i<Nsite; i++)	{
			int state = currentseq[i];
			double t = 0;

			while (t < l)	{

				double dt = Q[i]->DrawWaitingTime(state);
				t += dt;
				if (t < l)	{
					int newstate = Q[i]->DrawOneStep(state);
					count++;
					if (codonstatespace->Synonymous(state,newstate))	{
						dscount++;
					}
					else	{
						dncount++;
					}
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
		CodonSequenceAlignment* codonali = new CodonSequenceAlignment(data,names,Nsite,codonstatespace,taxonset);

		cerr << "write ali\n";
		ofstream callos((basename + ".ali").c_str());
		codonali->ToStream(callos);
	
		ProteinSequenceAlignment* protali = new ProteinSequenceAlignment(codonali);
		if (mask)	{
			protali->Mask(protdata);
		}
		ofstream protos((basename + "_prot.ali").c_str());
		protali->ToStream(protos);

		ofstream sos((basename + ".summary").c_str());

		sos << "tot number of subs per site : " << ((double) count) / Nsite << '\n';
		sos << "tot number of reps per site : " << ((double) dncount) / Nsite << '\n';
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
	int dscount;
	int dncount;

	int mask;

	double pseudocount;
	int focus;
	double** aafreq;
	HBAACodonMutSelProfileSubMatrix** Q;

	int Ntaxa;
	int Nsite;
	int* currentseq;
	int Nrr;

	int Ngene;
	int* genealloc;
	int* genesize;
	int* genefirst;

	double* nucrr;
	double rralpha;
	double* nucstat;
	// gc
	double gc;
	double gcalpha;

	double** genenucrr;
	double** genenucstat;

	double nonsynrate;
	double mu;
	double omega;
	double Ne;

	Tree* tree;

	// with branch lengths, measured in time
	map<const Branch*, double> bl;
	double totlength;

	map<const Node*, int*> nodeseq;

	SequenceAlignment* protdata;
	SequenceAlignment** genedata;
	const TaxonSet* taxonset;
	CodonStateSpace* codonstatespace;

	string basename;

};

int main(int argc, char* argv[])	{

	string datafile = "";
	string treefile = "";
	string partitionfile = "";
	string profilefile = "None";
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
			else if ((s == "-f") || (s == "-freq"))	{
				i++;
				profilefile = argv[i];
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
		cerr << "simucodon -p <paramfile> -t <treefile> -d <datafile> [-f <profilefile> -m] <basename>\n";
		cerr << '\n';
		exit(1);
	}

	cerr << "new sim\n";
	Simulator* sim = new Simulator(datafile,treefile,partitionfile,paramfile,profilefile,mask,basename);

	cerr << "simulate\n";
	sim->Simulate();

	cerr << "write simu\n";
	sim->WriteSimu();

	cerr << "simu written in " << basename << '\n';
	cerr << '\n';

}

