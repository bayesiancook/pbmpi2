
#include "Tree.h"
#include "SequenceAlignment.h"
#include "Random.h"
#include <map>
#include <iostream>
#include <sstream>
#include <vector>

#include "StringStreamUtils.h"
#include "GTRSubMatrix.h"

#include "Parallel.h"
MPI_Datatype Propagate_arg;

class Simulator : public NewickTree {

	public:

	Simulator(string datafile, string treefile, string partitionfile, string rrfile, string paramfile, string profilefile, int inmask, int inrandomprofiles, int innfold, int ngene, int targetnsite, string inbasename)	{

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
		statespace = protdata->GetStateSpace();

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
		int mixprior = 0;
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

			if (profilefile != "None")	{
                ifstream is(profilefile.c_str());
                int nprofile;
                is >> nprofile;
                double** profile = new double*[nprofile];
                for (int k=0; k<nprofile; k++)	{
                    profile[k] = new double[Naa];
                    int tmp;
                    is >> tmp;
                    for (int i=0; i<Naa; i++)	{
                        is >> profile[k][i];
                    }
                }

				if (randomprofiles)	{
					for (int k=0; k<Nsite; k++)	{
						int cat = (int) (nprofile * rnd::GetRandom().Uniform());
						for (int i=0; i<Naa; i++)	{
							stat[k][i] = profile[cat][i];
						}
					}
				}
				else	{
                    if (Nsite > nprofile)   {
                        cerr << "error: not enough profiles in file\n";
                        exit(1);
                    }
					for (int k=0; k<Nsite; k++)	{
						for (int i=0; i<Naa; i++)	{
							stat[k][i] = profile[k][i];
						}
					}
				}
                for (int k=0; k<nprofile; k++)	{
                    delete[] profile[k];
                }
                delete[] profile;
			}
			else	{
				if (Ngene > 1)	{
					for (int gene=0; gene<Ngene; gene++)	{
						genedata[gene]->GetSiteEmpiricalFreq(stat + genefirst[gene],pseudocount,focus);
					}
				}
				else	{
					protdata->GetSiteEmpiricalFreq(stat,pseudocount,focus);
				}
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
			else if (tmp == "+P")	{
				mixprior = 1;
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
			
			if (mixprior)	{
				Ncat0 = Ncat;
				Ncat = Nsite;
				prmis >> tmp;
				if (tmp != "Concentrations")	{
					cerr << "error: missing Concentrations keyword\n";
					cerr << tmp << '\n';
					exit(1);
				}
				double concfactor = 1.0;
				prmis >> concfactor;
				conc = new double[Ncat0];
				for (int k=0; k<Ncat0; k++)	{
					prmis >> conc[k];
					conc[k] *= concfactor;
				}
				prmis >> tmp;
				if (tmp != "Stats")	{
					cerr << "error: missing Weights keyword\n";
					cerr << tmp << '\n';
					exit(1);
				}
				stat0 = new double*[Ncat0];
				for (int k=0; k<Ncat0; k++)	{
					stat0[k] = new double[Naa];
				}
				for (int k=0; k<Ncat0; k++)	{
					double tot = 0;
					for (int i=0; i<Naa; i++)	{
						prmis >> stat0[k][i];
						tot += stat0[k][i];
					}
					for (int i=0; i<Naa; i++)	{
						stat0[k][i] /= tot;
					}
				}
				alloc = new int[Nsite];
				for (int i=0; i<Nsite; i++)	{
					alloc[i] = i;
				}
				stat = new double*[Nsite];
				for (int k=0; k<Nsite; k++)	{
					stat[k] = new double[Naa];
				}
				for (int k=0; k<Nsite; k++)	{
					int tmp = rnd::GetRandom().FiniteDiscrete(Ncat0,mixweight);
					double tot = 0;
					for (int i=0; i<Naa; i++)	{
						stat[k][i] = rnd::GetRandom().sGamma(stat0[tmp][i] * conc[tmp]);
						tot += stat[k][i];
					}
					for (int i=0; i<Naa; i++)	{
						stat[k][i] /= tot;
					}
				}
			}
			else	{
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
		}

		// relative exchangeabilites

		Nrrcat = 0;
		rr = 0;
		Nrr = Naa * (Naa-1) / 2;

		prmis >> tmp;
		if (tmp != "RR")	{
			cerr << "error: missing RR keyword\n";
			cerr << tmp << '\n';
			exit(1);
		}
		prmis >> tmp;
		if (tmp == "Poisson")	{
			Nrrcat = 1;
		}
		else if (tmp == "LG")	{
			Nrrcat = 1;
			rr = new double*[Nrrcat];
			rr[0] = new double[Nrr];
			for (int i=0; i<Nrr; i++)	{
				rr[0][i] = LG_RR[i];
			}
		}
		else if (tmp == "WAG")	{
			Nrrcat = 1;
			rr = new double*[Nrrcat];
			rr[0] = new double[Nrr];
			for (int i=0; i<Nrr; i++)	{
				rr[0][i] = WAG_RR[i];
			}
		}
		else if (tmp == "Partition")	{
			Nrrcat = 5;
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
		}

		generralloc = new int[Ngene];
		if (Ngene == 1)	{
			generralloc[0] = 0;
		}
		else	{
			if (rrfile != "None")	{
				int fromngene;
				ifstream is(rrfile.c_str());
				is >> fromngene;
				string* rralloc = new string[fromngene];
				for (int gene=0; gene<fromngene; gene++)	{
					is >> rralloc[gene];
				}
				int genecount = 0;
				for (int gene=0; gene<Ngene; gene++)	{
					string tmp = rralloc[genecount];
					genecount++;
					if (genecount == fromngene)	{
						genecount = 0;
					}
					if (tmp == "WAG")	{
						generralloc[gene] = 0;
					}
					else if (tmp == "WAGF")	{
						generralloc[gene] = 0;
					}
					else if (tmp == "BLOSUM62")	{
						generralloc[gene] = 1;
					}
					else if (tmp == "BLOSUM62F")	{
						generralloc[gene] = 1;
					}
					else if (tmp == "JTT")	{
						generralloc[gene] = 2;
					}
					else if (tmp == "JTTF")	{
						generralloc[gene] = 2;
					}
					else if (tmp == "DCMUT")	{
						generralloc[gene] = 3;
					}
					else if (tmp == "DCMUTF")	{
						generralloc[gene] = 3;
					}
					else if (tmp == "MTREV")	{
						generralloc[gene] = 4;
					}
					else if (tmp == "MTREVF")	{
						generralloc[gene] = 4;
					}
					else if (tmp == "DAYHOFF")	{
						generralloc[gene] = 3;
					}
					else if (tmp == "DAYHOFFF")	{
						generralloc[gene] = 3;
					}
					else	{
						cerr << "error: did not recognize rr type: " << tmp << '\n';
						exit(1);
					}
				}
				delete[] rralloc;
			}
			else	{
                cerr << "random allocation of rr to genes\n";
				for (int gene=0; gene<Ngene; gene++)	{
					generralloc[gene] = (int) (Nrrcat * rnd::GetRandom().Uniform());
				}
			}
		}

		ofstream pos((basename + ".param").c_str());

		if (Nrrcat == 5)	{
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
		}

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
		
		ofstream pros((basename + ".trueprofiles").c_str());
		for (int i=0; i<Nsite; i++)	{
			for (int k=0; k<Naa; k++)	{
				pros << stat[alloc[i]][k] << '\t';
			}
			pros << '\n';
		}

		rate = new double[Nsite];
		pos << "site-specific rates\n";
		for (int i=0; i<Nsite; i++)	{
			rate[i] = rnd::GetRandom().Gamma(alpha,alpha);
			pos << rate[i] << '\t';
		}
		pos << '\n';
		
		if (rr)	{
			CreateMatrices();
			UpdateMatrices();
		}

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
		if (rr)	{
			for (int i=0; i<Nsite; i++)	{
				currentseq[i] = Q[i]->DrawFromStationary();
			}
		}
		else	{
			for (int i=0; i<Nsite; i++)	{
				currentseq[i] = rnd::GetRandom().DrawFromDiscreteDistribution(stat[alloc[i]],Naa);
			}
		}
	}

	void BranchSimulate(const Link* from)	{


		int* ancseq = nodeseq[from->Out()->GetNode()];
		for (int i=0; i<Nsite; i++)	{
			currentseq[i] = ancseq[i];
		}

		for (int i=0; i<Nsite; i++)	{
			double l = bl[from->GetBranch()] * mu * rate[i];

			if (rr)	{
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
			else	{
				double u = rnd::GetRandom().Uniform();
				if (u > exp(-l))	{
					currentseq[i] = rnd::GetRandom().DrawFromDiscreteDistribution(stat[alloc[i]],Naa);
				}
			}

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

		SequenceAlignment* protali = new SequenceAlignment(data,names,Nsite,statespace,taxonset);
		cerr << "ok\n";
		if (mask)	{
			protali->Mask(protdata);
		}

		ofstream callos((basename + ".ali").c_str());
		protali->ToStream(callos);
	
		ofstream sos((basename + ".summary").c_str());

		sos << "tot number of subs per site : " << ((double) count) / Nsite << '\n';
		sos << '\n';
		sos << "fraction inv colums: " << ((double) protali->GetNumberConstantColumns()) / Nsite << '\n';
		sos << "ref frac inv colums: " << ((double) protdata->GetNumberConstantColumns()) / protdata->GetNsite() << '\n';
		sos << '\n';
		sos << "mean diversity     : " << protali->GetMeanDiversity() << '\n';
		sos << "ref  diversity     : " << protdata->GetMeanDiversity() << '\n';
		sos << '\n';
		sos << "mean pairwise diff : " << protali->GetMeanPairwiseDiff() << '\n';
		sos << "ref  pairwise diff : " << protdata->GetMeanPairwiseDiff() << '\n';
		sos << '\n';

		ofstream pfos((basename + ".txt.pf").c_str());

		pfos << "## BRANCHLENGTHS: linked | unlinked ##\n";
		pfos << "branchlengths = linked;\n";
		pfos << "\n";
		pfos << "## MODELS OF EVOLUTION for PartitionFinder: all | raxml | mrbayes | <list> ##\n";
		pfos << "##              for PartitionFinderProtein: all_protein | <list> ##\n";
		pfos << "models = all_protein;\n";
		pfos << "\n";
		pfos << "# MODEL SELECCTION: AIC | AICc | BIC #\n";
		pfos << "model_selection = BIC;\n";
		pfos << "\n";
		pfos << "## DATA BLOCKS: see manual for how to define ##\n";
		pfos << "[data_blocks]\n";
		for (int gene=0; gene<Ngene; gene++)	{
			pfos << "gene" << gene << " = " << genefirst[gene]+1 << "-" << genefirst[gene]+genesize[gene] << ";\n";
		}
		pfos << "\n";
		pfos << "## SCHEMES, search: all | user | greedy ##\n";
		pfos << "[schemes]\n";
		pfos << "search = greedy;\n";
		pfos << "#search = hcluster;\n";
		pfos << "\n";
		pfos << "#user schemes go here if search=user. See manual for how to define.#\n";

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
	int randomprofiles;
	int nfold;

	int Ntaxa;
	int Nsite;
	int* currentseq;
	int Nrr;

	int Ngene;
	int* genealloc;
    vector<int> genesize;
    vector<int> genefirst;
    vector<SequenceAlignment*> genedata;

    /*
	int* genesize;
	int* genefirst;
	SequenceAlignment** genedata;
    */

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
	double** stat0;
	int Ncat0;
	double* conc;

	double mu;

	Tree* tree;

	// with branch lengths, measured in time
	map<const Branch*, double> bl;
	double totlength;

	map<const Node*, int*> nodeseq;

	SequenceAlignment* protdata;
	const TaxonSet* taxonset;
	StateSpace* statespace;

	string basename;

};

int main(int argc, char* argv[])	{

	string datafile = "";
	string treefile = "";
	string partitionfile = "None";
	string profilefile = "None";
	string paramfile = "";
	string rrfile = "None";
	string basename = "";
	int mask = 1;
	int randomprofiles = 0;
	int nfold = 1;
	int ngene = -1;
    int targetnsite = -1;

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
			else if (s == "-nomask")	{
				mask = 0;
			}
			else if ((s == "-t") || (s == "-tree"))	{
				i++;
				treefile = argv[i];
			}
			else if (s == "-part")	{
				i++;
				partitionfile = argv[i];
			}
			else if (s == "-rr")	{
				i++;
				rrfile = argv[i];
			}
			else if ((s == "-p") || (s == "-param"))	{
				i++;
				paramfile = argv[i];
			}
			else if ((s == "-f") || (s == "-freq"))	{
				i++;
				profilefile = argv[i];
			}
			else if (s == "-randprof")	{
				randomprofiles = 1;
			}
			else if (s == "-nfold")	{
				i++;
				nfold = atoi(argv[i]);
			}
			else if (s == "-ngene")	{
				i++;
				ngene = atoi(argv[i]);
			}
            else if (s == "-targetnsite") {
                i++;
                targetnsite = atoi(argv[i]);
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
		cerr << "simucodon -p <paramfile> -t <treefile> -d <datafile> [-f <aaprofilefile> -part <partitionfile> -m] <basename>\n";
		cerr << '\n';
		exit(1);
	}

	cerr << "new sim\n";
	Simulator* sim = new Simulator(datafile,treefile,partitionfile,rrfile,paramfile,profilefile,mask,randomprofiles,nfold,ngene,targetnsite,basename);

	cerr << "simulate\n";
	sim->Simulate();

	cerr << "write simu\n";
	sim->WriteSimu();

	cerr << "simu written in " << basename << '\n';
	cerr << '\n';

}

