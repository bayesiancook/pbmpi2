
#ifndef CPPH
#define CPPH

#include "Tree.h"
#include "Random.h"
#include <vector>
#include <algorithm>

typedef pair<vector<int>,double> bpentry;

class BpList	{

	public:


	static bool compentry(bpentry e1, bpentry e2)	{
		return (e1.second > e2.second);
	}

	BpList(const TaxonSet* intaxset)	{
		taxset = intaxset;
		totweight = 0;

		// number of unrooted trees of size n is equal to number of rooted trees of size n-1
		if (taxset)	{
			logNTrees = GetLogNPlantedTrees(GetNtaxa()-1);
		}
	}

	void SetLogNTrees()	{
		logNTrees = GetLogNPlantedTrees(GetNtaxa()-1);
	}
		
	void Read(string name, int burnin, int every, int until)	{

		ifstream is(name.c_str());
		Tree* tree = new Tree(is);
		if (! taxset)	{
			taxset = tree->GetTaxonSet();
		}
		int count = 1;
		string line;
		while (! is.eof())	{
			is >> line;
			count++;
		}
		count--;

		if (until == -1)	{
			until = count;
		}

		if (burnin > until)	{
			burnin = until;
		}

		ifstream is2(name.c_str());
		for (int i=0; i<burnin; i++)	{
			is2 >> line;
		}

		totweight = 0;
		int cycle = 0;
		count = burnin;
		while (count < until)	{
			tree->ReadFromStream(is2);
			if (!cycle)	{
				Add(tree,1);
			}
			count++;
			cycle++;
			if (cycle == every)	{
				cycle = 0;
			}
		}
	}

	const TaxonSet* GetTaxonSet() {return taxset;}

	int GetTaxonIndex(string taxon) {
		int tmp = GetTaxonSet()->GetTaxonIndex(taxon);
		if (tmp == -1)	{
			cerr << "error: taxon not found\n";
			exit(1);
		}
		return tmp;
	}

	int GetNtaxa() {return GetTaxonSet()->GetNtaxa();}

	double GetWeight(vector<int> bp)	{
		if (! isObserved(bp))	{
			return 0;
		}
		return weight[bp] / totweight;
	}

	bool isObserved(vector<int> bp)	{
		map<vector<int>,double>::iterator i = weight.find(bp);
		return (i != weight.end());
	}

	virtual double GetLogWeight(vector<int> bp, double cutoff)	{
		double weight = (1 - cutoff) * GetWeight(bp) + cutoff * exp(GetLogPriorBpMass(bp));
		return log(weight);
	}

	virtual double GetLength(vector<int> bp)	{
		if (isObserved(bp))	{
			return length[bp] / weight[bp];
		}
		return 0;
	}

	double GetLogNPlantedTrees(int n)	{

		double tot = 0;
		for (int i=2; i<n; i++)	{
			tot += log(2*i-3);
		}
		return tot;
	}

	double GetLogPriorBpMass(vector<int> bp)	{

		int n = 0;
		for (int i=0; i<GetNtaxa(); i++)	{
			if (bp[i] == 1)	{
				n++;
			}
		}
		if ((n < 2) || (n > GetNtaxa() - 2))	{
			cerr << "error: bp is not informative\n";
			exit(1);
		}
			
		return GetLogNPlantedTrees(n) + GetLogNPlantedTrees(GetNtaxa() - n) - logNTrees;
	}

	void Add(Tree* tree, double inweight)	{

		RecursiveAdd(tree->GetRoot(), inweight);
		totweight += inweight;
	}

	virtual vector<int> RecursiveAdd(const Link* from, double inweight)	{

		vector<int> ret(GetNtaxa(),0);

		if (from->isLeaf())	{
			int index = GetTaxonIndex(from->GetNode()->GetName());
			if ((index < 0) || (index >= GetNtaxa()))	{
				cerr << "error in BpList::RecursiveAdd: taxon index out of bound\n";
				cerr << index << '\t' << from->GetNode()->GetName() << '\n';
				exit(1);
			}
			ret[index] = 1;
		}

		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			vector<int> tmp = RecursiveAdd(link->Out(), inweight);
			for (int i=0; i<GetNtaxa(); i++)	{
				ret[i] |= tmp[i];
			}
		}

		if ((! from->isRoot()) && (! from->isLeaf()))	{
			// polarize vector
			vector<int> cp = Polarize(ret);

			// and add it to bplist
			weight[cp] += inweight;
			length[cp] += inweight * from->GetBranch()->GetLength();
		}
		return ret;
	}

	virtual vector<int> Polarize(vector<int> in)	{

		vector<int> cp = in;
		if (cp[0])	{
			for (int i=0; i<GetNtaxa(); i++)	{
				cp[i] = 1-cp[i];
			}
		}
		return cp;
	}

	void Normalize()	{
		for (map<vector<int>,double>::iterator i=weight.begin(); i!=weight.end(); i++)	{
			i->second /= totweight;
			length[i->first] /= totweight;
		}
		totweight = 1.0;
	}

	void Merge(BpList* inlist)	{
		int count = 0;
		for (map<vector<int>,double>::iterator i=inlist->weight.begin(); i!=inlist->weight.end(); i++)	{
			weight[i->first] += i->second;
			length[i->first] += inlist->length[i->first];
			count++;
		}
		totweight += inlist->totweight;
	}

	virtual void ToStream(ostream& os)	{
		vector<bpentry> bpvec;
		for (map<vector<int>,double>::iterator i=weight.begin(); i!=weight.end(); i++)	{
			bpvec.push_back(bpentry(i->first,i->second));
		}
		sort(bpvec.begin(),bpvec.end(),&compentry);

		os << bpvec.size() << '\t' << totweight << '\n';
		for (unsigned int i=0; i<bpvec.size(); i++)	{
			vector<int> tmp = bpvec[i].first;
			for (int i=0; i<GetNtaxa(); i++)	{
				os << tmp[i];
			}
			os << '\t' << GetWeight(tmp) << '\t' << GetLength(tmp) << '\n';
		}
	}

	void FromStream(istream& is)	{

		weight.clear();
		length.clear();
		int size;
		is >> size >> totweight;
		for (int i=0; i<size; i++)	{
			vector<int> bp(GetNtaxa(),0);

			string tmp;
			is >> tmp;
			if (tmp.length() != GetNtaxa())	{
				cerr << "error in BpList::FromStream: " << tmp << '\n';
				exit(1);
			}
			for (int j=0; j<GetNtaxa(); j++)	{
				char c = tmp[j];
				if (c == '0')	{
					bp[j] = 0;
				}
				else if (c == '1')	{
					bp[j] = 1;
				}
				else if (c == '2')	{
					bp[j] = 2;
				}
				else	{
					cerr << "error in BpList::FromStream: invalid partition\n";
					exit(1);
				}
			}

			double w, l;
			is >> w >> l;
			weight[bp] = w;
			length[bp] = l;
		}
	}

	virtual void ToStreamSexy(ostream& os)	{

		os << "total weight : " << totweight << '\n';
		os << '\n';

		vector<bpentry> bpvec;
		for (map<vector<int>,double>::iterator i=weight.begin(); i!=weight.end(); i++)	{
			bpvec.push_back(bpentry(i->first,i->second));
		}
		sort(bpvec.begin(),bpvec.end(),&compentry);

		for (unsigned int i=0; i<bpvec.size(); i++)	{
			vector<int> tmp = bpvec[i].first;
			for (int i=0; i<GetNtaxa(); i++)	{
				if (tmp[i] == 1)	{
					os << '*';
				}
				else	{
					os << '.';
				}
			}
			os << '\t' << (int) (100 * GetWeight(tmp)) << '\t' << GetLength(tmp) << '\n';
		}
	}

	// protected:

	const TaxonSet* taxset;

	map<vector<int>,double> weight;
	map<vector<int>,double> length;
	double totweight;
	double logNTrees;
};

class TpList : public BpList {

	public:

	TpList(const TaxonSet* intaxset) : BpList(intaxset) {}

	virtual double GetLength(vector<int> bp)	{
		return 0;
	}

	virtual double GetLogWeight(vector<int> bp, double cutoff)	{
		double weight = (1 - cutoff) * GetWeight(bp) + cutoff * exp(GetLogPriorTpMass(bp));
		return log(weight);
	}

	double GetLogPriorTpMass(vector<int> bp)	{

		int n1 = 0;
		int n2 = 0;
		for (int i=0; i<GetNtaxa(); i++)	{
			if (bp[i] == 1)	{
				n1++;
			}
			if (bp[i] == 2)	{
				n2++;
			}
		}
		return GetLogNPlantedTrees(n1) + GetLogNPlantedTrees(n2) + GetLogNPlantedTrees(GetNtaxa() - n1 - n2) - logNTrees;
	}

	virtual vector<int> RecursiveAdd(const Link* from, double inweight)	{

		vector<int> ret(GetNtaxa(),0);

		if (from->isLeaf())	{
			int index = GetTaxonIndex(from->GetNode()->GetName());
			if ((index < 0) || (index >= GetNtaxa()))	{
				cerr << "error in BpList::RecursiveAdd: taxon index out of bound\n";
				cerr << index << '\t' << from->GetNode()->GetName() << '\n';
				exit(1);
			}
			ret[index] = 1;
		}
		else	{
			int index = 0;
			for (const Link* link=from->Next(); link!=from; link=link->Next())	{
				index++;
				if (index == 3)	{
					if (! from->isRoot())	{
						cerr << "error in TpList: tree is not fully bifurcating\n";
						exit(1);
					}
					index = 0;
				}
				vector<int> tmp = RecursiveAdd(link->Out(), inweight);
				for (int i=0; i<GetNtaxa(); i++)	{
					if (tmp[i])	{
						ret[i] = index;
					}
				}
			}
		}

		if (! from->isLeaf())	{
			// polarize vector
			vector<int> cp = Polarize(ret);

			// and add it to bplist
			weight[cp] += inweight;
		}
		return ret;
	}

	vector<int> Polarize(vector<int> in)	{

		// first taxon should be 0
		// first taxon not 0 should be 1
		vector<int> cp = in;
		if (cp[0])	{
			for (int i=0; i<GetNtaxa(); i++)	{
				if (cp[i] != 2)	{
					cp[i] = 1-cp[i];
				}
			}
		}
		int i = 0;
		while ((i<GetNtaxa()) && (!cp[i]))	{
			i++;
		}
		if (i < GetNtaxa())	{
			if (cp[i] == 2)	{
				for (int j=i; j<GetNtaxa(); j++)	{
					if (cp[j] == 2)	{
						cp[j] = 1;
					}
					else if (cp[j] == 1)	{
						cp[j] = 2;
					}
				}
			}
		}
		return cp;
	}

};

class MultiBpList	{

	public:

	MultiBpList(const TaxonSet* intaxset)	{
		bp0 = new BpList(intaxset);
	}

	void Add(BpList* inlist)	{
		bpl.push_back(inlist);
	}

	const TaxonSet* GetTaxonSet() {return bp0->GetTaxonSet();}

	int GetNtaxa() {return GetTaxonSet()->GetNtaxa();}

	void Merge()	{

		for (unsigned int i=0; i<bpl.size(); i++)	{
			bpl[i]->Normalize();
			bp0->Merge(bpl[i]);
		}
		bp0->Normalize();

		globalmaxdiff = 0;
		globalmeandiff = 0;
		double globalweight = 0;
		for (map<vector<int>,double>::iterator i=bp0->weight.begin(); i!=bp0->weight.end(); i++)	{

			double diff = 0;
			for (unsigned int j=0; j<bpl.size(); j++)	{
				for (unsigned int k=j+1; k<bpl.size(); k++)	{
					double tmp = fabs(bpl[j]->GetWeight(i->first) - bpl[k]->GetWeight(i->first));
					if (diff < tmp)	{
						diff = tmp;
					}
				}
			}
			maxdiff[i->first] = diff;
			if (globalmaxdiff < diff)	{
				globalmaxdiff = diff;
			}
			globalmeandiff += bp0->GetWeight(i->first) * diff;
			globalweight += bp0->GetWeight(i->first);
		}

		globalmeandiff /= globalweight;
	}

	void ToStream(ostream& os)	{
		bp0->ToStream(os);
	}
	
	void ToStreamSexy(ostream& os)	{

		os << "maxdiff : " << globalmaxdiff << '\n';
		os << "meandiff : " << globalmeandiff << '\n';
		os << '\n';

		vector<bpentry> bpvec;
		for (map<vector<int>,double>::iterator i=maxdiff.begin(); i!=maxdiff.end(); i++)	{
			bpvec.push_back(bpentry(i->first,i->second));
		}
		sort(bpvec.begin(),bpvec.end(),&BpList::compentry);
		cerr << "size : " << bpvec.size() << '\n';

		for (unsigned int i=0; i<bpvec.size(); i++)	{
			vector<int> tmp = bpvec[i].first;
			for (int i=0; i<GetNtaxa(); i++)	{
				if (tmp[i])	{
					os << '*';
				}
				else	{
					os << '.';
				}
			}
			os << '\t' << (int) (100 * maxdiff[tmp]) << '\t'; //  << bp0->GetWeight(tmp) << '\t';
			for (unsigned int i=0; i<bpl.size(); i++)	{
				os << '\t' << (int) (100 * bpl[i]->GetWeight(tmp));
			}
			os << '\t';
			for (unsigned int i=0; i<bpl.size(); i++)	{
				os << '\t' << bpl[i]->GetLength(tmp);
			}
			os << '\n';
		}
	}

	protected:

	TaxonSet* taxset;

	double globalmaxdiff;
	double globalmeandiff;
	BpList* bp0;
	vector<BpList*> bpl;
	map<vector<int>,double> maxdiff;

};

class UnrootedCCP	{

	protected:

	BpList* bplist;
	TpList* tplist;
	double cutoff;
	TaxonSet* taxset;
	map<const Link*, vector<double> > probmap;
	map<pair<Link*,Link*>,double> logprobmap;

	public:

	UnrootedCCP(string filename, double incutoff)	{

		cutoff = incutoff;
		ifstream is(filename.c_str());
		ReadFromStream(is);
	}

	void ReadFromStream(istream& is)	{

		int ntaxa;
		is >> ntaxa;
		string* names = new string[ntaxa];
		for (int i=0; i<ntaxa; i++)	{
			is >> names[i];
		}
		taxset = new TaxonSet(names,ntaxa);
		
		bplist = new BpList(taxset);
		bplist->FromStream(is);
		tplist = new TpList(taxset);
		tplist->FromStream(is);
	}

	int GetNtaxa()	{
		return taxset->GetNtaxa();
	}

	const TaxonSet* GetTaxonSet() {return taxset;}

	int GetTaxonIndex(string taxon) {
		int tmp = GetTaxonSet()->GetTaxonIndex(taxon);
		if (tmp == -1)	{
			cerr << "error: taxon not found\n";
			exit(1);
		}
		return tmp;
	}

	double GetLogProb(Tree* tree)	{

		double logprob = 0;
		RecursiveGetLogProb(tree->GetRoot(),logprob);
		return logprob;
	}

	vector<int> RecursiveGetLogProb(const Link* from, double& logprob)	{

		vector<int> ret(GetNtaxa(),0);

		if (from->isLeaf())	{
			int index = GetTaxonIndex(from->GetNode()->GetName());
			if ((index < 0) || (index >= GetNtaxa()))	{
				cerr << "error in CCP::RecursiveGetLogProb: taxon index out of bound\n";
				cerr << index << '\t' << from->GetNode()->GetName() << '\n';
				exit(1);
			}
			ret[index] = 1;
		}
		else	{
			int index = 0;
			for (const Link* link=from->Next(); link!=from; link=link->Next())	{
				index++;
				if (index == 3)	{
					if (! from->isRoot())	{
						cerr << "error in CCP::RecursiveGetLogProb: tree is not fully bifurcating\n";
						exit(1);
					}
					index = 0;
				}
				vector<int> tmp = RecursiveGetLogProb(link->Out(),logprob);
				for (int i=0; i<GetNtaxa(); i++)	{
					if (tmp[i])	{
						ret[i] = index;
					}
				}
			}
		}

		if (! from->isLeaf())	{

			// tripartition
			// non root: (A,B,\AB)
			// root    : (A,B,C)
			vector<int> tp = tplist->Polarize(ret);
			logprob += tplist->GetLogWeight(tp,cutoff);

			if (! from->isRoot())	{
				// marginalize into bipartition
				// (AB,\AB)
				for (int i=0; i<GetNtaxa(); i++)	{
					if (ret[i] == 2)	{
						ret[i] = 1;
					}
				}
				vector<int> bp = bplist->Polarize(ret);
				logprob -= bplist->GetLogWeight(bp,cutoff);
			}			
		}
		return ret;
	}

	int FullGibbsSPR(Tree* tree)	{

		logprobmap.clear();

		Link* trailer = tree->GetRoot();
		for (Link* link=tree->GetRoot()->Next(); link!=tree->GetRoot(); link=link->Next())	{
			RecursiveFullGibbsSPR(tree,link->Out(),trailer);
			trailer = trailer->Next();
		}

		double max = 0;
		int found = 1;
		for (map<pair<Link*,Link*>,double>::iterator i=logprobmap.begin(); i!=logprobmap.end(); i++)	{
			if ((!found) || (max < i->second))	{
				max = i->second;
			}
		}

		double total = 0;
		for (map<pair<Link*,Link*>,double>::iterator i=logprobmap.begin(); i!=logprobmap.end(); i++)	{
			total += exp(i->second-max);
		}

		double u = total * rnd::GetRandom().Uniform();
		map<pair<Link*,Link*>,double>::iterator i=logprobmap.begin();
		double cumul = exp(i->second-max);
		while ((i!=logprobmap.end()) && (u > cumul))	{
			i++;
			if (i == logprobmap.end())	{
				cerr << "error: overflow\n";
				exit(1);
			}
			cumul += exp(i->second-max);
		}

		Link* subtree = i->first.first;
		Link* subtreeup = i->first.second;

		Link* down = 0;
		Link* up = 0;

		Link* fromdown = tree->Detach(subtree,subtreeup);
		GetSummedLogProb(tree,subtree,down,up);
		int accepted = (fromdown != down);
		tree->Attach(subtree,subtreeup,down,up);
		return accepted;
	}

	void RecursiveFullGibbsSPR(Tree* tree, Link* from, Link* fromup)	{

		Link* trailer = from;
		for (Link* link=from->Next(); link!=from; link=link->Next())	{

			Link* subtree = link->Out();
			Link* subtreeup = trailer;
			Link* fromdown = tree->Detach(subtree,subtreeup);

			Link* down = 0;
			Link* up = 0;
			double tmp = GetSummedLogProb(tree,subtree,down,up);

			logprobmap[pair<Link*,Link*>(subtree,subtreeup)] = tmp;

			tree->Attach(subtree,subtreeup,fromdown,fromup);

			RecursiveFullGibbsSPR(tree,link->Out(),trailer);

			trailer = trailer->Next();
		}
	}

	double GetSummedLogProb(Tree* tree, const Link* subtree, Link*& down, Link*& up)	{

		double phisubtree = 0;
		vector<int> bpsubtree = RecursiveGetLogProb(subtree,phisubtree);
		probmap.clear();

		double phi = 0;
		double mu = 0;

		RecursiveBackwardRegraft(tree->GetRoot(),subtree,bpsubtree,phisubtree,phi,mu);

		down = 0;
		up = 0;

		int index = 0;
		double cumul = probmap[tree->GetRoot()][index];
		double u = rnd::GetRandom().Uniform();

		Link* trailer = tree->GetRoot();
		for (Link* link=tree->GetRoot()->Next(); link!=tree->GetRoot(); link=link->Next())	{
			if ((! down) && (u < cumul))	{
				RecursiveForwardRegraft(link->Out(),trailer,down,up);
			}
			index++;
			trailer = trailer->Next();
			cumul += probmap[tree->GetRoot()][index];
		}
	}

	vector<int> RecursiveBackwardRegraft(const Link* from, const Link* subtree, vector<int> bpsubtree, double phisubtree, double& phi, double& mu)	{

		vector<int> ret(GetNtaxa(),0);

		if (from->isLeaf())	{
			int index = GetTaxonIndex(from->GetNode()->GetName());
			if ((index < 0) || (index >= GetNtaxa()))	{
				cerr << "error in CCP::RecursiveGetLogProb: taxon index out of bound\n";
				cerr << index << '\t' << from->GetNode()->GetName() << '\n';
				exit(1);
			}
			ret[index] = 1;
			for (int i=0; i<GetNtaxa(); i++)	{
				if (bpsubtree[i])	{
					if (ret[i])	{
						cerr << "error in backward\n";
						exit(1);
					}
					ret[i] = 2;
				}
			}
			vector<int> tp = tplist->Polarize(ret);
			mu = phisubtree + tplist->GetLogWeight(tp,cutoff);
			for (int i=0; i<GetNtaxa(); i++)	{
				if (bpsubtree[i])	{
					ret[i] = 1;
				}
			}
			vector<int> bp = bplist->Polarize(ret);
			mu -= bplist->GetLogWeight(bp,cutoff);
			for (int i=0; i<GetNtaxa(); i++)	{
				if (bpsubtree[i])	{
					ret[i] = 0;
				}
			}
			phi = 0;
		}
		else	{

			if (from->isRoot())	{
				int index = 0;
				double tmpmu[3];
				double tmpphi[3];
				const Link* tmplink[3];
				for (const Link* link=from->Next(); link!=from; link=link->Next())	{
					vector<int> tmp = RecursiveBackwardRegraft(link->Out(),subtree,bpsubtree,phisubtree,tmpphi[index],tmpmu[index]);
					tmplink[index] = link->Out();

					for (int i=0; i<GetNtaxa(); i++)	{
						if (tmp[i])	{
							ret[i] = index;
						}
					}
					index++;
				}

				double submu[3];
				submu[0] = tmpmu[0] + tmpphi[1] + tmpphi[2];
				submu[1] = tmpphi[0] + tmpmu[1] + tmpphi[2];
				submu[2] = tmpphi[0] + tmpphi[1] + tmpmu[2];

				double max = 0;
				for (int index=0; index<3; index++)	{
					for (int i=0; i<GetNtaxa(); i++)	{
						if (bpsubtree[i])	{
							ret[i] = index;
						}
					}
					vector<int> tp = tplist->Polarize(ret);
					submu[index] += tplist->GetLogWeight(tp,cutoff);
					if ((!index) || (max < submu[index]))	{
						max = submu[index];
					}
				}

				double total = 0;
				for (int index=0; index<3; index++)	{
					double tmp = exp(submu[index] - max);
					probmap[from].push_back(tmp);
					total += tmp;
				}
				for (int index=0; index<3; index++)	{
					probmap[from][index] /= total;
				}

				for (int i=0; i<GetNtaxa(); i++)	{
					if (bpsubtree[i])	{
						ret[i] = 0;
					}
					if (ret[i] == 2)	{
						ret[i] = 1;
					}
				}

				vector<int> bp = bplist->Polarize(ret);
				double logw = bplist->GetLogWeight(bp,cutoff);
				mu = log(total) + max - logw;

				phi = tmpphi[0] + tmpphi[1] + tmpphi[2];
	
			}
			else	{

				double tmpphi[2];
				double tmpmu[2];
				int index = 0;
				for (const Link* link=from->Next(); link!=from; link=link->Next())	{
					vector<int> tmp = RecursiveBackwardRegraft(link->Out(),subtree,bpsubtree,phisubtree,tmpphi[index],tmpmu[index]);
					index++;
					for (int i=0; i<GetNtaxa(); i++)	{
						if (tmp[i])	{
							ret[i] = index;
						}
					}
				}

				double submu[3];

				submu[0] = tmpphi[0] + tmpphi[1] + phisubtree;
				phi = tmpphi[0] + tmpphi[1];

				// (A,B,\AB)				
				vector<int> tp = tplist->Polarize(ret);
				phi += tplist->GetLogWeight(tp,cutoff);

				// (AB,S,\(ASB))
				vector<int> s_ab(GetNtaxa(),0);
				for (int i=0; i<GetNtaxa(); i++)	{
					if (ret[i])	{
						s_ab[i] = 1;
					}
					else if (bpsubtree[i])	{
						s_ab[i] = 2;
					}
				}
				vector<int> tp2 = tplist->Polarize(s_ab);
				submu[0] += tplist->GetLogWeight(tp2,cutoff);

				// (AB,\AB)
				for (int i=0; i<GetNtaxa(); i++)	{
					if (bpsubtree[i])	{
						s_ab[i] = 0;
					}
				}
				vector<int> bp = bplist->Polarize(s_ab);
				double tmp = bplist->GetLogWeight(bp,cutoff);
				phi -= tmp;
				submu[0] -= tmp;
			
				submu[1] = tmpmu[0] + tmpphi[1];
				submu[2] = tmpphi[0] + tmpmu[1];

				double max = 0;
				for (int index=0; index<3; index++)	{

					// 0: (A,B,\AB)
					// 1: (AS,B,\ASB)
					// 2: (A,BS,\ASB)
					for (int i=0; i<GetNtaxa(); i++)	{
						if (bpsubtree[i])	{
							ret[i] = index;
						}
					}
					vector<int> tp = tplist->Polarize(ret);
					submu[index] += tplist->GetLogWeight(tp,cutoff);

					if ((!index) || (max < submu[index]))	{
						max = submu[index];
					}
				}
				
				double total = 0;
				for (int index=0; index<3; index++)	{
					double tmp = exp(submu[index] - max);
					probmap[from].push_back(tmp);
					total += tmp;
				}
				for (int index=0; index<3; index++)	{
					probmap[from][index] /= total;
				}

				// common to all three cases
				// (ASB,\ASB)
				for (int i=0; i<GetNtaxa(); i++)	{
					// already imp,lied by last interation above
					/*
					if (bpsubtree[i])	{
						ret[i] = 1;
					}
					*/
					if (ret[i] == 2)	{
						ret[i] = 1;
					}
				}
				vector<int> bp2 = bplist->Polarize(ret);
				double logw = bplist->GetLogWeight(bp2,cutoff);
				mu = log(total) + max - logw;

				// restore bipartition to: (AB,\AB)
				for (int i=0; i<GetNtaxa(); i++)	{
					if (bpsubtree[i])	{
						ret[i] = 0;
					}
				}
			}
		}

		return ret;
	}

	void RecursiveForwardRegraft(Link* from, Link* fromup, Link*& down, Link*& up)	{

		// choose
		if (from->isLeaf())	{
			down = from;
			up = fromup;
		}
		else	{
			double u = rnd::GetRandom().Uniform();
			if (u < probmap[from][0])	{
				down = from;
				up = fromup;
			}
			else if (u < (probmap[from][0] + probmap[from][1]))	{
				RecursiveForwardRegraft(from->Next()->Out(),from,down,up);
			}
			else	{
				RecursiveForwardRegraft(from->Next()->Next()->Out(),from->Next(),down,up);
			}
		}
	}
};

class UnrootedBPP	{

	protected:

	BpList* bplist;
	double cutoff;
	double beta;
	TaxonSet* taxset;
	map<const Link*, vector<double> > probmap;
	map<pair<Link*,Link*>,double> logprobmap;
	map<pair<Link*,Link*>,double> logprobmap2;

	public:

	UnrootedBPP(string filename, double incutoff, double inbeta = 1)	{

		cutoff = incutoff;
		beta = inbeta;
		ifstream is(filename.c_str());
		ReadFromStream(is);
	}

	void ReadFromStream(istream& is)	{

		int ntaxa;
		is >> ntaxa;
		string* names = new string[ntaxa];
		for (int i=0; i<ntaxa; i++)	{
			is >> names[i];
		}
		taxset = new TaxonSet(names,ntaxa);
		
		bplist = new BpList(taxset);
		bplist->FromStream(is);
	}

	int GetNtaxa()	{
		return taxset->GetNtaxa();
	}

	const TaxonSet* GetTaxonSet() {return taxset;}

	int GetTaxonIndex(string taxon) {
		int tmp = GetTaxonSet()->GetTaxonIndex(taxon);
		if (tmp == -1)	{
			cerr << "error: taxon not found\n";
			exit(1);
		}
		return tmp;
	}

	double GetLogProb(Tree* tree);

	vector<int> RecursiveGetLogProb(const Link* from, double& logprob);

	double ProposeFullGibbsSPR(Tree* tree, Link*& subtree, Link*& subtreeup, Link*& fromdown, Link*& fromup, Link*& down, Link*& up, double& finallogp);

	int FullGibbsSPR(Tree* tree);

	void RecursiveFullGibbsSPR(Tree* tree, Link* from, Link* fromup, double reflogp);

	double GetSummedLogProb(Tree* tree, const Link* subtree, Link*& down, Link*& up);
	double GetSummedLogProb2(Tree* tree, Link* subtree, Link* subtreeup, Link* fromdown, Link*& down, Link*& up);
	void RecursiveGetSummedLogProb2(Tree* tree, Link* subtree, Link* subtreeup, Link* fromdown, Link* from, Link* up);

	vector<int> RecursiveBackwardRegraft(const Link* from, const Link* subtree, vector<int> bpsubtree, double phisubtree, double& phi, double& mu);

	void RecursiveForwardRegraft(Link* from, Link* fromup, Link*& down, Link*& up);
};

#endif
