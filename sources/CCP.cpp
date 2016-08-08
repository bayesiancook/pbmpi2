#include "CCP.h"

void BpList::RegisterWithTaxonSet(const TaxonSet* intaxset)	{

	if (! taxset)	{
		taxset = intaxset;
	}
	else	{

		if (intaxset->GetNtaxa() != taxset->GetNtaxa())	{
			cerr << "error: taxon sets do not match\n";
			exit(1);
		}
		int* permut = new int[GetNtaxa()];
		for (int i=0; i<GetNtaxa(); i++)	{
			permut[i] = taxset->GetTaxonIndex(intaxset->GetTaxon(i));
		}

		map<vector<int>,double> weight2;
		map<vector<int>,double> length2;
		for (map<vector<int>,double>::iterator i=weight.begin(); i!=weight.end(); i++)	{
			vector<int> from = i->first;
			vector<int> to(GetNtaxa(),0);
			for (int i=0; i<GetNtaxa(); i++)	{
				to[i] = from[permut[i]];
			}
			vector<int> bp = Polarize(to);
			weight2[bp] = weight[from];
			length2[bp] = length[from];
		}

		weight = weight2;
		length = length2;
		taxset = intaxset;
		delete[] permut;
	}
}




double UnrootedBPP::GetLogProb(Tree* tree)	{

	double logprob = 0;
	RecursiveGetLogProb(tree->GetRoot(),logprob);
	return logprob;
}


double UnrootedCCP::GetLogProb(Tree* tree)	{

	double logprob = 0;
	if (tpactivated)	{
		RecursiveGetLogProb(tree->GetRoot(),logprob);
	}
	else	{
		UnrootedBPP::RecursiveGetLogProb(tree->GetRoot(),logprob);
	}
	return logprob;
}

vector<int> UnrootedBPP::RecursiveGetLogProb(const Link* from, double& logprob)	{

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
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			vector<int> tmp = RecursiveGetLogProb(link->Out(),logprob);
			for (int i=0; i<GetNtaxa(); i++)	{
				if (tmp[i])	{
					ret[i] = 1;
				}
			}
		}
	}

	if ((! from->isLeaf()) && (! from->isRoot()))	{
		vector<int> bp = bplist->Polarize(ret);
		logprob += beta * bplist->GetLogWeight(bp,cutoff);
	}
	return ret;
}

int UnrootedBPP::FullGibbsSPR(Tree* tree)	{

	Link* subtree = 0;
	Link* subtreeup = 0;
	Link* fromdown = 0;
	Link* fromup = 0;
	Link* down = 0;
	Link* up = 0;

	double logp1 = GetLogProb(tree);
	double logh = ProposeFullGibbsSPR(tree,subtree,subtreeup,fromdown,fromup,down,up);

	tree->Detach(subtree,subtreeup);

	tree->Attach(subtree,subtreeup,down,up);

	double logp2 = GetLogProb(tree);

	double logratio = logp2 - logp1 + logh;

	int accepted = (log(rnd::GetRandom().Uniform()) < logratio);
	if (! accepted)	{

		tree->Detach(subtree,subtreeup);
		tree->Attach(subtree,subtreeup,fromdown,fromup);
	}
	return accepted;
}

double UnrootedBPP::ProposeFullGibbsSPR(Tree* tree, Link*& subtree, Link*& subtreeup, Link*& fromdown, Link*& fromup, Link*& down, Link*& up)	{

	double logp1 = GetLogProb(tree);

	logprobmap.clear();

	for (Link* link=tree->GetRoot()->Next(); link!=tree->GetRoot(); link=link->Next())	{
		RecursiveFullGibbsSPR(tree,link->Out(),logp1);
	}

	if (logprobmap.size() != (2*GetNtaxa() - 6))	{
		cerr << "error in map size : " << logprobmap.size() << '\t' << 2*GetNtaxa() - 6 << '\n';
		exit(1);
	}

	double max = 0;
	int found = 1;
	for (map<pair<Link*,Link*>,double>::iterator i=logprobmap.begin(); i!=logprobmap.end(); i++)	{
		if ((!found) || (max < i->second))	{
			max = i->second;
			found = 0;
		}
	}

	double total = 0;
	for (map<pair<Link*,Link*>,double>::iterator i=logprobmap.begin(); i!=logprobmap.end(); i++)	{
		total += exp(i->second-max);
	}
	double logZ1 = log(total) + max;


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

	subtree = i->first.first;
	subtreeup = i->first.second;

	/*
	cerr << '(';
	tree->ToStream(cerr,subtree);
	cerr << ");\n";
	cerr << '\n';
	*/
	
	down = 0;
	up = 0;

	fromdown = tree->Detach(subtree,subtreeup);
	fromup = tree->GetAncestor(fromdown);
	double tmp = GetSummedLogProb2(tree,subtree,subtreeup,fromdown,down,up);
	if (fromdown == down)	{
		cerr << "error: proposed tree is identical to initial tree\n";
		exit(1);
	}
	if (fabs(tmp - i->second) > 1e-3)	{
		cerr << "error in nested choice: non matching prob\n";
		cerr << tmp << '\t' << i->second << '\n';
		exit(1);
	}
	tree->Attach(subtree,subtreeup,down,up);
	double logp2 = GetLogProb(tree);

	double logpfwd = logp2 - logZ1;

	// reverse move

	logprobmap.clear();

	for (Link* link=tree->GetRoot()->Next(); link!=tree->GetRoot(); link=link->Next())	{
		RecursiveFullGibbsSPR(tree,link->Out(),logp2);
	}

	if (logprobmap.size() != (2*GetNtaxa() - 6))	{
		cerr << "error in map size : " << logprobmap.size() << '\t' << 2*GetNtaxa() - 6 << '\n';
		exit(1);
	}
	max = 0;
	found = 1;
	for (map<pair<Link*,Link*>,double>::iterator i=logprobmap.begin(); i!=logprobmap.end(); i++)	{
		if ((!found) || (max < i->second))	{
			max = i->second;
			found = 0;
		}
	}

	total = 0;
	for (map<pair<Link*,Link*>,double>::iterator i=logprobmap.begin(); i!=logprobmap.end(); i++)	{
		total += exp(i->second-max);
	}
	double logZ2 = log(total) + max;

	double logprev = logp1 - logZ2;

	double logh = logprev - logpfwd;

	// restore tree
	tree->Detach(subtree,subtreeup);
	tree->Attach(subtree,subtreeup,fromdown,fromup);

	double logp3 = GetLogProb(tree);
	if (fabs(logp3 - logp1) > 1e-6)	{
		cerr << "error when restoring tree\n";
		cerr << logp1 << '\t' << logp3 << '\n';
		exit(1);
	}

	return logh;
}

void UnrootedBPP::RecursiveFullGibbsSPR(Tree* tree, Link* from, double reflogp)	{

	for (Link* link=from->Next(); link!=from; link=link->Next())	{

		Link* subtree = link->Out();
		Link* subtreeup = from;

		double tmp = GetSummedLogProb(tree,subtree);

		double p = 1 - exp(reflogp - tmp);
		if (p < -1e-9)	{
			cerr << "error in RecursiveFullGibbsSPR: negative prob\n";
			cerr << "reflogp : " << reflogp << '\n';
			cerr << "check   : " << GetLogProb(tree) << '\n';
			cerr << "totlogp : " << tmp << '\n';
			cerr << "diff : " << reflogp - tmp << '\n';
			cerr << tree->GetSize() << '\t' << tree->GetSize(subtree) << '\n';
			tree->ToStream(cerr);
			double phi = 0;
			vector<int> bpsubtree = RecursiveGetLogProb(subtree,phi);
			vector<int> bp = bplist->Polarize(bpsubtree);
			int size = 0;
			for (int i=0; i<GetNtaxa(); i++)	{
				if (bpsubtree[i] == 1)	{
					size++;
				}
				cerr << bpsubtree[i];
			}
			cerr << '\n';
			cerr << size << '\n';
			cerr << phi << '\n';
			exit(1);
		}
		if (p < 0)	{
			p = 0;
		}
		tmp += log(p);

		logprobmap[pair<Link*,Link*>(subtree,subtreeup)] = tmp;

		RecursiveFullGibbsSPR(tree,link->Out(),reflogp);
	}
}

double UnrootedBPP::ProposeFragileSPR(Tree* tree, Link*& subtree, Link*& subtreeup)	{

	proposemap.clear();

	for (Link* link=tree->GetRoot()->Next(); link!=tree->GetRoot(); link=link->Next())	{
		RecursiveProposeFragileSPR(tree,link->Out());
	}

	if (proposemap.size() != (2*GetNtaxa() - 6))	{
		cerr << "error in map size : " << proposemap.size() << '\t' << 2*GetNtaxa() - 6 << '\n';
		cerr << "propose fragile\n";
		exit(1);
	}

	double total = 0;
	for (map<pair<Link*,Link*>,double>::iterator i=proposemap.begin(); i!=proposemap.end(); i++)	{
		total += i->second;
	}

	double u = total * rnd::GetRandom().Uniform();
	map<pair<Link*,Link*>,double>::iterator i=proposemap.begin();
	double cumul = i->second;
	while ((i!=proposemap.end()) && (u > cumul))	{
		i++;
		if (i == proposemap.end())	{
			cerr << "error: overflow\n";
			exit(1);
		}
		cumul += i->second;
	}

	subtree = i->first.first;
	subtreeup = i->first.second;

	double ret = i->second / total;

	cerr << '(';
	tree->ToStream(cerr,subtree);
	cerr << ");\n";
	cerr << ret << '\t' << total << '\t' << 1.0 / total << '\t' << proposemap.size() << '\n';
	cerr << '\n';
	
	return ret;
}

double UnrootedBPP::GetFragileSPRProposalProb(Tree* tree, Link* subtree, Link* subtreeup)	{

	proposemap.clear();

	for (Link* link=tree->GetRoot()->Next(); link!=tree->GetRoot(); link=link->Next())	{
		RecursiveProposeFragileSPR(tree,link->Out());
	}

	if (proposemap.size() != (2*GetNtaxa() - 6))	{
		cerr << "error in map size : " << proposemap.size() << '\t' << 2*GetNtaxa() - 6 << '\n';
		cerr << "get proposal prob\n";
		exit(1);
	}

	double total = 0;
	for (map<pair<Link*,Link*>,double>::iterator i=proposemap.begin(); i!=proposemap.end(); i++)	{
		total += i->second;
	}

	double ret = proposemap[pair<Link*,Link*>(subtree,subtreeup)] / total;

	return ret;
}

vector<int> UnrootedBPP::RecursiveProposeFragileSPR(Tree* tree, Link* from)	{

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

		int index = 1;
		for (Link* link=from->Next(); link!=from; link=link->Next())	{

			vector<int> bp = RecursiveProposeFragileSPR(tree,link->Out());

			for (int i=0; i<GetNtaxa(); i++)	{
				if (bp[i])	{
					ret[i] = index;
				}
			}
			index++;
		}

		double tmp = GetFragilityScore(ret);

		for (Link* link=from->Next(); link!=from; link=link->Next())	{
			proposemap[pair<Link*,Link*>(link->Out(),from)] = tmp;
		}

	}
	return ret;
}

double UnrootedBPP::GetFragilityScore(vector<int> tp)	{

	cerr << "in bpp fragility score\n";
	exit(1);
	return 0;
}

double UnrootedCCP::GetFragilityScore(vector<int> tp)	{

	vector<int> poltp = tplist->Polarize(tp);
	double ret = 0;
	double w = tplist->GetWeight(poltp);
	if (w > cutoff)	{
		ret = beta;
	}
	else	{
		ret = 1;
	}
	return ret;
}

double UnrootedBPP::GetSummedLogProb(Tree* tree, const Link* subtree)	{

	double phisubtree = 0;
	vector<int> bpsubtree = RecursiveGetLogProb(subtree,phisubtree);
	// probmap.clear();

	double phi = 0;
	double mu = 0;

	RecursiveBackwardRegraft(tree->GetRoot(),subtree,bpsubtree,phisubtree,phi,mu);

	return mu;
}

double UnrootedBPP::GetSummedLogProb2(Tree* tree, Link* subtree, Link* subtreeup, Link* fromdown, Link*& down, Link*& up)	{

	logprobmap2.clear();

	RecursiveGetSummedLogProb2(tree,subtree,subtreeup,fromdown,tree->GetRoot(),tree->GetRoot());

	double max = 0;
	int notfound = 1;
	for (map<pair<Link*,Link*>,double>::iterator i=logprobmap2.begin(); i!=logprobmap2.end(); i++)	{
		if (notfound || (max < i->second))	{
			max = i->second;
			notfound = 0;
		}
	}

	double total = 0;
	for (map<pair<Link*,Link*>,double>::iterator i=logprobmap2.begin(); i!=logprobmap2.end(); i++)	{
		total += exp(i->second - max);
	}

	double u = total * rnd::GetRandom().Uniform();
	map<pair<Link*,Link*>,double>::iterator i=logprobmap2.begin();
	double cumul = exp(i->second - max);
	while ((i!=logprobmap2.end()) && (u > cumul))	{
		i++;
		if (i == logprobmap2.end())	{
			cerr << "error in get summed log prob2: oveeflow\n";
			exit(1);
		}
		cumul += exp(i->second - max);
	}
	
	down = i->first.first;
	up = i->first.second;

	return log(total) + max;
}

void UnrootedBPP::RecursiveGetSummedLogProb2(Tree* tree, Link* subtree, Link* subtreeup, Link* fromdown, Link* from, Link* up)	{

	if ((! from->isRoot()) && (from != fromdown))	{

		tree->Attach(subtree,subtreeup,from,up);
		double tmp = GetLogProb(tree);
		logprobmap2[pair<Link*,Link*>(from,up)] = tmp;
		tree->Detach(subtree,subtreeup);
	}

	for (Link* link = from->Next(); link!=from; link=link->Next())	{
		RecursiveGetSummedLogProb2(tree,subtree,subtreeup,fromdown,link->Out(),from);
	}
}

vector<int> UnrootedBPP::RecursiveBackwardRegraft(const Link* from, const Link* subtree, vector<int> bpsubtree, double phisubtree, double& phi, double& mu)	{

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
				if (ret[i] != 0)	{
					cerr << "error in backward when adding subtree: leaf\n";
					exit(1);
				}
				ret[i] = 1;
			}
		}
		vector<int> bp = bplist->Polarize(ret);
		mu = phisubtree + beta * bplist->GetLogWeight(bp,cutoff);
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
			for (const Link* link=from->Next(); link!=from; link=link->Next())	{
				vector<int> tmp = RecursiveBackwardRegraft(link->Out(),subtree,bpsubtree,phisubtree,tmpphi[index],tmpmu[index]);

				for (int i=0; i<GetNtaxa(); i++)	{
					if (tmp[i])	{
						if (ret[i] != 0)	{
							cerr << "error in backward when combining subclades: root\n";
							exit(1);
						}
						ret[i] = 1;
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
				if ((!index) || (max < submu[index]))	{
					max = submu[index];
				}
			}

			double total = 0;
			for (int index=0; index<3; index++)	{
				double tmp = exp(submu[index] - max);
				total += tmp;
			}

			mu = log(total) + max;
			phi = tmpphi[0] + tmpphi[1] + tmpphi[2];

		}
		else	{

			int sub = 0;
			double tmpphi[2];
			double tmpmu[2];
			int index = 0;
			for (const Link* link=from->Next(); link!=from; link=link->Next())	{
				if (link->Out() == subtree)	{
					sub = 1;
				}
				else	{
					vector<int> tmp = RecursiveBackwardRegraft(link->Out(),subtree,bpsubtree,phisubtree,tmpphi[index],tmpmu[index]);
					index++;
					for (int i=0; i<GetNtaxa(); i++)	{
						if (tmp[i])	{
							if (ret[i])	{
								cerr << "error in backward when combining subclades: internal\n";
								exit(1);
							}
							ret[i] = 1;
						}
					}
				}
			}

			if (sub)	{
				phi = tmpphi[0];
				mu = tmpmu[0];
			}
			else	{
				double submu[3];

				submu[0] = tmpphi[0] + tmpphi[1] + phisubtree;
				submu[1] = tmpmu[0] + tmpphi[1];
				submu[2] = tmpphi[0] + tmpmu[1];

				phi = tmpphi[0] + tmpphi[1];

				// (AB,\AB)				
				vector<int> bp = bplist->Polarize(ret);
				double tmp = beta * bplist->GetLogWeight(bp,cutoff);
				phi += tmp;
				submu[0] += tmp;

				double max = 0;
				for (int index=0; index<3; index++)	{
					if ((!index) || (max < submu[index]))	{
						max = submu[index];
					}
				}
				
				double total = 0;
				for (int index=0; index<3; index++)	{
					double tmp = exp(submu[index] - max);
					total += tmp;
				}

				// (ASB,\(ASB))
				for (int i=0; i<GetNtaxa(); i++)	{
					if (bpsubtree[i])	{
						if (ret[i])	{
							cerr << "error in backward when adding subtree: internal\n";
							exit(1);
						}
						ret[i] = 1;
					}
				}
				vector<int> bp2 = bplist->Polarize(ret);
				double logw = beta * bplist->GetLogWeight(bp2,cutoff);

				mu = log(total) + max + logw;

				// restore bipartition to: (AB,\AB)
				for (int i=0; i<GetNtaxa(); i++)	{
					if (bpsubtree[i])	{
						ret[i] = 0;
					}
				}
			}
		}
	}

	return ret;
}

vector<int> UnrootedCCP::RecursiveGetLogProb(const Link* from, double& logprob)	{

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
		// non root: (\AB.A,B)
		// root    : (C,A,B)
		vector<int> tp = tplist->Polarize(ret);
		logprob += beta * tplist->GetLogWeight(tp,cutoff);

		if (! from->isRoot())	{
			// marginalize into bipartition
			// (AB,\AB)
			for (int i=0; i<GetNtaxa(); i++)	{
				if (ret[i] == 2)	{
					ret[i] = 1;
				}
			}
			vector<int> bp = bplist->Polarize(ret);
			logprob -= beta * bplist->GetLogWeight(bp,cutoff);
		}			
	}
	return ret;
}

vector<int> UnrootedCCP::RecursiveBackwardRegraft(const Link* from, const Link* subtree, vector<int> bpsubtree, double phisubtree, double& phi, double& mu)	{

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
		mu = phisubtree + beta * tplist->GetLogWeight(tp,cutoff);
		for (int i=0; i<GetNtaxa(); i++)	{
			if (bpsubtree[i])	{
				ret[i] = 1;
			}
		}
		vector<int> bp = bplist->Polarize(ret);
		mu -= beta * bplist->GetLogWeight(bp,cutoff);
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
			for (const Link* link=from->Next(); link!=from; link=link->Next())	{
				vector<int> tmp = RecursiveBackwardRegraft(link->Out(),subtree,bpsubtree,phisubtree,tmpphi[index],tmpmu[index]);

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

			for (int index=0; index<3; index++)	{
				for (int i=0; i<GetNtaxa(); i++)	{
					if (bpsubtree[i])	{
						ret[i] = index;
					}
				}
				vector<int> tp = tplist->Polarize(ret);
				submu[index] += beta * tplist->GetLogWeight(tp,cutoff);
			}

			double max = 0;
			for (int index=0; index<3; index++)	{
				if ((!index) || (max < submu[index]))	{
					max = submu[index];
				}
			}

			double total = 0;
			for (int index=0; index<3; index++)	{
				double tmp = exp(submu[index] - max);
				total += tmp;
			}

			mu = log(total) + max;

			// phi does not really have any meaning at the root
			// but not used afterwards anyway
			phi = tmpphi[0] + tmpphi[1] + tmpphi[2];

			// restoring bipartition, but not used by the calling function
			for (int i=0; i<GetNtaxa(); i++)	{
				if (bpsubtree[i])	{
					ret[i] = 0;
				}
				else	{
					ret[i] = 1;
				}
			}
		}
		else	{

			int sub = 0;
			double tmpphi[2];
			double tmpmu[2];
			int index = 0;
			for (const Link* link=from->Next(); link!=from; link=link->Next())	{
				if (link->Out() == subtree)	{
					sub = 1;
				}
				else	{
					vector<int> tmp = RecursiveBackwardRegraft(link->Out(),subtree,bpsubtree,phisubtree,tmpphi[index],tmpmu[index]);
					index++;
					for (int i=0; i<GetNtaxa(); i++)	{
						if (tmp[i])	{
							ret[i] = index;
						}
					}
				}
			}

			if (sub)	{
				if (index != 1)	{
					cerr << "error in backward: incorrect number of subtrees\n";
					exit(1);
				}
				phi = tmpphi[0];
				mu = tmpmu[0];
			}
			else	{

				double submu[3];

				submu[0] = tmpphi[0] + tmpphi[1] + phisubtree;
				phi = tmpphi[0] + tmpphi[1];

				// (A,B,\AB)				
				vector<int> tp = tplist->Polarize(ret);
				phi += beta * tplist->GetLogWeight(tp,cutoff);

				// (AB,S,\(ASB))
				vector<int> s_ab(GetNtaxa(),0);
				for (int i=0; i<GetNtaxa(); i++)	{
					if (bpsubtree[i])	{
						if (ret[i])	{
							cerr << "error: ret and bpsubtree both non null\n";
							exit(1);
						}
						s_ab[i] = 2;
					}
					if (ret[i])	{
						s_ab[i] = 1;
					}
				}
				vector<int> tp2 = tplist->Polarize(s_ab);
				submu[0] += beta * tplist->GetLogWeight(tp2,cutoff);

				// (AB,\AB)
				for (int i=0; i<GetNtaxa(); i++)	{
					if (bpsubtree[i])	{
						s_ab[i] = 0;
					}
				}
				vector<int> bp = bplist->Polarize(s_ab);
				double tmp = beta * bplist->GetLogWeight(bp,cutoff);
				phi -= tmp;
				submu[0] -= tmp;
			
				submu[1] = tmpmu[0] + tmpphi[1];
				submu[2] = tmpphi[0] + tmpmu[1];

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
					double tmp = beta * tplist->GetLogWeight(tp,cutoff);
					submu[index] += tmp;
				}

				double max = 0;
				for (int index=0; index<3; index++)	{
					if ((!index) || (max < submu[index]))	{
						max = submu[index];
					}
				}
				
				double total = 0;
				for (int index=0; index<3; index++)	{
					double tmp = exp(submu[index] - max);
					total += tmp;
				}

				// common to all three cases
				// (ASB,\ASB)
				for (int i=0; i<GetNtaxa(); i++)	{
					if (bpsubtree[i])	{
						ret[i] = 1;
					}
					if (ret[i] == 2)	{
						ret[i] = 1;
					}
				}
				vector<int> bp2 = bplist->Polarize(ret);
				double logw = beta * bplist->GetLogWeight(bp2,cutoff);
				mu = log(total) + max - logw;

				// restore bipartition to: (AB,\AB)
				for (int i=0; i<GetNtaxa(); i++)	{
					if (bpsubtree[i])	{
						ret[i] = 0;
					}
				}
			}
		}
	}

	return ret;
}

/*
double UnrootedBPP::GetSummedLogProb(Tree* tree, const Link* subtree, Link*& down, Link*& up)	{

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

	return mu;
}

vector<int> UnrootedBPP::RecursiveBackwardRegraft(const Link* from, const Link* subtree, vector<int> bpsubtree, double phisubtree, double& phi, double& mu)	{

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
				if (ret[i] != 0)	{
					cerr << "error in backward when adding subtree: leaf\n";
					exit(1);
				}
				ret[i] = 1;
			}
		}
		vector<int> bp = bplist->Polarize(ret);
		mu = phisubtree + beta * bplist->GetLogWeight(bp,cutoff);
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
			for (const Link* link=from->Next(); link!=from; link=link->Next())	{
				vector<int> tmp = RecursiveBackwardRegraft(link->Out(),subtree,bpsubtree,phisubtree,tmpphi[index],tmpmu[index]);

				for (int i=0; i<GetNtaxa(); i++)	{
					if (tmp[i])	{
						if (ret[i] != 0)	{
							cerr << "error in backward when combining subclades: root\n";
							exit(1);
						}
						ret[i] = 1;
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

			mu = log(total) + max;
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
						if (ret[i])	{
							cerr << "error in backward when combining subclades: internal\n";
							exit(1);
						}
						ret[i] = 1;
					}
				}
			}

			double submu[3];

			submu[0] = tmpphi[0] + tmpphi[1] + phisubtree;
			submu[1] = tmpmu[0] + tmpphi[1];
			submu[2] = tmpphi[0] + tmpmu[1];

			phi = tmpphi[0] + tmpphi[1];

			// (AB,\AB)				
			vector<int> bp = bplist->Polarize(ret);
			double tmp = beta * bplist->GetLogWeight(bp,cutoff);
			phi += tmp;
			submu[0] += tmp;

			double max = 0;
			for (int index=0; index<3; index++)	{
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

			// (ASB,\(ASB))
			for (int i=0; i<GetNtaxa(); i++)	{
				if (bpsubtree[i])	{
					if (ret[i])	{
						cerr << "error in backward when adding subtree: internal\n";
						exit(1);
					}
					ret[i] = 1;
				}
			}
			vector<int> bp2 = bplist->Polarize(ret);
			double logw = beta * bplist->GetLogWeight(bp2,cutoff);

			mu = log(total) + max + logw;

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

void UnrootedBPP::RecursiveForwardRegraft(Link* from, Link* fromup, Link*& down, Link*& up)	{

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
*/

