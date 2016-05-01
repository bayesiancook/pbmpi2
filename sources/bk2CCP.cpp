#include "CCP.h"

double UnrootedBPP::GetLogProb(Tree* tree)	{

	double logprob = 0;
	RecursiveGetLogProb(tree->GetRoot(),logprob);
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
		logprob += bplist->GetLogWeight(bp,cutoff);
	}
	return ret;
}

double UnrootedBPP::ProposeFullGibbsSPR(Tree* tree, Link*& subtree, Link*& subtreeup, Link*& fromdown, Link*& fromup, Link*& down, Link*& up, double& finallogp)	{

	/*
	ostringstream s1;
	tree->ToStream(s1);
	*/

	double logp1 = GetLogProb(tree);

	logprobmap.clear();

	Link* trailer = tree->GetRoot();
	for (Link* link=tree->GetRoot()->Next(); link!=tree->GetRoot(); link=link->Next())	{
		RecursiveFullGibbsSPR(tree,link->Out(),trailer,logp1);
		trailer = trailer->Next();
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

	down = 0;
	up = 0;

	fromdown = 0;
	fromup = 0;

	tree->Detach2(subtree,subtreeup,fromdown,fromup);
	double tmp = GetSummedLogProb2(tree,subtree,subtreeup,fromdown,down,up);
	// double tmp = GetSummedLogProb(tree,subtree,down,up);
	if (fromdown == down)	{
		cerr << "error: proposed tree is identical to initial tree\n";
		exit(1);
	}
	// cerr << fromdown << '\t' << down << '\n';
	/*
	double tmp2 = GetSummedLogProb2(tree,subtree,subtreeup,down,up);
	if (fabs(tmp - tmp2)> 1e-6)	{
		cerr << "error: summed log probs don't match\n";
		exit(1);
	}
	*/
	/*
	double p = 1 - exp(logp1 - tmp);
	if (p <= 0)	{
		cerr << "error in RecursiveFullGibbsSPR: negative prob\n";
	}
	tmp -= log(p);
	*/
	if (fabs(tmp - i->second) > 1e-3)	{
		cerr << "error in nested choice: non matching prob\n";
		cerr << tmp << '\t' << i->second << '\n';
		exit(1);
	}
	tree->Attach(subtree,subtreeup,down,up);
	double logp2 = GetLogProb(tree);

	finallogp = logp2;

	double logpfwd = logp2 - logZ1;

	// reverse move

	logprobmap.clear();

	trailer = tree->GetRoot();
	for (Link* link=tree->GetRoot()->Next(); link!=tree->GetRoot(); link=link->Next())	{
		RecursiveFullGibbsSPR(tree,link->Out(),trailer,logp2);
		trailer = trailer->Next();
	}

	unsigned int size = logprobmap.size();
	if (size != (2*GetNtaxa() - 6))	{
		cerr << "error in map size : " << size << '\t' << 2*GetNtaxa() - 6 << '\n';
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

	/*
	ostringstream s2;
	tree->ToStream(s2);
	if (s1.str() != s2.str())	{
		cerr << "error: non matching trees\n";
		cerr << s1.str() << '\n';
		cerr << s2.str() << '\n';
		exit(1);
	}
	*/
	double logp3 = GetLogProb(tree);
	if (fabs(logp3 - logp1) > 1e-6)	{
		cerr << "error when restoring tree\n";
		cerr << logp1 << '\t' << logp3 << '\n';
		exit(1);
	}

	return logh;
}

int UnrootedBPP::FullGibbsSPR(Tree* tree)	{

	Link* subtree = 0;
	Link* subtreeup = 0;
	Link* fromdown = 0;
	Link* fromup = 0;
	Link* down = 0;
	Link* up = 0;

	double logp1 = GetLogProb(tree);
	double finallogp = 0;
	double logh = ProposeFullGibbsSPR(tree,subtree,subtreeup,fromdown,fromup,down,up,finallogp);

	Link* fromdown2 = 0;
	Link* fromup2 = 0;
	tree->Detach2(subtree,subtreeup,fromdown2,fromup2);
	/*
	if ((fromup != fromup2) || (fromdown != fromdown2))	{
		cerr << "error in full gibbs spr: non matching froms\n";
		exit(1);
	}
	*/

	tree->Attach(subtree,subtreeup,down,up);

	double logp2 = GetLogProb(tree);
	if (fabs(finallogp - logp2) > 1e-6)	{
		cerr << "error in gibbs: probs of final tree do not match\n";
		exit(1);
	}

	double logratio = logp2 - logp1 + logh;
	// cerr << logratio << '\t' << logh << '\t' << logp1 << '\t' << logp2 << '\n';
	
	int accepted = (log(rnd::GetRandom().Uniform()) < logratio);
	if (! accepted)	{

		Link* down2 = 0;
		Link* up2 = 0;
		tree->Detach2(subtree,subtreeup,down2,up2);
		/*
		if ((up != up2) || (down != down2))	{
			cerr << "error in full gibbs spr: non matching to's\n";
			exit(1);
		}
		*/
		tree->Attach(subtree,subtreeup,fromdown,fromup);
	}
	return accepted;
}

void UnrootedBPP::RecursiveFullGibbsSPR(Tree* tree, Link* from, Link* fromup, double reflogp)	{

	Link* trailer = from;
	for (Link* link=from->Next(); link!=from; link=link->Next())	{

		double logpbefore = GetLogProb(tree);
		Link* subtree = link->Out();
		Link* subtreeup = trailer;
		/*
		Link* fromdown = 0;
		tree->Detach2(subtree,subtreeup,fromdown,fromup);
		*/
		Link* fromdown = tree->Detach(subtree,subtreeup);
		/*
		if (fromup->Next()->Out() != fromdown)	{
			cerr << "eror in recursive full gibbs spr\n";
			exit(1);
		}
		*/

		Link* down = 0;
		Link* up = 0;
		// double tmp = GetSummedLogProb(tree,subtree,down,up);
		double tmp = GetSummedLogProb2(tree,subtree,subtreeup,fromdown,down,up);
		/*
		if (fabs(tmp - tmp2)> 1e-6)	{
			cerr << "error: summed log probs don't match\n";
			exit(1);
		}
		*/
		/*
		double p = 1 - exp(reflogp - tmp);
		if (p <= 0)	{
			cerr << "error in RecursiveFullGibbsSPR: negative prob\n";
		}
		tmp -= log(p);
		*/

		if (logprobmap[pair<Link*,Link*>(subtree,subtreeup)])	{
			cerr << "already done\n";
			exit(1);
		}
		logprobmap[pair<Link*,Link*>(subtree,subtreeup)] = tmp;

		tree->Attach(subtree,subtreeup,fromdown,fromup);

		double logpafter = GetLogProb(tree);

		if (fabs(logpbefore - logpafter) > 1e-6)	{
			cerr << "error in recursive gibbs: tree not correctly restored\n";
			exit(1);
		}

		trailer = trailer->Next();
	}
	trailer = from;
	for (Link* link=from->Next(); link!=from; link=link->Next())	{
		RecursiveFullGibbsSPR(tree,link->Out(),trailer,reflogp);

		trailer = trailer->Next();
	}

}

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

	Link* trailer = from;
	for (Link* link = from->Next(); link!=from; link=link->Next())	{
		RecursiveGetSummedLogProb2(tree,subtree,subtreeup,fromdown,link->Out(),trailer);
		trailer = trailer->Next();
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
		mu = phisubtree + bplist->GetLogWeight(bp,cutoff);
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
			double tmp = bplist->GetLogWeight(bp,cutoff);
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
			double logw = bplist->GetLogWeight(bp2,cutoff);

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

