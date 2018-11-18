#ifndef PARS_H
#define PARS_H

#include "PoissonSubstitutionProcess.h"
#include "PhyloProcess.h"

int PhyloProcess::Parsimony(int* sitescore, double** profile)   {

    int score = BackwardParsimony(GetRoot(),sitescore);
    cerr << "score: " << score << '\n';
    ForwardParsimony(GetRoot(),profile);
    return score;
}

int PhyloProcess::BackwardParsimony(const Link* from, int* sitescore)  {

    int n = 0;
    if (from->isLeaf()) {
        Initialize(GetConditionalLikelihoodVector(from),GetData(from),true);
    }
    else    {

		Reset(GetConditionalLikelihoodVector(from),true,0);
        int nchildren = 0;
        for (const Link* link=from->Next(); link!=from; link=link->Next())  {
            n += BackwardParsimony(link->Out(),sitescore);
			Add(GetConditionalLikelihoodVector(link->Out()),GetConditionalLikelihoodVector(from),true);
            nchildren++;
        }
        n += GetParsimonyBranchScore(nchildren,GetConditionalLikelihoodVector(from),sitescore);
    }
    return n;
}

int SubstitutionProcess::GetParsimonyBranchScore(int nchildren, double*** t, int* sitescores) {

    int total = 0;
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i))	{
            int j = ratealloc[i];
            double* to = t[i][j];
            int nstate = GetNstate(i);
            int max = 0;
            for (int k=0; k<nstate; k++)	{
                if (max < to[k])    {
                    max = to[k];
                }
            }
            int n = nchildren - max;
            if (max > 1)  {
                for (int k=0; k<nstate; k++)	{
                    if (to[k] == max)    {
                        to[k] = 1;
                    }
                    else    {
                        to[k] = 0;
                    }
                }
            }
            if (sitescores)  {
                sitescores[i] += n;
            }
            total += n;
        }
	}
    return total;
}

void SubstitutionProcess::ChooseParsimonyStates(const int* upstates, int* downstates, double*** t)    {

	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i))	{
            int j = ratealloc[i];
            double* to = t[i][j];
            int nstate = GetNstate(i);

            if ((upstates != downstates) && (to[upstates[i]])) {
                downstates[i] = upstates[i];
            }
            else    {
                double n = 0;
                double cumul[nstate];
                for (int k=0; k<nstate; k++)    {
                    n += to[k];
                    cumul[k] = n;
                }
                double u = rnd::GetRandom().Uniform() * n;
                int k = 0;
                while ((k<nstate) && (u>cumul[k]))  {
                    k++;
                }
                if (k == nstate)    {
                    cerr << "error in choose pars state: overflow\n";
                    exit(1);
                }
                downstates[i] = k;
            }
        }
    }
}

void SubstitutionProcess::ParsimonyProfiles(const int* upstates, const int* nodestates, double** profile)   {

    for (int i=GetSiteMin(); i<GetSiteMax(); i++)   {
        if (ActiveSite(i))  {
            if ((nodestates == upstates) || (nodestates[i] != upstates[i]))   {
                profile[i][nodestates[i]]++;
            }
        }
    }
}

void PoissonSubstitutionProcess::ParsimonyProfiles(const int* upstates, const int* nodestates, double** profile)    {

    for (int i=GetSiteMin(); i<GetSiteMax(); i++)   {
        if (ActiveSite(i))  {
            if (! GetOrbitSize(i)) {
                for (int k=0; k<GetDim(); k++)  {
                    profile[i][k] = 1.0 / GetDim();
                }
            }
            else    {
            if ((nodestates == upstates) || (nodestates[i] != upstates[i]))   {
                if (nodestates[i] >= GetOrbitSize(i))   {
                    cerr << "error in parsimony get state from zip: " << nodestates[i] << '\t' << GetOrbitSize(i) << '\n';
                    cerr << "site : " << i << '\n';
                    exit(1);
                }
                int truestate = GetStateFromZip(i,nodestates[i]);
                profile[i][truestate]++;
            }
            }
        }
    }
}

void PhyloProcess::ForwardParsimony(const Link* from, double** profile)    {
	
    int* nodestates = GetStates(from->GetNode());
    int* upstates = GetStates(from->Out()->GetNode());
    ChooseParsimonyStates(nodestates,upstates,GetConditionalLikelihoodVector(from));
    ParsimonyProfiles(upstates,nodestates,profile);

	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
        ForwardParsimony(link->Out(),profile);
    }
}

#endif
