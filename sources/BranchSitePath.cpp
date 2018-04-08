
#include "BranchSitePath.h"

void Plink::SetRelativeTime(double inrel_time) {
	rel_time = inrel_time;
	if (isnan(rel_time))	{
		cerr << "in Plink::SetRelativeTime: nan\n";
		exit(1);
	}
}

void BranchSitePath::AddRateSuffStat(double& count, double& beta, double factor, const double* rr, const double* stat, int nstate)	{
	Plink* link = Init();
	while (link)	{
		int state = link->GetState();
		double tmp = 0;
		for (int k=0; k<nstate; k++)	{
			if (k != state)	{
				tmp += rr[rrindex(state,k,nstate)] * stat[k];
			}
		}
		beta += GetRelativeTime(link) * factor * tmp;
		if (link != last)	{
			count ++;
		}
		link = link->Next();
	}
}

void BranchSitePath::AddProfileSuffStat(double* count, double* beta, double factor, const double* rr, int nstate)	{
	Plink* link = Init();
	while (link)	{
		int state = link->GetState();
		for (int k=0; k<nstate; k++)	{
			if (k!=state)	{
				beta[k] += GetRelativeTime(link) * factor * rr[rrindex(state,k,nstate)];
			}
		}
		if (link != last)	{
			count[link->Next()->GetState()] ++;
		}
		link = link->Next();
	}
}

void BranchSitePath::AddRRSuffStat(double* count, double* beta, double factor, SubMatrix* mat, const double* rr, int nstate)	{
	Plink* link = Init();
	while (link)	{
		int state = link->GetState();
		for (int k=0; k<nstate; k++)	{
			if (k!=state)	{
				beta[rrindex(state,k,nstate)] += GetRelativeTime(link) * factor * (*mat)(state,k) / rr[rrindex(state,k,nstate)];
			}
		}
		if (link != last)	{
			count[rrindex(state,link->Next()->GetState(),nstate)] ++;
		}
		link = link->Next();
	}
}

/*
void BranchSitePath::AddRRSuffStat(double* count, double* beta, double factor, const double* stat, int nstate)	{
	Plink* link = Init();
	while (link)	{
		int state = link->GetState();
		for (int k=0; k<nstate; k++)	{
			if (k!=state)	{
				beta[rrindex(state,k,nstate)] += GetRelativeTime(link) * factor * stat[k];
			}
		}
		if (link != last)	{
			count[rrindex(state,link->Next()->GetState(),nstate)] ++;
		}
		link = link->Next();
	}
}
*/

void BranchSitePath::AddGeneralPathRateSuffStat(double& count, double& beta, double factor, SubMatrix* mat)	{
	Plink* link = Init();
	while (link)	{
		int state = link->GetState();
		beta -= GetRelativeTime(link) * factor * (*mat)(state,state);
		if (link != last)	{
			count ++;
		}
		link = link->Next();
	}
}

void BranchSitePath::AddGeneralPathSuffStat(map<pair<int,int>,int>& paircount, map<int,double>& waitingtime, double factor)	{
	Plink* link = Init();
	while (link)	{
		int state = link->GetState();
		waitingtime[state] += GetRelativeTime(link) * factor;
		if (link != last)	{
			paircount[pair<int,int>(state,link->Next()->GetState())]++;
		}
		link = link->Next();
	}
}

