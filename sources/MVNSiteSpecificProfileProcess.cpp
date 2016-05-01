
#include "MVNSiteSpecificProfileProcess.h"

void MVNSiteSpecificProfileProcess::Create()	{

	if (! logprofile)	{
		MVNProfileProcess::Create();
		SiteSpecificProfileProcess::Create();
		alloclogprofile = new double[GetNsite() * GetDim()];
		logprofile = new double*[GetNsite()];
		for (int i=0; i<GetNsite(); i++)	{
			logprofile[i] = alloclogprofile + i*GetDim();
		}
	}
}

void MVNSiteSpecificProfileProcess::Delete()	{

	if (logprofile)	{
		delete[] alloclogprofile;
		delete[] logprofile;
		logprofile = 0;
		SiteSpecificProfileProcess::Delete();
		MVNProfileProcess::Delete();
	}
}

void MVNSiteSpecificProfileProcess::SampleStat(int i)	{
	SampleFrequencyStat(logprofile[i]);
	UpdateSite(i);
}

void MVNSiteSpecificProfileProcess::ResampleCovMatrix()	{

	if (GetNsite() != GetNactiveSite())	{
		cerr << "error in MVNSiteSpecificProfileProcess::ResampleCovMatrix: should implement active / inactive sites\n";
		exit(1);
	}
	covmatrix->GibbsResample(logprofile,GetNsite());
}

double MVNSiteSpecificProfileProcess::MoveProfile(int site, double tuning, int n, int nrep)	{

	int naccepted = 0;
	double* bk = new double[GetDim()];
	for (int k=0; k<GetDim(); k++)	{
		bk[k] = logprofile[site][k];
	} 
	for (int rep=0; rep<nrep; rep++)	{
		double deltalogprob = - LogStatPrior(site) - ProfileSuffStatLogProb(site);
		double loghastings = LogProfileProposeMove(logprofile[site],tuning,n);
		UpdateSite(site);
		deltalogprob += LogStatPrior(site) + ProfileSuffStatLogProb(site);
		deltalogprob += loghastings;
		int accepted = (log(rnd::GetRandom().Uniform()) < deltalogprob);
		if (accepted)	{
			naccepted ++;
			for (int k=0; k<GetDim(); k++)	{
				bk[k] = logprofile[site][k];
			} 
		}
		else	{
			for (int k=0; k<GetDim(); k++)	{
				logprofile[site][k] = bk[k];
			} 
			UpdateSite(site);
		}
	}
	delete[] bk;
	return naccepted / nrep;
}
