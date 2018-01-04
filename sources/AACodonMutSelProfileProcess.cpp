
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "AACodonMutSelProfileProcess.h"
// include things here...


//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
//	* AACodonMutSelProfileProcess
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

void AACodonMutSelProfileProcess::Create()	{
	
	if ( (!nucrr) && (!nucstat) && (!codonprofile) )	{
		ProfileProcess::Create();
		Nnucrr = Nnuc * (Nnuc-1) / 2;
		nucrr = new double[Nnucrr];
		nucstat = new double[Nnuc];
		codonprofile = new double[GetNcodon()];

		SampleNucRR();
		SampleNucStat();
		SampleCodonProfile();
		SampleOmega();
	}
	else	{
		cerr << "Create of AAMutSelProfileProcess, nucrr and/or nucstat and/or codonprofile are/is not 0.\n";
		exit(1);
	}
}


void AACodonMutSelProfileProcess::Delete()	{

	if ( (nucrr) && (nucstat) )	{
		delete[] nucrr;
		delete[] nucstat;
		delete[] codonprofile;
		nucrr = 0;
		nucstat = 0;
		codonprofile = 0;

		/*
		delete[] nucstatsuffstatcount;
		delete[] nucstatsuffstatbeta;
		delete[] nucrrsuffstatcount;
		delete[] nucrrsuffstatbeta;
		*/

		ProfileProcess::Delete();
	}
	else	{
		cerr << "Delete of AAMutSelProfileProcess, nucrr and/or are/is 0.\n";
		exit(1);
	}
}


double AACodonMutSelProfileProcess::LogNucRRPrior()	{
	double total = 0;
	//for (int i=0; i<GetNnucrr(); i++)  {
	//	total -= nucrr[i];
	//}
	return total;
}


void AACodonMutSelProfileProcess::SampleNucRR()	{
	//cerr << "nucrr set for testing...\n";
	double total = 0;
	for (int i=0; i<GetNnucrr(); i++)  {
		nucrr[i] = rnd::GetRandom().sExpo();
		/*
		//check loglikelihood against paml with kappa = 2
		if (i==0) nucrr[i] = 1.0; // AC
		if (i==1) nucrr[i] = 2.0; // AG
		if (i==2) nucrr[i] = 1.0; // AT
		if (i==3) nucrr[i] = 1.0; // CG
		if (i==4) nucrr[i] = 2.0; // CT
		if (i==5) nucrr[i] = 1.0; // GT
		*/
		total += nucrr[i];
	}
	for (int i=0; i<GetNnucrr(); i++)  {
		nucrr[i] /= total;
	}
}


double AACodonMutSelProfileProcess::LogNucStatPrior()	{
	return 0;
}


void AACodonMutSelProfileProcess::SampleNucStat()	{

	//cerr << "nucstat set to specific values for testing...\n";
	double total = 0;
	for (int i=0; i<Nnuc; i++)  {
		//nucstat[i] = rnd::GetRandom().sExpo();
		//*
		//check loglikelihood against paml with specific values
		if (i==0) nucstat[i] = 0.20792;
		if (i==1) nucstat[i] = 0.26089;
		if (i==2) nucstat[i] = 0.29248;
		if (i==3) nucstat[i] = 0.23870;
		//*/
		total += nucstat[i];
	}
	for (int i=0; i<Nnuc; i++)  {
		nucstat[i] /= total;
	}
}


double AACodonMutSelProfileProcess::LogCodonProfilePrior()	{
	return 0;
}


void AACodonMutSelProfileProcess::SampleCodonProfile()	{
	double total = 0;
	for (int i=0; i<GetNcodon(); i++)	{
		//codonprofile[i] = rnd::GetRandom().sExpo();
		codonprofile[i] = 1.0/GetNcodon();
		total += codonprofile[i];
	}
	for (int i=0; i<GetNcodon(); i++)	{
		codonprofile[i] /= total;
	}
}

double AACodonMutSelProfileProcess::GlobalParametersMove()	{

	GlobalUpdateParameters();
	GlobalUpdateSiteProfileSuffStat();
	GlobalUpdateModeProfileSuffStat();
	double tuning = 1.0;
	int n = 1;
	if (! fixcodonprofile)	{
		MoveCodonProfile(tuning,30,10);
		MoveNucStatCodonProfile(tuning,30,10);
		MoveNucRR(tuning,2);
		MoveCodonProfile(tuning*0.5,30,10);
		MoveNucStatCodonProfile(tuning*0.5,30,10);
		MoveNucRR(tuning*0.5,2);
		MoveCodonProfile(tuning*0.1,30,10);
		MoveNucStatCodonProfile(tuning*0.1,30,10);
		MoveNucRR(tuning*0.1,2);
		MoveCodonProfile(tuning*0.01,30,10);
		MoveNucStatCodonProfile(tuning*0.01,30,10);
		MoveNucRR(tuning*0.01,2);
		MoveCodonProfile(tuning*0.001,30,10);
		MoveNucStatCodonProfile(tuning*0.001,30,10);
		MoveNucRR(tuning*0.001,2);
	}
	else	{
		MoveNucRR(tuning,2);
		MoveNucStat(tuning,2);
		MoveNucRR(tuning*0.1,2);
		MoveNucStat(tuning*0.1,2);
	}

	if (! fixomega)	{
		//cerr << "###############################\n";
		//cerr << "#   In !fixomega\n";
		//cerr << "###############################\n";
		//exit(1);

		GlobalUpdateParameters();
		GlobalUpdateSiteOmegaSuffStat();
		// CheckSuffStatLogProb();
		MoveOmega(tuning);
		// MoveOmega(tuning*0.3);
	}

	GlobalUpdateParameters();
}

double AACodonMutSelProfileProcess::MoveNucRR(double tuning)	{

	//exponental prior

	int naccepted = 0;
        for (int i=0; i<GetNnucrr(); i++)  {
                double deltalogratio = - LogNucRRPrior() - ProfileSuffStatLogProb();
                double bk = nucrr[i];
                double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
                // double m = tuning * (Random::Uniform() - 0.5);
                double e = exp(m);
                nucrr[i] *= e;
                UpdateMatrices();
                deltalogratio += LogNucRRPrior() + ProfileSuffStatLogProb();
                deltalogratio += m;
                int accepted = (rnd::GetRandom().Uniform() < exp(deltalogratio));
                // int accepted = (Random::Uniform() < exp(deltalogratio));
                if (accepted)   {
                        naccepted++;
                }
                else    {
                        nucrr[i] = bk;
                        UpdateMatrices();
                }
        }
        return ((double) naccepted) / GetNnucrr();
}

double AACodonMutSelProfileProcess::MoveNucRR(double tuning, int n)	{

	//dirichlet prior

	int naccepted = 0;
	double* bk = new double[GetNnucrr()];
	for (int k=0; k<GetNnucrr(); k++)  {
		bk[k] = nucrr[k];
	}
	double deltalogprob = -ProfileSuffStatLogProb();
	double loghastings = ProfileProposeMove(nucrr,tuning,n,GetNnucrr(),0,0);
	UpdateMatrices();
	deltalogprob += ProfileSuffStatLogProb();
	deltalogprob += loghastings;
	int accepted = (rnd::GetRandom().Uniform() < exp(deltalogprob));
	// int accepted = (Random::Uniform() < exp(deltalogprob));
	if (accepted)   {
		naccepted ++;
	}
	else    {
		for (int k=0; k<GetNnucrr(); k++)  {
			nucrr[k] = bk[k];
		}
		UpdateMatrices();
	}
	delete[] bk;
	return naccepted; 
}


double AACodonMutSelProfileProcess::MoveNucStat(double tuning, int n)	{

	int naccepted = 0;
	double* bk = new double[Nnuc];
	for (int k=0; k<Nnuc; k++)  {
		bk[k] = nucstat[k];
	}
	double deltalogprob = -ProfileSuffStatLogProb();
	double loghastings = ProfileProposeMove(nucstat,tuning,n,Nnuc,0,0);
	UpdateMatrices();
	deltalogprob += ProfileSuffStatLogProb();
	deltalogprob += loghastings;
	int accepted = (rnd::GetRandom().Uniform() < exp(deltalogprob));
	// int accepted = (Random::Uniform() < exp(deltalogprob));
	if (accepted)   {
		naccepted ++;
	}
	else    {
		for (int k=0; k<Nnuc; k++)  {
			nucstat[k] = bk[k];
		}
		UpdateMatrices();
	}
	delete[] bk;
	return naccepted; 
}

double AACodonMutSelProfileProcess::MoveCodonProfile(double tuning, int n, int nrep)	{

	int naccepted = 0;
	double* bk = new double[GetNcodon()];
	for (int rep=0; rep<nrep; rep++)	{
		for (int k=0; k<GetNcodon(); k++)  {
			bk[k] = codonprofile[k];
		}
		double deltalogprob = -ProfileSuffStatLogProb();
		double loghastings = ProfileProposeMove(codonprofile,tuning,n,GetNcodon(),0,0);
		UpdateMatrices();
		deltalogprob += ProfileSuffStatLogProb();
		deltalogprob += loghastings;
		int accepted = (rnd::GetRandom().Uniform() < exp(deltalogprob));
		if (accepted)   {
			naccepted ++;
		}
		else    {
			for (int k=0; k<GetNcodon(); k++)  {
				codonprofile[k] = bk[k];
			}
			UpdateMatrices();
		}
	}
	delete[] bk;
	return naccepted; 
}

double AACodonMutSelProfileProcess::MoveNucStatCodonProfile(double tuning, int n, int nrep)	{

	int naccepted = 0;
	double* bkcodon = new double[GetNcodon()];
	double* bknuc = new double[Nnuc];
	for (int rep=0; rep<nrep; rep++)	{
		for (int k=0; k<GetNcodon(); k++)  {
			bkcodon[k] = codonprofile[k];
		}
		for (int k=0; k<Nnuc; k++)  {
			bknuc[k] = nucstat[k];
		}
		double deltalogprob = -ProfileSuffStatLogProb();
		ProfileProposeMove(codonprofile,tuning,n,GetNcodon(),0,0);
		ProfileProposeMove(nucstat,tuning,n,Nnuc,0,0);
		UpdateMatrices();
		deltalogprob += ProfileSuffStatLogProb();
		int accepted = (rnd::GetRandom().Uniform() < exp(deltalogprob));
		if (accepted)   {
			naccepted ++;
		}
		else    {
			for (int k=0; k<GetNcodon(); k++)  {
				codonprofile[k] = bkcodon[k];
			}
			for (int k=0; k<Nnuc; k++)  {
				nucstat[k] = bknuc[k];
			}
			UpdateMatrices();
		}
	}
	delete[] bkcodon;
	delete[] bknuc;
	return naccepted; 
}

double AACodonMutSelProfileProcess::GetCodonProfileEntropy()	{
	double ent = 0;
	double tmp;
	for (int i=0; i<GetNcodon(); i++)	{
		tmp = codonprofile[i];
		if (tmp > 1e-6)	{
			ent -= tmp * log(tmp);
		}
	}
	return ent;
}

void AACodonMutSelProfileProcess::UpdateSiteOmegaSuffStat()	{
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i))	{

			siteomegasuffstatcount[i] = 0;
			siteomegasuffstatbeta[i] = 0;
			map<pair<int,int>, int>& paircount = GetSitePairCount(i);
			map<int,double>& waitingtime = GetSiteWaitingTime(i);
			AACodonMutSelProfileSubMatrix* codonmatrix = GetCodonMatrix(i);
			/*
			if (fabs(codonmatrix->GetOmega() - GetSiteOmega(i)) > 1e-6)	{
				cerr << "error: non matching omega : \n";
				double tmp = 0;
				cerr << codonmatrix->GetOmega() << '\t' << GetSiteOmega(i) << '\n';
				cerr << GetSiteOmega(i) << '\n';
				// cerr << GetSiteOmega(i) << '\t' << GetSiteOmega(0) << '\t' << GetSiteOmega(i) << '\t' << GetSiteOmega(0) << '\n';
				cerr << "site : " << i << '\t' << GetNsite() << '\n';
				cerr << GetSiteOmega(0) << '\n';
				exit(1);
			}
			*/
			for (map<int,double>::iterator j = waitingtime.begin(); j!= waitingtime.end(); j++)	{
				siteomegasuffstatbeta[i] += j->second * codonmatrix->RateAwayNonsyn(j->first) / codonmatrix->GetOmega();
				//siteomegasuffstatbeta[i] += j->second * codonmatrix->RateAwayNonsyn(j->first) / GetSiteOmega(i);
			}

			for (map<pair<int,int>, int>::iterator j = paircount.begin(); j!= paircount.end(); j++)	{
				if (! codonmatrix->Synonymous(j->first.first,j->first.second) )	{
					siteomegasuffstatcount[i] += j->second;
				}
			}
		}
	}
}

