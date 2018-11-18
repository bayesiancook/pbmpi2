
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "StringStreamUtils.h"
#include "RASIIDDirichletGammaPhyloProcess.h"
#include "Parallel.h"
#include <string>

void RASIIDDirichletGammaPhyloProcess::EM(double cutoff, int nrep)   {

    if ((!cutoff) && (!nrep))   {
        cerr << "error in RASCATFiniteGammaPhyloProcess::EM: either cutoff or nrep should be strictly positive\n";
        exit(1);
    }

	modesitelogL = new double*[GetNsite()];
	modesitepostprob = new double*[GetNsite()];
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i))	{
            modesitelogL[i] = new double[GetNcat()];
            modesitepostprob[i] = new double[GetNcat()];
		}
	}

    if (sitefreq != "free") {
        cerr << "set site profiles\n";
        SetSiteProfiles();
    }

    // max parsimony : not working yet
    if (initprofile == 1)    {

        cerr << "max pars\n";

        int* sitescore = new int[GetNsite()];

        double* tmp = GetEmpiricalFreq();
        for (int i=0; i<GetNsite(); i++)    {
            sitescore[i] = 0;
            for (int k=0; k<GetDim(); k++)  {
                profile[i][k] = pseudocount * tmp[k];
            }
        }

        int totscore = Parsimony(sitescore,profile);

        cerr << "parsimony score: " << totscore << '\n';
        cerr << "per site : " << ((double) totscore) / GetNsite() << '\n';

        RenormalizeProfiles();

        ofstream os("sitescores");
        for (int i=0; i<GetNsite(); i++)    {
            os << sitescore[i] << '\t';
            for (int k=0; k<GetDim(); k++)  {
                os << profile[i][k] << '\t';
            }
            os << '\n';
        }

        delete[] sitescore;
    }

    int rep = 0;
    double diff = 2*cutoff;
    double currentlogl = 0;
    while ((nrep && (rep<nrep)) || (cutoff && (diff > cutoff))) {

        double logl1 = EMUpdateMeanSuffStat();
        EM_UpdateBranchLengths();

        if (! fixprofile)   {
            // double logl3 = EMUpdateMeanSuffStat();
            EM_UpdateProfiles(pseudocount);
        }

        double logl2 = EMUpdateMeanSuffStat();
        EM_UpdateAlpha(0.1,10.0,0.01);

        cout << logl2 << '\t';
        cout << GetRenormTotalLength() << '\t' << GetAlpha();
        // if (! fixprofile)   {
           cout << '\t' << GetStatEnt() << '\t' << GetMeanDirWeight();
        // }
        cout << '\n';

        if (rep)    {
            diff = logl2 - currentlogl;
        }
        rep++;
        currentlogl = logl2;
    }

	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i))	{
            delete[] modesitelogL[i];
            delete[] modesitepostprob[i];
		}
	}

	delete[] modesitelogL;
	delete[] modesitepostprob;

}

void RASIIDDirichletGammaPhyloProcess::EM_UpdateProfiles(double pseudocount)  {

    double* tmp = GetEmpiricalFreq();
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i))	{
            double total = 0;
            for (int k=0; k<GetDim(); k++)	{
                profile[i][k] = siteprofilesuffstatcount[i][k] + pseudocount * tmp[k];
                total += profile[i][k];
            }
            for (int k=0; k<GetDim(); k++)	{
                profile[i][k] /= total;
            }
		}
	}
    UpdateZip();
}

double RASIIDDirichletGammaPhyloProcess::EMUpdateMeanSuffStat()  {

    FillMissingMap(0);
    InactivateSumOverRateAllocations();

    UpdateZip();

    for (int l=0; l<Ncat; l++)  {

        for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
            if (ActiveSite(i))  {
                ratealloc[i] = l;
            }
        }

        UpdateConditionalLikelihoods();

        for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
            if (ActiveSite(i))	{
                modesitelogL[i][l] = sitelogL[i];
            }
        }
    }

    // calculate marginal site log likelihoods and allocation post probs
	double totlogL = 0;
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i))	{

			double max = modesitelogL[i][0];
            for (int l=0; l<Ncat; l++)  {
                if (max < modesitelogL[i][l])	{
                    max = modesitelogL[i][l];
                }
            }

			double total = 0;
            for (int l=0; l<Ncat; l++)  {
                double tmp = 1.0 / Ncat * exp(modesitelogL[i][l] - max);
                modesitepostprob[i][l] = tmp;
                total += tmp;
			}

            for (int l=0; l<Ncat; l++)  {
                modesitepostprob[i][l] /= total;
            }

			double sitetotlogL = log(total) + max;
			totlogL += sitetotlogL;
		}
	}

    for (int k=0; k<Ncat; k++)  {
        ratesuffstatcount[k] = 0;
        ratesuffstatbeta[k] = 0;
    }

    // reset suffstats
    for (int j=0; j<GetNbranch(); j++)	{
        branchlengthsuffstatcount[j] = 0;
        branchlengthsuffstatbeta[j] = 0;
    }

    for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
        if (ActiveSite(i))	{
            for (int k=0; k<GetDim(); k++)	{
                siteprofilesuffstatcount[i][k] = 0;
            }
        }
    }

    double* sitepostprob = new double[GetNsite()];

    // now, redo a pass over all components and for all sites
    // this time, sum sufficient satistics along the way

    for (int l=0; l<Ncat; l++)  {

        for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
            if (ActiveSite(i))  {
                ratealloc[i] = l;
            }
        }

        UpdateConditionalLikelihoods();

        for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
            if (ActiveSite(i))  {
                sitepostprob[i] = modesitepostprob[i][l];
            }
        }

        for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
            if (ActiveSite(i))	{
                siteratesuffstatcount[i] = 0;
                siteratesuffstatbeta[i] = 0;
            }
        }

        RecursiveUpdateMeanSuffStat(GetRoot(),condlmap[0],sitepostprob);

        for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
            if (ActiveSite(i))	{
                ratesuffstatcount[l] += siteratesuffstatcount[i];
                ratesuffstatbeta[l] += siteratesuffstatbeta[i];
            }
        }
    }

    delete[] sitepostprob;
	return totlogL;
}

void RASIIDDirichletGammaPhyloProcess::GlobalUpdateParameters()	{

	if (GetNprocs() > 1)	{

        int nd = 2 + GetNbranch();
        if (! fixprofile)   {
           nd += GetNsite()*GetDim() + GetDim();
        }
        double* dvector = new double[nd];

        MESSAGE signal = PARAMETER_DIFFUSION;
        MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

        // First we assemble the vector of doubles for distribution
        int index = 0;
        dvector[index] = GetAlpha();
        index ++;
        dvector[index] = GetPinv();
        index ++;

        for(int i=0; i<GetNbranch(); i++) {
            dvector[index] = blarray[i];
            index++;
        }
        
        if (! fixprofile)   {
            for(int i=0; i<GetNsite(); ++i) {
                for(int j=0; j<GetDim(); ++j) {
                    dvector[index] = profile[i][j];
                    index++;
                }
            }
            for (int i=0; i<GetDim(); i++)	{
                dvector[index] = dirweight[0][i];
                index++;
            }
        }

        // Now send out the doubles and ints over the wire...
        MPI_Bcast(dvector,nd,MPI_DOUBLE,0,MPI_COMM_WORLD);

        delete[] dvector;

	}
	else	{
		UpdateZip();
	}
}

void RASIIDDirichletGammaPhyloProcess::SlaveUpdateParameters()	{

	int nd = 2 + GetNbranch();
    if (! fixprofile)   {
       nd += GetNsite()*GetDim() + GetDim();
    }
    double* dvector = new double[nd];

	MPI_Bcast(dvector,nd,MPI_DOUBLE,0,MPI_COMM_WORLD);

	int index = 0;
	SetRateParams(dvector[index],dvector[index+1]);
    index+=2;

	for(int i=0; i<GetNbranch(); ++i) {
		blarray[i] = dvector[index];
		index++;
	}

    if (! fixprofile)   {
        for(int i=0; i<GetNsite(); ++i) {
            for(int j=0; j<GetDim(); ++j) {
                profile[i][j] = dvector[index];
                index++;
            }
        }
        for (int i=0; i<GetDim(); i++)	{
            dirweight[0][i] = dvector[index];
            index++;
        }
    }

	delete[] dvector;

    if (! fixprofile)   {
        UpdateZip();
    }
}

void RASIIDDirichletGammaPhyloProcess::SlaveExecute(MESSAGE signal)	{

	switch(signal) {

	case UPDATE_RATE:
		SlaveUpdateRateSuffStat();
		break;
    case MOVE_SPROFILE:
        ResampleSiteProfiles();
        break;
    case COLLECT_SPROFILE:
        SlaveCollectSiteProfiles();
        break;
	default:
		PhyloProcess::SlaveExecute(signal);
	}
}


