
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
#include "IIDDirichletIIDGammaPhyloProcess.h"
#include "Parallel.h"
#include <string>

void IIDDirichletIIDGammaPhyloProcess::VarBayes()    {

    double* meanrate = new double[GetNsite()];
    for (int i=0; i<GetNsite(); i++)    {
        meanrate[i] = 0;
    }
    double* meanblarray = new double[GetNbranch()];
    for (int j=0; j<GetNbranch(); j++)  {
        meanblarray[j] = 0;
    }
    double** meanprofile = new double*[GetNsite()];
    for (int i=0; i<GetNsite(); i++)    {
        meanprofile[i] = new double[GetDim()];
        for (int k=0; k<GetDim(); k++)  {
            meanprofile[i][k] = 0;
        }
    }

    int burnin = 10;
    int nrep = 100;

    for (int rep=0; rep<burnin; rep++)  {
        Move(1.0);
        Trace(cerr);
    }
    cerr << '\n';

    for (int rep=0; rep<nrep; rep++)    {
        Move(1.0);
        for (int i=0; i<GetNsite(); i++)    {
            meanrate[i] += rate[i];
        }
        for (int j=0; j<GetNbranch(); j++)  {
            meanblarray[j] += blarray[j];
        }
        for (int i=0; i<GetNsite(); i++)    {
            for (int k=0; k<GetDim(); k++)  {
                meanprofile[i][k] += profile[i][k];
            }
        }
        Trace(cerr);
    }

    cerr << '\n';

    double tot = 0;
    for (int i=0; i<GetNsite(); i++)    {
        meanrate[i] /= nrep;
        rate[i] = meanrate[i];
        tot += rate[i];
    }
    tot /= GetNsite();
    for (int i=0; i<GetNsite(); i++)    {
        meanrate[i] /= tot;
    }

    for (int j=0; j<GetNbranch(); j++)  {
        meanblarray[j] /= nrep;
        blarray[j] = meanblarray[j];
    }
    for (int i=0; i<GetNsite(); i++)    {
        for (int k=0; k<GetDim(); k++)  {
            meanprofile[i][k] /= nrep;
            profile[i][k] = meanprofile[i][k];
        }
    }

    UpdateZip();

    while(1)    {

        UpdateMeanSuffStat();
        int nrep = 200;
        for (int rep=0; rep<nrep; rep++)    {
            UpdateVarLengths();
            UpdateVarRates();
        }
        UpdateVarProfiles();
        UpdateZip();
        UpdateConditionalLikelihoods();
        double logp = GetVarLogMarginalLikelihood();
        double logl = GetLogLikelihood();
        cerr << logl << '\t' << (logl-logp)/GetNsite() << '\t' << logp << '\t';
        Trace(cerr);
    }
}

double IIDDirichletIIDGammaPhyloProcess::GetVarLogMarginalLikelihood()  {

    double logprior = 0;
    double logpost = 0;
    double logl = GetLogLikelihood();

    for (int j=1; j<GetNbranch(); j++)    {
        double a = branchalpha;
        double b = branchbeta;
        double x = blarray[j];
        logprior += a*log(b) - rnd::GetRandom().logGamma(a) + (a-1)*log(x) - b*x;
        a += branchlengthsuffstatcount[j];
        b += branchlengthsuffstatbeta[j];
        logpost += a*log(b) - rnd::GetRandom().logGamma(a) + (a-1)*log(x) - b*x;
    }

    for (int i=0; i<GetNsite(); i++)    {
        double a = alpha;
        double b = alpha;
        double x = rate[i];
        logprior += a*log(b) - rnd::GetRandom().logGamma(a) + (a-1)*log(x) - b*x;
        a += siteratesuffstatcount[i];
        b += siteratesuffstatbeta[i];
        logpost += a*log(b) - rnd::GetRandom().logGamma(a) + (a-1)*log(x) - b*x;
    }

    for (int i=0; i<GetNsite(); i++)    {

        double priorw = 0;
        double postw = 0;
        for (int k=0; k<GetDim(); k++)  {

            logprior += (dirweight[k] - 1)*log(profile[i][k]);
            logprior -= rnd::GetRandom().logGamma(dirweight[k]);
            priorw += dirweight[k];

            logpost += (dirweight[k] + siteprofilesuffstatcount[i][k] -1)*log(profile[i][k]);
            logpost -= rnd::GetRandom().logGamma(dirweight[k] + siteprofilesuffstatcount[i][k]);
            postw += dirweight[k] + siteprofilesuffstatcount[i][k];

        }

        logprior += rnd::GetRandom().logGamma(priorw);
        logpost += rnd::GetRandom().logGamma(postw);
    }

    double logp = logl + logprior - logpost;

    return logp;
}

void IIDDirichletIIDGammaPhyloProcess::UpdateVarLengths() {

    double A = GetMeanRate() * GetNsite();
	for (int j=0; j<GetNbranch(); j++)	{
		branchlengthsuffstatbeta[j] = A;
    }

    if (blvarmode == 0)   { // EM
        for (int j=1; j<GetNbranch(); j++)  {
            blarray[j] = branchlengthsuffstatcount[j] / branchlengthsuffstatbeta[j];
        }
    }
    else if (blvarmode == 1)   { // Mean Var Bayes
        for (int j=1; j<GetNbranch(); j++)  {
            blarray[j] = (branchalpha + branchlengthsuffstatcount[j]) / (branchbeta + branchlengthsuffstatbeta[j]);
        }
    }
    else if (blvarmode == 2) {
        for (int j=1; j<GetNbranch(); j++)  {
            blarray[j] = exp(rnd::GetRandom().Psi(branchalpha + branchlengthsuffstatcount[j]) - log(branchbeta + branchlengthsuffstatbeta[j])) + 1e-8;
        }
    }
    for (int j=1; j<GetNbranch(); j++)  {
        if (isnan(blarray[j]))  {
            cerr << "in varlength: nan\n";
            exit(1);
        }
        if (isinf(blarray[j]))  {
            cerr << "in var length: inf\n";
            exit(1);
        }
    }
}

void IIDDirichletIIDGammaPhyloProcess::UpdateVarRates() {

    double B = GetTotalLength();
    double totcount = 0;
    for (int i=0; i<GetNsite(); i++)    {
        siteratesuffstatbeta[i] = B;
        totcount += siteratesuffstatcount[i];
    }

    if (ratevarmode == 0)   {
        double tot = 0;
        for (int i=0; i<GetNsite(); i++)    {
            rate[i] = siteratesuffstatcount[i] / siteratesuffstatbeta[i];
            tot += rate[i];
        }
    }
    else if (ratevarmode == 1)   {
        double tot = 0;
        for (int i=0; i<GetNsite(); i++)    {
            rate[i] = (alpha + siteratesuffstatcount[i]) / (alpha + siteratesuffstatbeta[i]);
            tot += rate[i];
        }
        tot /= GetNsite();
        for (int i=0; i<GetNsite(); i++)    {
            rate[i] /= tot;
        }
    }
    else if (ratevarmode == 2)   {
        for (int i=0; i<GetNsite(); i++)    {
            rate[i] = exp(rnd::GetRandom().Psi(alpha + siteratesuffstatcount[i]) - log(alpha + siteratesuffstatbeta[i])) + 1e-8;
            if (isnan(rate[i])) {
                cerr << "nan\n";
                cerr << alpha << '\t' << siteratesuffstatcount[i] << '\t' << siteratesuffstatbeta[i] << '\n';
                cerr << rnd::GetRandom().Psi(alpha + siteratesuffstatcount[i]) << '\n';
                cerr << log(alpha + siteratesuffstatbeta[i]) << '\n';
                exit(1);
            }
        }
    }
    else    {
        cerr << "rate var mode not recognized\n";
        exit(1);
    }
    for (int i=0; i<GetNsite(); i++)    {
        if (isnan(rate[i])) {
            cerr << "in varrates: nan\n";
            exit(1);
        }
        if (isinf(rate[i])) {
            cerr << "in varrates: inf\n";
            exit(1);
        }
    }
}

void IIDDirichletIIDGammaPhyloProcess::UpdateVarProfiles() {

    if (profilevarmode == 0)    {
        for (int i=0; i<GetNsite(); i++)    {
            double total = 0;
            for (int k=0; k<GetDim(); k++)   {
                profile[i][k] = siteprofilesuffstatcount[i][k] + 1e-8;
                total += profile[i][k];
            }
            for (int k=0; k<GetDim(); k++)  {
                profile[i][k] /= total;
            }
        }
    }
    else if (profilevarmode == 1)   {
        for (int i=0; i<GetNsite(); i++)    {
            double total = 0;
            for (int k=0; k<GetDim(); k++)  {
                profile[i][k] = dirweight[k] + siteprofilesuffstatcount[i][k];
                total += profile[i][k];
            }
            for (int k=0; k<GetDim(); k++)  {
                profile[i][k] /= total;
            }
        }
    }
    else if (profilevarmode == 2)   {
        for (int i=0; i<GetNsite(); i++)    {
            double norm = 0;
            double tot = 0;
            for (int k=0; k<GetDim(); k++)  {
                profile[i][k] = exp(rnd::GetRandom().Psi(dirweight[k] + siteprofilesuffstatcount[i][k]));
                norm+= profile[i][k];
                // tot += dirweight[k] + siteprofilesuffstatcount[i][k];
            }
            /*
            double z = exp(rnd::GetRandom().Psi(tot));
            cerr << z << '\t' << norm << '\n';
            cerr << log(z) << '\t' << log(norm) << '\n';
            for (int k=0; k<GetDim(); k++)  {
                cerr << dirweight[k] << '\t' << siteprofilesuffstatcount[i][k] << '\n';
            }
            exit(1);
            if (fabs(z - norm) > 1e-6)  {
                cerr << "in var profiles: normalization error\n";
                exit(1);
            }
            */
            for (int k=0; k<GetDim(); k++)  {
                profile[i][k] /= norm;
            }
        }
    }
}


void IIDDirichletIIDGammaPhyloProcess::GlobalUpdateParameters()	{

	if (GetNprocs() > 1)	{

        RenormalizeProfiles();

        int nd = 1 + GetNsite() + GetNbranch() + GetNsite()*GetDim() + GetDim();
        double* dvector = new double[nd];

        MESSAGE signal = PARAMETER_DIFFUSION;
        MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

        // First we assemble the vector of doubles for distribution
        int index = 0;
        dvector[index] = GetAlpha();
        index++;
        for (int i=0; i<GetNsite(); i++)    {
            dvector[index] = rate[i];
            index++;
        }

        for(int i=0; i<GetNbranch(); i++) {
            dvector[index] = blarray[i];
            index++;
        }
        
        for(int i=0; i<GetNsite(); ++i) {
            for(int j=0; j<GetDim(); ++j) {
                dvector[index] = profile[i][j];
                index++;
            }
        }
        for (int i=0; i<GetDim(); i++)	{
            dvector[index] = dirweight[i];
            index++;
        }

        // Now send out the doubles and ints over the wire...
        MPI_Bcast(dvector,nd,MPI_DOUBLE,0,MPI_COMM_WORLD);

        delete[] dvector;

	}
	else	{
		UpdateZip();
	}
}

void IIDDirichletIIDGammaPhyloProcess::SlaveUpdateParameters()	{

	int nd = 1 + GetNsite() + GetNbranch() + GetNsite()*GetDim() + GetDim();
    double* dvector = new double[nd];

	MPI_Bcast(dvector,nd,MPI_DOUBLE,0,MPI_COMM_WORLD);

	int index = 0;
    alpha = dvector[index];
    index++;

    for (int i=0; i<GetNsite(); i++)    {
        rate[i] = dvector[index];
        index++;
    }

	for(int i=0; i<GetNbranch(); ++i) {
		blarray[i] = dvector[index];
		index++;
	}

	for(int i=0; i<GetNsite(); ++i) {
		for(int j=0; j<GetDim(); ++j) {
			profile[i][j] = dvector[index];
			index++;
		}
	}
	for (int i=0; i<GetDim(); i++)	{
		dirweight[i] = dvector[index];
		index++;
	}

	delete[] dvector;

	UpdateZip();
}

void IIDDirichletIIDGammaPhyloProcess::SlaveExecute(MESSAGE signal)	{

	switch(signal) {

	case UPDATE_HYPERRATE:
		SlaveUpdateRateHyperSuffStat();
		break;
	case UPDATE_HYPERPROFILE:
		SlaveUpdateProfileHyperSuffStat();
		break;
    case MOVE_SRATE:
        ResampleSiteRates();
        break;
    case MOVE_SPROFILE:
        ResampleSiteProfiles();
        break;
    case COLLECT_SRATE:
        SlaveCollectSiteRates();
        break;
    case COLLECT_SPROFILE:
        SlaveCollectSiteProfiles();
        break;
	default:
		PhyloProcess::SlaveExecute(signal);
	}
}


