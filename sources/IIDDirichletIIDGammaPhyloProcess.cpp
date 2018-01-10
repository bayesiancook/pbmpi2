
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

void IIDDirichletIIDGammaPhyloProcess::ComputeMeanRates()   {

    totmeanrate = 0;
    for (int i=0; i<GetNsite(); i++)    {
        meanrate[i] = alphastar[i] / betastar[i];
        totmeanrate += meanrate[i];
        meanlograte[i] = rnd::GetRandom().Psi(alphastar[i]) - log(betastar[i]);
    }
}

void IIDDirichletIIDGammaPhyloProcess::ComputeMeanLengths() {

    totmeanlength = 0;
    for (int j=1; j<GetNbranch(); j++)  {
        meanbl[j] = lambdastar[j] / mustar[j];
        totmeanlength += meanbl[j];
        meanlogbl[j] = rnd::GetRandom().Psi(lambdastar[j]) - log(mustar[j]);
    }
}

void IIDDirichletIIDGammaPhyloProcess::ComputeMeanProfiles()    {

    totmeanlogprofilenorm = 0;
    for (int i=0; i<GetNsite(); i++)    {
        double totgamma = 0;
        for (int k=0; k<GetDim(); k++)  {
            meanlogprofile[i][k] = rnd::GetRandom().Psi(gammastar[i][k]);
            totgamma += gammastar[i][k];
        }
        double total = 0;
        for (int k=0; k<GetDim(); k++)  {
            meanlogprofile[i][k] -= rnd::GetRandom().Psi(totgamma);
            total += exp(meanlogprofile[i][k]);
        }
        meanlogprofilenorm[i] = log(total);
        totmeanlogprofilenorm += meanlogprofilenorm[i];
    }
}

void IIDDirichletIIDGammaPhyloProcess::PosteriorMean(int burnin, int nrep) {

    for (int i=0; i<GetNsite(); i++)    {
        meanlograte[i] = 0;
        meanrate[i] = 0;
    }
    for (int j=0; j<GetNbranch(); j++)  {
        meanlogbl[j] = 0;
        meanbl[j] = 0;
    }
    for (int i=0; i<GetNsite(); i++)    {
        for (int k=0; k<GetDim(); k++)  {
            meanlogprofile[i][k] = 0;
        }
    }

    for (int rep=0; rep<burnin; rep++)  {
        Move(1.0);
        Trace(cerr);
    }
    cerr << '\n';

    for (int rep=0; rep<nrep; rep++)    {
        Move(1.0);
        for (int i=0; i<GetNsite(); i++)    {
            meanlograte[i] += log(rate[i]);
            meanrate[i] += rate[i];
        }
        for (int j=0; j<GetNbranch(); j++)  {
            meanlogbl[j] += log(blarray[j]);
            meanbl[j] += blarray[j];
        }
        for (int i=0; i<GetNsite(); i++)    {
            for (int k=0; k<GetDim(); k++)  {
                meanlogprofile[i][k] += log(profile[i][k]);
            }
        }
        Trace(cerr);
    }

    cerr << '\n';

    totmeanrate = 0;
    for (int i=0; i<GetNsite(); i++)    {
        meanlograte[i] /= nrep;
        meanrate[i] /= nrep;
        totmeanrate += meanrate[i];
    }

    totmeanlength = 0;
    for (int j=0; j<GetNbranch(); j++)  {
        meanlogbl[j] /= nrep;
        meanbl[j] /= nrep;
        totmeanlength += meanbl[j];
    }
    totmeanlogprofilenorm = 0;
    for (int i=0; i<GetNsite(); i++)    {
        double tot = 0;
        for (int k=0; k<GetDim(); k++)  {
            meanlogprofile[i][k] /= nrep;
            tot += exp(meanlogprofile[i][k]);
        }
        meanlogprofilenorm[i] = log(tot);
        totmeanlogprofilenorm += log(tot);
    }
}

void IIDDirichletIIDGammaPhyloProcess::AddRateCorrection(int sign)  {

    for (int i=0; i<GetNsite(); i++)    {
        if (sign > 0)   {
            rate[i] *= bkrate[i];
        }
        else    {
            rate[i] /= bkrate[i];
        }
    }
}

void IIDDirichletIIDGammaPhyloProcess::VarBayes()    {

    varfreebl = 1;
    varfreerate = 1;
    varfreeprofile = 1;

    int nrep = 20;
    ofstream os("finalscore");


    for (int rep=0; rep<10; rep++)  {
        if (varfreebl)  {
            SampleLength();
        }
        else    {
            // get length from file
            ifstream is("bl");
            for (int j=0; j<GetNbranch(); j++)  {
                is >> blarray[j];
            }
        }

        for (int j=0; j<GetNbranch(); j++)  {
            meanbl[j] = blarray[j];
        }

        if (varfreerate)    {
            SampleRate();
        }
        else    {
            for (int i=0; i<GetNsite(); i++)    {
                rate[i] = 1.0;
            }
        }

        for (int i=0; i<GetNsite(); i++)    {
            meanrate[i] = rate[i];
        }

        if (varfreeprofile) {
            SampleStat();
        }
        else    {
            for (int i=0; i<GetNsite(); i++)    {
                for (int k=0; k<GetDim(); k++)  {
                    profile[i][k] = 1.0/GetDim();
                }
            }
        }
        UpdateZip();

        PosteriorMean(20,100);
        SetNewParameters();

        totmeanrate = 0;
        for (int i=0; i<GetNsite(); i++)    {
            totmeanrate += meanrate[i];
        }
        totmeanlength = 0;
        for (int j=1; j<GetNbranch(); j++)  {
            totmeanlength += meanbl[j];
        }

        double logl = GetLogLikelihood();
        double logp = 0;
        double diff = 1.0;
        while (diff > 1e-4) {

            UpdateMeanSuffStat();

            if (varfreerate && varfreebl)   {
                for (int rep=0; rep<500; rep++)    {
                    UpdateVarLengths();
                    UpdateVarRates();
                }
            }
            else if (varfreerate)   {
                UpdateVarRates();
            }
            else if (varfreebl) {
                UpdateVarLengths();
            }
            if (varfreeprofile)    {
                UpdateVarProfiles();
            }
            SetNewParameters();
            UpdateZip();
            UpdateConditionalLikelihoods();
            double tmp = GetVarLogMarginalLikelihood();
            diff = tmp - logp;
            if (! logp) {
                diff= 1;
            }
            logp = tmp;
            logl = GetLogLikelihood();
            cerr << diff << '\t' << logl << '\t' << (logl-logp)/GetNsite() << '\t' << logp << '\t';
            Trace(cerr);
        }
        cerr << '\n';
        os << logl << '\t' << (logl-logp)/GetNsite() << '\t' << logp << '\t';
        Trace(os);
        os.flush();
    }
}

void IIDDirichletIIDGammaPhyloProcess::SetNewParameters()   {

    if (varfreebl)  {
        for (int j=1; j<GetNbranch(); j++)  {
            blarray[j] = exp(meanlogbl[j]);
        }
    }

    if (varfreerate && varfreeprofile)    {
        for (int i=0; i<GetNsite(); i++)    {
            rate[i] = exp(meanlograte[i] + meanlogprofilenorm[i]);
        }
    }
    else if (varfreerate)   {
        for (int i=0; i<GetNsite(); i++)    {
            rate[i] = exp(meanlograte[i]);
        }
    }
    else if (varfreeprofile)    {
        for (int i=0; i<GetNsite(); i++)    {
            rate[i] = exp(meanlogprofilenorm[i]);
        }
    }

    if (varfreeprofile) {
        for (int i=0; i<GetNsite(); i++)    {
            double tot = 0;
            for (int k=0; k<GetDim(); k++)  {
                profile[i][k] = exp(meanlogprofile[i][k] - meanlogprofilenorm[i]);
            }
        }
    }
}

double IIDDirichletIIDGammaPhyloProcess::GetVarLogMarginalLikelihood()  {

    double logprior = 0;
    double logpost = 0;
    double logl = GetLogLikelihood();

    if (varfreebl)  {
        for (int j=1; j<GetNbranch(); j++)    {
            double a = branchalpha;
            double b = branchbeta;
            logprior += a*log(b) - rnd::GetRandom().logGamma(a) + (a-1)*meanlogbl[j] - b*meanbl[j];
            a = lambdastar[j];
            b = mustar[j];
            logpost += a*log(b) - rnd::GetRandom().logGamma(a) + (a-1)*meanlogbl[j] - b*meanbl[j];
        }
    }

    if (varfreerate)    {
        for (int i=0; i<GetNsite(); i++)    {
            double a = alpha;
            double b = alpha;
            logprior += a*log(b) - rnd::GetRandom().logGamma(a) + (a-1)*meanlograte[i] - b*meanrate[i];
            a = alphastar[i];
            b = betastar[i];
            logpost += a*log(b) - rnd::GetRandom().logGamma(a) + (a-1)*meanlograte[i] - b*meanrate[i];
        }
    }

    if (varfreeprofile) {
        for (int i=0; i<GetNsite(); i++)    {

            double priorw = 0;
            double postw = 0;
            for (int k=0; k<GetDim(); k++)  {

                logprior += (dirweight[k] - 1)*meanlogprofile[i][k];
                logprior -= rnd::GetRandom().logGamma(dirweight[k]);
                priorw += dirweight[k];

                logpost += (gammastar[i][k] - 1)*meanlogprofile[i][k];
                logpost -= rnd::GetRandom().logGamma(gammastar[i][k]);
                postw += gammastar[i][k];
            }

            logprior += rnd::GetRandom().logGamma(priorw);
            logpost += rnd::GetRandom().logGamma(postw);
        }
    }

    double logp = logl + logprior - logpost;

    if (varfreeprofile) {
        logp += totmeanlogprofilenorm;
    }
    double tmp = GetRateLengthCorrection();
    logp += tmp;

    return logp;
}

/*
double IIDDirichletIIDGammaPhyloProcess::GetVarLogMarginalLikelihood()  {

    double logprior = 0;
    double logpost = 0;
    double logl = GetLogLikelihood();

    if (varfreebl)  {
        for (int j=1; j<GetNbranch(); j++)    {
            double a = branchalpha;
            double b = branchbeta;
            double x = blarray[j];
            logprior += a*log(b) - rnd::GetRandom().logGamma(a) + (a-1)*log(x) - b*x;
            a = lambdastar[j];
            b = mustar[j];
            logpost += a*log(b) - rnd::GetRandom().logGamma(a) + (a-1)*log(x) - b*x;
        }
    }

    if (varfreerate)    {
        for (int i=0; i<GetNsite(); i++)    {
            double a = alpha;
            double b = alpha;
            double x = rate[i];
            logprior += a*log(b) - rnd::GetRandom().logGamma(a) + (a-1)*log(x) - b*x;
            a = alphastar[i];
            b = betastar[i];
            logpost += a*log(b) - rnd::GetRandom().logGamma(a) + (a-1)*log(x) - b*x;
        }
    }

    if (varfreeprofile) {
        for (int i=0; i<GetNsite(); i++)    {

            double priorw = 0;
            double postw = 0;
            for (int k=0; k<GetDim(); k++)  {

                logprior += (dirweight[k] - 1)*log(profile[i][k]);
                logprior -= rnd::GetRandom().logGamma(dirweight[k]);
                priorw += dirweight[k];

                logpost += (gammastar[i][k] - 1)*log(profile[i][k]);
                logpost -= rnd::GetRandom().logGamma(gammastar[i][k]);
                postw += gammastar[i][k];
            }

            logprior += rnd::GetRandom().logGamma(priorw);
            logpost += rnd::GetRandom().logGamma(postw);
        }
    }

    double logp = logl + logprior - logpost;

    if (varfreeprofile) {
        logp += totmeanlogprofilenorm;
    }
    logp += GetRateLengthCorrection();

    return logp;
}
*/

double IIDDirichletIIDGammaPhyloProcess::GetRateLengthCorrection()    {

    double tot = 0;
    for (int i=0; i<GetNsite(); i++)    {
        for (int j=1; j<GetNbranch(); j++)  {
            tot += rate[i]*blarray[j] - meanrate[i]*meanbl[j];
            /*
            if (varfreerate && varfreebl)   {
                tot += rate[i]*blarray[j] + meanrate[i]*meanbl[j] - meanrate[i]*blarray[j] - rate[i]*meanbl[j];
            }
            else if (!varfreerate || varfreebl)  {
            }
            else    {
                tot += rate[i]*blarray[j] - meanrate[i]*meanbl[j];
            }
            */
        }
    }
    return tot;
}


void IIDDirichletIIDGammaPhyloProcess::UpdateVarLengths() {

	for (int j=1; j<GetNbranch(); j++)	{
        lambdastar[j] = branchalpha + branchlengthsuffstatcount[j];
        mustar[j] = branchbeta + totmeanrate;
    }
    ComputeMeanLengths();

    /*
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
    */
}

void IIDDirichletIIDGammaPhyloProcess::UpdateVarRates() {

    for (int i=0; i<GetNsite(); i++)    {
        alphastar[i] = alpha + siteratesuffstatcount[i];
        betastar[i] = alpha + totmeanlength;
    }
    ComputeMeanRates();

    /*
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
    */
}

void IIDDirichletIIDGammaPhyloProcess::UpdateVarProfiles() {

    for (int i=0; i<GetNsite(); i++)    {
        for (int k=0; k<GetDim(); k++)  {
            gammastar[i][k] = dirweight[k] + siteprofilesuffstatcount[i][k];
        }
    }
    ComputeMeanProfiles();

    /*
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
            double tot = 0;
            for (int k=0; k<GetDim(); k++)  {
                profile[i][k] = exp(rnd::GetRandom().Psi(dirweight[k] + siteprofilesuffstatcount[i][k]));
                tot += dirweight[k] + siteprofilesuffstatcount[i][k];
            }
            double z = exp(rnd::GetRandom().Psi(tot));
            double norm = 0;
            for (int k=0; k<GetDim(); k++)  {
                profile[i][k] /= z;
                norm+= profile[i][k];
            }
            for (int k=0; k<GetDim(); k++)  {
                profile[i][k] /= norm;
            }
            bkrate[i] = norm;
        }
    }
    */
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


