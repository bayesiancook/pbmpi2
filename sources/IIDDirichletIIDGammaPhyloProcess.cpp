
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
    for (int i=GetSiteMin(); i<GetSiteMax(); i++)    {
        if (ActiveSite(i))  {
            meanrate[i] = alphastar[i] / betastar[i];
            totmeanrate += meanrate[i];
            meanlograte[i] = rnd::GetRandom().Psi(alphastar[i]) - log(betastar[i]);
        }
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
    for (int i=GetSiteMin(); i<GetSiteMax(); i++)    {
        if (ActiveSite(i))  {
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
}

void IIDDirichletIIDGammaPhyloProcess::PosteriorMean(int burnin, int nrep) {

    for (int i=GetSiteMin(); i<GetSiteMax(); i++)    {
        meanlograte[i] = 0;
        meanrate[i] = 0;
    }
    for (int j=0; j<GetNbranch(); j++)  {
        meanlogbl[j] = 0;
        meanbl[j] = 0;
    }
    for (int i=GetSiteMin(); i<GetSiteMax(); i++)    {
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
        // PushParam(rep);
        for (int i=GetSiteMin(); i<GetSiteMax(); i++)    {
            meanlograte[i] += log(rate[i]);
            meanrate[i] += rate[i];
        }
        for (int j=0; j<GetNbranch(); j++)  {
            meanlogbl[j] += log(blarray[j]);
            meanbl[j] += blarray[j];
        }
        for (int i=GetSiteMin(); i<GetSiteMax(); i++)    {
            for (int k=0; k<GetDim(); k++)  {
                meanlogprofile[i][k] += log(profile[i][k]);
            }
        }
        Trace(cerr);
    }

    cerr << '\n';

    totmeanrate = 0;
    for (int i=GetSiteMin(); i<GetSiteMax(); i++)    {
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
    for (int i=GetSiteMin(); i<GetSiteMax(); i++)    {
        double tot = 0;
        for (int k=0; k<GetDim(); k++)  {
            meanlogprofile[i][k] /= nrep;
            tot += exp(meanlogprofile[i][k]);
        }
        meanlogprofilenorm[i] = log(tot);
        totmeanlogprofilenorm += log(tot);
    }
}

void IIDDirichletIIDGammaPhyloProcess::InitializeState()    {

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
            for (int i=GetSiteMin(); i<GetSiteMax(); i++)    {
                rate[i] = 1.0;
            }
        }

        for (int i=GetSiteMin(); i<GetSiteMax(); i++)    {
            meanrate[i] = rate[i];
        }

        if (varfreeprofile) {
            SampleStat();
        }
        else    {
            for (int i=GetSiteMin(); i<GetSiteMax(); i++)    {
                for (int k=0; k<GetDim(); k++)  {
                    profile[i][k] = 1.0/GetDim();
                }
            }
        }
        UpdateZip();
}

void IIDDirichletIIDGammaPhyloProcess::VarBayes()    {


    cerr << name << '\n';

    varfreebl = 1;
    varfreerate = 1;
    varfreeprofile = 1;

    int nrep = 10;
    int nsub = 3;
    ofstream os((name + ".finalscore").c_str());
    os.precision(12);
    CreateBKArrays(nsub);

    // CreateParamArrays(nrep);
    // PosteriorMean(10,nrep);
    
    double* bksitescore = new double[GetNsite()];
    double* sitescore = new double[GetNsite()];
    double* sitevar = new double[GetNsite()];
    double* sitemean = new double[GetNsite()];
    for (int i=0; i<GetNsite(); i++)    {
        sitemean[i] = 0;
        sitevar[i] = 0;
    }

    for (int rep=0; rep<nrep; rep++)  {

        // PopParam(rep);
        InitializeState();

        PosteriorMean(20,20);
        SetNewParameters();

        totmeanrate = 0;
        for (int i=0; i<GetNsite(); i++)    {
            totmeanrate += meanrate[i];
        }
        totmeanlength = 0;
        for (int j=1; j<GetNbranch(); j++)  {
            totmeanlength += meanbl[j];
        }

        double score = VBEM(0,1e-4*GetNsite());
        os << score << '\t';
        Trace(os);
        os.flush();

        cerr << score << '\t';
        Trace(cerr);
        Push(score);

        GetVarLogMarginalLikelihood(sitescore);
        for (int i=0; i<GetNsite(); i++)    {
            sitemean[i] += sitescore[i];
            sitevar[i] += sitescore[i]*sitescore[i];
        }
    }

    cerr << '\n';

    double varcutoff = 0.01;

    double c1 = 0;
    double c01 = 0;
    double c001 = 0;
    double c0001 = 0;
    int index = 0;
    int* permut = new int[GetNsite()];
    for (int i=0; i<GetNsite(); i++)    {
        permut[i] = -1;
    }
    for (int i=0; i<GetNsite(); i++)    {
        sitemean[i] /= nrep;
        sitevar[i] /= nrep;
        sitevar[i] -= sitemean[i]*sitemean[i];
        if (sitevar[i] > 1.0)   {
            c1++;
        }
        if (sitevar[i] > 0.1)   {
            c01++;
        }
        if (sitevar[i] > 0.01)   {
            permut[i] = index++;
            c001++;
        }
        if (sitevar[i] > 0.001)   {
            c0001++;
        }
    }
    for (int i=0; i<GetNsite(); i++)    {
        if (permut[i] == -1)    {
            permut[i] = index++;
        }
    }

    NonMPIReshuffleSites(permut);

    os << " > 1.000 : " << c1 / GetNsite() << '\n';
    os << " > 0.100 : " << c01 / GetNsite() << '\n';
    os << " > 0.010 : " << c001 / GetNsite() << '\n';
    os << " > 0.001 : " << c0001 / GetNsite() << '\n';

    for (int sub=0; sub<nsub; sub++)    {

        double scorebefore = Pop(sub,0);
        os << scorebefore << '\t';
        Trace(os);
        os.flush();
        double score = VBEM(0,1e-6*GetNsite());
        os << score << '\t';
        Trace(os);
        os.flush();

        double* bkalphastar = new double[GetNsite()];
        double* bkbetastar = new double[GetNsite()];
        double** bkgammastar = new double*[GetNsite()];
        for (int i=0; i<GetNsite(); i++)    {
            bkgammastar[i] = new double[GetDim()];
        }

        for (int r=0; r<10; r++)    {

            varfreebl = 0;

            SetMinMax(0,0.1);

            UpdateConditionalLikelihoods();
            // double currentscore = GetVarLogMarginalLikelihood();
            double currentscore = VBEM(0,1e-6*GetNsite(),bksitescore);

            for (int l=0; l<10; l++)    {

                for (int i=0; i<GetNsite(); i++)    {
                    if (ActiveSite(i))  {
                    // if (sitevar[i] > varcutoff) {
                        bkalphastar[i] = alphastar[i];
                        bkbetastar[i] = betastar[i];
                        for (int k=0; k<GetDim(); k++)  {
                            bkgammastar[i][k] = gammastar[i][k];
                        }

                        double tun = 4;
                        alphastar[i] *= rnd::GetRandom().Gamma(tun,tun);
                        betastar[i] *= rnd::GetRandom().Gamma(tun,tun);
                        for (int k=0; k<GetDim(); k++)  {
                            gammastar[i][k] *= rnd::GetRandom().Gamma(tun,tun);
                        }
                    }
                }

                ComputeMeanRates();
                ComputeMeanProfiles();
                SetNewParameters();

                UpdateConditionalLikelihoods();

                double score = VBEM(0,1e-6*GetNsite(),bksitescore);

                for (int i=0; i<GetNsite(); i++)    {
                    if (ActiveSite(i))  {
                    // if (sitevar[i] > varcutoff) {
                        if (sitescore[i] < bksitescore[i])   {
                            alphastar[i] = bkalphastar[i];
                            betastar[i] = bkbetastar[i];
                            for (int k=0; k<GetDim(); k++)  {
                                gammastar[i][k] = bkgammastar[i][k];
                            }
                            ComputeMeanLengths();
                            ComputeMeanRates();
                            ComputeMeanProfiles();

                            SetNewParameters();

                            UpdateConditionalLikelihoods();
                            GetVarLogMarginalLikelihood();

                            double rescore = VBEM(0,1e-4);
                        }
                        else    {
                            bksitescore[i] = sitescore[i];
                        }
                    }
                }
            }

            varfreebl = 1;
            SetMinMax(0,1);

            UpdateConditionalLikelihoods();
            double score = VBEM(0,1e-6*GetNsite());

            os << score << '\t';
            Trace(os);
            os.flush();
        }
        os << '\n';
        os.flush();
    }
}

/*
void IIDDirichletIIDGammaPhyloProcess::VarBayes()    {


    cerr << name << '\n';

    varfreebl = 1;
    varfreerate = 1;
    varfreeprofile = 1;

    int nrep = 10;
    int nsub = 3;
    ofstream os((name + ".finalscore").c_str());
    os.precision(12);
    CreateBKArrays(nsub);

    // CreateParamArrays(nrep);
    // PosteriorMean(10,nrep);
    
    double* sitescore = new double[GetNsite()];
    double* sitevar = new double[GetNsite()];
    double* sitemean = new double[GetNsite()];
    for (int i=0; i<GetNsite(); i++)    {
        sitemean[i] = 0;
        sitevar[i] = 0;
    }

    for (int rep=0; rep<nrep; rep++)  {

        // PopParam(rep);
        InitializeState();

        PosteriorMean(20,20);
        SetNewParameters();

        totmeanrate = 0;
        for (int i=0; i<GetNsite(); i++)    {
            totmeanrate += meanrate[i];
        }
        totmeanlength = 0;
        for (int j=1; j<GetNbranch(); j++)  {
            totmeanlength += meanbl[j];
        }

        double score = VBEM(0,1e-4*GetNsite());
        os << score << '\t';
        Trace(os);
        os.flush();

        cerr << score << '\t';
        Trace(cerr);
        Push(score);

        GetVarLogMarginalLikelihood(sitescore);
        for (int i=0; i<GetNsite(); i++)    {
            sitemean[i] += sitescore[i];
            sitevar[i] += sitescore[i]*sitescore[i];
        }
    }

    cerr << '\n';

    double varcutoff = 0.01;

    double c1 = 0;
    double c01 = 0;
    double c001 = 0;
    double c0001 = 0;
    for (int i=0; i<GetNsite(); i++)    {
        sitemean[i] /= nrep;
        sitevar[i] /= nrep;
        sitevar[i] -= sitemean[i]*sitemean[i];
        if (sitevar[i] > 1.0)   {
            c1++;
        }
        if (sitevar[i] > 0.1)   {
            c01++;
        }
        if (sitevar[i] > 0.01)   {
            c001++;
        }
        if (sitevar[i] > 0.001)   {
            c0001++;
        }
    }

    os << " > 1.000 : " << c1 / GetNsite() << '\n';
    os << " > 0.100 : " << c01 / GetNsite() << '\n';
    os << " > 0.010 : " << c001 / GetNsite() << '\n';
    os << " > 0.001 : " << c0001 / GetNsite() << '\n';

    for (int sub=0; sub<nsub; sub++)    {

        double scorebefore = Pop(sub,0);
        os << scorebefore << '\t';
        Trace(os);
        os.flush();
        double score = VBEM(0,1e-6*GetNsite());
        os << score << '\t';
        Trace(os);
        os.flush();

        for (int r=0; r<10; r++)    {
            int smin = sitemin[0];
            int smax = sitemax[0];

            varfreebl = 0;

            for (int i=0; i<GetNsite(); i++)    {

                if (sitevar[i] > varcutoff) {

                    sitemin[0] = i;
                    sitemax[0] = i+1;

                    UpdateConditionalLikelihoods();
                    // double currentscore = GetVarLogMarginalLikelihood();
                    double currentscore = VBEM(0,1e-4);

                    for (int l=0; l<10; l++)    {

                        double bkalphastar = alphastar[i];
                        double bkbetastar = betastar[i];
                        double bkgammastar[GetDim()];
                        for (int k=0; k<GetDim(); k++)  {
                            bkgammastar[k] = gammastar[i][k];
                        }

                        double tun = 4;
                        alphastar[i] *= rnd::GetRandom().Gamma(tun,tun);
                        betastar[i] *= rnd::GetRandom().Gamma(tun,tun);
                        for (int k=0; k<GetDim(); k++)  {
                            gammastar[i][k] *= rnd::GetRandom().Gamma(tun,tun);
                        }

                        ComputeMeanRates();
                        ComputeMeanProfiles();
                        SetNewParameters();

                        UpdateConditionalLikelihoods();

                        double score = VBEM(0,1e-4);
                        if (score < currentscore)   {
                            alphastar[i] = bkalphastar;
                            betastar[i] = bkbetastar;
                            for (int k=0; k<GetDim(); k++)  {
                                gammastar[i][k] = bkgammastar[k];
                            }
                            ComputeMeanLengths();
                            ComputeMeanRates();
                            ComputeMeanProfiles();

                            SetNewParameters();

                            UpdateConditionalLikelihoods();
                            GetVarLogMarginalLikelihood();

                            double rescore = VBEM(0,1e-4);
                        }
                    }
                }
            }

            sitemin[0] = 0;
            sitemax[0] = GetNsite();
            varfreebl = 1;

            UpdateConditionalLikelihoods();
            double score = VBEM(0,1e-6*GetNsite());

            os << score << '\t';
            Trace(os);
            os.flush();
        }
        os << '\n';
        os.flush();
    }
}
*/

double IIDDirichletIIDGammaPhyloProcess::VBEM(int nrep, double diffmin, double* sitescore)  {

    double logp = 0;
    int rep = 0;

    double diff = 1000.0;
    while ((diff > diffmin) && ((nrep == 0) || (rep<nrep))) {

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
        double tmp = GetVarLogMarginalLikelihood(sitescore);
        diff = tmp - logp;
        if (! logp) {
            diff= 1;
        }
        logp = tmp;
        if (GetSiteMax() - GetSiteMin() == GetNsite())    {
            cerr << diff << '\t' << logL << '\t' << (logL-logp)/GetNsite() << '\t' << logp << '\t';
            Trace(cerr);
        }
    }
    if (GetSiteMax() - GetSiteMin() == GetNsite())    {
        cerr << '\n';
    }
    return logp;
}

void IIDDirichletIIDGammaPhyloProcess::SetNewParameters()   {

    if (varfreebl)  {
        for (int j=1; j<GetNbranch(); j++)  {
            blarray[j] = exp(meanlogbl[j]);
        }
    }

    if (varfreerate && varfreeprofile)    {
        for (int i=GetSiteMin(); i<GetSiteMax(); i++)    {
            if (ActiveSite(i))  {
                rate[i] = exp(meanlograte[i] + meanlogprofilenorm[i]);
            }
        }
    }
    else if (varfreerate)   {
        for (int i=GetSiteMin(); i<GetSiteMax(); i++)    {
            if (ActiveSite(i))  {
                rate[i] = exp(meanlograte[i]);
            }
        }
    }
    else if (varfreeprofile)    {
        for (int i=GetSiteMin(); i<GetSiteMax(); i++)    {
            if (ActiveSite(i))  {
                rate[i] = exp(meanlogprofilenorm[i]);
            }
        }
    }

    if (varfreeprofile) {
        for (int i=GetSiteMin(); i<GetSiteMax(); i++)    {
            if (ActiveSite(i))  {
                double tot = 0;
                for (int k=0; k<GetDim(); k++)  {
                    profile[i][k] = exp(meanlogprofile[i][k] - meanlogprofilenorm[i]);
                }
            }
        }
    }
    UpdateZip();
}

double IIDDirichletIIDGammaPhyloProcess::GetVarLogMarginalLikelihood(double* sitescore)  {

    double logprior = 0;
    double logpost = 0;

    if (sitescore)  {
        double check = 0;
        for (int i=GetSiteMin(); i<GetSiteMax(); i++)    {
            if (ActiveSite(i))  {
                check += sitelogL[i];
                sitescore[i] = sitelogL[i];
            }
        }
        if (fabs(check - logL) > 1e-7)    {
            cerr << "error in get var log marginal likelihood: site logl not correctly updated\n";
            exit(1);
        }
    }

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
        for (int i=GetSiteMin(); i<GetSiteMax(); i++)    {
            if (ActiveSite(i))  {
                double a = alpha;
                double b = alpha;
                logprior += a*log(b) - rnd::GetRandom().logGamma(a) + (a-1)*meanlograte[i] - b*meanrate[i];
                if (sitescore)  {
                    sitescore[i] += a*log(b) - rnd::GetRandom().logGamma(a) + (a-1)*meanlograte[i] - b*meanrate[i];
                }
                a = alphastar[i];
                b = betastar[i];
                logpost += a*log(b) - rnd::GetRandom().logGamma(a) + (a-1)*meanlograte[i] - b*meanrate[i];
                if (sitescore)  {
                    sitescore[i] -= a*log(b) - rnd::GetRandom().logGamma(a) + (a-1)*meanlograte[i] - b*meanrate[i];
                }
            }
        }
    }

    if (varfreeprofile) {
        for (int i=GetSiteMin(); i<GetSiteMax(); i++)    {

            if (ActiveSite(i))  {
                double priorw = 0;
                double postw = 0;
                for (int k=0; k<GetDim(); k++)  {

                    logprior += (dirweight[0][k] - 1)*meanlogprofile[i][k];
                    logprior -= rnd::GetRandom().logGamma(dirweight[0][k]);
                    priorw += dirweight[0][k];

                    if (sitescore)  {
                        sitescore[i] += (dirweight[0][k] - 1)*meanlogprofile[i][k];
                        sitescore[i] -= rnd::GetRandom().logGamma(dirweight[0][k]);
                    }

                    logpost += (gammastar[i][k] - 1)*meanlogprofile[i][k];
                    logpost -= rnd::GetRandom().logGamma(gammastar[i][k]);
                    postw += gammastar[i][k];

                    if (sitescore)  {
                        sitescore[i] -= (gammastar[i][k] - 1)*meanlogprofile[i][k];
                        sitescore[i] += rnd::GetRandom().logGamma(gammastar[i][k]);
                    }
                }

                logprior += rnd::GetRandom().logGamma(priorw);
                logpost += rnd::GetRandom().logGamma(postw);

                if (sitescore)  {
                    sitescore[i] += rnd::GetRandom().logGamma(priorw);
                    sitescore[i] -= rnd::GetRandom().logGamma(postw);
                }
            }
        }
    }

    double logp = logL + logprior - logpost;

    if (varfreeprofile) {
        logp += totmeanlogprofilenorm;
        if (sitescore)  {
            for (int i=GetSiteMin(); i<GetSiteMax(); i++)    {
                if (ActiveSite(i))  {
                    sitescore[i] += meanlogprofilenorm[i];
                }
            }
        }
    }
    double tmp = GetRateLengthCorrection();
    logp += tmp;

    return logp;
}

/*
double IIDDirichletIIDGammaPhyloProcess::GetVarLogMarginalLikelihood()  {

    double logprior = 0;
    double logpost = 0;

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

                logprior += (dirweight[0][k] - 1)*log(profile[i][k]);
                logprior -= rnd::GetRandom().logGamma(dirweight[0][k]);
                priorw += dirweight[0][k];

                logpost += (gammastar[i][k] - 1)*log(profile[i][k]);
                logpost -= rnd::GetRandom().logGamma(gammastar[i][k]);
                postw += gammastar[i][k];
            }

            logprior += rnd::GetRandom().logGamma(priorw);
            logpost += rnd::GetRandom().logGamma(postw);
        }
    }

    double logp = logL + logprior - logpost;

    if (varfreeprofile) {
        logp += totmeanlogprofilenorm;
    }
    logp += GetRateLengthCorrection();

    return logp;
}
*/

double IIDDirichletIIDGammaPhyloProcess::GetRateLengthCorrection(double* sitescore)    {

    double tot = 0;
    for (int i=GetSiteMin(); i<GetSiteMax(); i++)    {
        if (ActiveSite(i))  {
            double totsite = 0;
            for (int j=1; j<GetNbranch(); j++)  {
                totsite += rate[i]*blarray[j] - meanrate[i]*meanbl[j];
            }
            tot += totsite;
            if (sitescore)  {
                sitescore[i] += totsite;
            }
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

    for (int i=GetSiteMin(); i<GetSiteMax(); i++)    {
        if (ActiveSite(i))  {
            alphastar[i] = alpha + siteratesuffstatcount[i];
            betastar[i] = alpha + totmeanlength;
        }
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

    for (int i=GetSiteMin(); i<GetSiteMax(); i++)    {
        if (ActiveSite(i))  {
            for (int k=0; k<GetDim(); k++)  {
                gammastar[i][k] = dirweight[0][k] + siteprofilesuffstatcount[i][k];
            }
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
                profile[i][k] = dirweight[0][k] + siteprofilesuffstatcount[i][k];
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
                profile[i][k] = exp(rnd::GetRandom().Psi(dirweight[0][k] + siteprofilesuffstatcount[i][k]));
                tot += dirweight[0][k] + siteprofilesuffstatcount[i][k];
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
            dvector[index] = dirweight[0][i];
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
		dirweight[0][i] = dvector[index];
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


