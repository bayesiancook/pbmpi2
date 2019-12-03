
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef IIDDIRGAM_H
#define IIDDIRGAM_H

#include "PoissonSubstitutionProcess.h"
#include "PoissonPhyloProcess.h"
#include "GammaRateProcess.h"
#include "PoissonSiteSpecificProfileProcess.h"
#include "GammaBranchProcess.h"

class IIDDirichletIIDGammaSubstitutionProcess : public virtual PoissonSubstitutionProcess, public virtual GammaRateProcess, public virtual PoissonSiteSpecificProfileProcess {

	public:

	using PoissonSubstitutionProcess::UpdateZip;

	IIDDirichletIIDGammaSubstitutionProcess() {}
	virtual ~IIDDirichletIIDGammaSubstitutionProcess() {}

	protected:

	virtual void UpdateZip(int site)	{
		PoissonSubstitutionProcess::UpdateZip(site);
	}

	virtual void Create()	{
		PoissonSubstitutionProcess::Create();
		GammaRateProcess::Create();
		PoissonSiteSpecificProfileProcess::Create();
	}

	virtual void Delete()	{
		PoissonSiteSpecificProfileProcess::Delete();
		GammaRateProcess::Delete();
		PoissonSubstitutionProcess::Delete();
	}

};

class IIDDirichletIIDGammaPhyloProcess : public virtual PoissonPhyloProcess, public virtual IIDDirichletIIDGammaSubstitutionProcess, public virtual GammaBranchProcess	{

	public:

    virtual void SlaveExecute(MESSAGE);
	virtual void GlobalUpdateParameters();
	virtual void SlaveUpdateParameters();

	IIDDirichletIIDGammaPhyloProcess() {}

	IIDDirichletIIDGammaPhyloProcess(istream& is, int inmyid, int innprocs)	{

		// generic
		FromStreamHeader(is);
		SetMPI(inmyid,innprocs);
		Open(is);
	}

	void ToStreamHeader(ostream& os)	{
		PhyloProcess::ToStreamHeader(os);
	}

	~IIDDirichletIIDGammaPhyloProcess() {
		Delete();
	}

	void TraceHeader(ostream& os)	{
		os << "#iter\ttime\ttopo\tloglik\tlength\trate\talpha";
		os << "\tstatent\tstatalpha";
		os << '\n'; 
	}

	void Trace(ostream& os)	{

		os << GetSize();
		if (chronototal.GetTime())	{
			os << '\t' << chronototal.GetTime() / 1000;
			os << '\t' << ((int) (propchrono.GetTime() / chronototal.GetTime() * 100));
			chronototal.Reset();
			propchrono.Reset();
		}
		else	{
			os << '\t' << 0;
			os << '\t' << 0;
		}

		os << '\t' << GetLogLikelihood() << '\t' << GetTotalLength() << '\t' << GetMeanRate() << '\t' << GetAlpha();
		// os << '\t' << GetLogLikelihood() << '\t' << GetRenormTotalLength() << '\t' << GetMeanRate() << '\t' << GetAlpha();
		os << '\t' << GetStatEnt();
		os << '\t' << GetMeanDirWeight();
		os << '\n';
	}

	virtual double Move(double tuning = 1.0)	{

		chronototal.Start();
		propchrono.Start();
		// BranchLengthMove(tuning);
		// BranchLengthMove(0.1 * tuning);
		if (! fixtopo)	{
			MoveTopo();
		}
		propchrono.Stop();

		for (int rep=0; rep<5; rep++)	{
            GlobalCollapse();

            GammaBranchProcess::Move(tuning,10);
            GlobalUpdateParameters();

            GammaRateProcess::Move();
            GlobalUpdateParameters();

            PoissonSiteSpecificProfileProcess::Move();
            GlobalUpdateParameters();

            GlobalUnfold();
        }

        GlobalCollectSiteRates();
        GlobalCollectSiteProfiles();

        chronototal.Stop();
		return 1;
	
	}

    virtual void VarBayes();
    double GetVarLogMarginalLikelihood(double* sitescore = 0);
    void UpdateVarLengths();
    void UpdateVarRates();
    void UpdateVarProfiles();

    void PosteriorMean(int burnin, int nrep);

    void UpdateSite(int i)  {
        UpdateZip(i);
    }

	// virtual void ReadPB(int argc, char* argv[]);

	void ToStream(ostream& os)	{
		GammaBranchProcess::ToStream(os);
		GammaRateProcess::ToStream(os);
		PoissonSiteSpecificProfileProcess::ToStream(os);
	}

	void FromStream(istream& is)	{
		GammaBranchProcess::FromStream(is);
		GammaRateProcess::FromStream(is);
		PoissonSiteSpecificProfileProcess::FromStream(is);
		GlobalUpdateParameters();
	}


	virtual void Create()	{
		PoissonPhyloProcess::Create();
		IIDDirichletIIDGammaSubstitutionProcess::Create();
		GammaBranchProcess::Create();
        CreateArrays();
    }

    void CreateArrays() {
        meanrate = new double[GetNsite()];
        meanlograte = new double[GetNsite()];
        alphastar = new double[GetNsite()];
        betastar = new double[GetNsite()];

        meanbl = new double[GetNbranch()];
        meanlogbl = new double[GetNbranch()];
        lambdastar = new double[GetNbranch()];
        mustar = new double[GetNbranch()];

        meanlogprofile = new double*[GetNsite()];
        for (int i=0; i<GetNsite(); i++)    {
            meanlogprofile[i] = new double[GetDim()];
        }
        meanlogprofilenorm = new double[GetNsite()];
        gammastar = new double*[GetNsite()];
        for (int i=0; i<GetNsite(); i++)    {
            gammastar[i] = new double[GetDim()];
        }

	}
		
    void DeleteArrays() {

    }

	virtual void Delete()	{
		GammaBranchProcess::Delete();
		IIDDirichletIIDGammaSubstitutionProcess::Delete();
		PoissonPhyloProcess::Delete();
	}

    void CreateParamArrays(int nrep)    {

        nparam = nrep;
        bllist = new double*[nparam];
        ratelist = new double*[nparam];
        profilelist = new double**[nparam];
        for (int rep=0; rep<nparam; rep++)  {
            bllist[rep] = new double[GetNbranch()];
            ratelist[rep] = new double[GetNsite()];
            profilelist[rep] = new double*[GetNsite()];
            for (int i=0; i<GetNsite(); i++)    {
                profilelist[rep][i] = new double[GetDim()];
            }
        }
    }

    void PushParam(int rep)   {
        for (int j=1; j<GetNbranch(); j++)  {
            bllist[rep][j] = blarray[j];
        }
        for (int i=0; i<GetNsite(); i++)    {
            ratelist[rep][i] = rate[i];
        }
        for (int i=0; i<GetNsite(); i++)    {
            for (int k=0; k<GetDim(); k++)  {
                profilelist[rep][i][k] = profile[i][k];
            }
        }
    }

    void PopParam(int rep)   {
        for (int j=1; j<GetNbranch(); j++)  {
            blarray[j] = bllist[rep][j];
        }
        for (int i=0; i<GetNsite(); i++)    {
            rate[i] = ratelist[rep][i];
        }
        for (int i=0; i<GetNsite(); i++)    {
            for (int k=0; k<GetDim(); k++)  {
                profile[i][k] = profilelist[rep][i][k];
            }
        }

        totmeanrate = 0;
        for (int i=0; i<GetNsite(); i++)    {
            totmeanrate += rate[i];
        }
        totmeanlength = 0;
        for (int j=1; j<GetNbranch(); j++)  {
            totmeanlength += blarray[j];
        }

        UpdateZip();
        UpdateConditionalLikelihoods();
    }

    void CreateBKArrays(int nrep)   {


        alphalist = new double*[nrep];
        betalist = new double*[nrep];
        lambdalist = new double*[nrep];
        mulist = new double*[nrep];
        gammalist = new double**[nrep];

        scorelist = new double[nrep];
        sitescorelist = new double*[nrep];

        for (int rep=0; rep<nrep; rep++)    {
            lambdalist[rep] = new double[GetNbranch()];
            mulist[rep] = new double[GetNbranch()];

            alphalist[rep] = new double[GetNsite()];
            betalist[rep] = new double[GetNsite()];
            gammalist[rep] = new double*[GetNsite()];
            for (int i=0; i<GetNsite(); i++)    {
                gammalist[rep][i] = new double[GetDim()];
            }
            sitescorelist[rep] = new double[GetNsite()];
        }

        nstore = nrep;
        bkindex = 0;
    }

    void Push(double score)    {
        int index = bkindex;
        if (bkindex == nstore)  {
            double min = 0;
            index = 0;
            for (int rep=0; rep<nstore; rep++)  {
                if ((!rep) || (min > scorelist[rep]))   {
                    min = scorelist[rep];
                    index = rep;
                }
            }
            if (score < min)    {
                index = -1;
            }
        }

        if (index != -1)    {
            scorelist[index] = GetVarLogMarginalLikelihood(sitescorelist[index]);

            for (int i=1; i<GetNbranch(); i++)  {
                lambdalist[index][i] = lambdastar[i];
                mulist[index][i] = mustar[i];
            }

            for (int i=0; i<GetNsite(); i++)    {
                alphalist[index][i] = alphastar[i];
                betalist[index][i] = betastar[i];

                for (int k=0; k<GetDim(); k++)  {
                    gammalist[index][i][k] = gammastar[i][k];
                }
            }
        }

        if (bkindex < nstore)   {
            bkindex++;
        }
    }

    double Pop(int index, int bestpersite)  {

        if (index == -1)    {
            index = 0;
            double max = scorelist[0];
            for (int rep=1; rep<nstore; rep++)  {
                if (scorelist[rep] > max)   {
                    max = scorelist[rep];
                    index = rep;
                }
            }
        }

        for (int i=1; i<GetNbranch(); i++)  {
            lambdastar[i] = lambdalist[index][i];
            mustar[i] = mulist[index][i];
        }

        if (bestpersite)    {
            for (int i=0; i<GetNsite(); i++)    {
                int lmax = 0;
                int max = 0;
                for (int l=0; l<nstore; l++)    {
                    if ((!l) || (max < sitescorelist[l][i]))    {
                        max = sitescorelist[l][i];
                        lmax = l;
                    }
                }
                alphastar[i] = alphalist[lmax][i];
                betastar[i] = betalist[lmax][i];
                for (int k=0; k<GetDim(); k++)  {
                    gammastar[i][k] = gammalist[lmax][i][k];
                }
            }
        }
        else    {
            for (int i=0; i<GetNsite(); i++)    {
                alphastar[i] = alphalist[index][i];
                betastar[i] = betalist[index][i];
                for (int k=0; k<GetDim(); k++)  {
                    gammastar[i][k] = gammalist[index][i][k];
                }
            }
        }

        ComputeMeanLengths();
        ComputeMeanRates();
        ComputeMeanProfiles();

        SetNewParameters();

        UpdateConditionalLikelihoods();
        GetVarLogMarginalLikelihood();
        return scorelist[index];
    }

    void ComputeMeanRates();
    void ComputeMeanLengths();
    void ComputeMeanProfiles();

    void SetNewParameters();

    double GetRateLengthCorrection(double* sitescore = 0);

    double VBEM(int nrep, double diff,double* sitescore = 0);
    void InitializeState();

    double* meanrate;
    double* meanlograte;
    double totmeanrate;
    double* meanbl;
    double* meanlogbl;
    double totmeanlength;
    double** meanlogprofile;
    double* meanlogprofilenorm;
    double totmeanlogprofilenorm;

    double* alphastar;
    double* betastar;
    double* lambdastar;
    double* mustar;
    double** gammastar;

    int varfreebl;
    int varfreerate;
    int varfreeprofile;

    double** alphalist;
    double** betalist;
    double** lambdalist;
    double** mulist;
    double*** gammalist;
    double** sitescorelist;
    double* scorelist;
    int bkindex;
    int nstore;

    int nparam;
    double** bllist;
    double** ratelist;
    double*** profilelist;
};

#endif

