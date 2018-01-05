
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef POISSONSSPROFILE_H
#define POISSONSSPROFILE_H

#include "PoissonProfileProcess.h"
#include "DirichletSiteSpecificProfileProcess.h"

// superclass for Poisson (F81) implementations
class PoissonSiteSpecificProfileProcess: public virtual PoissonProfileProcess, public virtual DirichletSiteSpecificProfileProcess	{

	public:

	PoissonSiteSpecificProfileProcess() : logsum(0) {}
	virtual ~PoissonSiteSpecificProfileProcess() {}

	protected:

	virtual void Create();
	virtual void Delete();

    double Move()   {
        if (GetNprocs() > 1)    {
            return MPIMove();
        }
        return NonMPIMove();
    }

	double MPIMove()    {
		GlobalUpdateSiteProfileSuffStat();
        GlobalResampleSiteProfiles();
        GlobalUpdateProfileHyperSuffStat();
		MoveDirWeights(1,100);
        MoveDirWeights(0.3,100);
        MoveDirWeights(0.1,100);
        return 1.0;
	}

	double NonMPIMove()	{
		UpdateSiteProfileSuffStat();
        ResampleSiteProfiles();
        UpdateProfileHyperSuffStat();
        /*
		MoveDirWeights(1,100);
        MoveDirWeights(0.3,100);
        MoveDirWeights(0.1,100);
        */
		return 1.0;
	}

	double MoveDirWeights(double tuning, int nrep);

    void GlobalResampleSiteProfiles();
	void ResampleSiteProfiles();
    void ResampleSiteProfile(int site);

	virtual double ProfileSuffStatLogProb(int site);

    void GlobalUpdateProfileHyperSuffStat();
    void SlaveUpdateProfileHyperSuffStat();
    void UpdateProfileHyperSuffStat();

    double ProfileHyperSuffStatLogProb();

    /*
    double LogStatIntPrior();
    double LogStatIntPrior(int site);
    */

    void CreateSite(int i) {}
    void DeleteSite(int i) {}

	double GetNormRate(int k)	{

		double tot = 0;
		for (int i=0; i<GetDim(); i++)	{
			for (int j=i+1; j<GetDim(); j++)	{
				tot += profile[k][i] * profile[k][j];
			}
		}
		return 2*tot;
	}

	virtual double GetNormalizationFactor()	{
		double norm = 0;
		for (int k=0; k<GetNsite(); k++)	{
            norm += GetNormRate(k);
		}
		norm /= GetNsite();
		return norm;
	}

    double* logsum;
};

#endif

