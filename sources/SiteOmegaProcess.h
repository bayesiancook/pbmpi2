/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef SITEOMEGA_H
#define SITEOMEGA_H

#include "OmegaProcess.h"

class SiteOmegaProcess : public virtual OmegaProcess	{

	public:

	SiteOmegaProcess() : omegaarray(0) {}
	virtual ~SiteOmegaProcess() {}

	protected:

	double GetSiteOmega(int site)	{
		return omegaarray[site];
	}

	double* GetSiteOmegaPtr(int site)	{
		return &(omegaarray[site]);
	}
	
	double GetMeanOmega()	{
		double mean = 0;
		for (int i=0; i<GetNsite(); i++)	{
			mean += omegaarray[i];
		}
		mean /= GetNsite();
		return mean;
	}

	double GetVarOmega()	{
		double mean = 0;
		double var = 0;
		for (int i=0; i<GetNsite(); i++)	{
			mean += omegaarray[i];
			mean += omegaarray[i] * omegaarray[i];
		}
		mean /= GetNsite();
		var /= GetNsite();
		var -= mean*mean;
		return var;
	}

	double OmegaSuffStatLogProb(int site)	{
		return siteomegasuffstatcount[site] * log(omegaarray[site]) - siteomegasuffstatbeta[site] * omegaarray[site];
	}

	virtual double OmegaSuffStatLogProb()	{
		double total = 0;
		for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
			total += OmegaSuffStatLogProb(i);	
		}
		return total;
	}

	// omega
	double LogOmegaPrior(int site)	{
		return omegaalpha*log(omegabeta) - rnd::GetRandom().logGamma(omegaalpha) + (omegaalpha-1)*log(omegaarray[site]) - omegabeta*omegaarray[site];
	}

	virtual double LogOmegaPrior()	{
		double total = 0;
		for (int i=0; i<GetNsite(); i++)	{
			total += LogOmegaPrior(i);
		}
		return total;
	}

	double LogOmegaHyperPrior()	{
		return -omegaalpha - omegabeta;
	}

	virtual void SampleOmega()	{
		omegaalpha = omegabeta = 1.0;
		for (int i=0; i<GetNsite(); i++)	{
			omegaarray[i] = rnd::GetRandom().Gamma(omegaalpha,omegabeta);
		}
	}

	double MoveOmega(double tuning) {}

	protected:

	virtual void Create()	{
		if (! omegaarray)	{
			OmegaProcess::Create();
			omegaarray = new double[GetNsite()];
		}
	}

	virtual void Delete()	{
		if (omegaarray)	{
			delete omegaarray;
			omegaarray = 0;
			OmegaProcess::Delete();
		}
	}
	
	void UpdateOmegaSuffStat() {}
	void GlobalUpdateOmegaSuffStat() {}
	void SlaveUpdateOmegaSuffStat() {}
	//void UpdateSiteOmegaSuffStat() {}

	double* omegaarray;
	double omegaalpha;
	double omegabeta;
};

#endif

