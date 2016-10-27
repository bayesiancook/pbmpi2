
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef OMEGA_H
#define OMEGA_H

#include <cstdlib>
#include <iostream>
#include <cmath>
using namespace std;

#include "MPIModule.h"
#include "Random.h"

class OmegaProcess : public virtual MPIModule	{

	public:

	OmegaProcess() : fixomega(0), omegaprior(0), siteomegasuffstatcount(0) {}
	virtual ~OmegaProcess() {}

	double GetMeanOmega()	{
		double tot = 0;
		for (int i=0; i<GetNsite(); i++)	{
			tot += GetSiteOmega(i);
		}
		return tot/GetNsite();
	}

	double GetMinOmega()	{
		double min = 10;
		for (int i=0; i<GetNsite(); i++)	{
			if (min > GetSiteOmega(i))	{
				min = GetSiteOmega(i);
			}
		}
		return min;
	}
	
	double GetRelVarOmega()	{
		double mean = 0;
		double var = 0;
		for (int i=0; i<GetNsite(); i++)	{
			double tmp = GetSiteOmega(i);
			mean += tmp;
			var += tmp * tmp;
		}
		mean /= GetNsite();
		var /= GetNsite();
		var -= mean*mean;
		var /= mean*mean;
		return var;
	}

	double GetProportionOmegaGreaterThan(double c)	{
		double tot = 0;
		for (int i=0; i<GetNsite(); i++)	{
			if (GetSiteOmega(i) > c)	{
				tot++;
			}
		}
		tot /= GetNsite();
		return tot;
	}

	protected:

	// omega

	virtual void UpdateOmegaSuffStat()	{
		cerr << "in OmegaProcess::UpdateOmegaSuffStat\n";
		exit(1);
	}

	virtual void GlobalUpdateOmegaSuffStat()	{
		cerr << "in OmegaProcess::GlobalUpdateOmegaSuffStat\n";
		exit(1);
	}

	virtual void SlaveUpdateOmegaSuffStat()	{
		cerr << "in OmegaProcess::SlaveUpdateOmegaSuffStat\n";
		exit(1);
	}	

	virtual double LogOmegaPrior()	{
		cerr << "in OmegaProcess::LogOmegaPrior\n";
		exit(1);
	}

	virtual double OmegaSuffStatLogProb()	{
		cerr << "in OmegaProcess::OmegaSuffStatLogProb\n";
		exit(1);
	}

	virtual double CountOmegaSuffStatLogProb()	{
		cerr << "in OmegaProcess::OmegaSuffStatLogProb\n";
		exit(1);
	}

	virtual double BetaOmegaSuffStatLogProb()	{
		cerr << "in OmegaProcess::OmegaSuffStatLogProb\n";
		exit(1);
	}

	virtual void SampleOmega()	{
		cerr << "in OmegaProcess::SampleOmega\n";
		exit(1);
	}

	virtual double MoveOmega(double tuning)	{
		cerr << "in OmegaProcess::MoveOmega\n";
		exit(1);
	}

	virtual void UpdateOmega()	{
		cerr << "in OmegaProcess::UpdateOmega\n";
		exit(1);
	}

	virtual double GetSiteOmega(int site)	{
		cerr << "in OmegaProcess::GetSiteOmega\n";
		exit(1);
	}

	virtual double* GetSiteOmegaPtr(int site)	{
		cerr << "in OmegaProcess::GetSiteOmega\n";
		exit(1);
	}

	virtual void UpdateSiteOmegaSuffStat()	{
		cerr << "in OmegaProcess::UpdateSiteOmegaSuffStat\n";
	}
	virtual void GlobalUpdateSiteOmegaSuffStat();
	virtual void SlaveUpdateSiteOmegaSuffStat();


	virtual void Create()	{
		if (! siteomegasuffstatcount)	{
			siteomegasuffstatbeta = new double[GetNsite()];
			siteomegasuffstatcount = new int[GetNsite()];
		}
	}

	virtual void Delete()	{
		if (siteomegasuffstatcount)	{
			delete[] siteomegasuffstatcount;
			delete[] siteomegasuffstatbeta;
		}
		siteomegasuffstatcount = 0;
	}

	int fixomega;
	int omegaprior;
	double* siteomegasuffstatbeta;
	int* siteomegasuffstatcount;
};


#endif



