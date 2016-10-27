
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef SINGLEOMEGA_H
#define SINGLEOMEGA_H

#include "OmegaProcess.h"

class SingleOmegaProcess : public virtual OmegaProcess	{

	public:

	SingleOmegaProcess() : omega(0)	{}
	virtual ~SingleOmegaProcess() {}

	double GetOmega()	{
		return *omega;
	}

	double GetSiteOmega(int site)	{
		return *omega;
	}

	double* GetSiteOmegaPtr(int site)	{
		return omega;
	}
	
	virtual double OmegaSuffStatLogProb()	{
		return omegasuffstatcount * log(*omega) - omegasuffstatbeta * *omega;
	}

	virtual double CountOmegaSuffStatLogProb()	{
		return omegasuffstatcount * log(*omega);
	}

	virtual double BetaOmegaSuffStatLogProb()	{
		return - omegasuffstatbeta * *omega;
	}

	// omega
	virtual double LogOmegaPrior();
	virtual void SampleOmega();
	double SimpleMoveOmega(double tuning);
	double MoveOmega(double tuning); 

	protected:

	virtual void Create()	{
		if (! omega)	{
			OmegaProcess::Create();
			omega = new double;
			*omega = 1.0;
		}
	}

	virtual void Delete()	{
		if (omega)	{
			delete omega;
			omega = 0;
			OmegaProcess::Delete();
		}
	}
	
	void UpdateOmegaSuffStat();
	void GlobalUpdateOmegaSuffStat();
	void SlaveUpdateOmegaSuffStat();
	//void UpdateSiteOmegaSuffStat();

	double* omega;
	int omegasuffstatcount;
	double omegasuffstatbeta;
	
};

#endif


