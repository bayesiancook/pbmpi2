
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef MULTIGENERATEROCESS_H
#define MULTIGENERATEPROCESS_H

#include "MultiGeneMPIModule.h"
#include "DGamRateProcess.h"

class MultiGeneRateProcess : public virtual DGamRateProcess, public virtual MultiGeneMPIModule	{

	public:

	MultiGeneRateProcess() : globalalpha(1), genealpha(0) {}
	virtual ~MultiGeneRateProcess() {}

	virtual void Create();
	virtual void Delete();

	virtual void SlaveUpdateSiteRateSuffStat();
	virtual void UpdateRateSuffStat();

	void SetGlobalAlpha(int inmode)	{
		globalalpha = inmode;
	}
		
	int GlobalAlpha()	{
		return globalalpha;
	}

	double GetMeanAlpha()	{
		if (globalalpha)	{
			return GetAlpha();
		}
		return GlobalGetMeanAlpha();
	}

	double GlobalGetMeanAlpha();
	void SlaveGetMeanAlpha();

	void GlobalCollectGeneAlphas();
	void SlaveCollectGeneAlphas();

	virtual void SampleRate();
	virtual void PriorSampleRate();

	virtual double Move(double tuning = 1, int nrep = 1);
	double MoveGeneRateHyperParams(double tuning = 1, int nrep = 1);

	virtual double LogRatePrior();

	DGamRateProcess* GetRateProcess(int gene)	{

		DGamRateProcess* tmp = dynamic_cast<DGamRateProcess*>(process[gene]);
		if (!tmp)	{
			cerr << "error in GetRateProcess\n";
			exit(1);
		}
		return tmp;
	}

	int globalalpha;
	double* genealpha;
	double* tmpgenealpha;
};

#endif
