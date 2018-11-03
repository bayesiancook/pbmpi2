
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef MULTIGENESBDPPROFILEPROCESS_H
#define MULTIGENESBDPPROFILEPROCESS_H

#include "MultiGeneProfileProcess.h"
#include "SBDPProfileProcess.h"

class MultiGeneSBDPProfileProcess : public virtual SBDPProfileProcess, public virtual MultiGeneProfileProcess	{

	public:

	MultiGeneSBDPProfileProcess() : kappaarray(0) {}
	virtual ~MultiGeneSBDPProfileProcess() {}

	virtual void Create();
	virtual void Delete();

	void GlobalCollectKappas();
	void SlaveCollectKappas();

	// override SBDPProfileProcess functions
	double LogHyperPrior();
	void SampleHyper();
	void PriorSampleHyper();

	double MoveHyper(double tuning, int nrep);
	double MoveKappaMean(double tuning);
	double MoveKappaRelVar(double tuning);

	SBDPProfileProcess* GetSBDPProfileProcess(int gene)	{

		SBDPProfileProcess* tmp = dynamic_cast<SBDPProfileProcess*>(process[gene]);
		if (!tmp)	{
			cerr << "error in MultiGeneSBDPProfileProcess::GetSBDPProfileProcess\n";
			exit(1);
		}
		return tmp;
	}

	double GlobalGetMeanNcomponent();
	double GlobalGetMeanStatEnt();
	double GlobalGetMeanStatAlpha();
	double GlobalGetMeanKappa();

	void SlaveGetMeanNcomponent();
	void SlaveGetMeanStatEnt();
	void SlaveGetMeanStatAlpha();
	void SlaveGetMeanKappa();

	double* kappaarray;
    double* statalphaarray;
    double* nmodearray;
};

#endif
