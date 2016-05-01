
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef POISSONSBDPPROFILE_H
#define POISSONSBDPPROFILE_H

#include "PoissonDPProfileProcess.h"
#include "SBDPProfileProcess.h"

// superclass for Poisson (F81) implementations
class PoissonSBDPProfileProcess: public virtual PoissonDPProfileProcess, public virtual SBDPProfileProcess	{

	public:

	PoissonSBDPProfileProcess() : InitIncremental(0) {}
	virtual ~PoissonSBDPProfileProcess() {}

	virtual double LogProxy(int site, int cat)	{
		return DiffLogSampling(cat,site);
	}

	virtual double Move(double tuning = 1, int n = 1, int nrep = 1)	{
		return SBDPProfileProcess::Move(tuning,n,nrep);
	}

	protected:

	virtual void Create()	{
		PoissonDPProfileProcess::Create();
		SBDPProfileProcess::Create();
	}

	virtual void Delete()	{
		SBDPProfileProcess::Delete();
		PoissonDPProfileProcess::Delete();
	}

	double GlobalMixMove(int nrep, int nallocrep, double epsilon, int nprofilerep);
	void SlaveMixMove();
	double MixMove(int nrep, int nallocrep, double epsilon, int nprofilerep);

	double IncrementalDPMove(int nrep, double c)	{
		cerr << "error : in poisson sbdp incremental\n";
		exit(1);
		return 1;
	}
	double IncrementalDPMove(int nrep)	{
		cerr << "error : in poisson sbdp incremental\n";
		exit(1);
		return 1;
	}

	virtual void SwapComponents(int cat1, int cat2)	{
		SBDPProfileProcess::SwapComponents(cat1,cat2);
	}

	// virtual void ToStream(ostream& os);
	virtual void FromStream(istream& is)	{
		PoissonDPProfileProcess::FromStream(is);
		ResampleWeights();
	}

	int InitIncremental;

};

#endif

