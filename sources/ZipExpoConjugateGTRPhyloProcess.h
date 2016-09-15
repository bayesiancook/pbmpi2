
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef ZIPEXPCONGTRPHYLO_H
#define ZIPEXPCONGTRPHYLO_H

#include "ExpoConjugateGTRPhyloProcess.h"
#include "ZipExpoConjugateGTRSubstitutionProcess.h"
#include "ZipMatrixPhyloProcess.h"

class ZipExpoConjugateGTRPhyloProcess : public virtual ExpoConjugateGTRPhyloProcess, public virtual ZipExpoConjugateGTRSubstitutionProcess, public virtual ZipMatrixPhyloProcess {

	public:

	ZipExpoConjugateGTRPhyloProcess() : zipmatcreated(0) {}
	virtual ~ZipExpoConjugateGTRPhyloProcess() {}

	virtual int GetNstate() {return ZipExpoConjugateGTRSubstitutionProcess::GetNstate();}
	virtual int GetNstate(int i) {return ZipExpoConjugateGTRSubstitutionProcess::GetNstate(i);}

	virtual SequenceAlignment* GetData()	{
		if (zipmode)	{
			return zipdata;
		}
		if (! truedata)	{
			return data;
		}
		return truedata;
	}

	virtual const SequenceAlignment* GetData() const	{
		if (zipmode)	{
			return zipdata;
		}
		if (! truedata)	{
			return data;
		}
		return truedata;
	}

	protected:

	virtual void Create()	{
		ZipMatrixPhyloProcess::Create();
		ExpoConjugateGTRPhyloProcess::Create();
		ZipExpoConjugateGTRSubstitutionProcess::Create();
	}

	virtual void Delete()	{
		ZipExpoConjugateGTRSubstitutionProcess::Delete();
		ExpoConjugateGTRPhyloProcess::Delete();
		ZipMatrixPhyloProcess::Delete();
	}

	void Unfold();
	void Collapse();

	// double ZipTopoMoveCycle(int nrep, double tuning);

	virtual void GlobalActivateFastTopo();
	virtual void GlobalInactivateFastTopo();

	int zipmatcreated;

};

#endif

