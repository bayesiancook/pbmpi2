
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef ZIPMATRIXPHYLO_H
#define ZIPMATRIXPHYLO_H

#include "MatrixPhyloProcess.h"
#include "ZipMatrixSubstitutionProcess.h"

class ZipMatrixPhyloProcess : public virtual MatrixPhyloProcess, public virtual ZipMatrixSubstitutionProcess	{

	public:

	ZipMatrixPhyloProcess() : zipdata(0), truedata(0) {}
	virtual ~ZipMatrixPhyloProcess() {}

	// protected:

	// true data here !
	virtual void Create();
	virtual void Delete();

	// in fact, same object as GetData, but now with its true type
	ZippedSequenceAlignment* GetZipData()	{
		if (! zipdata)	{
			cerr << "null zip\n";
			exit(1);
		}
		return zipdata;
	}

	virtual int GetNstate() {return truedata->GetNstate();}
	int GetZipSize(int site) {return GetZipData()->GetZipSize(site);}
	int GetOrbitSize(int site) {return GetZipData()->GetOrbitSize(site);}
	int GetStateFromZip(int site, int state) {return GetZipData()->GetStateFromZip(site,state);}
	bool InOrbit(int site, int state) {return GetZipData()->InOrbit(site,state);}
	int* GetZipIndices(int site) {return zipdata->GetZipIndices(site);}
	
	void SetDataFromLeaves()	{
		SampleTrueNodeStates(GetRoot());
		PhyloProcess::SetDataFromLeaves();
	}

	void SampleTrueNodeStates(const Link* from) {
		cerr << "error: sample true node states not yet implemented in zip phyloprocess\n";
		exit(1);
	}

	void GlobalSetTestData() {}
	void SlaveSetTestData() {}

	// private:

	ZippedSequenceAlignment* zipdata;
	SequenceAlignment* truedata;
};


#endif
