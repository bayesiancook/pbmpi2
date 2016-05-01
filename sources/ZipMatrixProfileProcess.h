
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef ZIPMATRIXPROFILE_H
#define ZIPMATRIXPROFILE_H

#include "MatrixProfileProcess.h"
#include "ZipSubMatrix.h"

// superclass for all Matrix implementations
class ZipMatrixProfileProcess : public virtual MatrixProfileProcess	{

	public:

	ZipMatrixProfileProcess() : zipmatrixarray(0) {}
	virtual ~ZipMatrixProfileProcess() {}

	// access to matrices
	virtual SubMatrix* GetMatrix(int site)	{
		return zipmatrixarray[site];
	}

	protected:

	virtual void Create();
	virtual void Delete();

	virtual int GetZipSize(int site) = 0;
	virtual int* GetZipIndices(int site) = 0;
	virtual int GetOrbitSize(int site) = 0;
	virtual int GetStateFromZip(int site, int state) = 0;

	// create/delete all matrices
	// Create called when deactivating sufficient statistics and activating pruning-based computation (Unfold  in PhyloProcess)
	virtual void CreateMatrices()	{
		if (! GetMyid())	{
			cerr << "error: master in ZipMatrixProfileProcess::CreateMatrices\n";
			exit(1);
		}	
		for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
			if (ActiveSite(i))	{
				if (! zipmatrixarray[i])	{
					SubMatrix* temp = GetOriginalMatrix(i);
					zipmatrixarray[i] = new ZipSubMatrix(temp,GetZipSize(i),GetZipIndices(i));
				}
			}
		}
	}

	virtual void DeleteMatrices()	{
		for (int i=0; i<GetNsite(); i++)	{
			if (zipmatrixarray[i])	{
				delete zipmatrixarray[i];
				zipmatrixarray[i] = 0;
			}
		}
	}

	// updates all matrices
	// (should be called, e.g. when performing a Metropolis on relative exchangeabilities or global mutation parameters)
	// virtual void UpdateMatrices() = 0;
	// virtual void DiagonaliseMatrices() = 0;

	ZipSubMatrix** zipmatrixarray;
};

#endif

