
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef MATRIXSSPROFILE_H
#define MATRIXSSPROFILE_H

#include "MatrixProfileProcess.h"
#include "SiteSpecificProfileProcess.h"

class MatrixSiteSpecificProfileProcess : public virtual MatrixProfileProcess, public virtual SiteSpecificProfileProcess {

	public:

	MatrixSiteSpecificProfileProcess() : matrixarray(0) {}
	virtual ~MatrixSiteSpecificProfileProcess() {}

	SubMatrix* GetOriginalMatrix(int site)	{
		if (! matrixarray[site]) 	{
			cerr << "error in get matrix : null matrix , site " << site << '\n';
			exit(1);
		}
		return matrixarray[site];
	}

	protected:

	Chrono sschrono;

	// called at the beginning and the end of the run
	virtual void Create();
	virtual void Delete();

	// should be called each time global parameters are modified
	virtual void UpdateMatrices()	{
		for (int k=0; k<GetNsite(); k++)	{
			if (ActiveSite(k))	{
				UpdateMatrix(k);
			}
		}
	}

	virtual void CreateMatrices()	{
		for (int k=0; k<GetNsite(); k++)	{
			if (ActiveSite(k))	{
				if (! matrixarray[k])	{
					CreateMatrix(k);
				}
			}
		}
	}

	virtual void DeleteMatrices()	{
		for (int k=0; k<GetNsite(); k++)	{
			if (matrixarray[k])	{
				DeleteMatrix(k);
			}
		}
	}

	virtual void CreateMatrix(int k) = 0;

	virtual void DeleteMatrix(int k)	{
		delete matrixarray[k];
		matrixarray[k] = 0;
	}

	virtual void UpdateMatrix(int k)	{
		if (matrixarray[k])	{
			matrixarray[k]->CorruptMatrix();
		}
	}

	// update component k (in particular, will be used for updating the matrix
	// typically, after the profile has changed)
	virtual void UpdateSite(int k)	{
		UpdateMatrix(k);
	}

	SubMatrix** matrixarray;
};

#endif
