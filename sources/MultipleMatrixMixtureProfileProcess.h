
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef MULTIPLEMATRIXMIXTUREPROFILE_H
#define MULTIPLEMATRIXMIXTUREPROFILE_H

#include "MatrixProfileProcess.h"
#include "MixtureProfileProcess.h"

class MultipleMatrixMixtureProfileProcess : public virtual MatrixProfileProcess, public virtual MixtureProfileProcess {

	public:

	MultipleMatrixMixtureProfileProcess() : matrixarray(0) {}
	virtual ~MultipleMatrixMixtureProfileProcess() {}

	// gives the allocation of site to one of the components of the 'sub-mixture' (e.g. mixture of omega values for mut sel models)
	virtual int GetSubAlloc(int site) = 0;
	virtual int GetNSubAlloc() = 0;

	SubMatrix* GetOriginalMatrix(int site)	{
		if (! matrixarray[alloc[site]][GetSubAlloc(site)]) 	{
			cerr << "error in get matrix : null matrix , site " << site << " and alloc : " << alloc[site] << '\t' << GetSubAlloc(site) << '\n';
			exit(1);
		}
		return matrixarray[alloc[site]][GetSubAlloc(site)];
	}

	// should be called each time global parameters are modified
	virtual void UpdateMatrices()	{
		for (int k=0; k<GetNcomponent(); k++)	{
			for (int l=0; l<GetNSubAlloc(); l++)	{
				UpdateMatrix(k,l);
			}
		}
	}

	protected:

	Chrono sschrono;

	// called at the beginning and the end of the run
	virtual void Create();
	virtual void Delete();

	virtual void CreateMatrices()	{
		for (int k=0; k<GetNmodeMax(); k++)	{
			for (int l=0; l<GetNSubAlloc(); l++)	{
				CreateMatrix(k,l);
			}
		}
	}

	virtual void DeleteMatrices()	{
		for (int k=0; k<GetNmodeMax(); k++)	{
			for (int l=0; l<GetNSubAlloc(); l++)	{
				DeleteMatrix(k,l);
			}
		}
	}

	virtual void CreateMatrix(int alloc, int suballoc) = 0;

	virtual void DeleteMatrix(int alloc, int suballoc)	{
		delete matrixarray[alloc][suballoc];
		matrixarray[alloc][suballoc] = 0;
	}

	virtual void UpdateMatrix(int alloc, int suballoc)	{
		if (matrixarray[alloc][suballoc])	{
			matrixarray[alloc][suballoc]->CorruptMatrix();
		}
	}

	// update component k (in particular, will be used for updating the matrix
	// typically, after the profile has changed)
	virtual void UpdateComponent(int k)	{
		for (int l=0; l<GetNSubAlloc(); l++)	{
			UpdateMatrix(k,l);
		}
	}

	SubMatrix*** matrixarray;
};

#endif
