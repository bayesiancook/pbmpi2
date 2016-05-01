
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef MATRIXMIXTUREPROFILE_H
#define MATRIXMIXTUREPROFILE_H

#include "MatrixProfileProcess.h"
#include "MixtureProfileProcess.h"

class MatrixMixtureProfileProcess : public virtual MatrixProfileProcess, public virtual MixtureProfileProcess {

	public:

	MatrixMixtureProfileProcess() : matrixarray(0) {}
	virtual ~MatrixMixtureProfileProcess() {}

	SubMatrix* GetOriginalMatrix(int site)	{
		if (! matrixarray[alloc[site]]) 	{
			cerr << "error in get matrix : null matrix , site " << site << " and alloc : " << alloc[site] << '\n';
			exit(1);
		}
		return matrixarray[alloc[site]];
	}

	protected:

	Chrono sschrono;

	// called at the beginning and the end of the run
	virtual void Create();
	virtual void Delete();

	// should be called each time global parameters are modified
	virtual void UpdateMatrices()	{
		for (int k=0; k<GetNcomponent(); k++)	{
			UpdateMatrix(k);
		}
	}

	virtual void CreateMatrices()	{
		for (int k=0; k<GetNmodeMax(); k++)	{
			if (! matrixarray[k])	{
				CreateMatrix(k);
			}
		}
		/*
		for (int k=0; k<GetNcomponent(); k++)	{
			if (! matrixarray[k])	{
				CreateMatrix(k);
			}
		}
		for (int k=GetNcomponent(); k<GetNmodeMax(); k++)	{
			if (matrixarray[k])	{
				DeleteMatrix(k);
			}
			matrixarray[k] = 0;
		}
		*/
	}

	virtual void DeleteMatrices()	{
		for (int k=0; k<GetNmodeMax(); k++)	{
		// for (int k=0; k<GetNcomponent(); k++)	{
			DeleteMatrix(k);
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
	virtual void UpdateComponent(int k)	{
		UpdateMatrix(k);
	}

	SubMatrix** matrixarray;
};

#endif
