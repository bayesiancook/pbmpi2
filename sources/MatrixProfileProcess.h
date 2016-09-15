
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef MATRIXPROFILE_H
#define MATRIXPROFILE_H

#include "SubMatrix.h"
#include "ProfileProcess.h"

// superclass for all Matrix implementations
class MatrixProfileProcess : public virtual ProfileProcess	{

	public:

	MatrixProfileProcess() {}
	virtual ~MatrixProfileProcess() {}

	// access to matrices
	virtual SubMatrix* GetMatrix(int site) {return GetOriginalMatrix(site);}

	virtual SubMatrix* GetOriginalMatrix(int site) = 0;

	protected:

	// create/delete all matrices
	virtual void CreateMatrices() = 0;
	virtual void DeleteMatrices() = 0;

	// updates all matrices
	virtual void UpdateMatrices() = 0;
};

#endif

