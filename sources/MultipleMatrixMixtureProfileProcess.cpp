
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "MultipleMatrixMixtureProfileProcess.h"
#include "Random.h"
#include "Parallel.h"


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* MatrixMixtureProfileProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void MultipleMatrixMixtureProfileProcess::Create()	{
	if (! matrixarray)	{
		MixtureProfileProcess::Create();
		matrixarray = new SubMatrix**[GetNmodeMax()];
		for (int k=0; k<GetNmodeMax(); k++)	{
			matrixarray[k] = new SubMatrix*[GetNSubAlloc()];
			for (int l=0; l<GetNSubAlloc(); l++)	{
				matrixarray[k][l] = 0;
			}
		}
		CreateMatrices();
	}
}

void MultipleMatrixMixtureProfileProcess::Delete() {
	if (matrixarray)	{
		DeleteMatrices();
		for (int k=0; k<GetNmodeMax(); k++)	{
			delete[] matrixarray[k];
		}
		delete[] matrixarray;
		matrixarray = 0;
		MixtureProfileProcess::Delete();
	}
}

