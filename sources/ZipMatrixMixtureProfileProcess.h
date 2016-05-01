
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef ZIPMATMIXTUREPROFILE_H
#define ZIPMATMIXTUREPROFILE_H

#include "MatrixMixtureProfileProcess.h"
#include "ZipMatrixProfileProcess.h"

// general subclass for all Matrix Mixture mixtures using generic sufficient statistics
class ZipMatrixMixtureProfileProcess : public virtual ZipMatrixProfileProcess, public virtual MatrixMixtureProfileProcess  {

	public:

	ZipMatrixMixtureProfileProcess() {}
	virtual ~ZipMatrixMixtureProfileProcess() {}

	protected:

	virtual void Create()	{
		ZipMatrixProfileProcess::Create();
		MatrixMixtureProfileProcess::Create();
	}

	virtual void Delete()	{
		MatrixMixtureProfileProcess::Delete();
		ZipMatrixProfileProcess::Delete();
	}

	virtual void CreateMatrices()	{
		MatrixMixtureProfileProcess::CreateMatrices();
		if (GetMyid())	{
			ZipMatrixProfileProcess::CreateMatrices();
		}
	}

	virtual void DeleteMatrices()	{
		if (GetMyid())	{
			ZipMatrixProfileProcess::DeleteMatrices();
		}
		MatrixMixtureProfileProcess::DeleteMatrices();
	}

	virtual void UpdateMatrix(int k)	{
		MatrixMixtureProfileProcess::UpdateMatrix(k);
		if (GetMyid())	{
			for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
				if (ActiveSite(i))	{
					if ((alloc[i] == k) && (zipmatrixarray[i]))	{
						if (zipmatrixarray[i]->GetFrom() != matrixarray[k])	{
							zipmatrixarray[i]->SetFrom(matrixarray[k]);
						}
						zipmatrixarray[i]->CorruptMatrix();
					}
				}
			}
		}
	}	
};

#endif

