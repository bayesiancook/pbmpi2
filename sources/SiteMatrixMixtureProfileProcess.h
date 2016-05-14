
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef SITEMATRIXMIXTUREPROFILE_H
#define SITEMATRIXMIXTUREPROFILE_H

#include "MatrixProfileProcess.h"
#include "MixtureProfileProcess.h"

class SiteMatrixMixtureProfileProcess : public virtual MatrixProfileProcess, public virtual MixtureProfileProcess {

	public:

	SiteMatrixMixtureProfileProcess() : matrixarray(0) {}
	virtual ~SiteMatrixMixtureProfileProcess() {}

	SubMatrix* GetOriginalMatrix(int site)	{
		if (! matrixarray[site]) 	{
			cerr << "error in get matrix : null matrix , site " << site << " and alloc : " << alloc[site] << '\n';
			exit(1);
		}
		return matrixarray[site];
	}

	virtual void UpdateMatrices()	{
		// for (int i=0; i<GetNsite(); i++)	{
		for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
			if (ActiveSite(i))	{
				UpdateMatrix(i);
			}
		}
	}

	protected:

	Chrono sschrono;

	// called at the beginning and the end of the run
	virtual void Create()	{
		if (! matrixarray)	{
			MixtureProfileProcess::Create();
			matrixarray = new SubMatrix*[GetNsite()];
			for (int i=0; i<GetNsite(); i++)	{
				matrixarray[i] = 0;
			}
			CreateMatrices();
		}
	}

	virtual void Delete()	{
		if (matrixarray)	{
			DeleteMatrices();
			delete[] matrixarray;
			matrixarray = 0;
			MixtureProfileProcess::Delete();
		}
	}

	virtual void CreateMatrices()	{
		// for (int i=0; i<GetNsite(); i++)	{
		for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
			if (ActiveSite(i))	{
				CreateMatrix(i);
			}
		}
	}

	virtual void DeleteMatrices()	{
		for (int i=0; i<GetNsite(); i++)	{
			DeleteMatrix(i);
		}
	}

	virtual void CreateMatrix(int site) = 0;

	virtual void DeleteMatrix(int site)	{
		delete matrixarray[site];
		matrixarray[site] = 0;
	}

	virtual void UpdateMatrix(int site)	{
		if (matrixarray[site])	{
			matrixarray[site]->CorruptMatrix();
		}
	}

	// update component k (in particular, will be used for updating the matrix
	// typically, after the profile has changed)
	virtual void UpdateComponent(int k)	{
		// for (int i=0; i<GetNsite(); i++)	{
		for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
			if (ActiveSite(i) && (alloc[i] == k))	{
				UpdateMatrix(i);
			}
		}
	}

	SubMatrix** matrixarray;
};

#endif
