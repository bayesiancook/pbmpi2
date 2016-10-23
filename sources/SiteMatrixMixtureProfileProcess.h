// in this version, SiteMatrixMixture now combines two systems of matrices: per component and per site

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

class SiteMatrixMixtureProfileProcess : public virtual MatrixMixtureProfileProcess	{

	public:

	SiteMatrixMixtureProfileProcess() : sitematrixarray(0) {}
	virtual ~SiteMatrixMixtureProfileProcess() {}

	SubMatrix* GetOriginalMatrix(int site)	{
		if (! sitematrixarray[site]) 	{
			cerr << "error in get matrix : null matrix , site " << site << " and alloc : " << alloc[site] << '\n';
			exit(1);
		}
		return sitematrixarray[site];
	}

	virtual void UpdateMatrices()	{
		cerr << "in SiteMatrixMixtureProfileProcess::UpdateMatrices: should we update site matrices?\n";
		exit(1);
	}

	virtual void UpdateSiteMatrices()	{
		// for (int i=0; i<GetNsite(); i++)	{
		for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
			if (ActiveSite(i))	{
				UpdateSiteMatrix(i);
			}
		}
	}

	protected:

	Chrono sschrono;

	// called at the beginning and the end of the run
	virtual void Create()	{
		if (! sitematrixarray)	{
			MatrixMixtureProfileProcess::Create();
			sitematrixarray = new SubMatrix*[GetNsite()];
			for (int i=0; i<GetNsite(); i++)	{
				sitematrixarray[i] = 0;
			}
		}
	}

	virtual void Delete()	{
		if (sitematrixarray)	{
			for (int i=0; i<GetNsite(); i++)	{
				delete sitematrixarray[i];
			}
			delete[] sitematrixarray;
			sitematrixarray = 0;
			MatrixMixtureProfileProcess::Delete();
		}
	}

	virtual void CreateSiteMatrices()	{
		for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
			if (ActiveSite(i))	{
				CreateSiteMatrix(i);
			}
		}
	}

	virtual void DeleteSiteMatrices()	{
		for (int i=0; i<GetNsite(); i++)	{
			DeleteSiteMatrix(i);
		}
	}

	virtual void CreateSiteMatrix(int site) = 0;

	virtual void DeleteSiteMatrix(int site)	{
		delete sitematrixarray[site];
		sitematrixarray[site] = 0;
	}

	virtual void UpdateSiteMatrix(int site)	{
		if (sitematrixarray[site])	{
			sitematrixarray[site]->CorruptMatrix();
		}
	}

	/*
	virtual void UpdateComponent(int k)	{
		MatrixMixtureProfileProcess::UpdateComponent(k);
		for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
			if (ActiveSite(i) && (alloc[i] == k))	{
				UpdateSiteMatrix(i);
			}
		}
	}
	*/

	SubMatrix** sitematrixarray;
};

#endif
