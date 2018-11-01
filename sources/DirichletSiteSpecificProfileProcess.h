
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef DIRSSPROFILE_H
#define DIRSSPROFILE_H

#include <cmath>

#include "SiteSpecificProfileProcess.h"
#include "DirichletProfileProcess.h"

class DirichletSiteSpecificProfileProcess : public virtual DirichletProfileProcess, public virtual SiteSpecificProfileProcess	{

	public:

	DirichletSiteSpecificProfileProcess() {}
	virtual ~DirichletSiteSpecificProfileProcess() {}

	virtual void Create()	{
		DirichletProfileProcess::Create();
		SiteSpecificProfileProcess::Create();
	}

	virtual void Delete()	{
		SiteSpecificProfileProcess::Delete();
		DirichletProfileProcess::Delete();
	}

    virtual void ToStream(ostream& os)  {
        if (Nstatcomp > 1)  {
            cerr << "error in DirichletSiteSpecificProfileProcess ToStream\n";
            exit(1);
        }
        for (int k=0; k<GetDim(); k++)  {
            os << dirweight[0][k] << '\t';
        }
        os << '\n';
        for (int i=0; i<GetNsite(); i++)  {
            for (int k=0; k<GetDim(); k++)  {
                os << profile[i][k] << '\t';
            }
            os << '\n';
        }
    }

    virtual void FromStream(istream& is)    {
        if (Nstatcomp > 1)  {
            cerr << "error in DirichletSiteSpecificProfileProcess FromStream\n";
            exit(1);
        }
        for (int k=0; k<GetDim(); k++)  {
            is >> dirweight[0][k];
        }
        for (int i=0; i<GetNsite(); i++)  {
            for (int k=0; k<GetDim(); k++)  {
                is >> profile[i][k];
            }
        }
    }
};


#endif

