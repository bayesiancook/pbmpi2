
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef AACODONMUTSELMVNSSPROFILE_H
#define AACODONMUTSELMVNSSPROFILE_H

#include "SingleOmegaProcess.h"
#include "AACodonMutSelProfileProcess.h"
#include "GeneralPathSuffStatMatrixMVNSiteSpecificProfileProcess.h"

class AACodonMutSelMVNSiteSpecificProfileProcess : public virtual SingleOmegaProcess, public virtual AACodonMutSelProfileProcess, public virtual GeneralPathSuffStatMatrixMVNSiteSpecificProfileProcess	{

	public:

	AACodonMutSelMVNSiteSpecificProfileProcess() {}
	virtual ~AACodonMutSelMVNSiteSpecificProfileProcess() {}

	protected:

	void Create()	{
		GeneralPathSuffStatMatrixMVNSiteSpecificProfileProcess::Create();
		SingleOmegaProcess::Create();
		AACodonMutSelProfileProcess::Create();
	}
	
	void Delete()	{
		AACodonMutSelProfileProcess::Delete();
		SingleOmegaProcess::Delete();
		GeneralPathSuffStatMatrixMVNSiteSpecificProfileProcess::Delete();
	}

	void ToStream(ostream& os) {
		for (int i=0; i<Nnuc; i++)	{
			os << GetNucStat(i) << '\t';
		}
		os << '\n';
		os << '\n';
		for (int i=0; i<GetNnucrr(); i++)	{
			os << GetNucRR(i) << '\t';
		}
		os << '\n';
		os << '\n';		
		for (int i=0; i<GetNcodon(); i++)	{
			os << codonprofile[i] << '\t';
		}
		os << '\n';
		os << '\n';
		os << *omega << '\n';

		for (int j=0; j<GetDim(); j++)	{
			os << kappa[j] << '\t';
		}
		os << '\n';
		os << df << '\n';
		for (int j=0; j<GetDim(); j++)	{
			for (int k=0; k<GetDim(); k++)	{
				os << (*covmatrix)[j][k] << '\t';
			}
			os << '\n';
		}
		os << '\n';
		for (int i=0; i<GetNsite(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				os << logprofile[i][j] << '\t';
			}
			os << '\n';
		}
		os << '\n';
	}

	void FromStream(istream& is) {
		for (int i=0; i<Nnuc; i++)	{
			is >> nucstat[i];
		}
		for (int i=0; i<GetNnucrr(); i++)	{
			is >> nucrr[i];
		}

		for (int i=0; i<GetNcodon(); i++)	{
			is >> codonprofile[i];
		}

		is >> *omega;
		for (int j=0; j<GetDim(); j++)	{
			is >> kappa[j];
		}
		is >> df;
		for (int j=0; j<GetDim(); j++)	{
			for (int k=0; k<GetDim(); k++)	{
				is >> (*covmatrix)[j][k];
			}
		}
		for (int i=0; i<GetNsite(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				is >> logprofile[i][j];
			}
		}
	}

	void CreateMatrix(int k)	{
		if (matrixarray[k])	{
			cerr << "error in AACodonMutSelMVNSiteSpecificProfileProcess: matrixarray is not 0\n";
			exit(1);
		}
		matrixarray[k] = new AACodonMutSelProfileSubMatrix(GetCodonStateSpace(),nucrr,nucstat,codonprofile,profile[k],omega,true);
	}

};

#endif

