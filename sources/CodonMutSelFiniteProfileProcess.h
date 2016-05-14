
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef CODONMUTSELFINITEPROFILE_H
#define CODONMUTSELFINITEPROFILE_H

#include "FiniteProfileProcess.h"
#include "CodonMutSelProfileProcess.h"
#include "GeneralPathSuffStatMatrixMixtureProfileProcess.h"

class CodonMutSelFiniteProfileProcess : public virtual FiniteProfileProcess, public virtual CodonMutSelProfileProcess, public virtual GeneralPathSuffStatMatrixMixtureProfileProcess	{

	// implementer les fonctions create matrix et delete matrix
	// ainsi que CreateComponent(int k) and DeleteComponent(k)

	// s'inspirer de GeneralPathSuffStatGTRFiniteProfileProcess

	public:

	CodonMutSelFiniteProfileProcess() {}
	virtual ~CodonMutSelFiniteProfileProcess() {}

	protected:

	void Create()	{
		CodonMutSelProfileProcess::Create();
		FiniteProfileProcess::Create();
		GeneralPathSuffStatMatrixMixtureProfileProcess::Create();
	}
	
	void Delete()	{
		CodonMutSelProfileProcess::Delete();
		GeneralPathSuffStatMatrixMixtureProfileProcess::Delete();
		FiniteProfileProcess::Delete();
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

		os << Ncomponent << '\n';
		for (int j=0; j<GetDim(); j++)	{
			os << dirweight[j] << '\t';
		}
		os << '\n';
		os << '\n';
		for (int i=0; i<Ncomponent; i++)	{
			for (int j=0; j<GetDim(); j++)	{
				os << profile[i][j] << '\t';
			}
			os << '\n';
		}
		for (int i=0; i<GetNsite(); i++)	{
			if (ActiveSite(i))	{
				os << alloc[i] << '\t';
			}
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

		is >> Ncomponent;
		for (int j=0; j<GetDim(); j++)	{
			is >> dirweight[j];
		}
		for (int i=0; i<Ncomponent; i++)	{
			for (int j=0; j<GetDim(); j++)	{
				is >> profile[i][j];
			}
		}
		for (int i=0; i<GetNsite(); i++)	{
			if (ActiveSite(i))	{
				is >> alloc[i];
			}
		}
	}

	void CreateMatrix(int k)	{
		if (matrixarray[k])	{
			cerr << "error in AACodonMutSelFiniteProfileProcess: matrixarray is not 0\n";
			exit(1);
		}
		matrixarray[k] = new CodonMutSelProfileSubMatrix(GetCodonStateSpace(),nucrr,nucstat,profile[k],true);
	}
};

#endif

