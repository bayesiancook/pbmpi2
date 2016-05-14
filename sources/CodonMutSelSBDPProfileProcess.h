
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef CODONMUTSELSBDPPROFILE_H
#define CODONMUTSELSBDPPROFILE_H

#include "SBDPProfileProcess.h"
#include "CodonMutSelProfileProcess.h"
#include "GeneralPathSuffStatMatrixMixtureProfileProcess.h"

class CodonMutSelSBDPProfileProcess : public virtual SBDPProfileProcess, public virtual CodonMutSelProfileProcess, public virtual GeneralPathSuffStatMatrixMixtureProfileProcess	{

	// implementer les fonctions create matrix et delete matrix
	// ainsi que CreateComponent(int k) and DeleteComponent(k)

	// s'inspirer de GeneralPathSuffStatGTRFiniteProfileProcess

	public:

	CodonMutSelSBDPProfileProcess() {}
	virtual ~CodonMutSelSBDPProfileProcess() {}

	protected:

	void Create()	{
		CodonMutSelProfileProcess::Create();
		SBDPProfileProcess::Create();
		GeneralPathSuffStatMatrixMixtureProfileProcess::Create();
	}
	
	void Delete()	{
		CodonMutSelProfileProcess::Delete();
		GeneralPathSuffStatMatrixMixtureProfileProcess::Delete();
		SBDPProfileProcess::Delete();
	}

	void ToStream(ostream& os)	{
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

		os << kappa << '\n';
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
	void FromStream(istream& is)	{
		for (int i=0; i<Nnuc; i++)	{
			is >> nucstat[i];
		}
		for (int i=0; i<GetNnucrr(); i++)	{
			is >> nucrr[i];
		}

		is >> kappa;
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
			cerr << "error in CodonMutSelSBDPProfileProcess: matrixarray is not 0\n";
			exit(1);
		}
		matrixarray[k] = new CodonMutSelProfileSubMatrix(GetCodonStateSpace(),nucrr,nucstat,profile[k],true);
	}
};

#endif

