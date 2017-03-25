
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef AACODONMUTSELFINITEOMEGASBDPPROFILE_H
#define AACODONMUTSELFINITEOMEGASBDPPROFILE_H

#include "SBDPProfileProcess.h"
#include "FiniteOmegaProcess.h"
#include "AACodonMutSelProfileProcess.h"
#include "GeneralPathSuffStatMultipleMatrixMixtureProfileProcess.h"

class AACodonMutSelFiniteOmegaSBDPProfileProcess : public virtual SBDPProfileProcess, public virtual FiniteOmegaProcess, public virtual AACodonMutSelProfileProcess, public virtual GeneralPathSuffStatMultipleMatrixMixtureProfileProcess	{

	public:

	AACodonMutSelFiniteOmegaSBDPProfileProcess() {}
	virtual ~AACodonMutSelFiniteOmegaSBDPProfileProcess() {}

	int GetNSubAlloc()	{
		return GetNomega();
	}

	int GetSubAlloc(int site)	{
		return GetOmegaSiteAlloc(site);
	}
	
	protected:

	void Create()	{
		FiniteOmegaProcess::Create();
		AACodonMutSelProfileProcess::Create();
		SBDPProfileProcess::Create();
		GeneralPathSuffStatMultipleMatrixMixtureProfileProcess::Create();
	}
	
	void Delete()	{
		GeneralPathSuffStatMultipleMatrixMixtureProfileProcess::Delete();
		SBDPProfileProcess::Delete();
		AACodonMutSelProfileProcess::Delete();
		FiniteOmegaProcess::Delete();
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

		for (int i=0; i<GetNcodon(); i++)	{
			os << codonprofile[i] << '\t';
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
		os << Nomega << '\n';
		for (int i=0; i<Nomega; i++)	{
			os << omega[i] << '\t'; 	
		}

		os << '\n';
		os << '\n';
		for (int i=0; i<GetNsite(); i++)	{
			if (ActiveSite(i))	{
				os << alloc[i] << '\t';
				os << omegaalloc[i] << '\t';
			}
		}
		os << '\n';
		os << '\n';
		
	}
	void FromStream(istream& is)	{
		for (int i=0; i<Nnuc; i++)	{
			is >> nucstat[i];
		}
		for (int i=0; i<GetNnucrr(); i++)	{
			is >> nucrr[i];
		}
		for (int i=0; i<GetNcodon(); i++)	{
			is >> codonprofile[i];
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
		is >> Nomega;
		for (int i=0; i<Nomega; i++)	{
			is >> omega[i]; 	
		}
		for (int i=0; i<GetNsite(); i++)	{
			if (ActiveSite(i))	{
				is >> alloc[i];
				is >> omegaalloc[i];
			}
		}
		ResampleWeights();
		ResampleOmegaWeights();
	}


	void CreateMatrix(int alloc, int suballoc)	{
		if (matrixarray[alloc][suballoc])	{
			cerr << "error in AACodonMutSelSBDPProfileProcess: matrixarray is not 0\n";
			exit(1);
		}
		matrixarray[alloc][suballoc] = new AACodonMutSelProfileSubMatrix(GetCodonStateSpace(),nucrr,nucstat,codonprofile,profile[alloc],GetOmegaPointer(suballoc),true);
	}
};

#endif

