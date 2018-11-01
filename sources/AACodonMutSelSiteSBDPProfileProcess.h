
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef AACODONMUTSELSITESBDPPROFILE_H
#define AACODONMUTSELSITESBDPPROFILE_H

#include "SBDPProfileProcess.h"
#include "AACodonMutSelSiteMatrixMixtureProfileProcess.h"
#include "SingleOmegaProcess.h"

class AACodonMutSelSiteSBDPProfileProcess : public virtual SBDPProfileProcess, public virtual SingleOmegaProcess, public virtual AACodonMutSelSiteMatrixMixtureProfileProcess	{

	public:

	AACodonMutSelSiteSBDPProfileProcess() {}
	virtual ~AACodonMutSelSiteSBDPProfileProcess() {}

	protected:

	void Create()	{
		SingleOmegaProcess::Create();
		AACodonMutSelSiteMatrixMixtureProfileProcess::Create();
		SBDPProfileProcess::Create();
	}
	
	void Delete()	{
		SBDPProfileProcess::Delete();
		AACodonMutSelSiteMatrixMixtureProfileProcess::Delete();
		SingleOmegaProcess::Delete();
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
		os << *omega << '\n';

		os << kappa << '\n';
		os << Ncomponent << '\n';
        for (int k=0; k<Nstatcomp; k++) {
            os << statweight[k];
            for (int j=0; j<GetDim(); j++)	{
                os << dirweight[k][j] << '\t';
            }
            os << '\n';
        }
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
		for (int i=0; i<GetNcodon(); i++)	{
			is >> codonprofile[i];
		}

		is >> *omega;
		is >> kappa;
		is >> Ncomponent;
        for (int k=0; k<Nstatcomp; k++) {
            is >> statweight[k];
            for (int j=0; j<GetDim(); j++)	{
                is >> dirweight[k][j];
            }
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
		ResampleWeights();
	}
};

#endif

