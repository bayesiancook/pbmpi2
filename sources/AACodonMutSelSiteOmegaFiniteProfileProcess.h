
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef AACODONMUTSELSITEOMFINITEPROFILE_H
#define AACODONMUTSELSITEOMFINITEPROFILE_H

#include "FiniteProfileProcess.h"
#include "SiteOmegaProcess.h"
#include "AACodonMutSelSiteMatrixMixtureProfileProcess.h"

class AACodonMutSelSiteOmegaFiniteProfileProcess : public virtual FiniteProfileProcess, public virtual SiteOmegaProcess, public virtual AACodonMutSelSiteMatrixMixtureProfileProcess	{

	public:

	AACodonMutSelSiteOmegaFiniteProfileProcess() {}
	virtual ~AACodonMutSelSiteOmegaFiniteProfileProcess() {}

	protected:

	void Create()	{
		SiteOmegaProcess::Create();
		AACodonMutSelSiteMatrixMixtureProfileProcess::Create();
		FiniteProfileProcess::Create();
	}
	
	void Delete()	{
		FiniteProfileProcess::Delete();
		AACodonMutSelSiteMatrixMixtureProfileProcess::Delete();
		SiteOmegaProcess::Delete();
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
		for (int i=0; i<GetNsite(); i++)	{
			if (ActiveSite(i))	{
				os << omegaarray[i] << '\t';
			}
		}
		os << '\n';
		os << omegaalpha << '\t' << omegabeta << '\n';
		os << '\n';

		os << Ncomponent << '\n';
        for (int k=0; k<Nstatcomp; k++) {
            os << statweight[k] << '\t';
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

		for (int i=0; i<GetNsite(); i++)	{
			is >> omegaarray[i];
		}
		is >> omegaalpha >> omegabeta;

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

