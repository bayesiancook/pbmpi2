
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
		for (int i=0; i<GetNcodon(); i++)	{
			is >> codonprofile[i];
		}

		is >> *omega;
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
		ResampleWeights();
	}


	/*
	void CreateMatrix(int site)	{
		if (matrixarray[site])	{
			cerr << "error in AACodonMutSelSSiteBDPProfileProcess: matrixarray is not 0\n";
			exit(1);
		}
		matrixarray[site] = new AACodonMutSelProfileSubMatrix(GetCodonStateSpace(),nucrr,nucstat,codonprofile,GetProfile(site),omega,true);
	}
	*/
	virtual void CreateMatrix(int k)	{
		if (matrixarray[k])	{
			cerr << "error in AACodonMutSelSiteMatrixMixtureProfileProcess: matrixarray is not 0\n";
			exit(1);
		}
		matrixarray[k] = new AACodonMutSelProfileSubMatrix(GetCodonStateSpace(),nucrr,nucstat,codonprofile,profile[k],omega0,true);
	}

	virtual void CreateSiteMatrix(int i)	{
		if (sitematrixarray[i])	{
			cerr << "error in AACodonMutSelSiteMatrixMixtureProfileProcess: sitematrixarray is not 0\n";
			exit(1);
		}
		if (*omega < 1e-6)	{
			cerr << "in create site matrix:small omega : " << *omega << '\n';
			exit(1);
		}
		if (GetProfile(i) != profile[alloc[i]])	{
			cerr << "error1\n";
			exit(1);
		}
		if (alloc[i] == -1)	{
			cerr << "error 2 \n";
			exit(1);
		}
		sitematrixarray[i] = new AACodonMutSelProfileSubMatrix(GetCodonStateSpace(),nucrr,nucstat,codonprofile,GetProfile(i),omega,true);
	}

};

#endif

