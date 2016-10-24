
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef AACODONMUTSELSBDPPROFILE_H
#define AACODONMUTSELSBDPPROFILE_H

#include "SBDPProfileProcess.h"
#include "OmegaProcess.h"
#include "AACodonMutSelProfileProcess.h"
#include "GeneralPathSuffStatMatrixMixtureProfileProcess.h"

class AACodonMutSelSBDPProfileProcess : public virtual SBDPProfileProcess, public virtual SingleOmegaProcess, public virtual AACodonMutSelProfileProcess, public virtual GeneralPathSuffStatMatrixMixtureProfileProcess	{

	// implementer les fonctions create matrix et delete matrix
	// ainsi que CreateComponent(int k) and DeleteComponent(k)

	// s'inspirer de GeneralPathSuffStatGTRFiniteProfileProcess

	public:

	AACodonMutSelSBDPProfileProcess() {}
	virtual ~AACodonMutSelSBDPProfileProcess() {}

	/*
	// still to be implemented
	virtual void UpdateNucStatSuffStat();
	virtual void UpdateNucRRSuffStat();
	virtual void UpdateOmegaSuffStat(); done
	*/

	///*
	using ProfileProcess::ProfileSuffStatLogProb;
	using MixtureProfileProcess::BetaProfileSuffStatLogProb;
	using MixtureProfileProcess::CountProfileSuffStatLogProb;
	double OmegaSuffStatLogProb()	{
		return ProfileSuffStatLogProb();
	}
	//*/

	/*
	void CheckSuffStatLogProb()	{

		double diff1 = OmegaSuffStatLogProb();
		UpdateMatrices();
		double diff2 = ProfileSuffStatLogProb();
		*omega /= 10;
		diff1 -= OmegaSuffStatLogProb();
		UpdateMatrices();
		diff2 -= ProfileSuffStatLogProb();
		if (fabs(diff1 - diff2) > 1e-6)	{
			cerr << "error in check suff stat log prob\n";
			cerr << diff1 - diff2 << '\t' << diff1 << '\t' << diff2 << '\t' << omegasuffstatcount << '\t' << omegasuffstatbeta << '\n';
			exit(1);
		}
		*omega *= 10;
		UpdateMatrices();
	}
	*/

	protected:

	void Create()	{
		SingleOmegaProcess::Create();
		AACodonMutSelProfileProcess::Create();
		SBDPProfileProcess::Create();
		GeneralPathSuffStatMatrixMixtureProfileProcess::Create();
	}
	
	void Delete()	{
		GeneralPathSuffStatMatrixMixtureProfileProcess::Delete();
		SBDPProfileProcess::Delete();
		AACodonMutSelProfileProcess::Delete();
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


	void CreateMatrix(int k)	{
		if (matrixarray[k])	{
			cerr << "error in AACodonMutSelSBDPProfileProcess: matrixarray is not 0\n";
			exit(1);
		}
		matrixarray[k] = new AACodonMutSelProfileSubMatrix(GetCodonStateSpace(),nucrr,nucstat,codonprofile,profile[k],omega,true);
	}

};

#endif

