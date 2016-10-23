
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef AACODONMUTSELSITEMATMIX_H
#define AACODONMUTSELSITEMATMIX_H

#include "AACodonMutSelProfileProcess.h"
#include "GeneralPathSuffStatSiteMatrixMixtureProfileProcess.h"

class AACodonMutSelSiteMatrixMixtureProfileProcess : public virtual AACodonMutSelProfileProcess, public virtual GeneralPathSuffStatSiteMatrixMixtureProfileProcess	{

	public:

	AACodonMutSelSiteMatrixMixtureProfileProcess() : omega0(0) {}
	virtual ~AACodonMutSelSiteMatrixMixtureProfileProcess() {}

	AACodonMutSelProfileSubMatrix* GetComponentCodonMatrix(int k)	{

		AACodonMutSelProfileSubMatrix* codonmatrix = dynamic_cast<AACodonMutSelProfileSubMatrix*>(matrixarray[k]);
		if (! codonmatrix)	{
			cerr << "error in AACodonMutSelProfileSubMatrix::GetComponentCodonMatrix: null matrix\n";
			cerr << matrixarray[k] << '\n';
			exit(1);
		}
		return codonmatrix;
	}

	protected:

	virtual void Create();
	virtual void Delete();

	virtual void CreateComponent(int k)	{
		occupancy[k] = 0;
		SampleStat(k);
		// useful?
		if (activesuffstat)	{
			profilepaircount[k].clear();
			profilerootcount[k].clear();
			profilewaitingtime[k].clear();
			profilenonsynwaitingtime[k].clear();
		}
		// CreateMatrix(k);
		UpdateMatrix(k);
	}

	virtual void DeleteComponent(int k)	{
		// DeleteMatrix(k);
	}

	virtual void UpdateModeProfileSuffStat();
	virtual void GlobalUpdateModeProfileSuffStat();
	virtual void SlaveUpdateModeProfileSuffStat();
	virtual double LogStatProb(int site, int cat);

	virtual void CreateMatrix(int k)	{
		if (matrixarray[k])	{
			cerr << "error in AACodonMutSelSiteMatrixMixtureProfileProcess: matrixarray is not 0\n";
			exit(1);
		}
		matrixarray[k] = new AACodonMutSelProfileSubMatrix(GetCodonStateSpace(),nucrr,nucstat,codonprofile,profile[k],omega0,true);
	}

	virtual void CreateSiteMatrix(int i)	{
		cerr << "in base create site matrix\n";
		if (sitematrixarray[i])	{
			cerr << "error in AACodonMutSelSiteMatrixMixtureProfileProcess: sitematrixarray is not 0\n";
			exit(1);
		}
		sitematrixarray[i] = new AACodonMutSelProfileSubMatrix(GetCodonStateSpace(),nucrr,nucstat,codonprofile,GetProfile(i),GetSiteOmegaPtr(i),true);
	}

	double* omega0;
	// componentwise
	map<int,double>* profilenonsynwaitingtime;

};

#endif

