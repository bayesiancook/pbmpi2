
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef GENPATHSSSITEMATMIXTUREPROFILE_H
#define GENPATHSSSITEMATMIXTUREPROFILE_H

#include "GeneralPathSuffStatMatrixProfileProcess.h"
#include "SiteMatrixMixtureProfileProcess.h"

// general subclass for all Matrix Mixture mixtures using generic sufficient statistics
class GeneralPathSuffStatSiteMatrixMixtureProfileProcess : public virtual SiteMatrixMixtureProfileProcess, public virtual GeneralPathSuffStatMatrixProfileProcess  {

	public:

	GeneralPathSuffStatSiteMatrixMixtureProfileProcess() : profilerootcount(0), profilepaircount(0) , profilewaitingtime(0) {}
	virtual ~GeneralPathSuffStatSiteMatrixMixtureProfileProcess() {}

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
		}
		// CreateMatrix(k);
		UpdateMatrix(k);
	}

	virtual void DeleteComponent(int k)	{
		// DeleteMatrix(k);
	}

	// necessary to keep track of componentwise sufficient statistics 
	// when updating the structure of the mixture
	void AddSite(int site, int cat);
	void RemoveSite(int site, int cat);

	virtual void UpdateModeProfileSuffStat() {}
	virtual void GlobalUpdateModeProfileSuffStat() {}
	virtual void SlaveUpdateModeProfileSuffStat() {}

	virtual double ProfileSuffStatLogProb(int cat);
	virtual double ProfileSuffStatLogProb();
	virtual void GlobalProfileSuffStatLogProb();
	virtual void SlaveProfileSuffStatLogProb();

	virtual double LogStatProb(int site, int cat);

	virtual double GlobalMoveProfile(double tuning = 1, int n = 1, int nrep = 1);
	virtual void SlaveMoveProfile();

	// componentwise
	map<int,int>* profilerootcount;
	map< pair<int,int>, int>* profilepaircount;
	map<int,double>* profilewaitingtime;

};

#endif

