
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef EXPCONPARTGTRPHYLO_H
#define EXPCONPARTGTRPHYLO_H

#include "PhyloProcess.h"
#include "ExpoConjugatePartitionGTRSubstitutionProcess.h"
#include "RateProcess.h"

// this is the Exponential Conjugate version
class ExpoConjugatePartitionGTRPhyloProcess : public virtual PhyloProcess, public virtual ExpoConjugatePartitionGTRSubstitutionProcess, public virtual RateProcess {

	public:

	ExpoConjugatePartitionGTRPhyloProcess() : siteprofilesuffstatcount(0), siteprofilesuffstatbeta(0), allocsiteprofilesuffstatcount(0), allocsiteprofilesuffstatbeta(0) {}
	virtual ~ExpoConjugatePartitionGTRPhyloProcess() {}

	protected:

	virtual void Create()	{
		RateProcess::Create();
		ExpoConjugatePartitionGTRSubstitutionProcess::Create();
		PhyloProcess::Create();
	}

	virtual void Delete()	{
		PhyloProcess::Delete();
		ExpoConjugatePartitionGTRSubstitutionProcess::Delete();
		RateProcess::Delete();
	}

	const int* GetSiteProfileSuffStatCount(int site) {return siteprofilesuffstatcount[site];}
	const double* GetSiteProfileSuffStatBeta(int site) {return siteprofilesuffstatbeta[site];}

	// protected:
	public:


	void CreateSuffStat();
	void DeleteSuffStat();

	void GlobalUpdateSiteProfileSuffStat();
	virtual void SlaveUpdateSiteProfileSuffStat();

	void UpdateSiteRateSuffStat();
	void UpdateBranchLengthSuffStat();
	void UpdateRRSuffStat();
	void UpdateSiteProfileSuffStat();

	int GlobalCountMapping();
	int CountMapping();
	int CountMapping(int site);

	int** siteprofilesuffstatcount;
	double** siteprofilesuffstatbeta;

	int* allocsiteprofilesuffstatcount;
	double* allocsiteprofilesuffstatbeta;

	double LengthRelRateMove(double tuning, int nrep);
};

#endif

