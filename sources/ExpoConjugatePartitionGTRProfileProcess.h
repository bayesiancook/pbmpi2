
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef EXPCONPARTGTRPROFILE_H
#define EXPCONPARTGTRPROFILE_H

#include "PartitionGTRProfileProcess.h"

// superclass for GTR-like models using the exponential conjugate relation between relative rates and substitution processes
class ExpoConjugatePartitionGTRProfileProcess : public virtual PartitionGTRProfileProcess {

	public:

	ExpoConjugatePartitionGTRProfileProcess() : rrsuffstatcount(0), rrsuffstatbeta(0) {}
	virtual ~ExpoConjugatePartitionGTRProfileProcess() {}

	const int* GetRRSuffStatCount(int part) {return rrsuffstatcount[part];}
	const double* GetRRSuffStatBeta(int part) {return rrsuffstatbeta[part];}
	// protected:

	virtual void Create();
	virtual void Delete();

	// profiles
	virtual const int* GetSiteProfileSuffStatCount(int site) = 0;
	virtual const double* GetSiteProfileSuffStatBeta(int site) = 0;

	// update of relative rates
	// conjugate Gibbs resampling
	virtual double MoveRR();

	virtual void GlobalUpdateRRSuffStat();
	virtual void SlaveUpdateRRSuffStat();

	int* allocrrsuffstatcount;
	double* allocrrsuffstatbeta;
	int** rrsuffstatcount;
	double** rrsuffstatbeta;
};

#endif

