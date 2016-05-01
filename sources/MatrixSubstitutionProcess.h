
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef MATRIXSUB_H
#define MATRIXSUB_H

#include "SubstitutionProcess.h"
#include "MatrixProfileProcess.h"

class MatrixSubstitutionProcess : public virtual SubstitutionProcess, public virtual MatrixProfileProcess	{

	using SubstitutionProcess::GetNstate;

	public:

	MatrixSubstitutionProcess() {}
	virtual ~MatrixSubstitutionProcess() {}

	// virtual int GetNstate(int site)	{return GetMatrix(site)->GetNstate();}

	virtual const double* GetStationary(int site)	{
		return GetMatrix(site)->GetStationary();
	}

	protected:

	// CPU Level 3: implementations of likelihood propagation and substitution mapping methods
	void SitePropagate(int site, double** from, double** to, double time, bool condalloc = false);
	void Propagate(double*** from, double*** to, double time, bool condalloc);
	BranchSitePath* SampleSitePath(int site, int stateup, int statedown, double time);
	BranchSitePath* SampleRootSitePath(int site, int rootstate);
	BranchSitePath* ResampleAcceptReject(int maxtrial, int stateup, int statedown, double rate, double totaltime, SubMatrix* matrix);
	BranchSitePath* ResampleUniformized(int stateup, int statedown, double rate, double totaltime, SubMatrix* matrix);

};

#endif

