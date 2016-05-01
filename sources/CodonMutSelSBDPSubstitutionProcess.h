
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef CODONMUTSELSBDPSUB_H
#define CODONMUTSELSBDPSUB_H

#include "CodonMutSelSBDPProfileProcess.h"
#include "UniformRateProcess.h"
#include "GeneralPathSuffStatMatrixSubstitutionProcess.h"


class CodonMutSelSBDPSubstitutionProcess : public virtual CodonMutSelSBDPProfileProcess, public virtual UniformRateProcess, public virtual GeneralPathSuffStatMatrixSubstitutionProcess {

	// s'inspirer de GeneralPathSuffStatGTRSubstitutionProcess
	// et GeneralPathSuffStatRASCATGTRSubstitutionProcess

	public:

	CodonMutSelSBDPSubstitutionProcess() {}
	virtual ~CodonMutSelSBDPSubstitutionProcess() {}

	virtual int GetNstate() {return GetNcodon();}

	protected:

	void Create()	{
		CodonMutSelSBDPProfileProcess::Create();
		UniformRateProcess::Create();
		GeneralPathSuffStatMatrixSubstitutionProcess::Create();
	}

	void Delete()	{
		GeneralPathSuffStatMatrixSubstitutionProcess::Delete();
		UniformRateProcess::Delete();
		CodonMutSelSBDPProfileProcess::Delete();
	}

};

#endif

