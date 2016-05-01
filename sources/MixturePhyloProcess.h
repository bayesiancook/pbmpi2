
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef MIXPHYLOPROCESS_H
#define MIXPHYLOPROCESS_H

#include "PhyloProcess.h"
#include "MixtureProfileProcess.h"

class MixturePhyloProcess : public virtual PhyloProcess, public virtual MixtureProfileProcess	{

	public:

	MixturePhyloProcess() {}
	virtual ~MixturePhyloProcess() {}

	// protected:

	virtual void SlaveComputeCVScore();
	virtual void SlaveComputeSiteLogL();

	virtual void Create()	{
		PhyloProcess::Create();
		MixtureProfileProcess::Create();
	}

	virtual void Delete()	{
		MixtureProfileProcess::Delete();
		PhyloProcess::Delete();
	}
};


#endif

