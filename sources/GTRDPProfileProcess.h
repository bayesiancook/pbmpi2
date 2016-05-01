
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef GTRDPPROFILE_H
#define GTRDPPROFILE_H

#include "DPProfileProcess.h"
#include "GTRMixtureProfileProcess.h"

class GTRDPProfileProcess : public virtual DPProfileProcess, public virtual GTRMixtureProfileProcess {

	public:

	GTRDPProfileProcess() {}
	virtual ~GTRDPProfileProcess() {}

	void ToStream(ostream& os);
	void FromStream(istream& is);

	protected:

	virtual void Create()	{
		DPProfileProcess::Create();
		GTRMixtureProfileProcess::Create();
	}

	virtual void Delete()	{
		DPProfileProcess::Delete();
		GTRMixtureProfileProcess::Delete();
	}

};

#endif

