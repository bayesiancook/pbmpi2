
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef ZIPSUB_H
#define ZIPSUB_H

#include "SubMatrix.h"

class ZipSubMatrix : public virtual SubMatrix	{

	public:

	ZipSubMatrix(SubMatrix* infrom, int inzipsize, int* inzipindices) : SubMatrix(inzipsize, infrom->isNormalised()) {
		from = infrom;
		zipindices = inzipindices;

	}
	~ZipSubMatrix() {}

	void SetFrom(SubMatrix* newfrom)	{
		from = newfrom;
	}

	SubMatrix* GetFrom()	{
		return from;
	}

	protected:

	// virtual void ComputeFullArray();
	virtual void ComputeArray(int state);
	virtual void ComputeStationary();

	SubMatrix* from;
	int* zipindices;
};

#endif

