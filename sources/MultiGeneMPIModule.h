

/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef MULTIGENEMPIMODULE_H
#define MULTIGENEMPIMODULE_H

#include "PhyloProcess.h"

class MultiGeneMPIModule	{

	public:

	MultiGeneMPIModule() : Ngene(0), genelnL(0)	{}
	virtual ~MultiGeneMPIModule() {}

	virtual void Create();
	virtual void Delete();

	protected:

	PhyloProcess** process;
	double* genelnL;
	double* tmpgenelnL;
	int Ngene;
	int* genealloc;
	int* genesize;
	string* genename;

	int GlobalNsite;
	int* globalnsite;

};

#endif
