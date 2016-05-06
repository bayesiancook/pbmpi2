
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include <cassert>
#include "Parallel.h"
#include <string.h>

#include "MultiGeneMPIModule.h"
#include "PhyloProcess.h"

void MultiGeneMPIModule::Create()	{

	if (! genelnL)	{
		genelnL = new double[Ngene];
		tmpgenelnL = new double[Ngene];
		process = new PhyloProcess*[Ngene];
		for (int gene=0; gene<Ngene; gene++)	{
			process[gene] = 0;
		}
	}
}

void MultiGeneMPIModule::Delete()	{

	if (genelnL)	{
		delete[] genelnL;
		delete[] tmpgenelnL;
		genelnL = 0;
		for (int gene=0; gene<Ngene; gene++)	{
			delete process[gene];
			process[gene] = 0;
		}
		delete[] process;
	}
}
