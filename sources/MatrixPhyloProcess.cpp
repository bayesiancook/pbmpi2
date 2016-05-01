
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "MatrixPhyloProcess.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* Matrix PhyloProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------


void MatrixPhyloProcess::Unfold()	{

	DeleteMappings();
	ActivateSumOverRateAllocations();
	// CreateMatrices();
	UpdateMatrices();
	// CreateConditionalLikelihoods();
	UpdateConditionalLikelihoods();
	if (!sumratealloc)	{
		DrawAllocations(0);
		InactivateSumOverRateAllocations(ratealloc);
	}
	activesuffstat = false;
}

void MatrixPhyloProcess::Collapse()	{

	// UpdateConditionalLikelihoods();
	if (sumratealloc)	{
		DrawAllocations(0);
		InactivateSumOverRateAllocations(ratealloc);
	}
	SampleNodeStates();
	// DeleteConditionalLikelihoods();
	SampleSubstitutionMappings(GetRoot());
	// DeleteMatrices();
	activesuffstat = true;
}

