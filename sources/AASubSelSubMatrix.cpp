
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "AASubSelSubMatrix.h"
#include <iostream>
using namespace std;

// ---------------------------------------------------------------------------
//		 AASubSelSubMatrix
// ---------------------------------------------------------------------------


AASubSelSubMatrix::AASubSelSubMatrix(int inNstate, const double* rr, const double* stat, bool innormalise) : SubMatrix(inNstate, innormalise)	{

	Nrr = Nstate * (Nstate-1) / 2;
	mRelativeRate = rr;
	ExternalStat = stat;
	CorruptMatrix();
}

void AASubSelSubMatrix::ComputeStationary()	{
	double total = 0;
	for (int k=0; k<Nstate; k++)	{
		mStationary[k] = ExternalStat[k];
	}
}

// ---------------------------------------------------------------------------
//		 ComputeArray
// ---------------------------------------------------------------------------

void	AASubSelSubMatrix::ComputeArray(int i)	{

	if (mRelativeRate)	{
		double total = 0;
		for (int j=0; j<Nstate; j++)	{
			if (i!=j)	{
				Q[i][j] = RelativeRate(i,j);;
				double deltaF = log(mStationary[j] / mStationary[i]);
				if (fabs(deltaF) < TOOSMALL)        {
					Q[i][j] /= ( 1.0 - (deltaF / 2) );
				}
				else    {
					Q[i][j] *=  (deltaF)/(1.0 - exp(-deltaF));
				}
				total += Q[i][j];
			}
		}

		// should always ensure that the diagonal entry of the matrix Q[i][i] is such that 
		// the sum over all entries of the row is equal to 0
		Q[i][i] = - total;
	}
	else	{
		for (int j=0; j<Nstate; j++)	{
			Q[i][j] = mStationary[j];
		}
		Q[i][i] -= 1;
	}
}

