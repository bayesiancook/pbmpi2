
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef AASSMIXTUREPROFILE_H
#define AASSMIXTUREPROFILE_H


#include "GTRProfileProcess.h"
#include "AASubSelSubMatrix.h"
#include "MatrixMixtureProfileProcess.h"

// general superclass for GTR-like Dirichlet-process mixture on profiles
class AASubSelMixtureProfileProcess : public virtual GTRMixtureProfileProcess	{

	static CONST double TOOSMALL = 1e-30;

	public:

	AASubSelMixtureProfileProcess() {}
	virtual ~AASubSelMixtureProfileProcess() {}

	protected:

	// simply creates/deletes GTR matrices for all currently existing components
	void CreateMatrix(int k)	{
		if (matrixarray[k])	{
			cerr << "error in gen path suff stat gtr dp profile process: matrixarray is not 0\n";
			cerr << matrixarray[k]->GetNstate() << '\n';
			exit(1);
		}
		matrixarray[k] = new AASubSelSubMatrix(GetDim(),rr,profile[k],false);
	}

	AASubSelSubMatrix* GetAASubSelMatrix(int k)	{
		AASubSelSubMatrix* tmp = dynamic_cast<AASubSelSubMatrix*>(matrixarray[k]);
		if (!tmp)	{
			cerr << "error in GetAASubSelMatrix: null matrix \n";
			exit(1);
		}
		return tmp;
	}

	double GetNormRate(int k)	{
		double tot = 0;
		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				if (i != j)	{
					double tmp = profile[k][i] * rr[rrindex(i,j,GetDim())];
					double deltaF = log(profile[k][j] / profile[k][i]);
					if (fabs(deltaF) < TOOSMALL)        {
						tmp /= ( 1.0 - (deltaF / 2) );
					}
					else    {
						tmp *=  (deltaF)/(1.0 - exp(-deltaF));
					}
					tot += tmp;
				}
			}
		}
		return tot;
	}

	virtual double GetNormalizationFactor()	{
		UpdateOccupancyNumbers();
		double norm = 0;
		int tot = 0;
		for (int k=0; k<GetNcomponent(); k++)	{
			if (occupancy[k])	{
				norm += (occupancy[k] + 1) * GetNormRate(k);
				tot += occupancy[k] + 1;
			}
		}
		if (!tot)	{
			return 1;
		}
		norm /= tot;
		return norm;
	}

};

#endif

