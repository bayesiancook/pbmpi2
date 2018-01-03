
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "GeneralPathSuffStatMatrixSiteSpecificProfileProcess.h"
#include "Random.h"
#include "Parallel.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* GeneralPathSuffStatMatrixSiteSpecificProfileProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

double GeneralPathSuffStatMatrixSiteSpecificProfileProcess::ProfileSuffStatLogProb(int site)	{

	double total = 0;
	SubMatrix* mat = matrixarray[site];
	if (! mat)	{
		cerr << "error : null matrix\n";
		exit(1);
	}
	map<pair<int,int>, int>& paircount = GetSitePairCount(site);
	map<int,double>& waitingtime = GetSiteWaitingTime(site);
	int rootstate = GetSiteRootState(site);
	const double* stat = matrixarray[site]->GetStationary();
	if (rootstate != -1)	{
		total += log(stat[rootstate]);
		for (map<int,double>::iterator i = waitingtime.begin(); i!= waitingtime.end(); i++)	{
			total += i->second * (*mat)(i->first,i->first);
		}
		for (map<pair<int,int>, int>::iterator i = paircount.begin(); i!= paircount.end(); i++)	{
			total += i->second * log((*mat)(i->first.first, i->first.second));
		}
	}
	return total;
}

