
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef MPIMODULE_H
#define MPIMODULE_H

#include <string>
using namespace std;

class MPIModule {

	public:

	void SetName(string inname)	{
		name = inname;
	}

	// protected:

	MPIModule() : sitemin(0), sitemax(0), myid(0), nprocs(0), nsite(0), fmin(0), fmax(1), bkfmin(0), bkfmax(1) {
		version = "2.0";
		name = "noname";
		sitesuffstat = 0;
	}

	virtual ~MPIModule()  {
		delete[] sitemin;
		sitemin = 0;
		delete[] sitemax;
		sitemax = 0;
		delete[] globalrank;
		// delete[] localrank;
	}

	void SetSiteSuffStat(int in)	{
		sitesuffstat = in;
	}

	string GetVersion() {return version;}

	string GetName() {return name;}

	// total number of sites, active or not, across processes
	int GetNsite()	{
		return nsite;
	}

	// total number of active sites across processes
	int GetNactiveSite();

	// min and max sites for this process
	int GetSiteMin()	{
		return sitemin[myid];
	}

	int GetSiteMax()	{
		return sitemax[myid];
	}

	// number of active sites for this process
	// int GetActiveSiteNumber();

	// total number of sites for this process
	int GetSiteNumber()	{
		return sitemax[myid] - sitemin[myid];
	}

	// min max and number of sites for a given process
	int GetProcSiteMin(int proc)	{
		return sitemin[proc];
	}

	int GetProcSiteMax(int proc)	{
		return sitemax[proc];
	}

	int GetProcSiteNumber(int proc)	{
		return sitemax[proc] - sitemin[proc];
	}

	// int GetProcActiveSiteNumber(int proc);

	int GetMaxSiteNumber()	{
		return maxwidth;
	}

	int GetNprocs() {
		return nprocs;
	}

	int GetMyid() {
		return myid;
	}

	// returns whether a given site is active 
	bool ActiveSite(int site)	{
		return ((globalrank[site] >= (fmin * nsite)) && (globalrank[site] < (fmax * nsite)));
	}

	bool NewlyActivated(int site)	{
		bool currentlyactive= ((globalrank[site] >= (fmin * nsite)) && (globalrank[site] < (fmax * nsite)));
		bool previouslyactive= ((globalrank[site] >= (bkfmin * nsite)) && (globalrank[site] < (bkfmax * nsite)));
		return (currentlyactive && (! previouslyactive));
	}

	void SetMinMax(double inmin, double inmax)	{
		bkfmin = fmin;
		bkfmax = fmax;
		fmin = inmin;
		fmax = inmax;
	}

	void RestoreMinMax()	{
		fmin = bkfmin;
		fmax = bkfmax;
	}

	void BackupMinMax()	{
		bkfmin = fmin;
		bkfmax = fmax;
	}

	// SMC
	void ResetNsite()	{
		fmin = fmax = 0;
		bkfmin = bkfmax = 0;
	}

	void IncrementNsite(int delta)	{
		if (fmin != 0)	{
			cerr << "error in IncrementNsite\n";
			exit(1);
		}
		double d = ((double) delta) / nsite;
		if (fmax + d > 1)	{
			SetMinMax(0,1);
		}
		else	{
			SetMinMax(0,fmax+d);
		}
	}

	void CreateMPI(int innsite);
	virtual void NonMPIReshuffleSites();
	virtual void GlobalReshuffleSites();
	virtual void SlaveReshuffleSites();
	virtual void GlobalWriteSiteRankToStream(ostream& os);
	virtual void GlobalReadSiteRankFromStream(istream& is);

	int myid;
	int nprocs;

	int* globalrank;
	// int* localrank;
	int* sitemin;
	int* sitemax;
	int nsite;

	int maxwidth;

	double fmin;
	double fmax;
	double bkfmin;
	double bkfmax;

	int sitesuffstat;

	string version;
	string name;

};


#endif

