
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef PARTGTRPROFILE_H
#define PARTGTRPROFILE_H

#include "GTRSubMatrix.h"
#include "MatrixProfileProcess.h"

// superclass for all GTR-like models
class PartitionGTRProfileProcess : public virtual MatrixProfileProcess {

	public:

	PartitionGTRProfileProcess() : rr(0), fixrr(false), rrtype("None") {}
	virtual ~PartitionGTRProfileProcess() {}

	int GetNrr()	{
		return Nrr;
	}

	const double* GetSiteRR(int site)	{
		return GetRR(partalloc[site]);
	}

	const double* GetRR(int part) {
		if (! rr)	{
			cerr << "error : getrr\n";
			exit(1);
		}
		return rr[part];
	}

	/*
	void SetRR(int part, const double* inrr)	{
		for (int i=0; i<Nrr; i++)	{
			rr[part][i] = inrr[i];
		}
	}
	*/

	double GetRRMean()	{
		double tot = 0;
		for (int part=0; part<Npart; part++)	{
			tot += GetRRMean(part);
		}
		tot /= Npart;
		return tot;
	}

	double GetRRMean(int part)	{
		double mean = 0;
		for (int i=0; i<Nrr; i++)	{
			mean += rr[part][i];
		}
		mean /= Nrr;
		return mean;
	}

	// file contains list of RR types for all partitions
	void SetRR(string file);
	void SetRR(string type, double* inrr);

	void SetFixRR(bool in)	{
		fixrr = in;
	}

	bool FixRR()	{
		return fixrr;
	}

	double GetRRVarCoeff()	{
		double tot = 0;
		for (int part=0; part<Npart; part++)	{
			tot += GetRRVarCoeff(part);
		}
		tot /= Npart;
		return tot;
	}

	double GetRRVarCoeff(int part)	{
		double mean = 0;
		double var = 0;
		for (int i=0; i<Nrr; i++)	{
			mean += rr[part][i];
			var += rr[part][i] * rr[part][i];
		}
		mean /= Nrr;
		var /= Nrr;
		var -= mean*mean;
		var /= mean * mean;
		return var;
	}

	double GetRREntropy()	{
		double tot = 0;
		for (int part=0; part<Npart; part++)	{
			tot += GetRREntropy(part);
		}
		tot /= Npart;
		return tot;
	}

	double GetRREntropy(int part)	{
		double total = 0;
		for (int i=0; i<Nrr; i++)	{
			total += rr[part][i];
		}
		double ent = 0;
		for (int i=0; i<Nrr; i++)	{
			double tmp = rr[part][i] / total;
			if (tmp > 1e-6)	{
				ent -= tmp * log(tmp);
			}
		}
		return ent;
	}

	static int rrindex(int i, int j, int nstate)	{
		return (i<j) ? (2 * nstate - i - 1) * i / 2 + j - i - 1 : (2 * nstate - j - 1) * j / 2 + i - j - 1 ;
	}

	// protected:

	virtual void Create();
	virtual void Delete();

	// relative rates
	virtual double LogRRPrior();
	virtual double LogRRPrior(int part);
	virtual void SampleRR();
	virtual void SampleGlobalParameters()	{
		SampleRR();
	}
	virtual void PriorSampleGlobalParameters()	{
		SampleRR();
	}

	// assumes that site-specific sufficient statistics are already updated
	// collect them into more compact sufficient statistics
	// (depending on what is permitted by the type of sufficient satistics used)
	virtual void UpdateRRSuffStat() = 0;
	virtual void GlobalUpdateRRSuffStat() = 0;
	virtual void SlaveUpdateRRSuffStat() = 0;

	// Metropolis or Gibbs Sampling algorithm,
	// (depending on what is permitted by the type of sufficient satistics used)
	virtual double MoveRR(double tuning, int n, int nrep);
	virtual double MoveRR(int part, double tuning, int n, int nrep);

	virtual double MoveRR()	{
		MoveRR(1.0,1,Nrr);
		MoveRR(0.1,1,Nrr);
		return 1;
	}

	virtual double GlobalParametersMove();

	int Nrr;
	double** rr;
	double* allocrr;

	bool fixrr;
	string rrtype;

};

#endif

