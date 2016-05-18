
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef AACODONMUTSEL_H
#define AACODONMUTSEL_H

#include "OmegaProcess.h"
#include "CodonSubMatrix.h"
#include "CodonSequenceAlignment.h"
#include "GeneralPathSuffStatMatrixProfileProcess.h"

class AACodonMutSelProfileProcess : public virtual GeneralPathSuffStatMatrixProfileProcess, public virtual OmegaProcess	{

	// y mettre les variables globales (taux de mutation essentiellement)

	// s'inspirer de GTRProfileProcess et GeneralPathSuffStatGTRProfileProcess

	public:

	AACodonMutSelProfileProcess() : nucrr(0), nucstat(0), codonprofile(0) {}
	virtual ~AACodonMutSelProfileProcess() {}

	int GetNnucrr()	{
		return Nnucrr;
	}	

	int GetNcodon()	{
		return GetCodonStateSpace()->GetNstate();
	}

	const double* GetNucRR()	{
		if (! nucrr)	{
			cerr << "error : getnucrr\n";
			exit(1);
		}
		return nucrr;
	}

	const double GetNucRR(int i)	{
		if ( (! nucrr) || (i>GetNnucrr()) || (i<0) )	{
			cerr << "error : getnucrr i\n";
			exit(1);
		}
		return nucrr[i];
	}

	const double* GetNucStat()	{
		if (! nucstat)	{
			cerr << "error : getnucstat\n";
			exit(1);
		}
		return nucstat;
	}

	const double GetNucStat(int i)	{
		if ( (! nucstat) || (i>Nnuc) || (i<0) )	{
			cerr << "error : getnucstat i\n";
			exit(1);
		}
		return nucstat[i];
	}

	const double GetCodonProfileEntry(int i)	{
		return codonprofile[i];	
	}

	double GetCodonProfileEntropy();

	/*
	virtual void UpdateNucStatSuffStat() = 0;
	virtual void UpdateNucRRSuffStat() = 0;

	double NucStatSuffStatLogProb();
	double NucRRSuffStatLogProb();
	*/

	virtual void UpdateSiteOmegaSuffStat()	{
		for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
			if (ActiveSite(i))	{
	
				map<pair<int,int>, int>& paircount = GetSitePairCount(i);
				map<int,double>& waitingtime = GetSiteWaitingTime(i);
				CodonSubMatrix* codonmatrix = dynamic_cast<CodonSubMatrix*>(GetMatrix(i));
				for (map<int,double>::iterator j = waitingtime.begin(); j!= waitingtime.end(); j++)	{
					siteomegasuffstatbeta[i] += j->second * codonmatrix->RateAwayNonsyn(j->first) / GetSiteOmega(i);
				}

				for (map<pair<int,int>, int>::iterator j = paircount.begin(); j!= paircount.end(); j++)	{
					if (! codonmatrix->Synonymous(j->first.first,j->first.second) )	{
						siteomegasuffstatcount[i] += j->second;
					}
				}
			}
		}
	}

	protected:

	CodonSequenceAlignment* GetCodonData()	{

		CodonSequenceAlignment* tmp = dynamic_cast<CodonSequenceAlignment*>(GetData());
		if (!tmp)	{
			if (GetData())	{
				cerr << "in AACodonMutSelProfileProcess: cast error on sequence alignment\n";
				exit(1);
			}
		}
		return tmp;
	}

	CodonStateSpace* GetCodonStateSpace()	{
		return GetCodonData()->GetCodonStateSpace();
	}

	virtual void Create();
	virtual void Delete();
	virtual double GetNormalizationFactor()	{return 1.0;}
	
	//double GetMinTotWeight() {return GetDim() / (GetDim()/2);}

	// nuc relative rates
	virtual double LogNucRRPrior();
	virtual void SampleNucRR();

	// nuc stationaries
	virtual double LogNucStatPrior();
	virtual void SampleNucStat();

	// codon profile
	virtual double LogCodonProfilePrior();
	virtual void SampleCodonProfile();

	virtual double GlobalParametersMove();

	// moves on global parameters
	double MoveNucRR(double tuning); 
	double MoveNucRR(double tuning, int n); 
	double MoveNucStat(double tuning, int n);
	double MoveCodonProfile(double tuning, int n, int nrep=1);
	double MoveNucStatCodonProfile(double tuning, int n, int nrep=1);
	
	int Nnucrr;
	double* nucrr;
	double* nucstat;
	double* codonprofile;

	/*
	int* nucstatsuffstatcount;
	double* nucstatsuffstatbeta;
	int* nucrrsuffstatcount;
	double* nucrrsuffstatbeta;
	*/

	int fixcodonprofile;

};

#endif
