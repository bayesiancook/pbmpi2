/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef SBDPOMEGA_H
#define SBDPOMEGA_H

#include "OmegaProcess.h"

class SBDPOmegaProcess : public virtual OmegaProcess	{

	public:
	SBDPOmegaProcess() : omegaarray(0), omegakappa(1), omegamovekappa(true), omegakappaprior(0), omegaV(0), maxweighterror(0) {}
	virtual ~SBDPOmegaProcess() {}

	protected:

	double GetOmegaKapp()	{
		return omegakappa;
	}
	
	int GetNDisplayedOmegaComponent()	{
		return GetNOccupiedOmegaComponent();
	}

	int GetNOccupiedOmegaComponent()	{
		int n = 0;
		for (int k=0; k<NcomponentOmega; k++)	{
			if (omegaoccupancy[k])	{
				n++;
			}
		}
		return n;
	}
		
	double GetSiteOmega(int site)	{
		return omegaarray[omegaalloc[site]];
	}

	int GetOmegaAlloc(int site)	{
		return omegaalloc[site];
	}

	double GetMeanOmega()	{
		double mean = 0;
		for (int i=0; i<GetNsite(); i++)	{
			mean += omegaarray[omegaalloc[i]];
		}
		mean /= GetNsite();
		return mean;
	}

	int GetNcomponentOmega()	{
		return NcomponentOmega;
	}

	double GetProportionOmegaGreaterThan(double c)	{
		double tot = 0;
		for (int i=0; i<GetNsite(); i++)	{
			if (omegaarray[omegaalloc[i]] > c)	{
				tot++;
			}
		}
		tot /= GetNsite();
		return tot;
	}

	double OmegaSiteSuffStatLogProb(int site)	{
		return siteomegasuffstatcount[site] * log(omegaarray[GetOmegaAlloc(site)]) - siteomegasuffstatbeta[site] * omegaarray[GetOmegaAlloc(site)];
	}

	double OmegaCompSuffStatLogProb(int k)	{
		double tot = 0;
		for (int i=0; i<GetNsite(); i++)	{
			if (ActiveSite(i))	{
				if (GetOmegaAlloc(i) == k)	{
					tot += 0; // not yet sure what to do here, or if relevant at all.
				}
			}	
		}	
	}

	virtual double OmegaSuffStatLogProb()	{
		double total = 0;
		for (int k=0; k<GetNcomponentOmega(); k++)	{
			if (omegaoccupancy[k])	{
				total += OmegaCompSuffStatLogProb(k);
			}
		}
		return total;
	}

	// omega
	double LogOmegaPrior(int component)	{
		return omegaalpha*log(omegabeta) - rnd::GetRandom().logGamma(omegaalpha) + (omegaalpha-1)*log(omegaarray[component]) - omegabeta*omegaarray[component];
	}

	virtual double LogOmegaPrior()	{
		double total = 0;
		for (int i=0; i<GetNcomponentOmega(); i++)	{
			if (omegaoccupancy[i])	{
				total += LogOmegaPrior(i);
			}
		}
		return total;
	}

	double LogOmegaHyperPrior()	{
		return -omegaalpha - omegabeta;
	}

	virtual void SampleOmega()	{
		omegaalpha = omegabeta = 1.0;
		for (int i=0; i<GetNcomponentOmega(); i++)	{
			omegaarray[i] = rnd::GetRandom().Gamma(omegaalpha,omegabeta);
		}
	}

	virtual int GetNomegaMax() {
		return GetNsite() > nomegamax ? nomegamax : GetNsite();
	}
	virtual void SetNomegaMax(int n) {nomegamax = n;}

	double MoveOmega(double tuning, int k);
	double MoveCompOmega(double tuning, int k);
	double MPIMoveOmega(double tuning);
	double NonMPIMoveOmega(double tuning);

	
	double MixMoveOmega(int nmix, double tuning, int nsitenrep, int nhyperrep);
	double GlobalMixMoveOmega(int nmix, double tuning, int nsitenrep, int nhyperrep);
	void SlaveMixMoveOmega();

	void GlobalUpdateOmegaSuffStat();
	void SlaveUpdateOmegaSuffStat();

	double MoveOmegas(double tuning, int nrep);
	double MoveOmegaHyper(double tuning, int nrep);
	double MoveOmegaAlpha(double tuning);
	double MoveOmegaBeta(double tuning);

	protected:

	virtual void Create()	{
		if (! omegaarray)	{
			OmegaProcess::Create();
			omegaarray = new double[GetNsite()];
		}
	}

	virtual void Delete()	{
		if (omegaarray)	{
			delete omegaarray;
			omegaarray = 0;
			OmegaProcess::Delete();
		}
	}
	
	void UpdateOmegaSuffStat();
	//void GlobalUpdateSiteOmegaSuffStat();
	//void SlaveUpdateSiteOmegaSuffStat();

	void ResampleOmegas();
	void GlobalCollectSiteOmegaSuffStats();
	void SlaveCollectSiteOmegaSuffStats();

	int NcomponentOmega;
	int* omegaoccupancy;
	int nomegamax;

	double* omegaarray;
	int* omegaalloc;
	double omegakappa;
	double omegaalpha;
	double omegabeta;
	double* omegaV;
	double* omegaweight;
	
	bool omegamovekappa;
	int omegakappaprior;

	double maxweighterror;
};

#endif

