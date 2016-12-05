/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "SiteOmegaProcess.h"
#include "Parallel.h"

void SiteOmegaProcess::UpdateOmegaSuffStat()	{
	// UpdateSiteOmegaSuffStat();
}

void SiteOmegaProcess::GlobalUpdateOmegaSuffStat() {

	if (GetNprocs() > 1)	{
		MPI_Status stat;
		MESSAGE signal = UPDATE_OMEGA;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	}
	else	{
		UpdateOmegaSuffStat();
	}
}

void SiteOmegaProcess::SlaveUpdateOmegaSuffStat() {
	UpdateOmegaSuffStat();
}

double SiteOmegaProcess::MoveSiteOmega(double tuning, int site)	{

	double bkomega = omegaarray[site];
	double deltalogprob = -LogOmegaPrior(site) - OmegaSuffStatLogProb(site);

	double h = tuning * (rnd::GetRandom().Uniform() -0.5);
	double e = exp(h);
	omegaarray[site] *= e;

	deltalogprob += h;
	deltalogprob += LogOmegaPrior(site) + OmegaSuffStatLogProb(site);

	int accepted = (rnd::GetRandom().Uniform() < exp(deltalogprob));
	if (! accepted)	{
		omegaarray[site] = bkomega;
	}
	return accepted;	
}

double SiteOmegaProcess::MoveSiteOmegas(double tuning, int nrep)	{

	double nacc = 0;
	double ntot = 0;
	for (int rep=0; rep<nrep; rep++)	{
		for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
			if (ActiveSite(i))	{
				nacc += MoveSiteOmega(tuning,i);
				ntot++;
			}
		}
	}
	return nacc/ntot;
}

double SiteOmegaProcess::MoveOmegaHyper(double tuning, int nrep)	{

	double nacc = 0;
	double ntot = 0;
	for (int rep=0; rep<nrep; rep++)	{
		nacc += MoveOmegaAlpha(tuning);
		ntot++;
		nacc += MoveOmegaBeta(tuning);
		ntot++;
	}
	return nacc/ntot;
}

double SiteOmegaProcess::MoveOmegaAlpha(double tuning)	{

	double bkomegaalpha = omegaalpha;
	double deltalogprob = - LogOmegaPrior() - LogOmegaHyperPrior();
	double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
	double e = exp(m);
	omegaalpha *= e;
	deltalogprob += LogOmegaPrior() + LogOmegaHyperPrior();
	deltalogprob += m;
	int accepted = (log(rnd::GetRandom().Uniform()) < deltalogprob);

	if (!accepted)	{
		omegaalpha = bkomegaalpha;
	}
	return accepted;
}

double SiteOmegaProcess::MoveOmegaBeta(double tuning)	{

	double bkomegabeta = omegabeta;
	double deltalogprob = - LogOmegaPrior() - LogOmegaHyperPrior();
	double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
	double e = exp(m);
	omegabeta *= e;
	deltalogprob += LogOmegaPrior() + LogOmegaHyperPrior();
	deltalogprob += m;
	int accepted = (log(rnd::GetRandom().Uniform()) < deltalogprob);

	if (!accepted)	{
		omegabeta = bkomegabeta;
	}
	return accepted;
}

double SiteOmegaProcess::MoveOmega(double tuning)	{

	double ret = 0;
	if (integrated)	{

		if (GetNprocs() > 1)	{
			GlobalCollectSiteOmegaSuffStats();
		}

		MoveOmegaHyperIntegrated(tuning,10);
		MoveOmegaHyperIntegrated(0.3*tuning,10);

		ResampleSiteOmegas();
	}
	else	{
		if (GetNprocs() > 1)	{
			ret = MPIMoveOmega(tuning);
		}
		else	{
			ret = NonMPIMoveOmega(tuning);
		}
	}
	return ret;
}

double SiteOmegaProcess::MPIMoveOmega(double tuning)	{

	return GlobalMixMoveOmega(5,tuning,1,10);	
}

double SiteOmegaProcess::NonMPIMoveOmega(double tuning)	{

	return MixMoveOmega(5,tuning,1,10);	
}

double SiteOmegaProcess::MixMoveOmega(int nmix, double tuning, int nsiterep, int nhyperrep)	{

	for (int mix=0; mix<nmix; mix++)	{
		MoveSiteOmegas(tuning, nsiterep);
		MoveOmegaHyper(tuning,nhyperrep);
		MoveOmegaHyper(0.1*tuning,nhyperrep);
	}
	UpdateOmega();
	return 1.0;
}

double SiteOmegaProcess::GlobalMixMoveOmega(int nmix, double tuning, int nsiterep, int nhyperrep)	{

	MESSAGE signal = MIXMOVEOMEGA;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&nmix,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&nsiterep,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&tuning,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

	double* tmp = new double[GetNsite()];

	for (int mix=0; mix<nmix; mix++)	{

		// slave move site omegas

		// collect site omegas
		MPI_Status stat;
		for(int i=1; i<GetNprocs(); i++) {
			MPI_Recv(tmp,GetNsite(),MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
			for(int j=GetProcSiteMin(i); j<GetProcSiteMax(i); ++j) {
				if (ActiveSite(j))	{
					omegaarray[j] = tmp[j];
				}
			}
		}

		// move hyperparam
		MoveOmegaHyper(tuning,nhyperrep);
		MoveOmegaHyper(0.1*tuning,nhyperrep);

		// send hyperparam
		MPI_Bcast(&omegaalpha,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(&omegabeta,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	}

	delete[] tmp;

	UpdateOmega();

	return 1.0;
}

void SiteOmegaProcess::SlaveMixMoveOmega()	{

	double tuning = 1;
	int nrep = 1;
	int nmix = 1;
	MPI_Bcast(&nmix,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&nrep,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&tuning,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

	for (int mix=0; mix<nmix; mix++)	{

		MoveSiteOmegas(tuning, nrep);
		MPI_Send(omegaarray,GetNsite(),MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);

		// get hyperparam
		MPI_Bcast(&omegaalpha,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(&omegabeta,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	}

	UpdateOmega();
}

double SiteOmegaProcess::OmegaSuffStatIntegratedLogProb(int site)	{

	double ret = rnd::GetRandom().logGamma(omegaalpha + siteomegasuffstatcount[site]) - rnd::GetRandom().logGamma(omegaalpha);
	ret += omegaalpha * log(omegabeta) - (omegaalpha + siteomegasuffstatcount[site]) * log(omegabeta + siteomegasuffstatbeta[site]);
	return ret;
}

double SiteOmegaProcess::OmegaSuffStatIntegratedLogProb()	{

	double tot = 0;
	for (int i=0; i<GetNsite(); i++)	{
		if (ActiveSite(i))	{
			tot += OmegaSuffStatIntegratedLogProb(i);
		}
	}
	return tot;
}

void SiteOmegaProcess::ResampleSiteOmegas()	{

	for (int i=0; i<GetNsite(); i++)	{
		if (ActiveSite(i))	{
			omegaarray[i] = rnd::GetRandom().Gamma(omegaalpha + siteomegasuffstatcount[i],omegabeta + siteomegasuffstatbeta[i]);
		}
	}
}

double SiteOmegaProcess::MoveOmegaHyperIntegrated(double tuning, int nrep)	{

	double nacc = 0;
	double ntot = 0;
	for (int rep=0; rep<nrep; rep++)	{
		nacc += MoveOmegaAlphaIntegrated(tuning);
		ntot++;
		nacc += MoveOmegaBetaIntegrated(tuning);
		ntot++;
	}
	return nacc/ntot;
}

double SiteOmegaProcess::MoveOmegaAlphaIntegrated(double tuning)	{

	double bkomegaalpha = omegaalpha;
	double deltalogprob = - OmegaSuffStatIntegratedLogProb() - LogOmegaHyperPrior();
	double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
	double e = exp(m);
	omegaalpha *= e;
	deltalogprob += OmegaSuffStatIntegratedLogProb() + LogOmegaHyperPrior();
	deltalogprob += m;
	int accepted = (log(rnd::GetRandom().Uniform()) < deltalogprob);

	if (!accepted)	{
		omegaalpha = bkomegaalpha;
	}
	return accepted;
}

double SiteOmegaProcess::MoveOmegaBetaIntegrated(double tuning)	{

	double bkomegabeta = omegabeta;
	double deltalogprob = - OmegaSuffStatIntegratedLogProb() - LogOmegaHyperPrior();
	double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
	double e = exp(m);
	omegabeta *= e;
	deltalogprob += OmegaSuffStatIntegratedLogProb() + LogOmegaHyperPrior();
	deltalogprob += m;
	int accepted = (log(rnd::GetRandom().Uniform()) < deltalogprob);

	if (!accepted)	{
		omegabeta = bkomegabeta;
	}
	return accepted;
}

void SiteOmegaProcess::GlobalCollectSiteOmegaSuffStats()	{

	MESSAGE signal = COLLECTOMEGASUFFSTAT;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	MPI_Status stat;
	int* itmp = new int[GetNsite()];
	double* tmp = new double[GetNsite()];

	for(int i=1; i<GetNprocs(); i++) {
		MPI_Recv(itmp,GetNsite(),MPI_INT,i,TAG1,MPI_COMM_WORLD,&stat);
		MPI_Recv(tmp,GetNsite(),MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
		for(int j=GetProcSiteMin(i); j<GetProcSiteMax(i); ++j) {
			siteomegasuffstatcount[j] = itmp[j];
			siteomegasuffstatbeta[j] = tmp[j];
		}
	}
	delete[] tmp;
	delete[] itmp;
}

void SiteOmegaProcess::SlaveCollectSiteOmegaSuffStats()	{

	MPI_Send(siteomegasuffstatcount,GetNsite(),MPI_INT,0,TAG1,MPI_COMM_WORLD);
	MPI_Send(siteomegasuffstatbeta,GetNsite(),MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);

}



