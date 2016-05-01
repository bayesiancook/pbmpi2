
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef GPSSZIPRASCATGTRSBDPPHYLO_H
#define GPSSZIPRASCATGTRSBDPPHYLO_H

#include "GeneralPathSuffStatRASCATGTRSBDPGammaPhyloProcess.h"
#include "ZipGeneralPathSuffStatMatrixMixtureProfileProcess.h"
#include "ZipMatrixPhyloProcess.h"

class ZipGeneralPathSuffStatGTRSBDPProfileProcess : public virtual GeneralPathSuffStatGTRSBDPProfileProcess, public virtual ZipGeneralPathSuffStatMatrixMixtureProfileProcess	{

	using MixtureProfileProcess::LogStatPrior;

	public:

	ZipGeneralPathSuffStatGTRSBDPProfileProcess() : rrsuffstatcount(0) {}
	virtual ~ZipGeneralPathSuffStatGTRSBDPProfileProcess() {}

	protected:

	virtual double Move(double tuning = 1, int n = 1, int nrep = 1, int nalloc = 1)	{
		totchrono.Start();

		for (int rep=0; rep<nrep; rep++)	{

			// relative rates
			if (! fixrr)	{
				GlobalUpdateParameters();
				GlobalUpdateSiteProfileSuffStat();
				if (proposemode)	{
					GlobalUpdateModeProfileSuffStat();
					ApproxSuffStatSampleRR();
				}
				else	{
					rracc += MoveRR(0.5,5,3);
					rracc += MoveRR(0.3,10,3);
					rracc += MoveRR(0.1,10,3);
					rrtry += 3;
				}
			}

			incchrono.Start();
			GlobalUpdateParameters();
			GlobalUpdateSiteProfileSuffStat();
			if (proposemode)	{
				GlobalUpdateModeProfileSuffStat();
				if (allocmode)	{
					profacc += GlobalMixMove(5,nalloc,0.1,40);
				}
				else	{
					profacc += GlobalMixMove(5,nalloc,0.001,40);
				}
			}
			else	{
				if (allocmode)	{
					profacc += ZipGlobalMixMove(5,nalloc,0.1,40);
				}
				else	{
					profacc += ZipGlobalMixMove(5,nalloc,0.001,40);
				}
			}
			/*
			double acc1 = MultipleTryGlobalMixMove(5,1,0.001,1,4,1,5);
			double acc2 = MultipleTryGlobalMixMove(5,1,0.001,1,4,1,10);
			double acc3 = MultipleTryGlobalMixMove(5,1,0.001,1,4,0.3,10);
			// cerr << acc1 << '\t' << acc2 << '\t' << acc3 << '\n';
			*/
			// MixMove(5,1,0.001,40);
			MoveOccupiedCompAlloc(5);
			MoveAdjacentCompAlloc(5);
			incchrono.Stop();

			MoveHyper(tuning,10);
			// UpdateMatrices();

		}
		totchrono.Stop();
		return 1;
	}

	virtual double LogStatProb(int site, int cat);
	virtual double ProfileSuffStatLogProb(int cat);
	virtual double ProfileSuffStatLogProb();

	virtual double MoveRR(double tuning, int n, int nrep);
	void SlaveMoveRR();

	// based on approximate suff stats: returns logh = log(approxprob(initial) / approxprob(final))
	double ApproxProfileSuffStatLogProb(int cat);
	double ApproxSuffStatProfileProposeMove(int cat, int nstep);
	double ApproxSuffStatSampleRR();

	virtual double ZipGlobalMixMove(int nrep, int nallocrep, double epsilon, int nprofilerep);
	virtual void ZipSlaveMixMove();
	virtual void SlaveMixMove();

	virtual double MultipleTryGlobalMixMove(int nrep, int nallocrep, double epsilon, int nprofilerep, int nmultryrep, double tuning, int nentry);
	virtual void MultipleTrySlaveMixMove();
	void SpecialProfileProposeMove(double* oldprofile, double* newprofile, double tuning, int nentry);
	double SpecialProfileProposeMoveLogProb(double* oldprofile, double* newprofile, double tuning, int nentry);

	virtual void Create();
	virtual void Delete();

	// approximate suff stats: based on suff stats already stored by slaves
	virtual void UpdateRRSuffStat();

	virtual void GlobalUpdateModeProfileSuffStat();
	virtual void SlaveUpdateModeProfileSuffStat();
	virtual void UpdateModeProfileSuffStat();

	// component-specific sufficient statistics
	int** profilesuffstatcount;
	double** profilesuffstatbeta;
	int* allocprofilesuffstatcount;
	double* allocprofilesuffstatbeta;
	int* rrsuffstatcount;
	double* rrsuffstatbeta;

};

class ZipGeneralPathSuffStatRASCATGTRSBDPSubstitutionProcess : public virtual GeneralPathSuffStatRASCATGTRSBDPSubstitutionProcess, public virtual ZipMatrixSubstitutionProcess, public virtual ZipGeneralPathSuffStatGTRSBDPProfileProcess	{

	public:

	ZipGeneralPathSuffStatRASCATGTRSBDPSubstitutionProcess() {}
	virtual ~ZipGeneralPathSuffStatRASCATGTRSBDPSubstitutionProcess() {}

	int GetNstate(int site) {return ZipMatrixSubstitutionProcess::GetNstate(site);}

	protected:

	virtual void Create()	{
		ZipMatrixSubstitutionProcess::Create();
		GeneralPathSuffStatRASCATGTRSBDPSubstitutionProcess::Create();
		ZipGeneralPathSuffStatGTRSBDPProfileProcess::Create();
	}

	virtual void Delete()	{
		ZipGeneralPathSuffStatGTRSBDPProfileProcess::Delete();
		GeneralPathSuffStatRASCATGTRSBDPSubstitutionProcess::Delete();
		ZipMatrixSubstitutionProcess::Delete();
	}
};

class ZipGeneralPathSuffStatRASCATGTRSBDPGammaPhyloProcess : public virtual GeneralPathSuffStatRASCATGTRSBDPGammaPhyloProcess, public virtual ZipGeneralPathSuffStatRASCATGTRSBDPSubstitutionProcess, public virtual ZipMatrixPhyloProcess	{

	public:

	ZipGeneralPathSuffStatRASCATGTRSBDPGammaPhyloProcess(int nratecat, string inrrtype, int inkappaprior)	{

		Ncat = nratecat;
		rrtype = inrrtype;
		kappaprior = inkappaprior;
	}

	ZipGeneralPathSuffStatRASCATGTRSBDPGammaPhyloProcess(istream& is, int inmyid, int innprocs)	{

		// generic
		FromStreamHeader(is);
		SetMPI(inmyid,innprocs);

		// specific
		is >> Ncat;
		is >> kappaprior;
		is >> rrtype;

		Open(is);

	}

	~ZipGeneralPathSuffStatRASCATGTRSBDPGammaPhyloProcess() {
		Delete();
	}

	void ToStreamHeader(ostream& os)	{
		PhyloProcess::ToStreamHeader(os);
		os << Ncat << '\n';
		os << kappaprior << '\n';
		os << rrtype << '\n';
	}

	void UpdateRRSuffStat()	{
		ZipGeneralPathSuffStatGTRSBDPProfileProcess::UpdateRRSuffStat();
	}

	virtual double Move(double tuning = 1.0)	{

		chronototal.Start();

		propchrono.Start();
		BranchLengthMove(tuning);
		BranchLengthMove(0.1 * tuning);
		if (! fixtopo)	{
			MoveTopo();
		}
		propchrono.Stop();

		GlobalCollapse();

		GammaBranchProcess::Move(tuning,10);

		GlobalUpdateParameters();
		DGamRateProcess::Move(0.3*tuning,10);
		DGamRateProcess::Move(0.03*tuning,10);

		ZipGeneralPathSuffStatGTRSBDPProfileProcess::Move(1,1,5);

		if (! fixrr)	{
			LengthRelRateMove(1,10);
			LengthRelRateMove(0.1,10);
			LengthRelRateMove(0.01,10);
		}


		GlobalUnfold();

		chronototal.Stop();

		return 1;
	}

	protected:

	virtual void SlaveExecute(MESSAGE signal);

	virtual void Create()	{
		ZipMatrixPhyloProcess::Create();
		ZipGeneralPathSuffStatRASCATGTRSBDPSubstitutionProcess::Create();
		GeneralPathSuffStatRASCATGTRSBDPGammaPhyloProcess::Create();
	}
		
	virtual void Delete()	{
		GeneralPathSuffStatRASCATGTRSBDPGammaPhyloProcess::Delete();
		ZipGeneralPathSuffStatRASCATGTRSBDPSubstitutionProcess::Delete();
		ZipMatrixPhyloProcess::Delete();
	}

	/*
	int** siteprofilesuffstatcount;
	double** siteprofilesuffstatbeta;

	int* allocsiteprofilesuffstatcount;
	double* allocsiteprofilesuffstatbeta;
	*/
};

#endif

