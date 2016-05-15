
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef AACODONMUTSELSBDPPROFILE_H
#define AACODONMUTSELSBDPPROFILE_H

#include "Parallel.h"
#include "SBDPProfileProcess.h"
#include "SingleOmegaAACodonMutSelProfileProcess.h"
#include "GeneralPathSuffStatMatrixMixtureProfileProcess.h"

class AACodonMutSelSBDPProfileProcess : public virtual SBDPProfileProcess, public virtual SingleOmegaAACodonMutSelProfileProcess, public virtual GeneralPathSuffStatMatrixMixtureProfileProcess	{

	// implementer les fonctions create matrix et delete matrix
	// ainsi que CreateComponent(int k) and DeleteComponent(k)

	// s'inspirer de GeneralPathSuffStatGTRFiniteProfileProcess

	public:

	AACodonMutSelSBDPProfileProcess() {}
	virtual ~AACodonMutSelSBDPProfileProcess() {}

	/*
	// still to be implemented
	virtual void UpdateNucStatSuffStat();
	virtual void UpdateNucRRSuffStat();
	virtual void UpdateOmegaSuffStat(); check
	*/
	
	virtual void UpdateOmegaSuffStat()	{
		omegasuffstatbeta = 0;
		omegasuffstatcount = 0;
	
		for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{

			if (ActiveSite(i))	{
		
				map<pair<int,int>, int>& paircount = GetSitePairCount(i);
				map<int,double>& waitingtime = GetSiteWaitingTime(i);
				int cat = alloc[i];
				CodonSubMatrix* codonmatrix = dynamic_cast<CodonSubMatrix*>(matrixarray[cat]);
				for (map<int,double>::iterator i = waitingtime.begin(); i!= waitingtime.end(); i++)	{
					omegasuffstatbeta += i->second * codonmatrix->RateAwayNonsyn(i->first) / *omega;
				}

				for (map<pair<int,int>, int>::iterator i = paircount.begin(); i!= paircount.end(); i++)	{
					if (! codonmatrix->Synonymous(i->first.first,i->first.second) )	{
						omegasuffstatcount += i->second;
					}
				}
			}
		}
	}

	void GlobalUpdateOmegaSuffStat()	{

		if (GetNprocs() > 1)	{
			// MPI2
			// should ask the slaves to call their UpdateRateSuffStat
			// and then gather the statistics;
			MPI_Status stat;
			MESSAGE signal = UPDATE_OMEGA;
			MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

			omegasuffstatcount = 0;
			omegasuffstatbeta = 0;	

			int ivalue;
			double dvalue;
			for(int i=1; i<GetNprocs(); i++) {
	                	MPI_Recv(&ivalue,1,MPI_INT,MPI_ANY_SOURCE,TAG1,MPI_COMM_WORLD,&stat);
				omegasuffstatcount += ivalue;
			}
        		MPI_Barrier(MPI_COMM_WORLD);
			for(int i=1; i<GetNprocs(); i++) {
				MPI_Recv(&dvalue,1,MPI_DOUBLE,MPI_ANY_SOURCE,TAG1,MPI_COMM_WORLD,&stat);
				omegasuffstatbeta += dvalue;
        		}
		}
		else	{
			UpdateOmegaSuffStat();	
		}
	}

	void SlaveUpdateOmegaSuffStat()	{
		UpdateOmegaSuffStat();
		
		MPI_Send(&omegasuffstatcount,1,MPI_INT,0,TAG1,MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Send(&omegasuffstatbeta,1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
	}
			

	protected:

	void Create()	{
		SingleOmegaAACodonMutSelProfileProcess::Create();
		SBDPProfileProcess::Create();
		GeneralPathSuffStatMatrixMixtureProfileProcess::Create();
	}
	
	void Delete()	{
		GeneralPathSuffStatMatrixMixtureProfileProcess::Delete();
		SBDPProfileProcess::Delete();
		SingleOmegaAACodonMutSelProfileProcess::Delete();
	}

	void ToStream(ostream& os)	{
		for (int i=0; i<Nnuc; i++)	{
			os << GetNucStat(i) << '\t';
		}
		os << '\n';
		os << '\n';
		for (int i=0; i<GetNnucrr(); i++)	{
			os << GetNucRR(i) << '\t';
		}
		os << '\n';
		os << '\n';

		for (int i=0; i<GetNcodon(); i++)	{
			os << codonprofile[i] << '\t';
		}
		os << '\n';
		os << '\n';
		os << *omega << '\n';

		os << kappa << '\n';
		os << Ncomponent << '\n';
		for (int j=0; j<GetDim(); j++)	{
			os << dirweight[j] << '\t';
		}
		os << '\n';
		os << '\n';
		for (int i=0; i<Ncomponent; i++)	{
			for (int j=0; j<GetDim(); j++)	{
				os << profile[i][j] << '\t';
			}
			os << '\n';
		}
		for (int i=0; i<GetNsite(); i++)	{
			if (ActiveSite(i))	{
				os << alloc[i] << '\t';
			}
		}
		os << '\n';

	}
	void FromStream(istream& is)	{
		for (int i=0; i<Nnuc; i++)	{
			is >> nucstat[i];
		}
		for (int i=0; i<GetNnucrr(); i++)	{
			is >> nucrr[i];
		}
		for (int i=0; i<GetNcodon(); i++)	{
			is >> codonprofile[i];
		}

		is >> *omega;
		is >> kappa;
		is >> Ncomponent;
		for (int j=0; j<GetDim(); j++)	{
			is >> dirweight[j];
		}
		for (int i=0; i<Ncomponent; i++)	{
			for (int j=0; j<GetDim(); j++)	{
				is >> profile[i][j];
			}
		}
		for (int i=0; i<GetNsite(); i++)	{
			if (ActiveSite(i))	{
				is >> alloc[i];
			}
		}
		ResampleWeights();
	}


	void CreateMatrix(int k)	{
		if (matrixarray[k])	{
			cerr << "error in AACodonMutSelSBDPProfileProcess: matrixarray is not 0\n";
			exit(1);
		}
		matrixarray[k] = new AACodonMutSelProfileSubMatrix(GetCodonStateSpace(),nucrr,nucstat,codonprofile,profile[k],omega,true);
	}

};

#endif

