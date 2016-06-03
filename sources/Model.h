
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "CodonSequenceAlignment.h"
#include "RASCATGTRFiniteGammaPhyloProcess.h"
#include "RASCATGTRSBDPGammaPhyloProcess.h"
#include "GeneralPathSuffStatRASCATGTRSBDPGammaPhyloProcess.h"
#include "AASubSelRASCATSBDPGammaPhyloProcess.h"
#include "RASCATGTRDPGammaPhyloProcess.h"
#include "RASCATGammaPhyloProcess.h"
#include "RASCATFiniteGammaPhyloProcess.h"
#include "RASCATSBDPGammaPhyloProcess.h"
#include "GeneralPathSuffStatRASCATGTRFiniteGammaPhyloProcess.h"
#include "AACodonMutSelFinitePhyloProcess.h"
#include "AACodonMutSelSBDPPhyloProcess.h"
#include "MultipleOmegaAACodonMutSelSBDPPhyloProcess.h"
#include "AACodonMutSelSiteSBDPPhyloProcess.h"
#include "AACodonMutSelMVNSiteSpecificPhyloProcess.h"
#include "CodonMutSelFinitePhyloProcess.h"
#include "CodonMutSelSBDPPhyloProcess.h"
#include "ZipRASCATGTRSBDPGammaPhyloProcess.h"
#include "ZipRASCATGTRFiniteGammaPhyloProcess.h"
#include "ZipGeneralPathSuffStatRASCATGTRSBDPGammaPhyloProcess.h"
#include "ZipGeneralPathSuffStatRASCATGTRFiniteGammaPhyloProcess.h"

#include "MultiGeneRASCATGTRSBDPGammaPhyloProcess.h"
#include "MultiGeneRASCATSBDPGammaPhyloProcess.h"

#include "Parallel.h"
#include <iostream>
#include <fstream>

using namespace std;

MPI_Datatype Propagate_arg;



class Model	{

	public:

	PhyloProcess* process;
	string type;
	string name;
	int every;
	int until;
	int saveall;

	Model(string datafile, string treefile, int multigene, int globalalpha, int globalbl, int modeltype, int nratecat, int mixturetype, int nmodemax, int ncat, GeneticCodeType codetype, int suffstat, int fixncomp, int empmix, string mixtype, string rrtype, int iscodon, int sis, int sisnfrac, int sisnrep, double siscutoff, int fixtopo, int fixroot, int topoburnin, int topobf, int bfburnin, int bfnfrac, int bfnrep, double blfactor, string blfile, int NSPR, int NMHSPR, int NTSPR, int temperedbl, int temperedgene, int temperedrate,double topolambda, double topomu, double toponstep, int NNNI, int nspec, int ntspec, string taxon1, string taxon2, string taxon3, string taxon4, int bpp, int nbpp, int ntbpp, int bppnstep, string bppname, double bppcutoff, double bppbeta, int fixcodonprofile, int fixomega, int Nomega, int fixbl, int sumovercomponents, int omegaprior, int kappaprior, int profilepriortype, int dc, int inevery, int inuntil, int insaveall, int zip, int proposemode, int allocmode, int fasttopo, double fasttopofracmin, int fasttoponstep, int fastcondrate, string inname, int myid, int nprocs, int sitesuffstat)	{

		every = inevery;
		until = inuntil;
		name = inname;
		saveall = insaveall;

		// 1 : CAT
		// 2 : CATGTR
		// 3 : MutSel

		// mixturetype
		// 0 : one for all
		// 1 : finite
		// 2 : dp
		// 3 : tdp : removed
		// 4 : sbdp
		// 5 : site specific
		
		if (multigene == 1)	{

			if (modeltype == 1)	{
				type = "MULTIGENECATSBDP";
				process = new MultiGeneRASCATSBDPGammaPhyloProcess(nratecat,kappaprior,nmodemax,globalalpha,globalbl);
			}
			else if (modeltype == 2)	{
				type = "MULTIGENECATGTRSBDP";
				process = new MultiGeneRASCATGTRSBDPGammaPhyloProcess(nratecat,rrtype,kappaprior,nmodemax,globalalpha,globalbl);
			}
			else	{
				cerr << "model not recognized\n";
				exit(1);
			}
		}
		else	{
		// CAT
		if (modeltype == 1)	{
			if (myid == 0) {
				// cerr << "cat model\n";
			}
			if (mixturetype == 1)	{
				type = "CATFINITE";
				process = new RASCATFiniteGammaPhyloProcess(nratecat,ncat,fixncomp,empmix,mixtype);
			}
			else if (mixturetype == 2)	{
				type = "CATDP";
				process = new RASCATGammaPhyloProcess(nratecat,kappaprior);
			}
			else	{
				type = "CATSBDP";
				process = new RASCATSBDPGammaPhyloProcess(nratecat,kappaprior);
			}
		}

		// CATGTR
		else if (modeltype == 2)	{
			if (myid == 0) {
				// cerr << "catgtr model\n";
			}
			if (mixturetype == 1)	{
				if (zip)	{
					if (suffstat)	{
						type = "ZIPCATGTRFINITE";
						process = new ZipRASCATGTRFiniteGammaPhyloProcess(nratecat,ncat,fixncomp,empmix,mixtype,rrtype);
					}
					else	{
						type = "ZIPGPSSCATGTRFINITE";
						process = new ZipGeneralPathSuffStatRASCATGTRFiniteGammaPhyloProcess(nratecat,ncat,fixncomp,empmix,mixtype,rrtype);
					}
				}
				else if (suffstat)	{
					type = "CATGTRFINITE";
					process = new RASCATGTRFiniteGammaPhyloProcess(nratecat,ncat,fixncomp,empmix,mixtype,rrtype);
				}
				else	{
					type = "GPSSCATGTRFINITE";
					process = new GeneralPathSuffStatRASCATGTRFiniteGammaPhyloProcess(nratecat,ncat,fixncomp,empmix,mixtype,rrtype);
				}
			}
			else if (mixturetype == 2)	{
				cerr << "simple dp deprecated\n";
				exit(1);
			}
			else if (mixturetype == 3)	{
				if (zip)	{
					if (suffstat)	{
						type = "ZIPCATGTRSBDP";
						process = new ZipRASCATGTRSBDPGammaPhyloProcess(nratecat,rrtype,kappaprior);
					}
					else	{
						type = "ZIPGPSSCATGTRSBDP";
						process = new ZipGeneralPathSuffStatRASCATGTRSBDPGammaPhyloProcess(nratecat,rrtype,kappaprior);
					}
				}
				else	{
					if (suffstat)	{
						type = "CATGTRSBDP";
						process = new RASCATGTRSBDPGammaPhyloProcess(nratecat,rrtype,kappaprior);
					}
					else	{
						type = "GPSSCATGTRSBDP";
						process = new GeneralPathSuffStatRASCATGTRSBDPGammaPhyloProcess(nratecat,rrtype,kappaprior);
					}
				}
			}
			else	{
				cerr << "mixture type " << mixturetype << " not recognized\n";
				exit(1);
			}
		}

		else if (modeltype == 3)	{
			cerr << "deprecated.\n";
			exit(1);
		}

		// CodonMutSel
		else if (modeltype == 4)	{
			if (mixturetype == 1)	{
				type = "CODONMUTSELFINITE";
				process = new CodonMutSelFinitePhyloProcess(ncat,fixncomp,empmix,mixtype);
			}
			else if (mixturetype == 3)	{
				type = "CODONMUTSELSBDP";
				process = new CodonMutSelSBDPPhyloProcess(kappaprior);
			}
			else	{
				cerr << "mixture type " << mixturetype << " not recognized or not yet implemented.\n";
				exit(1);
			}
		}

		// AACodonMutSel
		else if (modeltype == 5) {
			if (mixturetype == 1)	{
				type = "AACODONMUTSELFINITE";
				process = new AACodonMutSelFinitePhyloProcess(ncat,fixncomp,empmix,mixtype,fixcodonprofile,fixomega,omegaprior);
			}
			else if (mixturetype == 3)	{
				if (suffstat == 2)	{
					type = "AACODONMUTSELSITESBDP";
					process = new AACodonMutSelSiteSBDPPhyloProcess(fixcodonprofile,fixomega,omegaprior,kappaprior);
				}
				else	{
					type = "AACODONMUTSELSBDP";
					process = new AACodonMutSelSBDPPhyloProcess(fixcodonprofile,fixomega,omegaprior,kappaprior);
				}
			}
			else if (mixturetype == 5)	{
				type = "AACODONMUTSELMVNSS";
				process = new AACodonMutSelMVNSiteSpecificPhyloProcess(fixcodonprofile,fixomega,omegaprior);
			}
			else if (mixturetype == 6)	{
				type = "MULOMEGAAACODONMUTSELSBDP";
				process = new MultipleOmegaAACodonMutSelSBDPPhyloProcess(fixcodonprofile,Nomega,fixomega,omegaprior,kappaprior);
			}
			else	{
				cerr << "mixture type " << mixturetype << " not recognized or not yet implemented.\n";
				exit(1);
			}
			
		}

		// AA SubSel
		else if (modeltype == 6)	{
			if (myid == 0) {
				// cerr << "catgtr model\n";
			}
			if (mixturetype == 1)	{
				cerr << "finite aa sub sel not yet implemented\n";
				exit(1);
			}
			else if (mixturetype == 2)	{
				cerr << "simple dp deprecated\n";
				exit(1);
			}
			else if (mixturetype == 3)	{
				type = "AASUBSELSBDP";
				process = new AASubSelRASCATSBDPGammaPhyloProcess(nratecat,kappaprior);
			}
			else	{
				cerr << "mixture type " << mixturetype << " not recognized\n";
				exit(1);
			}
		}
		}

		process->SetParameters(datafile,treefile,iscodon,codetype,sis,sisnfrac,sisnrep,siscutoff,fixtopo,fixroot,topoburnin,topobf,bfburnin,bfnfrac,bfnrep,blfactor,blfile,NSPR,NMHSPR,NTSPR,temperedbl,temperedgene,temperedrate,topolambda,topomu,toponstep,NNNI,nspec,ntspec,taxon1,taxon2,taxon3,taxon4,bpp,nbpp,ntbpp,bppnstep,bppname,bppcutoff,bppbeta,profilepriortype,dc,fixbl,sumovercomponents,proposemode,allocmode,fasttopo,fasttopofracmin,fasttoponstep,fastcondrate);

		if (topobf == 1)	{
			until = bfburnin + bfnfrac*bfnrep;
		}
		if (topobf == 2)	{
			until = bfburnin + 2*bfnfrac*bfnrep;
		}
		process->SetName(name);

		// process->SetTrackTopo(1);
		process->SetSiteSuffStat(sitesuffstat);

		process->SetMPI(myid,nprocs);
		process->New();

	}

	Model(string inname, int myid, int nprocs)	{

		name = inname;

		ifstream is((name + ".param").c_str());
		if (! is)	{
			cerr << "error: cannot open " << name << ".param\n";
			exit(1);
		}

		is >> type;
		int size;
		is >> every >> until >> size;
		is >> saveall;
		
		if (type == "CATDP")	{
			process = new RASCATGammaPhyloProcess(is,myid,nprocs); 
		}
		else if (type == "CATSBDP")	{
			process = new RASCATSBDPGammaPhyloProcess(is,myid,nprocs); 
		}
		else if (type == "MULTIGENECATSBDP")	{
			process = new MultiGeneRASCATSBDPGammaPhyloProcess(is,myid,nprocs); 
		}
		else if (type == "CATFINITE")	{
			process = new RASCATFiniteGammaPhyloProcess(is,myid,nprocs); 
		}
		else if (type == "CATGTRDP")	{
			process = new RASCATGTRDPGammaPhyloProcess(is,myid,nprocs); 
		}
		else if (type == "CATGTRSBDP")	{
			process = new RASCATGTRSBDPGammaPhyloProcess(is,myid,nprocs); 
		}
		else if (type == "MULTIGENECATGTRSBDP")	{
			process = new MultiGeneRASCATGTRSBDPGammaPhyloProcess(is,myid,nprocs); 
		}
		else if (type == "ZIPCATGTRSBDP")	{
			process = new ZipRASCATGTRSBDPGammaPhyloProcess(is,myid,nprocs); 
		}
		else if (type == "ZIPGPSSCATGTRSBDP")	{
			process = new ZipGeneralPathSuffStatRASCATGTRSBDPGammaPhyloProcess(is,myid,nprocs); 
		}
		else if (type == "GPSSCATGTRSBDP")	{
			process = new GeneralPathSuffStatRASCATGTRSBDPGammaPhyloProcess(is,myid,nprocs); 
		}
		else if (type == "CATGTRFINITE")	{
			process = new RASCATGTRFiniteGammaPhyloProcess(is,myid,nprocs); 
		}
		else if (type == "GPSSCATGTRFINITE")	{
			process = new GeneralPathSuffStatRASCATGTRFiniteGammaPhyloProcess(is,myid,nprocs); 
		}
		else if (type == "ZIPCATGTRFINITE")	{
			process = new ZipGeneralPathSuffStatRASCATGTRFiniteGammaPhyloProcess(is,myid,nprocs); 
		}
		else if (type == "AACODONMUTSELFINITE")	{
			process = new AACodonMutSelFinitePhyloProcess(is,myid,nprocs);
		}
		else if (type == "AACODONMUTSELSBDP")	{
			process = new AACodonMutSelSBDPPhyloProcess(is,myid,nprocs);
		}
		else if (type == "AACODONMUTSELMVNSS")	{
			process = new AACodonMutSelMVNSiteSpecificPhyloProcess(is,myid,nprocs);
		}
		else if (type == "CODONMUTSELFINITE")	{
			process = new CodonMutSelFinitePhyloProcess(is,myid,nprocs);
		}
		else if (type == "CODONMUTSELSBDP")	{
			process = new CodonMutSelSBDPPhyloProcess(is,myid,nprocs);
		}
		else	{
			cerr << "error, does not recognize model type : " << type << '\n';
			exit(1);
		}

		process->SetSize(size);
		process->SetName(name);
		if ((!myid) && (process->topobf))	{
			ifstream mpis((name + ".mpi").c_str());
			process->GlobalReadSiteRankFromStream(mpis);
		}
	}

	void ToStream(ostream& os, bool header)	{
		if (header)	{
			os << type << '\n';
			os << every << '\t' << until << '\t' << GetSize() << '\n';
			os << saveall << '\n';
			process->ToStreamHeader(os);
		}
		process->ToStream(os);
	}

	~Model()	{
		// things here !!
	}

	void WaitLoop()	{
		process->WaitLoop();
	}

	double Move(double tuning, int nrep)	{
		double total = 0;
		for (int rep=0; rep<nrep; rep++)	{
			total += process->Move(tuning);
		}
		return total / nrep;
	}

	int RunningStatus()	{
		ifstream ris((name + ".run").c_str());
		int i;
		ris >> i;
		ris.close();
		return i;
	}

	int GetSize()	{
		return process->GetSize();
	}

	/*
	void SMC(int shortcycle, int mediumcycle, int longcycle, int mediumsize, int maxsize, int nrep)	{

		if (nrep)	{
			ofstream tos((name + ".trees").c_str());
			ofstream os((name + ".logweights").c_str());
			double logw[nrep];
			for (int rep=0; rep<nrep; rep++)	{
				logw[rep] = process->SMC(name,shortcycle,mediumcycle,longcycle,mediumsize,maxsize);
				os << logw[rep] << '\n';
				os.flush();
				process->RenormalizeBranchLengths();
				GetTree()->ToStream(tos);
				process->DenormalizeBranchLengths();
				tos.flush();
			}
		}
		else	{
			process->SMC(name,shortcycle,mediumcycle,longcycle,mediumsize,maxsize);
		}
	}
	*/

	void Run(int smc, int deltansite, int shortcycle, int longcycle, int cutoffsize, int nrep)	{

		ofstream ros((name + ".run").c_str());
		ros << 1 << '\n';
		ros.close();
	
		if (smc)	{
			cerr << "burnin\n";
			process->SMCBurnin(name,deltansite,shortcycle,longcycle,cutoffsize,nrep);
			cerr << "burnin done\n";
		}

		while (RunningStatus() && ((until == -1) || (GetSize() < until)))	{

			Move(1,every);
			
			process->IncSize();

			ofstream os((name + ".treelist").c_str(), ios_base::app);
			process->SetNamesFromLengths();
			process->RenormalizeBranchLengths();
			GetTree()->ToStream(os);
			process->DenormalizeBranchLengths();
			os.close();

			ofstream tos((name + ".trace").c_str(), ios_base::app);
			Trace(tos);
			tos.close();

			ofstream mos((name + ".monitor").c_str());
			process->Monitor(mos);
			mos.close();

			ofstream pos((name + ".param").c_str());
			pos.precision(12);
			ToStream(pos,true);
			pos.close();

			if (saveall)	{
				ofstream cos((name + ".chain").c_str(),ios_base::app);
				cos.precision(12);
				ToStream(cos,false);
				cos.close();
			}
		}	
		cerr << name << ": stopping after " << GetSize() << " points.\n";
		cerr << '\n';
	}

	NewickTree* GetTree() {return process->GetLengthTree();}

	void TraceHeader(ostream& os)	{
		process->TraceHeader(os);
	}

	void Trace(ostream& os)	{
		process->Trace(os);
	}

	void ReadPB(int argc, char* argv[])	{
		process->ReadPB(argc,argv);
	}
};
