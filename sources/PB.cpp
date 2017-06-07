
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "Model.h"

int main(int argc, char* argv[])	{

	int myid  = 0;
	int nprocs = 0;

	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

	int blockcounts[2] = {1,3};
	MPI_Datatype types[2] = {MPI_DOUBLE,MPI_INT};
	MPI_Aint dtex,displacements[2];
	
	displacements[0] = (MPI_Aint) 0;
	MPI_Type_extent(MPI_DOUBLE,&dtex);
	displacements[1] = dtex;
	MPI_Type_struct(2,blockcounts,displacements,types,&Propagate_arg);
	MPI_Type_commit(&Propagate_arg); 

	if (! myid)	{
		cerr << '\n';
	}

	string datafile = "";
	string treefile = "None";
	string name = "";
	GeneticCodeType type = Universal;

	string partitionfile = "None";

	int every = 1;
	int until = -1;

	int mixturetype = -1;
	// int mixturetype = 3;
	int modeltype = -1;
	// int modeltype = 2;
	int dgam = 4;
	int pinv = 0;
	int ncat = 100;
	int iscodon = 0;
	int omegaprior = 0;
	int omegamixturetype = -1; // -1 for  homogeneous omega, 1 for finite, 3 for sbdp, 4 for site-specific-gamma-distributed


	int dc = 0;
	int fixtopo = 0;
	int fixroot = 0;
	int fixcodonprofile = 1;
	int fixomega = 1;
	int nomega = -1;
	int NSPR = 10;
	int NMHSPR = 0;
	int NTSPR = 0;
	double topolambda = 1.0;
	double topomu = 0;
	int toponstep = 10;
	int NNNI = 0;
	string taxon1 = "None";
	string taxon2 = "None";
	string taxon3 = "None";
	string taxon4 = "None";
	string roottax1 = "None";
	string roottax2 = "None";
	int nspec = 0;
	int ntspec = 0;
	int fixbl = 0;
	int sumovercomponents = 0;
	int fixncomp = 0;
	int force = 0;
	int empmix = 0;
	string mixtype = "None";
	string rrtype = "None";

	int dirpriortype = 1;
	int nstatcomp = 1;
	int priorempmix = 0;
	string priormixtype = "None";
	int fixstatweight = 0;
	int fixstatalpha = 0;
	int fixstatcenter = 0;

	int kappaprior = 0;
	int profilepriortype =0;

	int suffstat = 1;

	int saveall = 1;

	int randfix = -1;

	int zip = 0;

	int smc = 0;
	int deltansite = 10;
	int shortcycle = 1;
	int longcycle = 10;
	int cutoffsize = 0;
	int nrep = 1;

	int bpp = 0;
	int nbpp = 0;
	int ntbpp = 0;
	int bppnstep = 0;
	string bppname = "";
	double bppcutoff = 0.1;
	double bppbeta = 1.0;

	int proposemode = 0;
	int allocmode = 0;
	// int sumratealloc = 1;

	int fasttopo = 0;
	double fasttopofracmin = 0.1;
	int fasttoponstep = 10;

	int fastcondrate = 0;

	int multigene = 0;
	int globalalpha = 1;
	int globalbl = 1;
	int mappsuffstat = 1;

	int topoburnin = 0;

	int topobf = 0;
	int bfburnin = 1000;
	int bfnfrac = 100;
	int bfnrep = 100;
	double bffrac = 0;
	double blfactor = 1.0;
	string blfile = "None";

	int sis = 0;
	double sisfrac = 0;
	int sisnfrac = 100;
	int sisnrep = 10;
	double siscutoff = 0.2;

	int nmodemax = 0;

	int temperedbl = 1;
	int temperedgene = 0;
	int temperedrate = 0;

	int reshuffle = 0;

	try	{

		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];
			if ((s == "-v") || (s == "--version"))	{
				if (! myid)	{
					cerr << "\n";
					cerr << "pb_mpi version 1.5\n";
					cerr << "\n";
				}
				MPI_Finalize();
				exit(1);
			}
			else if ((s == "-h") || (s == "--help"))	{
				throw(0);
			}
			else if (s == "-f")	{
				force = 1;
			}
			else if (s == "-rnd")	{
				i++;
				randfix = atoi(argv[i]);
			}
			else if (s == "-d")	{
				i++;
				datafile = argv[i];
			}
			else if (s == "-part")	{
				i++;
				partitionfile = argv[i];
			}
			else if (s == "-iscodon")	{
				iscodon = 1;
			}
			else if ((s == "-t") || (s == "-T"))	{
				i++;
				treefile = argv[i];
				if (s == "-T")	{
					fixtopo = 1;
				}
			}
			else if (s == "-fixroot")	{
				fixroot = 1;
			}
			else if (s == "-reroot")	{
				i++;
				roottax1 = argv[i];
				i++;
				roottax2 = argv[i];
			}
			else if (s == "-topoburnin")	{
				i++;
				topoburnin = atoi(argv[i]);
			}
			else if (s == "-sis")	{
				sis = 1;
				i++;
				sisnfrac = atoi(argv[i]);
				i++;
				sisnrep = atoi(argv[i]);
				i++;
				siscutoff = atof(argv[i]);
			}
			else if (s == "-fixsis")	{
				sis = 2;
				i++;
				int tmp = atoi(argv[i]);
				i++;
				sisnfrac = atoi(argv[i]);
				sisfrac = ((double) tmp) / sisnfrac;
				i++;
				sisnrep = atoi(argv[i]);
				i++;
				siscutoff = atof(argv[i]);
			}
			else if (s == "-fixtopobf")	{
				topobf = 3;
				i++;
				bfburnin = atoi(argv[i]);
				i++;
				int tmp = atoi(argv[i]);
				i++;
				bfnfrac = atoi(argv[i]);
				bffrac = ((double) tmp) / bfnfrac;
				i++;
				bfnrep = atoi(argv[i]);
				i++;
				taxon1 = argv[i];
				i++;
				taxon2 = argv[i];
				i++;
				taxon3 = argv[i];
				i++;
				taxon4 = argv[i];
			}
			else if (s == "-topobf")	{
				topobf = 1;
				i++;
				bfburnin = atoi(argv[i]);
				i++;
				bfnfrac = atoi(argv[i]);
				i++;
				bfnrep = atoi(argv[i]);
				i++;
				taxon1 = argv[i];
				i++;
				taxon2 = argv[i];
				i++;
				taxon3 = argv[i];
				i++;
				taxon4 = argv[i];
			}
			else if (s == "-blfile")	{
				i++;
				blfile = argv[i];
			}
			else if (s == "-topobl")	{
				topobf = 2;
				i++;
				bfburnin = atoi(argv[i]);
				i++;
				bfnfrac = atoi(argv[i]);
				i++;
				bfnrep = atoi(argv[i]);
				i++;
				blfactor = atof(argv[i]);
				i++;
				taxon1 = argv[i];
				i++;
				taxon2 = argv[i];
				i++;
				taxon3 = argv[i];
				i++;
				taxon4 = argv[i];
			}
			else if (s == "-fixtopobl")	{
				topobf = 4;
				i++;
				bfburnin = atoi(argv[i]);
				i++;
				bffrac = atof(argv[i]);
				i++;
				bfnfrac = atoi(argv[i]);
				i++;
				bfnrep = atoi(argv[i]);
				i++;
				blfactor = atof(argv[i]);
				i++;
				taxon1 = argv[i];
				i++;
				taxon2 = argv[i];
				i++;
				taxon3 = argv[i];
				i++;
				taxon4 = argv[i];
			}
			else if (s == "-smc")	{
				smc = 1;
				i++;
				deltansite = atoi(argv[i]);
				i++;
				cutoffsize = atoi(argv[i]);
			}
			else if (s == "-reshuffle")	{
				reshuffle = 1;
			}
			else if (s == "-Tbl")	{
				i++;
				treefile = argv[i];
				fixtopo = 1;
				fixbl = 1;
			}
			else if (s == "-fasttopo")	{
				fasttopo = 1;
				i++;
				fasttopofracmin = atof(argv[i]);
				i++;
				fasttoponstep = atoi(argv[i]);
				i++;
				topomu = atof(argv[i]);
			}
			else if (s == "-fastcondrate")	{
				fastcondrate = 1;
			}
			else if (s == "-fullsumcomp")	{
				sumovercomponents = -1;
			}
			else if (s == "-sumcomp")	{
				i++;
				sumovercomponents = atoi(argv[i]);
			}
			else if (s == "-nmodemax")	{
				i++;
				nmodemax = atoi(argv[i]);
			}
			else if (s == "-bpp")	{
				bpp = 1;
				i++;
				nbpp = atoi(argv[i]);
				i++;
				ntbpp = atoi(argv[i]);
				i++;
				bppnstep = atoi(argv[i]);
				i++;
				bppname = argv[i];
				i++;
				bppcutoff = atof(argv[i]);
				i++;
				bppbeta = atof(argv[i]);
			}
			else if (s == "-ccp")	{
				bpp = 2;
				i++;
				nbpp = atoi(argv[i]);
				i++;
				ntbpp = atoi(argv[i]);
				i++;
				bppnstep = atoi(argv[i]);
				i++;
				bppname = argv[i];
				i++;
				bppcutoff = atof(argv[i]);
				i++;
				bppbeta = atof(argv[i]);
			}
			else if (s == "-ftp")	{
				bpp = 3;
				i++;
				nbpp = atoi(argv[i]);
				i++;
				ntbpp = atoi(argv[i]);
				i++;
				bppnstep = atoi(argv[i]);
				i++;
				bppname = argv[i];
				i++;
				bppcutoff = atof(argv[i]);
				i++;
				bppbeta = atof(argv[i]);
			}
			else if (s == "-spr")	{
				i++;
				NSPR = atoi(argv[i]);
			}
			else if ((s == "-sprmh") || (s == "-mhspr"))	{
				i++;
				NMHSPR = atoi(argv[i]);
				i++;
				topolambda = atof(argv[i]);
			}

			else if (s == "-tspr")	{
				i++;
				NTSPR = atoi(argv[i]);
				i++;
				topolambda = atof(argv[i]);
				i++;
				topomu = atof(argv[i]);
				i++;
				toponstep = atoi(argv[i]);
			}
			else if (s == "-nni")	{
				i++;
				NNNI = atoi(argv[i]);
			}
			else if (s == "-specspr")	{
				i++;
				nspec = atoi(argv[i]);
				i++;
				taxon1 = argv[i];
				i++;
				taxon2 = argv[i];
				fixroot = 1;
			}
			else if (s == "-tspecspr")	{
				i++;
				ntspec = atoi(argv[i]);
				i++;
				topomu = atof(argv[i]);
				i++;
				toponstep = atoi(argv[i]);
				i++;
				taxon1 = argv[i];
				i++;
				taxon2 = argv[i];
				fixroot = 1;
			}
			else if (s == "-binspr")	{
				i++;
				nspec = atoi(argv[i]);
				i++;
				taxon1 = argv[i];
				i++;
				taxon2 = argv[i];
				i++;
				taxon3 = argv[i];
				i++;
				taxon4 = argv[i];
				fixroot = 1;
			}
			else if (s == "-tbinspr")	{
				i++;
				ntspec = atoi(argv[i]);
				i++;
				toponstep = atoi(argv[i]);
				i++;
				taxon1 = argv[i];
				i++;
				taxon2 = argv[i];
				i++;
				taxon3 = argv[i];
				i++;
				taxon4 = argv[i];
				fixroot = 1;
			}
			else if (s == "-temperedbl")	{
				i++;
				temperedbl = atoi(argv[i]);
			}
			else if (s == "-temperedgene")	{
				i++;
				temperedgene = atoi(argv[i]);
			}
			else if (s == "-temperedrate")	{
				i++;
				temperedrate = atoi(argv[i]);
			}
			else if (s == "-fixcodonprofile")	{
				fixcodonprofile = 1;
			}
			else if (s == "-freecodonprofile")	{
				fixcodonprofile = 0;
			}
			else if (s == "-fixomega")	{
				fixomega = 1;
			}
			else if (s == "-freeomega")	{
				fixomega = 0;
			}
			else if (s == "-Nomega")	{
				i++;
				nomega = atoi(argv[i]);
				// mixturetype = 3;
			}
			else if (s == "-nomega")	{
				i++;
				nomega = atoi(argv[i]);
				//fixomega = 0;
				//if (nomega > 1)	{
					mixturetype = 6;
				//}
			}
			else if (s == "-finiteomegafiniteaa")	{
				mixturetype = 7;
			}
			else if (s == "-dc")	{
				dc = 1;
			}
			else if (s == "-s")	{
				saveall = 1;
			}
			else if (s == "-S")	{
				saveall = 0;
			}
			else if (s == "-multigene")	{
				multigene = 1;
			}
			else if (s == "-alpha")	{
				i++;
				s = argv[i];
				if (s == "gene")	{
					globalalpha = 0;
				}
				else if (s == "global")	{
					globalalpha = 1;
				}
				else	{
					cerr << "error after -alpha: should be gene or global\n";
					throw(0);
				}
			}
			else if (s == "-bl")	{
				i++;
				s = argv[i];
				if (s == "gene")	{
					globalbl = 0;
				}
				else if (s == "global")	{
					globalbl = 1;
				}
				else	{
					cerr << "error after -bl: should be gene or global\n";
					throw(0);
				}
			}
			else if (s == "+mappsuffstat")	{
				mappsuffstat = 1;
			}
			else if (s == "-mappsuffstat")	{
				mappsuffstat = 0;
			}
			else if (s == "-zip")	{
				zip = 1;
			}
			else if ((s == "-poisson") || (s == "-f81"))	{
				modeltype = 1;
			}
			else if (s == "-gtr")	{
				modeltype = 2;
			}
			else if (s == "-globalomega")	{
				iscodon = 1;
				modeltype = 5;
				mixturetype = 1;
				ncat = 1;
				empmix = 1;
				fixncomp = 1;
				mixtype = "uniform";
				fixomega = 0;
			}
			else if (s == "-siteomega")	{
				iscodon = 1;
				modeltype = 5;
				mixturetype = 1;
				omegamixturetype = 4;
				ncat = 1;
				empmix = 1;
				fixncomp = 1;
				mixtype = "uniform";
				fixomega = 0;
			}
			else if (s == "-mutselc")	{
				iscodon = 1;
				modeltype = 4;
			}
			else if ((s == "-mutsel") || (s == "-mutselaa") || (s == "-mutselaac"))	{
				iscodon = 1;
				modeltype = 5;
			}
			else if (s == "-aasubsel")	{
				modeltype = 6;
			}
			else if (s == "-dgam")	{
				i++;
				dgam = atoi(argv[i]);
			}
			else if (s == "-igam")	{
				i++;
				dgam = atoi(argv[i]);
				pinv = 1;
			}
			else if (s == "-genpath")	{
				suffstat = 0;
			}
			else if (s == "-sitegenpath")	{
				suffstat = 2;
			}
			else if (s == "-olddp")	{
				mixturetype = 2;
			}
			else if ((s == "-finite") || (s == "-ncat"))	{
				mixturetype = 1;
				i++;
				ncat = atoi(argv[i]);
				fixncomp = 1;
			}
			else if (s == "-fixncomp")	{
				fixncomp = 1;
			}
			else if (s == "-freencomp")	{
				fixncomp = 0;
			}
			else if (s == "-conscatfix")	{
				mixturetype = 1;
				empmix = 3;
				fixncomp = 1;
				i++;
				mixtype = argv[i];
			}
			else if (s == "-catfix")	{
				mixturetype = 1;
				empmix = 1;
				fixncomp = 1;
				i++;
				mixtype = argv[i];
			}
			else if (s == "-priorcatfix")	{
				mixturetype = 1;
				empmix = 2;
				fixncomp = 1;
				i++;
				mixtype = argv[i];
			}
			else if (s == "-mixbaseprior")	{
				dirpriortype = 0;
				i++;
				nstatcomp = atoi(argv[i]);
			}
			else if (s == "-empmixbaseprior")	{
				dirpriortype = 0;
				i++;
				priorempmix = 1;
				priormixtype = argv[i];
				fixstatalpha = 1;
				fixstatcenter = 1;
			}
			else if (s == "-fixstatweight")	{
				fixstatweight = 1;
			}
			else if (s == "-movestatweight")	{
				fixstatweight = 0;
			}
			else if (s == "-fixstatalpha")	{
				fixstatalpha = 1;
			}
			else if (s == "-movestatalpha")	{
				fixstatalpha = 0;
			}
			else if (s == "-fixstatcenter")	{
				fixstatcenter = 1;
			}
			else if (s == "-movestatcenter")	{
				fixstatcenter = 0;
			}
			else if (s == "-rr")	{
				modeltype = 2;
				i++;
				rrtype = argv[i];
			}
			else if (s == "-lg")	{
				modeltype = 2;
				rrtype = "lg";
			}
			else if (s == "-wag")	{
				modeltype = 2;
				rrtype = "wag";
			}
			else if (s == "-jtt")	{
				modeltype = 2;
				rrtype = "jtt";
			}
			else if (s == "-mtzoa")	{
				modeltype = 2;
				rrtype = "mtzoa";
			}
			else if (s == "-mtrev")	{
				modeltype = 2;
				rrtype = "mtrev";
			}
			else if (s == "-mtart")	{
				modeltype = 2;
				rrtype = "mtart";
			}
			else if (s == "-mtzoa")	{
				modeltype = 2;
				rrtype = "mtzoa";
			}
			else if (s == "-cg6")	{
				mixturetype = 1;
				empmix = 1;
				fixncomp = 1;
				mixtype = "CG6";
			}
			else if (s == "-cg10")	{
				mixturetype = 1;
				empmix = 1;
				fixncomp = 1;
				mixtype = "CG10";
			}
			else if (s == "-cg20")	{
				mixturetype = 1;
				empmix = 1;
				fixncomp = 1;
				mixtype = "CG20";
			}
			else if (s == "-cg30")	{
				mixturetype = 1;
				empmix = 1;
				fixncomp = 1;
				mixtype = "CG30";
			}
			else if (s == "-cg40")	{
				mixturetype = 1;
				empmix = 1;
				fixncomp = 1;
				mixtype = "CG40";
			}
			else if (s == "-cg50")	{
				mixturetype = 1;
				empmix = 1;
				fixncomp = 1;
				mixtype = "CG50";
			}
			else if (s == "-cg60")	{
				mixturetype = 1;
				empmix = 1;
				fixncomp = 1;
				mixtype = "CG60";
			}
			else if (s == "-cgr6")	{
				modeltype = 2;
				rrtype = "CG6";
				mixturetype = 1;
				empmix = 1;
				fixncomp = 1;
				mixtype = "CG6";
			}
			else if (s == "-cgr10")	{
				modeltype = 2;
				rrtype = "CG10";
				mixturetype = 1;
				empmix = 1;
				fixncomp = 1;
				mixtype = "CG10";
			}
			else if (s == "-cgr20")	{
				modeltype = 2;
				rrtype = "CG20";
				mixturetype = 1;
				empmix = 1;
				fixncomp = 1;
				mixtype = "CG20";
			}
			else if (s == "-cgr30")	{
				modeltype = 2;
				rrtype = "CG30";
				mixturetype = 1;
				empmix = 1;
				fixncomp = 1;
				mixtype = "CG30";
			}
			else if (s == "-cgr40")	{
				modeltype = 2;
				rrtype = "CG40";
				mixturetype = 1;
				empmix = 1;
				fixncomp = 1;
				mixtype = "CG40";
			}
			else if (s == "-cgr50")	{
				modeltype = 2;
				rrtype = "CG50";
				mixturetype = 1;
				empmix = 1;
				fixncomp = 1;
				mixtype = "CG50";
			}
			else if (s == "-cgr60")	{
				modeltype = 2;
				rrtype = "CG60";
				mixturetype = 1;
				empmix = 1;
				fixncomp = 1;
				mixtype = "CG60";
			}
			
			else if ((s == "-dp") || (s == "-sbdp")	|| (s == "-cat")){
				mixturetype = 3;
			}
			/*
			else if (s == "-tdp")	{
				mixturetype = 4;
				i++;
				ncat = atoi(argv[i]);
			}
			*/
			else if (s == "-ss")	{
				mixturetype = 5;
			}
			else if (s == "-uni")	{
				type = Universal;
			}
			else if ((s == "-mtmam") || (s == "-MtMam") || (s == "mtvert") || (s == "MtVert"))	{
				type = MtMam;
			}
			else if (s == "-jeffkappa")	{
				kappaprior = 1;
			}
			else if (s == "-expkappa")	{
				kappaprior = 0;
			}
			else if (s == "-gamkappa")	{
				kappaprior = 2;
			}
			else if (s == "-rigidbaseprior")	{
				profilepriortype = 1;
			}
			else if (s == "-mvn")	{
				profilepriortype = 2;
			}
			else if (s == "-dirprofile")	{
				proposemode = 1;
			}
			else if (s == "-suballoc")	{
				allocmode = 1;
			}
			/*
			else if (s == "-condrate")	{
				sumratealloc = 0;
			}
			else if (s == "-sumrate")	{
				sumratealloc = 1;
			}
			*/
			else if (s == "-jeffomega")	{
				omegaprior = 1;
			}
			else if (s == "-omegafinite")	{
				omegamixturetype = 1;
			}
			else if (s == "-omegaspdp")	{
				omegamixturetype = 3;
			}
			else if (s == "-omegadp")	{
				omegamixturetype = 3;
			}
			else if (s == "-omegacat")	{
				omegamixturetype = 3;
			}
			else if (s == "-omegagamma")	{  // M5-like
				omegamixturetype = 4;
			}
			else if (s == "-nrep")	{
				i++;
				nrep = atoi(argv[i]);
			}
			else if ( (s == "-x") || (s == "-extract") )	{
				i++;
				if (i == argc) throw(0);
				every = atoi(argv[i]);
				i++;
				if (i == argc) throw(0);
				until = atoi(argv[i]);
			}
			else	{
				if (i != (argc -1))	{
					throw(0);
				}
				name = argv[i];
			}
			i++;
		}
		/*
		if ((datafile == "") && (argc != 2))	{
			throw(0);
		}
		*/
		/*
		if (nprocs <= 1)	{
			if (! myid)	{
				cerr << "error : pb_mpi requires at least 2 processes running in parallel (one master and at least one slave)\n";
			}
			MPI_Finalize();
			exit(1);
		}
		*/
	}
	catch(...)	{
		if (! myid)	{
			cerr << '\n';
			cerr << "mpirun -np <np> pb_mpi -d <datafile> [options] <chainname>\n";
			cerr << "\tcreates a new chain, sampling from the posterior distribution, conditional on specified data\n";
			cerr << "\n";
			cerr << "mpirun -np <np> pb_mpi <chainname>\n";
			cerr << "\tstarts an already existing chain\n";
			cerr << "\n";
			cerr << "\tmpirun -np <np>     : number of parallel processes (should be at least 2)\n";
			cerr << "\n";
			cerr << "\t-cat -dp            : infinite mixture (Dirichlet process) of equilibirium frequency profiles\n";
			cerr << "\t-ncat <ncat>        : finite mixture of equilibirium frequency profiles\n";
			cerr << "\t-catfix <pr>        : specifying a fixed pre-defined mixture of profiles\n";
			cerr << '\n';
			cerr << "\t-lg                 : Le and Gascuel 2008\n";
			cerr << "\t-wag                : Whelan and Goldman 2001\n";	
			cerr << "\t-jtt                : Jones, Taylor, Thornton 1992\n";	
			cerr << "\t-gtr                : general time reversible\n";
			cerr << "\t-poisson            : Poisson matrix, all relative exchangeabilities equal to 1 (Felsenstein 1981)\n";
			cerr << '\n';
			cerr << "\t-dgam <ncat>        : discrete gamma. ncat = number of categories (4 by default, 1 = uniform rates model)\n";
			cerr << "\t-igam <ncat>        : discrete gamma + invariable sites. ncat = number of categories (4 by default, 1 = uniform rates model)\n";
			cerr << '\n';
			cerr << "\t-dc                 : excludes constant columns\n";
			cerr << "\t-t <treefile>       : starts from specified tree\n"; 
			cerr << "\t-T <treefile>       : chain run under fixed, specified tree\n"; 
			cerr << '\n';
			cerr << "\t-x <every> <until>  : saving frequency, and chain length (until = -1 : forever)\n";
			cerr << "\t-f                  : forcing checks\n";
			cerr << "\t-s/-S               : -s : save all / -S : save only the trees\n";
			cerr << '\n';
			
			cerr << '\n';
			cerr << "see manual for details\n";
			cerr << '\n';

		}
		MPI_Finalize();
		exit(1);
	}

	if ((modeltype == -1) && (mixturetype == -1))	{
		modeltype = 2;
		mixturetype = 3;
	}
	else	{
		if (modeltype == -1)	{
			if (!myid)	{
			cerr << '\n';
			cerr << "error: incompletely specified model\n";
			cerr << "exchangeability parameters should be explicitly given\n";
			cerr << "-gtr -poisson (-f81) -lg -wag -jtt -mtrev -mtart -mtzoa or custom (-rr <filename>)\n";
			cerr << '\n';
			}
			MPI_Finalize();
			exit(1);
		}
		if (mixturetype == -1)	{
			if (!myid)	{
			cerr << '\n';
			cerr << "error: incompletely specified model\n";
			cerr << "mixture of equilibrium frequency profiles should be explicitly chosen among:\n";
			cerr << "-cat (or -dp) : infinite mixture (Dirichlet process)\n";
			cerr << "-ncat 1 : one matrix model\n";
			cerr << "-catfix <empmix>: empirical mixture (see manual for details)\n";
			cerr << '\n';
			}
			MPI_Finalize();
			exit(1);
		}
	}
	if (randfix != -1)	{
		rnd::init(1,randfix);
	}

	Model* model = 0;
	if (name == "")		{
		if (! myid)	{
			cerr << "error in command: no name was specified\n";
		}
		MPI_Finalize();
		exit(1);
	}
	if (datafile != "")	{
		if (myid == 0) {
			cerr << "model:\n";
			if (mixturetype == 1)	{
				if (empmix)	{
					cerr << "empirical mixture: " << mixtype << '\n';
				}
				else if (ncat == 1)	{
					cerr << "one-matrix model\n";
				}
				else	{
					cerr << "finite mixture of " << ncat << " components\n";
				}
			}
			else if (mixturetype == 3)	{
				cerr << "stick-breaking Dirichlet process mixture (cat)\n";
			}
			else if (mixturetype == 5)	{
				cerr << "site-specific profiles\n";
			}

			if (modeltype == 5)	{
				cerr << "codon mutation selection model\n";
			}
			else if (modeltype == 1)	{
				cerr << "exchangeabilities : f81 (Poisson)\n";
			}
			else if (modeltype == 2)	{
				if (rrtype == "None")	{
					cerr << "exchangeabilities estimated from data (gtr)\n";
				}
				else	{
					cerr << "exchangeabilities : " << rrtype << '\n';
				}
			}
			if (dgam == 1)	{
				cerr << "uniform rates across sites\n";
			}
			else	{
				if (modeltype != 3 && modeltype != 5)	{
					cerr << "discrete gamma distribution of rates across sites (" << dgam << " categories)\n";
				}
			}
			cerr << '\n';
		}
		if (! force)	{
			if (ifstream((name + ".param").c_str()))	{
				if (!myid)	{
					cerr << "a chain named " << name << " already exists; use -f to override\n";
					cerr << '\n';
				}
				MPI_Finalize();
				exit(1);
			}
		}
		model = new Model(datafile,treefile,partitionfile,multigene,globalalpha,globalbl,mappsuffstat,modeltype,dgam,pinv,mixturetype,nmodemax,ncat,type,suffstat,fixncomp,empmix,mixtype,dirpriortype,nstatcomp,priorempmix,priormixtype,fixstatweight,fixstatalpha,fixstatcenter,rrtype,iscodon,sis,sisfrac,sisnfrac,sisnrep,siscutoff,fixtopo,fixroot,roottax1,roottax2,topoburnin,topobf,bfburnin,bffrac,bfnfrac,bfnrep,blfactor,blfile,NSPR,NMHSPR,NTSPR,temperedbl,temperedgene,temperedrate,topolambda,topomu,toponstep,NNNI,nspec,ntspec,taxon1,taxon2,taxon3,taxon4,bpp,nbpp,ntbpp,bppnstep,bppname,bppcutoff,bppbeta,fixcodonprofile,fixomega,nomega,fixbl,sumovercomponents,omegaprior,omegamixturetype,kappaprior,profilepriortype,dc,every,until,saveall,zip,proposemode,allocmode,fasttopo,fasttopofracmin,fasttoponstep,fastcondrate,reshuffle,name,myid,nprocs);

		if (! myid)	{
			cerr << '\n';
			cerr << "chain name : " << name << '\n';
			// MPI master only
			ofstream os((name + ".treelist").c_str());
			if (NTSPR || fasttopo)	{
				ofstream tspros((name + ".temperedmove").c_str());
			}
			ofstream topos((name + ".topo").c_str());
			if (topobf)	{
				ofstream bos((name + ".bf").c_str());
			}
			if (sis)	{
				ofstream sos((name + ".sis").c_str());
			}
			ofstream tos((name + ".trace").c_str());
			model->TraceHeader(tos);
			tos.close();
			ofstream pos((name + ".param").c_str());
			model->ToStream(pos,true);
			pos.close();
			if (saveall)	{
				ofstream cos((name + ".chain").c_str());
			}
		}
	}
	else	{
		model = new Model(name,myid,nprocs);
		if (until != -1)	{
			model->until = until;
		}
	}

	if (myid == 0) {

		cerr << "run started\n";
		cerr << '\n';
		model->TraceHeader(cerr);
		model->Trace(cerr);
		cerr << '\n';
		model->Run(smc,deltansite,shortcycle,longcycle,cutoffsize,nrep);
		if (nprocs > 1)	{
			MESSAGE signal = KILL;
			MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
		}
	}
	else {
		// MPI slave
		model->WaitLoop();
	}
	MPI_Finalize();
}
