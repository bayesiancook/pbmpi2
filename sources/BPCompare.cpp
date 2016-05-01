
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#include "CCP.h"
#include "StringStreamUtils.h"

#include "Parallel.h"
MPI_Datatype Propagate_arg;

int main(int argc, char* argv[])	{

	string* name = new string[argc];
	for (int i=0; i<argc; i++)	{
		name[i] = "";
	}

	string outfile = "bpcomp";
	double cutoff = 0.05;
	int burnin = -1;
	double conscutoff = 0.5;
	int verbose = 1;

	bool bench = false;

	if (argc == 1)	{
		cerr << "bpcomp [-covx] ChainName1 ChainName2 ... \n";
		cerr << "\t-c <cutoff> : only partitions with max prob >  cutoff. (default 0.5)\n";
		cerr << "\t-o <output> : output into file\n"; 
		cerr << "\t-v          : verbose output\n";
		cerr << "\t-x <burnin> [<every> <until>]. default burnin = 1/10 of the chain\n";
		cerr << '\n';
		cerr << "\t compare bipartition frequencies between independent chains\n";
		cerr << "\t and build consensus based on merged lists of trees\n";
		cerr << '\n';
		exit(1);
	}


	int i = 1;
	int P = 0; // number of chains to be compared

	int every = 1;
	int until = -1;

	while (i < argc)	{
		string s = argv[i];
		if (s == "-c")	{
			i++;
			s = argv[i];
			conscutoff = atof(argv[i]);
		}
		else if ( (s == "-x") || (s == "-extract") )	{
			i++;
			if (i == argc) throw(0);
			s = argv[i];
			if (! IsInt(s))	{
				throw(0);
			}
			burnin = atoi(argv[i]);
			i++;
			if (i == argc) throw(0);
			s = argv[i];
			if (IsInt(s))	{
				every = atoi(argv[i]);
				i++;
				if (i == argc) throw(0);
				s = argv[i];
				if (IsInt(s))	{
					until = atoi(argv[i]);
				}
				else	{
					i--;
				}
			}
			else	{
				i--;
			}
		}
		else if (s == "-v")	{
			verbose = 1;
		}
		else if (s == "-o")	{
			i++;
			outfile = argv[i];
		}
		else	{
			name[P] = argv[i];
			P++;
		}
		i++;
	}

	ofstream os((outfile + ".bplist").c_str());

	if (P == 1)	{

		BpList* bplist = new BpList(0);
		bplist->Read(name[0],burnin,every,until);
		const TaxonSet* taxset = bplist->GetTaxonSet();
		TpList* tplist = new TpList(taxset);
		tplist->Read(name[0],burnin,every,until);

		taxset->ToStream(os);
		os << '\n';
		bplist->ToStream(os);
		os << '\n';
		tplist->ToStream(os);
	}
	else	{

		MultiBpList* globalbplist = 0;
		const TaxonSet* taxset = 0;	
		for (int i=0; i<P; i++)	{
			cerr << name[i] << '\n';
			BpList* bplist = new BpList(taxset);
			bplist->Read(name[i],burnin,every,until);
			cerr << bplist->totweight << '\n';
			if (!i)	{
				taxset = bplist->GetTaxonSet();
				globalbplist = new MultiBpList(taxset);
				
			}
			globalbplist->Add(bplist);
		}
		globalbplist->Merge();
		// globalbplist->ToStreamSexy(os);

		MultiBpList* globaltplist = new MultiBpList(taxset);
		for (int i=0; i<P; i++)	{
			cerr << name[i] << '\n';
			TpList* tplist = new TpList(taxset);
			tplist->Read(name[i],burnin,every,until);
			cerr << tplist->totweight << '\n';
			globaltplist->Add(tplist);
		}
		globaltplist->Merge();

		taxset->ToStream(os);
		os << '\n';
		globalbplist->ToStream(os);
		os << '\n';
		globaltplist->ToStream(os);
	}
}

