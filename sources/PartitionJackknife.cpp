
#include "Parallel.h"
#include "SequenceAlignment.h"
MPI_Datatype Propagate_arg;

int main(int argc, char* argv[])	{

	string datafile = argv[1];
	string partitionfile = argv[2];
	int N = atoi(argv[3]);
	int nrep = atoi(argv[4]);
	string name = argv[5];

	SequenceAlignment* data = new FileSequenceAlignment(datafile);
	int Nsite = data->GetNsite();

	ifstream is(partitionfile.c_str());
	int Ngene;
	is >> Ngene;
	int* genesize = new int[Ngene];
	int* genefirst = new int[Ngene];

	int count = 0;
	for (int gene=0; gene<Ngene; gene++)	{
		is >> genesize[gene];
		genefirst[gene] = count;
		count += genesize[gene];
	}
	if (count != Nsite)	{
		cerr << "error: non matching total size\n";
		cerr << count << '\t' << Nsite << '\n';
	}

	int* choose = new int[N];
	bool* genemask = new bool[Ngene];
	bool* mask = new bool[Nsite];
	for (int rep=0; rep<nrep; rep++)	{
		rnd::GetRandom().DrawFromUrn(choose,N,Ngene);
		for (int gene=0; gene<Ngene; gene++)	{
			genemask[gene] = false;
		}
		for (int k=0; k<N; k++)	{
			genemask[choose[k]] = true;
		}
		for (int gene=0; gene<Ngene; gene++)	{
			for (int i=0; i<genesize[gene]; i++)	{
				mask[genefirst[gene] + i] = genemask[gene];
			}
		}

		SequenceAlignment* replicate = new SequenceAlignment(data,mask);
		ostringstream s;
		s << name << rep;
		ofstream os((s.str() + ".ali").c_str());
		replicate->ToStream(os);
		ofstream os2((s.str() + ".partition").c_str());
		os2 << N << '\n';
		for (int gene=0; gene<Ngene; gene++)	{
			if (genemask[gene])	{
				os2 << genesize[gene] << '\t';
			}
		}
		os << '\n';
	}
}
		

