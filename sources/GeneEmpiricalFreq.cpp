
#include "SequenceAlignment.h"
#include "Random.h"

#include "Parallel.h"
MPI_Datatype Propagate_arg;


int main(int argc, char* argv[])	{

	string datafile = argv[1];
	string partitionfile = argv[2];
	string outfile = argv[3];

	SequenceAlignment* protdata = new FileSequenceAlignment(datafile);
	if (protdata->GetNstate() != Naa)	{
		cerr << "error: should be protein datafile\n";
		exit(1);
	}
	StateSpace* statespace = protdata->GetStateSpace();

	int Nsite = protdata->GetNsite();

	ifstream pis(partitionfile.c_str());
	int Ngene;
	pis >> Ngene;

	int* genesize = new int[Ngene];
	int* genefirst = new int[Ngene];
	SequenceAlignment** genedata = new SequenceAlignment*[Ngene];
	int count = 0;
	for (int gene=0; gene<Ngene; gene++)	{
		pis >> genesize[gene];
		genefirst[gene] = count;
		count += genesize[gene];
		genedata[gene] = new SequenceAlignment(protdata,genefirst[gene],genesize[gene]);
	}
	if (count != Nsite)	{
		cerr << "error: non matching total size\n";
		cerr << "sum over genes: " << count << '\n';
		cerr << "concatenation : " << Nsite << '\n';
		exit(1);
	}

	cerr << "number of genes: " << Ngene << '\n';
	ofstream os(argv[3]);

	double* freq = new double[Naa];
	for (int gene=0; gene<Ngene; gene++)	{
		genedata[gene]->GetEmpiricalFreq(freq);
		os << gene;
		for (int k=0; k<Naa; k++)	{
			os << '\t' << freq[k];
		}
		os << '\n';
	}
}

