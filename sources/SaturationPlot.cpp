
#include "Random.h"
#include "SequenceAlignment.h"
#include "Tree.h"


#include "Parallel.h"
MPI_Datatype Propagate_arg;

int main(int argc, char* argv[])	{

	SequenceAlignment* data = new FileSequenceAlignment(argv[1]);
	Tree* tree = new Tree(argv[2]);
	ofstream os(argv[3]);

	tree->RegisterWith(data->GetTaxonSet());

	map<pair<string,string>,double> treedist;
	tree->GetPairwiseDistances(treedist);

	map<pair<string,string>,double> seqdist;
	data->GetPairwiseDistances(seqdist);

	cerr << "size: " << treedist.size() << '\t' << seqdist.size() << '\n';

	for (map<pair<string,string>,double>::iterator i=treedist.begin(); i!=treedist.end(); i++)	{
		if (seqdist[i->first])	{
			os << i->first.first << '\t' << i->first.second << '\t' << i->second << '\t' << seqdist[i->first] << '\n';
		}
	}
}

