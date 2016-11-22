
#include "Random.h"
#include "Tree.h"
#include "Parallel.h"
MPI_Datatype Propagate_arg;

int main(int argc, char* argv[])	{

	string treefile = argv[1];
	string roottax1 = argv[2];
	string roottax2 = argv[3];
	ofstream os(argv[4]);

	Tree* tree = new Tree(treefile);
	
	Link* newroot = tree->GetLCA(roottax1,roottax2);
	if (!newroot)	{
		cerr << "error when rerooting\n";
		exit(1);
	}
	tree->RootAt(newroot);
	tree->ToStream(os);
}
