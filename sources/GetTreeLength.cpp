
#include "Random.h"
#include "Tree.h"
#include "Parallel.h"
MPI_Datatype Propagate_arg;

double RecursiveGetLength(const Link* from)	{
	double totlength = 0;
	if (! from->isRoot())	{
		double l = atof(from->GetBranch()->GetName().c_str());
		if (l<0)	{
			cerr << "negative branch length : " << from->GetBranch()->GetName() << '\n';
			exit(1);
		}
		totlength += l;
	}
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		totlength += RecursiveGetLength(link->Out());
	}
	return totlength;
}

int main(int argc, char* argv[])	{

	string treefile = argv[1];
	Tree* tree = new Tree(treefile);
	
	double length = RecursiveGetLength(tree->GetRoot());
	cout << "tree length : " << length << '\n';
}
