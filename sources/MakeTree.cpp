#include "Tree.h"
#include "Random.h"

int main(int argc, char* argv[])	{

	string name = argv[1];
	Tree* tree = new Tree(name);
	tree->ToStream(cerr);

	cerr << "\n2\n";
	ifstream is(argv[1]);
	tree->ReadFromStream(is);
	tree->ToStream(cerr);
}

