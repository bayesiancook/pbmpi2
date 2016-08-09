
#include "SequenceAlignment.h"

int main(int argc, char* argv[])	{

	string infile = argv[1];
	SequenceAlignment* data = new FileSequenceAlignment(infile);
	ofstream os(argv[2]);
	data->PrintWithoutAllMissingTaxa(os);
}
