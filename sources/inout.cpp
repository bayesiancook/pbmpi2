#include "Parallel.h"
#include "SequenceAlignment.h"
MPI_Datatype Propagate_arg;

int main(int argc, char* argv[])	{

	string datafile = argv[1];

	SequenceAlignment* data = new FileSequenceAlignment(datafile);
	ofstream os(argv[2]);
	data->ToStream(os);
}
