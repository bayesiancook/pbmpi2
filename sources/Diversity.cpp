#include "Parallel.h"
#include "SequenceAlignment.h"
MPI_Datatype Propagate_arg;

int main(int argc, char* argv[])	{

	string datafile = argv[1];

	SequenceAlignment* data = new FileSequenceAlignment(datafile);
	cout << "alphabet size : " << data->GetNstate() << '\n';
	cout << "mean diversity: " << data->GetMeanDiversity() << '\n';
	cout << "frac inv cols : " << ((double) data->GetNumberConstantColumns()) / data->GetNsite() << '\n';
}
