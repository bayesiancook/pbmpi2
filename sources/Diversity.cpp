#include "Parallel.h"
#include "SequenceAlignment.h"
MPI_Datatype Propagate_arg;

int main(int argc, char* argv[])	{

	string datafile = argv[1];

	SequenceAlignment* data = new FileSequenceAlignment(datafile);
	cout << "alphabet size : " << data->GetNstate() << '\n';
	cout << "mean diversity: " << data->GetMeanDiversity() << '\n';

	int n = 0;
	for (int i=0; i<data->GetNsite(); i++)	{
		if (data->ConstantColumn(i))	{
			n++;
		}
	}
	cout << "constant columns : " << n << " / " << data->GetNsite() << '\t' << 100 * ((double) n) / data->GetNsite() << "%\n";
}
