#include "Parallel.h"
#include "SequenceAlignment.h"
MPI_Datatype Propagate_arg;

int main(int argc, char* argv[])	{

	string datafile = argv[1];

	SequenceAlignment* data = new FileSequenceAlignment(datafile);
	cout << "alphabet size : " << data->GetNstate() << '\n';
	cout << "mean diversity: " << data->GetMeanDiversity() << '\n';
	cout << "frac inv cols : " << ((double) data->GetNumberConstantColumns()) / data->GetNsite() << '\n';

	double* freq = new double[data->GetNstate()];
	data->GetEmpiricalFreq(freq);
	for (int i=0; i<data->GetNstate(); i++)	{
		cout << freq[i] << '\t';
	}
	cout << '\n';
}
