#include "Random.h"
#include "BiologicalSequences.h"

#include "Parallel.h"
MPI_Datatype Propagate_arg;

int main(int argc, char* argv[])	{

	ifstream is(argv[1]);
	int N = atoi(argv[2]);
	ofstream os(argv[3]);
	double* pi = new double[Naa];
	for (int i=0; i<N; i++)	{
		
		double tot = 0;
		for (int k=0; k<Naa; k++)	{
			is >> pi[k];
			tot += pi[k];
		}
		if (fabs(tot-1)>1e-4)	{
			cerr << "error when reading profile: does not sum to 1\n";
			cerr << tot << '\t' << tot-1 << '\n';
			exit(1);
		}
		for (int k=0; k<Naa; k++)	{
			pi[k] /= tot;
		}
		os << GetHI7(pi) << '\t' << GetMolWeight(pi) << '\t' << GetEffSize(pi) << '\t' << GetTwoMaxRatio(pi) << '\t' << GetMinMaxRatio(pi);
		for (int k=0; k<Naa; k++)	{
			os << '\t' << pi[k];
		}
		os << '\n';
	}
}

