#include "Random.h"
#include "BiologicalSequences.h"

#include "Parallel.h"
MPI_Datatype Propagate_arg;

int main(int argc, char* argv[])	{

	ifstream is(argv[1]);
	int N = atoi(argv[2]);
	ofstream os(argv[3]);

	int ncat = 0;
	is >> ncat;
	double* weight = new double[ncat];
	double totweight = 0;
	for (int i=0; i<ncat; i++)	{
		is >> weight[i];
		totweight += weight[i];
	}
	for (int i=0; i<ncat; i++)	{
		weight[i] /= totweight;
	}

	double concfactor;
	is >> concfactor;
	double* conc = new double[ncat];
	for (int i=0; i<ncat; i++)	{
		is >> conc[i];
		conc[i] *= concfactor;
	}

	double** stat = new double*[ncat];
	for (int i=0; i<ncat; i++)	{
		stat[i] = new double[Naa];
		double tot = 0;
		for (int k=0; k<Naa; k++)	{
			is >> stat[i][k];
			tot += stat[i][k];
		}
		for (int k=0; k<Naa; k++)	{
			stat[i][k] /= tot;
		}
	}

	double* pi = new double[Naa];
	for (int i=0; i<N; i++)	{
		int cat = rnd::GetRandom().DrawFromDiscreteDistribution(weight,ncat);
		double tot = 0;
		for (int k=0; k<Naa; k++)	{
			pi[k] = rnd::GetRandom().sGamma(conc[cat] * stat[cat][k]);
			tot += pi[k];
		}
		double hi = 0;
		for (int k=0; k<Naa; k++)	{
			pi[k] /= tot;
			hi += pi[k] * HydrophobicityIndex_pH7[k];
		}
		os << hi;
		for (int k=0; k<Naa; k++)	{
			os << '\t' << pi[k];
		}
		os << '\n';
	}
}

