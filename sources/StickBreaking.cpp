#include "Random.h"
#include "BiologicalSequences.h"

#include "Parallel.h"
MPI_Datatype Propagate_arg;

int ncat;
double* weight;
double* conc;
double** stat;

double SampleProfile(double* pi)	{

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

	return hi;
}

int main(int argc, char* argv[])	{

	ifstream is(argv[1]);
	double alpha = atof(argv[2]);
	int N = atoi(argv[3]);
	int M = atoi(argv[4]);
	string name = argv[5];

	is >> ncat;
	weight = new double[ncat];
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
	conc = new double[ncat];
	for (int i=0; i<ncat; i++)	{
		is >> conc[i];
		conc[i] *= concfactor;
	}

	stat = new double*[ncat];
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

	double** pi = new double*[N];
	double* w = new double[N];
	double* h = new double[N];
	for (int i=0; i<N; i++)	{
		pi[i] = new double[Naa];
	}
	ofstream os((name + ".spikes").c_str());

	double w0 = 1.0;
	for (int i=0; i<N; i++)	{
		double q1 = rnd::GetRandom().sGamma(1.0);
		double q2 = rnd::GetRandom().sGamma(alpha);
		double r = q1 / (q1 + q2);
		w[i] = w0 * r;
		w0 *= (1-r);
		h[i] = SampleProfile(pi[i]);
		os << w[i] << '\t' << h[i];
		for (int k=0; k<Naa; k++)	{
			os << '\t' << pi[i][k];
		}
		os << '\n';
	}

	ofstream sos((name + ".sample").c_str());
	for (int j=0; j<M; j++)	{
		int i = rnd::GetRandom().DrawFromDiscreteDistribution(w,N);
		if ((i<0) || (i>=N))	{
			cerr << "error: draw out of range\n";
			exit(1);
		}
		sos << h[i];
		for (int k=0; k<Naa; k++)	{
			sos << '\t' << pi[i][k];
		}
		sos << '\n';
	}
}

