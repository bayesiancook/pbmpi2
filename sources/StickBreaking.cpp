#include "Random.h"
#include "BiologicalSequences.h"

#include "Parallel.h"
MPI_Datatype Propagate_arg;

// a mixture of ncat Dirichlet densities
// weight[k]: weights of component k
// conc[k]: concentration parameter of component k
// stat[k]: profile of compoment k

int ncat;
double* weight;
double* conc;
double** stat;

// sample a profile from the mixture
double SampleProfile(double* pi)	{

	// first, choose component
	int cat = rnd::GetRandom().DrawFromDiscreteDistribution(weight,ncat);

	// then draw Dirichlet random variate
	double tot = 0;
	for (int k=0; k<Naa; k++)	{
		pi[k] = rnd::GetRandom().sGamma(conc[cat] * stat[cat][k]);
		tot += pi[k];
	}

	// finally, calculate h-index
	double hi = 0;
	for (int k=0; k<Naa; k++)	{
		pi[k] /= tot;
		hi += pi[k] * HydrophobicityIndex_pH7[k];
	}

	return hi;
}

int main(int argc, char* argv[])	{

	if (argc == 1)	{
		cerr << "stickbreaking mixfile alpha N M name\n";
		cerr << "draws a truncated mixture from a Dirichlet process:\n";
		cerr << "alpha: stick-breaking parameter alpha\n";
		cerr << "base distribution is a mixture of Dirichlet densities specified in mixfile\n";
		cerr << "output:\n";
		cerr << "name.spikes: truncated (up to N) mixture\n";
		cerr << "name.sample: sample of M iid draws from this truncated mixture\n";
		cerr << "mixfile format:\n";
		cerr << "ncat\n";
		cerr << "weights\n";
		cerr << "concentration factor\n";
		cerr << "concentrations\n";
		cerr << "amino-acid frequency profiles\n";
		cerr << '\n';
		exit(1);
	}
	// read mixture of Dirichlet densities from file
	ifstream is(argv[1]);

	// concentration parameter of the stick-breaking process
	double alpha = atof(argv[2]);

	// 
	int N = atoi(argv[3]);
	int M = atoi(argv[4]);

	// output
	string name = argv[5];

	// read mixture of Dirichlet densities from file
	// number of compoments
	is >> ncat;

	// weights
	weight = new double[ncat];
	double totweight = 0;
	for (int i=0; i<ncat; i++)	{
		is >> weight[i];
		totweight += weight[i];
	}
	for (int i=0; i<ncat; i++)	{
		weight[i] /= totweight;
	}

	// concentrations
	double concfactor;
	is >> concfactor;
	conc = new double[ncat];
	for (int i=0; i<ncat; i++)	{
		is >> conc[i];
		conc[i] *= concfactor;
	}

	// mean profiles for each component
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

	// draw a truncated Dirichlet process (up to N components)
	// using stick-breaking process
	double** pi = new double*[N];
	double* w = new double[N];
	double* h = new double[N];
	for (int i=0; i<N; i++)	{
		pi[i] = new double[Naa];
	}

	// write infinite mixture configurations into .spikes file
	// each line:
	// weight, hi, profile
	// here, could also add effective size of state-space as well as ratio between two max amino-acids
	ofstream os((name + ".spikes").c_str());

	// draw from stick-breaking process
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

	// draw a sample of M independent draws from truncated Dirichlet process
	// can be used to plot an estimate of the cdf
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

