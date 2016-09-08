
#include "Random.h"

int main(int argc, char* argv[])	{

	int bfnfrac = atoi(argv[1]);
	int b = atoi(argv[2]);
	int n = atoi(argv[3]);
	string name = argv[4];

	double totvarlog = 0;
	double logbf = 0;

	double* w = new double[n];

	for (int frac=-bfnfrac; frac<bfnfrac-1; frac++)	{

		ostringstream s;
		s << name << frac << ".bf";
		ifstream is(s.str().c_str());

		double delta[n];

		for (int i=0; i<b; i++)	{
			int tmp1;
			double tmp2, tmp3;
			is >> tmp1 >> tmp2 >> tmp3;
			if (tmp1 != frac)	{
				cerr << "error in topo bf 2: read " << tmp1 << " instead of " << frac << '\n';
				exit(1);
			}
		}

		double max = 0;
		for (int i=0; i<n; i++)	{
			int tmp1;
			double tmp2, tmp3;
			is >> tmp1 >> tmp2 >> tmp3;
			if (tmp1 != frac)	{
				cerr << "error in topo bf 2: read " << tmp1 << " instead of " << frac << '\n';
				exit(1);
			}
			delta[i] = tmp2;
			if ((!i) || (max < tmp2))	{
				max = tmp2;
			}
		}

		double meanlog = 0;
		double varlog = 0;
		for (int i=0; i<n; i++)	{
			meanlog += delta[i];
			varlog += delta[i]*delta[i];
		}
		meanlog /= n;
		varlog /= n;
		varlog -= meanlog*meanlog;
		totvarlog += varlog;

		double tot = 0;
		for (int i=0; i<n; i++)	{
			w[i] = exp(delta[i] - max);
			tot += w[i];
		}
		double m2 = 0;
		for (int i=0; i<n; i++)	{
			w[i] /= tot;
			m2 += w[i] * w[i];
		}

		tot /= n;

		double logscore = log(tot) + max;

		logbf += logscore;

		cerr << frac << '\t' << 1.0 / m2 << '\n';
	}

	cout << '\n';
	cout << "log bf : " << logbf << '\n';
	cout << '\n';
	cout << "total log variance: " << totvarlog << '\n';
}
