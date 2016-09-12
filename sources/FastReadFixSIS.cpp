
#include "Random.h"

int main(int argc, char* argv[])	{

	int sisnfrac = atoi(argv[1]);
	int b = atoi(argv[2]);
	int n = atoi(argv[3]);
	string name = argv[4];

	double totvarlog = 0;
	double logbf = 0;

	double* w = new double[n];

	for (int frac=0; frac<sisnfrac; frac++)	{

		ostringstream s;
		s << name << frac << ".sis";
		ifstream is(s.str().c_str());
		double f = ((double) frac) / sisnfrac;

		double delta[n];

		for (int i=0; i<b; i++)	{
			double tmp1, tmp2;
			is >> tmp1 >> tmp2;
			if (fabs(tmp1 - f) > 1e-6)	{
				cerr << "error in topo bf 2: read " << tmp1 << " instead of " << f << '\n';
				exit(1);
			}
		}

		double max = 0;
		for (int i=0; i<n; i++)	{
			double tmp1, tmp2;
			is >> tmp1 >> tmp2;
			if (fabs(tmp1 - f) > 1e-6)	{
				cerr << "error in topo bf 2: read " << tmp1 << " instead of " << f << '\n';
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

		cerr << f << '\t' << 1.0 / m2 << '\n';
	}

	cout << '\n';
	cout << "log p : " << logbf << '\n';
	cout << '\n';
	cout << "total log variance: " << totvarlog << '\n';
}