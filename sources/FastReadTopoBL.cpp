
#include "Random.h"

int main(int argc, char* argv[])	{

	int bfnfrac = atoi(argv[1]);
	int bfnrep = atoi(argv[2]);
	double prop = atof(argv[3]);
	string name = argv[4];

	ifstream is((name + ".bf").c_str());

	int b = (1-prop) * bfnrep;
	int n = bfnrep - b;

	double totvarlog = 0;
	double logbf = 0;

	for (int frac=-bfnfrac; frac<bfnfrac-1; frac++)	{

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
			tot += exp(delta[i] - max);
		}
		tot /= n;

		double logscore = log(tot) + max;

		logbf += logscore;
	}

	cout << '\n';
	cout << "log bf : " << logbf << '\n';
	cout << '\n';
	cout << "total log variance: " << totvarlog << '\n';
}

