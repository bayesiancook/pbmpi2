
#include "Random.h"

int main(int argc, char* argv[])	{

	int sisnfrac = atoi(argv[1]);
	int sisnrep = atoi(argv[2]);
	double prop = atof(argv[3]);
	string name = argv[4];

	ifstream is((name + ".sis").c_str());

	int b = (1-prop) * sisnrep;
	int n = sisnrep - b;

	double totvarlog = 0;
	double logbf = 0;

	for (int frac=0; frac<sisnfrac; frac++)	{

		double f = ((double) frac) / sisnfrac;

		double delta[n];

		for (int i=0; i<b; i++)	{
			double tmp1, tmp2, tmp3;
			is >> tmp1 >> tmp2 >> tmp3;
			if (fabs(tmp1 - f) > 1e-6)	{
				cerr << "error in bf 2: read " << tmp1 << " instead of " << f << '\n';
				exit(1);
			}
		}
		double max = 0;
		for (int i=0; i<n; i++)	{
			double tmp1, tmp2, tmp3;
			is >> tmp1 >> tmp2 >> tmp3;
			if (fabs(tmp1 - f) > 1e-6)	{
				cerr << "error in bf 2: read " << tmp1 << " instead of " << f << '\n';
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

	ofstream os((name + ".logp").c_str());
	os << logbf << '\n';
	os.close();

	cout << '\n';
	cout << "log marginal likelihood : " << logbf << '\n';
	cout << '\n';
	cout << "total log variance: " << totvarlog << '\n';
}
