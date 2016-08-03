
#include "Random.h"

int main(int argc, char* argv[])	{

	ifstream is(argv[1]);
	int N = atoi(argv[2]);
	double* val = new double[N];
	for (int i=0; i<N; i++)	{
		is >> val[i];
	}

	double max = val[0];
	for (int i=0; i<N; i++)	{
		if (max < val[i])	{
			max = val[i];
		}
	}

	double meanexp = 0;
	double mean = 0;
	double var = 0;
	double* weight = new double[N];
	for (int i=0; i<N; i++)	{
		mean += val[i];
		meanexp += exp(val[i] - max);
		weight[i] = exp(val[i] - max);
		var += val[i] * val[i];
	}
	double sum2 = 0;
	for (int i=0; i<N; i++)	{
		weight[i] /= meanexp;
		sum2 += weight[i] * weight[i];
	}
	mean /= N;
	var /= N;
	var -= mean*mean;

	meanexp /= N;
	double logmeanexp = log(meanexp) + max;

	cout << "log mean exp   : " << logmeanexp << '\n';
	cout << "eff sample size: " << 1.0 / sum2 << '\n';
	cout << "log mean       : " << mean << '\t' <<  var << '\n';
	cout << "debiased mean  : " << mean + 0.5*var << '\n';
}
