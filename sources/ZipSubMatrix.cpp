
#include "ZipSubMatrix.h"

void ZipSubMatrix::ComputeStationary()	{

	const double* fromstat = from->GetStationary();
	int zipsize = GetNstate();
	int nstate = from->GetNstate();

	for (int k=0; k<zipsize; k++)	{
		mStationary[k] = 0;
	}
	for (int k=0; k<nstate; k++)	{
		mStationary[zipindices[k]] += fromstat[k];
		
	}
}

void ZipSubMatrix::ComputeArray(int state)	{

	const double* fromstat = from->GetStationary();
	double** fromq = from->GetQ();
	const double* stat = GetStationary();

	int zipsize = GetNstate();
	int nstate = from->GetNstate();

	double num[zipsize];
	for (int b=0; b<zipsize; b++)	{
		num[b] = 0;
	}

	for (int j=0; j<nstate; j++)	{
		int a = zipindices[j];
		if (a == state)	{
			for (int k=0; k<nstate; k++)	{
				int b = zipindices[k];
				num[b] += fromstat[j] * fromq[j][k];
			}
		}
	}

	double total = 0;
	for (int b=0; b<zipsize; b++)	{
		if (state != b)	{
			Q[state][b] = num[b] / stat[state];
			total += Q[state][b];
		}
	}
	Q[state][state] = -total;

}
