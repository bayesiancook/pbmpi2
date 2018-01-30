#include "Random.h"

int main(int argc, char* argv[])	{

	int Nsite = atoi(argv[1]);
	int Nind = atoi(argv[2]);
	int Nread = atof(argv[3]);
	double theta = atof(argv[4]);
	double epsilon = atof(argv[5]);
	string name = argv[6];

	const char dna[] = {'A','C','G','T'};

	double** freq = new double*[Nsite];
	for (int i=0; i<Nsite; i++)	{
		freq[i] = new double[4];
	}

	int** base1 = new int*[Nsite];
	int** base2 = new int*[Nsite];
	for (int i=0; i<Nsite; i++)	{
		base1[i] = new int[Nind];
		base2[i] = new int[Nind];
	}

	int*** data = new int**[Nsite];
	for (int i=0; i<Nsite; i++)	{
		data[i] = new int*[Nind];
		for (int j=0; j<Nind; j++)	{
			data[i][j] = new int[4];
		}
	}

	for (int i=0; i<Nsite; i++)	{
		double x = rnd::GetRandom().sGamma(theta);
		double y = rnd::GetRandom().sGamma(theta);
		double p = x /(x+y);
		int a = (int) (4 * rnd::GetRandom().Uniform());
		int b = (int) (3 * rnd::GetRandom().Uniform());
		if (b>=a)	{
			b++;
		}
		for (int k=0; k<4; k++)	{
			freq[i][k] = 0;
		}
		freq[i][a] = 1-p;
		freq[i][b] = p;
		for (int j=0; j<Nind; j++)	{
			int a1 = a;
			if (rnd::GetRandom().Uniform() < p)	{
				a1 = b;
			}
			int a2 = a;
			if (rnd::GetRandom().Uniform() < p)	{
				a2 = b;
			}
			if (a1 <= a2)	{
				base1[i][j] = a1;
				base2[i][j] = a2;
			}
			else	{
				base1[i][j] = a2;
				base2[i][j] = a1;
			}
			int* count = data[i][j];
			int nread = (int) (2*Nread*rnd::GetRandom().Uniform());
			for (int k=0; k<4; k++)	{
				count[k] = 0;
			}
			for (int k=0; k<nread; k++)	{
				int l = (int) (2*rnd::GetRandom().Uniform());
				int c = a1;
				if (l)	{
					c = a2;
				}
				if (rnd::GetRandom().Uniform() < epsilon)	{
					c = (int) (4 * rnd::GetRandom().Uniform());
				}
				count[c]++;
			}
		}
	}

	ofstream tos((name + ".true").c_str());
	for (int j=0; j<Nind; j++)	{
		tos << j << '\t';
		for (int i=0; i<Nsite; i++)	{
			tos << dna[base1[i][j]];
		}
		tos << '\n';
		tos << j << '\t';
		for (int i=0; i<Nsite; i++)	{
			tos << dna[base2[i][j]];
		}
		tos << '\n';
		tos << '\n';
	}
	
	ofstream os((name + ".counts").c_str());
	os << "ind\tbase";
	for (int i=0; i<Nsite; i++)	{
		os << '\t' << i;
	}
	os << '\n';
	for (int j=0; j<Nind; j++)	{
		for (int k=0; k<4; k++)	{
			os << j << '\t' << dna[k];
			for (int i=0; i<Nsite; i++)	{
				os << '\t' << data[i][j][k];
			}
			os << '\n';
		}
	}
}
