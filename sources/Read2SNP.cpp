
#include "Random.h"
#include "BiologicalSequences.h"

int main(int argc, char* argv[])	{

	ifstream is(argv[1]);
	int N = atoi(argv[2]);
	int P = atoi(argv[3]);
	double epsilon = atof(argv[4]);
	ofstream os(argv[5]);

	int** count = new int*[P];
	for (int j=0; j<P; j++)	{
		count[j] = new int[4];
	}

	int* totcount = new int[4];
	double* pi = new double[4];

	double** genotype = new double*[P];
	for (int j=0; j<P; j++)	{
		genotype[j] = new double[10];
	}

	int* totreliable = new int[P];
	int* tothetero = new int[P];
	double* meanposthetero = new double[P];
	int* totmeanpostreliable = new int[P];
	for (int j=0; j<P; j++)	{
		totreliable[j] = 0;
		tothetero[j] = 0;
		totmeanpostreliable[j] = 0;
		meanposthetero[j] = 0;
	}

	double totlogp = 0;

	for (int i=0; i<N; i++)	{
		int tmp;
		is >> tmp;
		if (tmp != i)	{
			cerr << "error when reading new line\n";
			exit(1);
		}

		for (int k=0; k<4; k++)	{
			totcount[k] = 0;
		}
		for (int j=0; j<P; j++)	{
			for (int k=0; k<4; k++)	{
				is >> count[j][k];
				totcount[k] += count[j][k];
			}
		}
		double tot = 0;
		for (int k=0; k<4; k++)	{
			pi[k] = totcount[k] + 0.01;
			tot += pi[k];
		}
		for (int k=0; k<4; k++)	{
			pi[k] /= tot;
		}

		os << i << '\t';
		for (int j=0; j<P; j++)	{
			int n = 0;
			for (int k=0; k<4; k++)	{
				n += count[j][k];
			}
			double logpost[10];
			double max = 0;
			int gmax = 0;
			int g = 0;
			for (int k=0; k<4; k++)	{
				for (int l=k; l<4; l++)	{
					double logp = 0;
					if (k == l)	{
						logp += 2*log(pi[k]);
						logp += count[j][k]*log(1-epsilon) + (n-count[j][k])*log(epsilon/3);
					}
					else	{
						logp += log(2.0) + log(pi[k]) + log(pi[l]);
						logp += (count[j][k]+count[j][l])*log(0.5*(1-2*epsilon/3)) + (n-count[j][k]-count[j][l])*log(epsilon/3);
					}
					if ((!g) || (max < logp))	{
						max = logp;
						gmax = g;
					}
					logpost[g] = logp;
					g++;
				}
			}
			double post[10];
			double tot = 0;
			for (int g=0; g<10; g++)	{
				post[g] = exp(logpost[g] - max);
				tot += post[g];
			}
			for (int g=0; g<10; g++)	{
				post[g] /= tot;
			}
			totlogp += log(tot) + max;
			
			g = 0;
			if (n >= 10)	{
				totmeanpostreliable[j]++;
			}
			for (int k=0; k<4; k++)	{
				for (int l=k; l<4; l++)	{
					if (n >= 10)	{
						if (k != l)	{
							meanposthetero[j] += post[g];
						}
					}
					if (g == gmax)	{
						os << DNAletters[k] << DNAletters[l] << "(" << (int) (100 * post[g]) << ")\t";
						if (post[g] > 0.95)	{
							totreliable[j]++;
							if (k != l)	{
								tothetero[j]++;
							}
						}
					}
					g++;
				}
			}
		}
		os << '\n';
	}
	for (int j=0; j<P; j++)	{
		cerr << j << '\t' << ((double) tothetero[j]) / totreliable[j] << '\t' << meanposthetero[j] / totmeanpostreliable[j] << '\n';
	}
	cerr << "totlogp : " << totlogp << '\n';
}

