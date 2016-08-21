
#include "DirichletMixtureProfileProcess.h"
#include "Random.h"

void DirichletMixtureProfileProcess::GetPostCount(int* count)	{

	if (dirpriortype)	{
		cerr << "error: in GetPostCount, under simple dirichlet prior\n";
		exit(1);
	}
	for (int k=0; k<Nstatcomp; k++)	{
		count[k] = 0;
	}

	double* post = new double[Nstatcomp];

	for (int i=0; i<GetNcomponent(); i++)	{
		GetPostStatProb(profile[i],post);
		int k = rnd::GetRandom().FiniteDiscrete(Nstatcomp,post);
		count[k]++;
	}

	delete[] post;
}
