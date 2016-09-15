
#include "ZipExpoConjugateGTRPhyloProcess.h"

void ZipExpoConjugateGTRPhyloProcess::Unfold()	{

	DeleteMappings();
	ActivateSumOverRateAllocations();
	ActivateZip();
	if (GetMyid() && (! zipmatcreated))	{
		CreateZipMatrices();
		zipmatcreated = 1;
	}
	UpdateMatrices();
	UpdateConditionalLikelihoods();
	/*
	if (!sumratealloc)	{
		DrawAllocations(0);
		InactivateSumOverRateAllocations();
	}
	*/
	activesuffstat = false;
}

void ZipExpoConjugateGTRPhyloProcess::Collapse()	{

	InactivateZip();
	MatrixPhyloProcess::Collapse();
}


void ZipExpoConjugateGTRPhyloProcess::GlobalActivateFastTopo()	{

	GlobalActivateZip();
	PhyloProcess::GlobalActivateFastTopo();
}

void ZipExpoConjugateGTRPhyloProcess::GlobalInactivateFastTopo()	{

	GlobalInactivateZip();
	PhyloProcess::GlobalInactivateFastTopo();
}

/*
double ZipExpoConjugateGTRPhyloProcess::ZipTopoMoveCycle(int nrep, double tuning)	{

	ziptopotry++;

	GlobalBackupTree();

	// GlobalUpdateConditionalLikelihoods();
	double logp1 = logL;

	GlobalActivateZip();
	GlobalUpdateConditionalLikelihoods();
	double sublogp1 = logL;
	double preacc = SimpleTopoMoveCycle(nrep,tuning);
	GlobalUpdateConditionalLikelihoods();
	double sublogp2 = logL;
	
	GlobalInactivateZip();
	GlobalUpdateConditionalLikelihoods();
	double logp2 = logL;


	double logratio = (logp2 - logp1) - (sublogp2 - sublogp1);
	int accepted = (log(rnd::GetRandom().Uniform()) < logratio);

	// cerr << accepted << '\t' << logp2 - logp1 - (sublogp2 - sublogp1) << '\t' << logp2 - logp1 << '\t' << sublogp2 - sublogp1 << '\n';

	if (accepted)	{
		ziptopoacc++;
	}
	else	{
		GlobalRestoreTree();
		// GlobalUpdateConditionalLikelihoods();
	}

	return ((double) accepted);
}
*/

