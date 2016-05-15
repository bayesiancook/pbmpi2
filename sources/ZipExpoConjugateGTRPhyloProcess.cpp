
#include "ZipExpoConjugateGTRPhyloProcess.h"


void ZipExpoConjugateGTRPhyloProcess::Unfold()	{

	DeleteMappings();
	ActivateSumOverRateAllocations();
	CreateMatrices();
	// ActivateZip();
	UpdateConditionalLikelihoods();
	/*
	if (!sumratealloc)	{
		DrawAllocations(0);
		InactivateSumOverRateAllocations(ratealloc);
	}
	*/
	activesuffstat = false;
}

void ZipExpoConjugateGTRPhyloProcess::Collapse()	{

	// InactivateZip();
	UpdateConditionalLikelihoods();
	// if (sumratealloc)	{
	DrawAllocations(0);
	InactivateSumOverRateAllocations(ratealloc);
	// }
	SampleNodeStates();
	FillMissingMap();
	SampleSubstitutionMappings(GetRoot());
	DeleteMatrices();
	activesuffstat = true;
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

