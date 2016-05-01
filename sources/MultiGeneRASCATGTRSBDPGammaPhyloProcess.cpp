
#include "MultiGeneRASCATGTRSBDPGammaPhyloProcess.h"
#include "Parallel.h"


void MultiGeneRASCATGTRSBDPPhyloProcess::Create()	{
	ExpoConjugateGTRProfileProcess::Create();
	RASCATGTRSBDPSubstitutionProcess::Create();
	GammaBranchProcess::Create();
	if (GetMyid())	{
		for (int gene=0; gene<Ngene; gene++)	{
			if (genealloc[gene] == myid)	{
			}
		}
	}
}
	
void MultiGeneRASCATGTRSBDPPhyloProcess::Delete()	{
	GammaBranchProcess::Delete();
	RASCATGTRSBDPSubstitutionProcess::Delete();
	ExpoConjugateGTRProfileProcess::Delete();
	if (GetMyid())	{
		for (int gene=0; gene<Ngene; gene++)	{
			if (genealloc[gene] == myid)	{
			}
		}
	}
}

void MultiGeneRASCATGTRSBDPPhyloProcess::SlaveGeneMove()	{

	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			process[gene]->ExpoConjugateGTRSBDPProfileProcess::Move(1,1,10);
		}
	}
}

double MultiGeneRASCATGTRSBDPPhyloProcess::GlobalGetMeanNcomponent()	{
	return 0;
}

double MultiGeneRASCATGTRSBDPPhyloProcess::GlobalGetMeanStatEnt()	{
	return 0;
}

double MultiGeneRASCATGTRSBDPPhyloProcess::GlobalGetMeanStatAlpha()	{
	return 0;
}

double MultiGeneRASCATGTRSBDPPhyloProcess::GlobalGetMeanStatKappa()	{
	return 0;
}

void MultiGeneRASCATGTRSBDPPhyloProcess::SlaveGetMeanKappa() {}
void MultiGeneRASCATGTRSBDPPhyloProcess::SlaveGetMeanStatEnt() {}
void MultiGeneRASCATGTRSBDPPhyloProcess::SlaveGetMeanStatAlpha() {}
void MultiGeneRASCATGTRSBDPPhyloProcess::SlaveGetMeanKappa() {}

void MultiGeneRASCATGTRSBDPPhyloProcess::UpdateRRSuffStat() {}
void MultiGeneRASCATGTRSBDPPhyloProcess::UpdateRateSuffStat() {}
void MultiGeneRASCATGTRSBDPPhyloProcess::UpdateBranchLengthSuffStat() {}

void MultiGeneRASCATGTRSBDPPhyloProcess::GlobalUpdateParameters() {}
void MultiGeneRASCATGTRSBDPPhyloProcess::SlaveUpdateParameters() {}

void MultiGeneRASCATGTRSBDPPhyloProcess::SlaveExecute(MESSAGE)	{}

