
#include "ZipGeneralPathSuffStatRASCATGTRFiniteGammaPhyloProcess.h"
#include "Parallel.h"

void ZipGeneralPathSuffStatRASCATGTRFiniteGammaPhyloProcess::SlaveExecute(MESSAGE signal)	{

	switch(signal) {

	case PROFILESSLOGPROB:
		SlaveProfileSuffStatLogProb();
		break;

	default:
		GeneralPathSuffStatRASCATGTRFiniteGammaPhyloProcess::SlaveExecute(signal);
	}
}

