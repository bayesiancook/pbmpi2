
#include "ZipMatrixPhyloProcess.h"

void ZipMatrixPhyloProcess::Create()	{
	if (! zipdata)	{
		truedata = data;
		zipdata = new ZippedSequenceAlignment(truedata);
		data = zipdata;
		PhyloProcess::Create();
		ZipMatrixProfileProcess::Create();
		truedata->GetEmpiricalFreq(empfreq);
	}
}

void ZipMatrixPhyloProcess::Delete()	{
	if (zipdata)	{
		ZipMatrixProfileProcess::Delete();
		PhyloProcess::Delete();
		delete zipdata;
		zipdata = 0;
	}
}

