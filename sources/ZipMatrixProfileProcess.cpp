
#include "ZipMatrixProfileProcess.h"

void ZipMatrixProfileProcess::Create()	{
	if (GetMyid() && (! zipmatrixarray))	{
		MatrixProfileProcess::Create();
		zipmatrixarray = new ZipSubMatrix*[GetNsite()];
		for (int i=0; i<GetNsite(); i++)	{
			zipmatrixarray[i] = 0;
		}
	}
}

void ZipMatrixProfileProcess::Delete()	{
	if (zipmatrixarray)	{
		DeleteZipMatrices();
		delete[] zipmatrixarray;
	}
	MatrixProfileProcess::Delete();
}


