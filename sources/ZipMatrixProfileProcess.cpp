
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
		for (int i=0; i<GetNsite(); i++)	{
			delete zipmatrixarray[i];
			zipmatrixarray[i] = 0;
		}
		delete[] zipmatrixarray;
	}
	MatrixProfileProcess::Delete();
}


