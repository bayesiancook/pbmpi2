
#include "GTRSBDPProfileProcess.h"

/*
void GTRSBDPProfileProcess::ToStream(ostream& os)	{

	os << Ncomponent << '\n';
	os << kappa << '\n';
	for (int j=0; j<GetDim(); j++)	{
		os << dirweight[j] << '\t';
	}
	os << '\n';
	os << '\n';
	for (int i=0; i<GetNrr(); i++)	{
		os << rr[i] << '\t';
	}
	os << '\n';
	os << '\n';
	for (int i=0; i<Ncomponent; i++)	{
		for (int j=0; j<GetDim(); j++)	{
			os << profile[i][j] << '\t';
		}
		os << '\n';
	}
	for (int i=0; i<GetNsite(); i++)	{
		if (ActiveSite(i))	{
			os << alloc[i] << '\t';
		}
		else	{
			os << -1 << '\t';
		}
	}
	os << '\n';
}

void GTRSBDPProfileProcess::FromStream(istream& is)	{

	is >> Ncomponent;
	is >> kappa;
	
	for (int i=0; i<GetDim(); i++)	{
		is >> dirweight[i];
	}

	for (int i=0; i<GetNrr(); i++)	{
		is >> rr[i];
	}

	for (int i=0; i<Ncomponent; i++)	{
		for (int j=0; j<GetDim(); j++)	{
			is >> profile[i][j];
		}
	}

	for (int i=0; i<GetNsite(); i++)	{
		is >> alloc[i];
	}

	ResampleWeights();
}
*/

