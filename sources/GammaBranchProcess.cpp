
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "GammaBranchProcess.h"
#include "Random.h"

void GammaBranchProcess::ToStream(ostream& os)	{

	SetNamesFromLengths();
	tree->ToStream(os);
	for (int j=0; j<GetNbranch(); j++)	{
		os << branchmean[j] << '\t';
	}
	os << '\n';
	for (int j=0; j<GetNbranch(); j++)	{
		os << branchrelvar[j] << '\t';
	}
	os << '\n';
}

void GammaBranchProcess::ToStreamWithLengths(ostream& os, const Link* from)	{

	if (from->isLeaf())	{
		os << from->GetNode()->GetName();
	}
	else	{
		os << "(";
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			ToStreamWithLengths(os, link->Out());
			if (link->Next() != from)	{
				os << ",";
			}
		}
		os << ")";
	}
	if (! from->isRoot())	{
		os << ":" << blarray[from->GetBranch()->GetIndex()];
	}
}


void GammaBranchProcess::FromStream(istream& is)	{

	tree->ReadFromStream(is);
	SetLengthsFromNames();
	for (int j=0; j<GetNbranch(); j++)	{
		is >> branchmean[j];
	}
	for (int j=0; j<GetNbranch(); j++)	{
		is >> branchrelvar[j];
	}
}
	
double GammaBranchProcess::LogBranchLengthPrior(const Branch* branch)	{
	int index = branch->GetIndex();
	double branchbeta = 1.0 / branchmean[index];
	double branchalpha = 1.0 / branchrelvar[index];
	return branchalpha * log(branchbeta) - rnd::GetRandom().logGamma(branchalpha) + (branchalpha-1) * log(blarray[index]) - branchbeta * blarray[index];
}

void GammaBranchProcess::SampleBranchLength(const Branch* branch)	{
	int index = branch->GetIndex();
	double branchbeta = 1.0 / branchmean[index];
	double branchalpha = 1.0 / branchrelvar[index];
	blarray[index] = rnd::GetRandom().Gamma(branchalpha,branchbeta);
}
	
void GammaBranchProcess::PriorSampleLength()	{
	if (! fixbl)	{
		double mean = 0.1 * rnd::GetRandom().sExpo();
		double relvar = rnd::GetRandom().sExpo();
		SetHyperParameters(mean,relvar);
		SampleLength();
	}
}

double GammaBranchProcess::LogHyperPrior()	{
	double branchbeta = 1.0 / branchmean[0];
	double branchalpha = 1.0 / branchrelvar[0];
	
	double total = -branchalpha;
	if (betaprior == 1)	{
		total -= log(branchbeta);
		if ((branchbeta < 1e-4) || (branchbeta > 1e4))	{
			total -= 1.0 / 0;
			// total = InfProb;
		}
	}
	else	{
		total -= 0.1 * branchbeta;
	}
	return total;
	// return -branchalpha - 0.1 * branchbeta ; // + 2 * log(branchbeta);
	// return -branchalpha - 10.0*branchbeta ; // + 2 * log(branchbeta);
}

double GammaBranchProcess::Move(double tuning, int nrep)	{
	double total = MoveLength();
	if (nrep)	{
		total += MoveBranchBeta(tuning,nrep);
	}
	return total;
}

double GammaBranchProcess::MoveBranchBeta(double tuning, int nrep)	{
	int Naccepted = 0;
	for (int rep=0; rep<nrep; rep++)	{
		double deltalogprob = - LogHyperPrior() - LogLengthPrior() - LengthSuffStatLogProb();
		double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
		double e = exp(m);
		double mean = branchmean[0];
		double relvar = branchrelvar[0];
		mean *= e;
		SetHyperParameters(mean,relvar);
		deltalogprob += LogHyperPrior() + LogLengthPrior() + LengthSuffStatLogProb();
		deltalogprob += m;
		int accepted = (log(rnd::GetRandom().Uniform()) < deltalogprob);
		if (accepted)	{
			Naccepted ++;
		}
		else	{
			mean /= e;
			SetHyperParameters(mean,relvar);
		}
	}
	return ((double) Naccepted) / nrep;
}

double GammaBranchProcess::MoveLength()	{

	GlobalUpdateBranchLengthSuffStat();
	for (int i=1; i<GetNbranch(); i++)	{
		double branchbeta = 1.0 / branchmean[i];
		double branchalpha = 1.0 / branchrelvar[i];
		blarray[i] = rnd::GetRandom().Gamma(branchalpha + GetBranchLengthSuffStatCount(i), branchbeta + GetBranchLengthSuffStatBeta(i));
	}
	return 1.0;
}

double GammaBranchProcess::NonMPIMoveLength()	{

	UpdateBranchLengthSuffStat();
	for (int i=1; i<GetNbranch(); i++)	{
		double branchbeta = 1.0 / branchmean[i];
		double branchalpha = 1.0 / branchrelvar[i];
		blarray[i] = rnd::GetRandom().Gamma(branchalpha + GetBranchLengthSuffStatCount(i), branchbeta + GetBranchLengthSuffStatBeta(i));
	}
	return 1.0;
}

double GammaBranchProcess::NonMPIMove(double tuning, int nrep)	{
	NonMPIMoveLength();
	MoveBranchBeta(tuning,nrep);
}

