
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

void GammaBranchProcess::Create()	{

	BranchProcess::Create();
	if ((!branchmean) && hierarchicallengthprior)	{
		branchmean = new double[GetNbranch()];
		branchrelvar = new double[GetNbranch()];
	}
}

void GammaBranchProcess::Delete()	{

	if (branchmean)	{
		delete[] branchmean;
		delete[] branchrelvar;
		branchmean = 0;
	}
	BranchProcess::Delete();
}

void GammaBranchProcess::ToStream(ostream& os)	{

	SetNamesFromLengths();
	tree->ToStream(os);
	// if hierarchical prior, the hyperparameters are already saved to stream by the multi gene branch process
	if (!hierarchicallengthprior)	{
		os << branchalpha << '\n';
		os << branchbeta << '\n';
	}
	/*
	if (hierarchicallengthprior)	{
		for (int j=0; j<GetNbranch(); j++)	{
			os << branchmean[j] << '\t';
		}
		os << '\n';
		for (int j=0; j<GetNbranch(); j++)	{
			os << branchrelvar[j] << '\t';
		}
		os << '\n';
	}
	else	{
		os << branchalpha << '\n';
		os << branchbeta << '\n';
	}
	*/
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
	if (! hierarchicallengthprior)	{
		is >> branchalpha;
		is >> branchbeta;
	}
	/*
	if (hierarchicallengthprior)	{
		for (int j=0; j<GetNbranch(); j++)	{
			is >> branchmean[j];
		}
		for (int j=0; j<GetNbranch(); j++)	{
			is >> branchrelvar[j];
		}
	}
	else	{
		is >> branchalpha;
		is >> branchbeta;
	}
	*/
}

void GammaBranchProcess::SampleBranchLength(const Branch* branch)	{
	int index = branch->GetIndex();
	if (hierarchicallengthprior)	{
		double branchbeta = 1.0 / branchmean[index];
		double branchalpha = 1.0 / branchrelvar[index];
		blarray[index] = rnd::GetRandom().Gamma(branchalpha,branchbeta);
	}
	blarray[index] = rnd::GetRandom().Gamma(branchalpha,branchbeta);
}
	
void GammaBranchProcess::PriorSampleLengthHyperParameters()	{

	branchalpha = rnd::GetRandom().sExpo();
	branchbeta = 10 * rnd::GetRandom().sExpo();
}

void GammaBranchProcess::SampleLengthHyperParameters()	{

	branchalpha = 1.0;
	branchbeta = 10.0;
}
	
double GammaBranchProcess::LogBranchLengthPrior(const Branch* branch)	{

	int index = branch->GetIndex();
	if (hierarchicallengthprior)	{
		double branchbeta = 1.0 / branchmean[index];
		double branchalpha = 1.0 / branchrelvar[index];
		return branchalpha * log(branchbeta) - rnd::GetRandom().logGamma(branchalpha) + (branchalpha-1) * log(blarray[index]) - branchbeta * blarray[index];
	}
	return branchalpha * log(branchbeta) - rnd::GetRandom().logGamma(branchalpha) + (branchalpha-1) * log(blarray[index]) - branchbeta * blarray[index];
}

double GammaBranchProcess::LogLengthHyperPrior()	{

	if (hierarchicallengthprior)	{
		cerr << "error: in default GammaBranchProcess::LogHyperPrior, with hierarchical prior\n";
		exit(1);
	}
	double total = -branchalpha;
	if (betaprior == 1)	{
		total -= log(branchbeta);
		if ((branchbeta < 1e-4) || (branchbeta > 1e4))	{
			total -= 1.0 / 0;
		}
	}
	else	{
		total -= 0.1 * branchbeta;
	}
	return total;
}

double GammaBranchProcess::Move(double tuning, int nrep)	{

	double total = 0;
	if (nrep && (! hierarchicallengthprior))	{
		total += MoveLengthHyperParameters(tuning,nrep);
		total += MoveLengthHyperParameters(tuning,0.1*nrep);
	}
	total += MPIMoveBranchLengths();
	if (nrep && (! hierarchicallengthprior))	{
		total += MoveLengthHyperParameters(tuning,0.1*nrep);
		total += MoveLengthHyperParameters(tuning,nrep);
	}
	return total;
}

double GammaBranchProcess::NonMPIMove(double tuning, int nrep)	{

	NonMPIMoveBranchLengths();
	if (! hierarchicallengthprior)	{
		MoveLengthHyperParameters(tuning,nrep);
	}
}

double GammaBranchProcess::MoveLengthHyperParameters(double tuning, int nrep)	{

	if (hierarchicallengthprior)	{
		cerr << "error: in GammaBranchProcess::MoveLengthHyperParams with hierarchical prior\n";
		exit(1);
	}
	int Naccepted = 0;
	for (int rep=0; rep<nrep; rep++)	{
		double deltalogprob = - LogLengthHyperPrior() - LogLengthPrior();
		// double deltalogprob = - LogLengthHyperPrior() - LogLengthPrior() - LengthSuffStatLogProb();
		double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
		double e = exp(m);
		branchbeta *= e;
		deltalogprob += LogLengthHyperPrior() + LogLengthPrior();
		// deltalogprob += LogLengthHyperPrior() + LogLengthPrior() + LengthSuffStatLogProb();
		deltalogprob += m;
		int accepted = (log(rnd::GetRandom().Uniform()) < deltalogprob);
		if (accepted)	{
			Naccepted ++;
		}
		else	{
			branchbeta /= e;
		}
	}
	return ((double) Naccepted) / nrep;
}

double GammaBranchProcess::MPIMoveBranchLengths()	{

	GlobalUpdateBranchLengthSuffStat();
	if (hierarchicallengthprior)	{
		for (int i=1; i<GetNbranch(); i++)	{
			double branchbeta = 1.0 / branchmean[i];
			double branchalpha = 1.0 / branchrelvar[i];
			blarray[i] = rnd::GetRandom().Gamma(branchalpha + GetBranchLengthSuffStatCount(i), branchbeta + GetBranchLengthSuffStatBeta(i));
			if (blarray[i] == 0)	{
				blarray[i] = 1e-6;
				/*
				cerr << "sampling a null branch length\n";
				cerr << branchalpha << '\t' << GetBranchLengthSuffStatCount(i) << '\t' << branchbeta << '\t' << GetBranchLengthSuffStatBeta(i) << '\n';
				exit(1);
				*/
			}
		}
	}
	else	{
		for (int i=1; i<GetNbranch(); i++)	{
			blarray[i] = rnd::GetRandom().Gamma(branchalpha + GetBranchLengthSuffStatCount(i), branchbeta + GetBranchLengthSuffStatBeta(i));
		}
	}
	return 1.0;
}

double GammaBranchProcess::NonMPIMoveBranchLengths()	{

	UpdateBranchLengthSuffStat();
	if (hierarchicallengthprior)	{
		for (int i=1; i<GetNbranch(); i++)	{
			double branchbeta = 1.0 / branchmean[i];
			double branchalpha = 1.0 / branchrelvar[i];
			blarray[i] = rnd::GetRandom().Gamma(branchalpha + GetBranchLengthSuffStatCount(i), branchbeta + GetBranchLengthSuffStatBeta(i));
		}
	}
	else	{
		for (int i=1; i<GetNbranch(); i++)	{
			blarray[i] = rnd::GetRandom().Gamma(branchalpha + GetBranchLengthSuffStatCount(i), branchbeta + GetBranchLengthSuffStatBeta(i));
		}
	}
	return 1.0;
}

