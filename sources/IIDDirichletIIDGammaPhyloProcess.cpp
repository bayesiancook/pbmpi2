
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "StringStreamUtils.h"
#include "IIDDirichletIIDGammaPhyloProcess.h"
#include "Parallel.h"
#include <string>

void IIDDirichletIIDGammaPhyloProcess::VarBayes()    {

    while(1)    {

        UpdateMeanSuffStat();
        UpdateVarLengths();
        UpdateVarRates();
        // UpdateVarProfiles();
        Trace(cerr);
    }
}

void IIDDirichletIIDGammaPhyloProcess::UpdateVarLengths() {

    cerr << "before length: " << GetMeanRate() << '\t' << GetTotalLength() << '\t' << GetMeanRate() * GetTotalLength() << '\n';
    double tot = 0;
    for (int j=1; j<GetNbranch(); j++)  {
        tot += branchlengthsuffstatcount[j];
    }
    cerr << "tot count: " << tot << '\t' << tot / GetNsite() << '\n';

    double A = GetMeanRate() * GetNsite();
	for (int j=0; j<GetNbranch(); j++)	{
		branchlengthsuffstatbeta[j] = A;
    }

    if (blvarmode == 0)   { // EM
        for (int j=1; j<GetNbranch(); j++)  {
            blarray[j] = branchlengthsuffstatcount[j] / branchlengthsuffstatbeta[j];
        }
    }
    else if (blvarmode == 1)   { // Mean Var Bayes
        for (int j=1; j<GetNbranch(); j++)  {
            blarray[j] = (branchalpha + branchlengthsuffstatcount[j]) / (branchbeta + branchlengthsuffstatbeta[j]);
        }
    }
    else if (blvarmode == 2) {
        for (int j=1; j<GetNbranch(); j++)  {
            blarray[j] = rnd::GetRandom().Psi(branchalpha + branchlengthsuffstatcount[j]) - log(branchbeta + branchlengthsuffstatbeta[j]);
        }
    }
    for (int j=1; j<GetNbranch(); j++)  {
        if (isnan(blarray[j]))  {
            cerr << "in varlength: nan\n";
            exit(1);
        }
        if (isinf(blarray[j]))  {
            cerr << "in var length: inf\n";
            exit(1);
        }
    }
    cerr << "after length: " << GetMeanRate() << '\t' << GetTotalLength() << '\t' << GetMeanRate() * GetTotalLength() << '\n';
}

void IIDDirichletIIDGammaPhyloProcess::UpdateVarRates() {

    double B = GetTotalLength();
    double totcount = 0;
    for (int i=0; i<GetNsite(); i++)    {
        siteratesuffstatbeta[i] = B;
        totcount += siteratesuffstatcount[i];
    }

    if (ratevarmode == 0)   {
        double tot = 0;
        for (int i=0; i<GetNsite(); i++)    {
            rate[i] = siteratesuffstatcount[i] / siteratesuffstatbeta[i];
            tot += rate[i];
        }
    }
    else if (ratevarmode == 1)   {
        for (int i=0; i<GetNsite(); i++)    {
            rate[i] = (alpha + siteratesuffstatcount[i]) / (alpha + siteratesuffstatbeta[i]);
        }
    }
    else if (ratevarmode == 2)   {
        for (int i=0; i<GetNsite(); i++)    {
            rate[i] = exp(rnd::GetRandom().Psi(alpha + siteratesuffstatcount[i]) - log(alpha + siteratesuffstatbeta[i]));
        }
    }
    else    {
        cerr << "rate var mode not recognized\n";
        exit(1);
    }
    for (int i=0; i<GetNsite(); i++)    {
        if (isnan(rate[i])) {
            cerr << "in varrates: nan\n";
            exit(1);
        }
        if (isinf(rate[i])) {
            cerr << "in varrates: inf\n";
            exit(1);
        }
    }
}

void IIDDirichletIIDGammaPhyloProcess::UpdateVarProfiles() {

    if (profilevarmode == 0)    {
        for (int i=0; i<GetNsite(); i++)    {
            double total = 0;
            for (int k=0; k<GetDim(); k++)   {
                profile[i][k] = siteprofilesuffstatcount[i][k];
                total += profile[i][k];
            }
            for (int k=0; k<GetDim(); k++)  {
                profile[i][k] /= total;
            }
        }
    }
    else if (profilevarmode == 1)   {
        for (int i=0; i<GetNsite(); i++)    {
            double total = 0;
            for (int k=0; k<GetDim(); k++)  {
                profile[i][k] = dirweight[k] + siteprofilesuffstatcount[i][k];
                total += profile[i][k];
            }
            for (int k=0; k<GetDim(); k++)  {
                profile[i][k] /= total;
            }
        }
    }
    else if (profilevarmode == 3)   {
        for (int i=0; i<GetNsite(); i++)    {
            double norm = 0;
            double tot = 0;
            for (int k=0; k<GetDim(); k++)  {
                profile[i][k] = exp(rnd::GetRandom().Psi(dirweight[k] + siteprofilesuffstatcount[i][k]));
                norm+= profile[i][k];
                tot += dirweight[k] + siteprofilesuffstatcount[i][k];
            }
            double z = exp(rnd::GetRandom().Psi(tot));
            if (fabs(z - norm) > 1e-6)  {
                cerr << "in var profiles: normalization error\n";
                exit(1);
            }
            for (int k=0; k<GetDim(); k++)  {
                profile[i][k] /= norm;
            }
        }
    }
}


void IIDDirichletIIDGammaPhyloProcess::GlobalUpdateParameters()	{

	if (GetNprocs() > 1)	{

        RenormalizeProfiles();

        int nd = 1 + GetNsite() + GetNbranch() + GetNsite()*GetDim() + GetDim();
        double* dvector = new double[nd];

        MESSAGE signal = PARAMETER_DIFFUSION;
        MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

        // First we assemble the vector of doubles for distribution
        int index = 0;
        dvector[index] = GetAlpha();
        index++;
        for (int i=0; i<GetNsite(); i++)    {
            dvector[index] = rate[i];
            index++;
        }

        for(int i=0; i<GetNbranch(); i++) {
            dvector[index] = blarray[i];
            index++;
        }
        
        for(int i=0; i<GetNsite(); ++i) {
            for(int j=0; j<GetDim(); ++j) {
                dvector[index] = profile[i][j];
                index++;
            }
        }
        for (int i=0; i<GetDim(); i++)	{
            dvector[index] = dirweight[i];
            index++;
        }

        // Now send out the doubles and ints over the wire...
        MPI_Bcast(dvector,nd,MPI_DOUBLE,0,MPI_COMM_WORLD);

        delete[] dvector;

	}
	else	{
		UpdateZip();
	}
}

void IIDDirichletIIDGammaPhyloProcess::SlaveUpdateParameters()	{

	int nd = 1 + GetNsite() + GetNbranch() + GetNsite()*GetDim() + GetDim();
    double* dvector = new double[nd];

	MPI_Bcast(dvector,nd,MPI_DOUBLE,0,MPI_COMM_WORLD);

	int index = 0;
    alpha = dvector[index];
    index++;

    for (int i=0; i<GetNsite(); i++)    {
        rate[i] = dvector[index];
        index++;
    }

	for(int i=0; i<GetNbranch(); ++i) {
		blarray[i] = dvector[index];
		index++;
	}

	for(int i=0; i<GetNsite(); ++i) {
		for(int j=0; j<GetDim(); ++j) {
			profile[i][j] = dvector[index];
			index++;
		}
	}
	for (int i=0; i<GetDim(); i++)	{
		dirweight[i] = dvector[index];
		index++;
	}

	delete[] dvector;

	UpdateZip();
}

void IIDDirichletIIDGammaPhyloProcess::SlaveExecute(MESSAGE signal)	{

	switch(signal) {

	case UPDATE_HYPERRATE:
		SlaveUpdateRateHyperSuffStat();
		break;
	case UPDATE_HYPERPROFILE:
		SlaveUpdateProfileHyperSuffStat();
		break;
    case MOVE_SRATE:
        ResampleSiteRates();
        break;
    case MOVE_SPROFILE:
        ResampleSiteProfiles();
        break;
    case COLLECT_SRATE:
        SlaveCollectSiteRates();
        break;
    case COLLECT_SPROFILE:
        SlaveCollectSiteProfiles();
        break;
	default:
		PhyloProcess::SlaveExecute(signal);
	}
}


