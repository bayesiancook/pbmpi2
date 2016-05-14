
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "OmegaProcess.h"
#include "Random.h"

double SingleOmegaProcess::LogOmegaPrior()        {

	if (omegaprior == 0)	{
		// ratio of exponential random variables
		return -2 * log(1 + *omega);
	}
	else	{	
		// jeffreys
		return -log(*omega);
	}
}

void SingleOmegaProcess::SampleOmega()        {

	*omega = 1.0;
}

double SingleOmegaProcess::MoveOmega(double tuning)        {

	int naccepted = 0;
	double bkomega = *omega;
	double deltalogprob = -LogOmegaPrior() - OmegaSuffStatLogProb();

	double h = tuning * (rnd::GetRandom().Uniform() -0.5);
	double e = exp(h);
	*omega *= e;

	UpdateOmega();
	deltalogprob += h;
	deltalogprob += LogOmegaPrior() + OmegaSuffStatLogProb();

	int accepted = (rnd::GetRandom().Uniform() < exp(deltalogprob));
	if (accepted)	{
		naccepted++;	
	}
	else	{
		*omega = bkomega;
		UpdateOmega();
	}
	return naccepted;	
}
