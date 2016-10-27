
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef MULTIPLEOMEGA_H
#define MULTIPLEOMEGA_H

#include "MixtureOmegaProcess.h"

class MultipleOmegaProcess : public virtual MixtureOmegaProcess	{

	public: 
	MultipleOmegaProcess() : omegaweightalpha(1.0) {}

	virtual ~MultipleOmegaProcess() {}

	// overall sampling of the omega part of the model
	void SampleOmega();

	// moves
	double MPIMoveOmega(double tuning, int nrep);
	double NonMPIMoveOmega(double tuning, int nrep);

	virtual void SampleOmegaWeights();
	void ResampleOmegaWeights();

	protected:
	virtual void Create()	{
		if (! omega)	{
			MixtureOmegaProcess::Create();
		}
	}

	virtual void Delete()	{
		if (omega)	{
			MixtureOmegaProcess::Delete();
		}
	}

	double omegaweightalpha;
};

#endif
