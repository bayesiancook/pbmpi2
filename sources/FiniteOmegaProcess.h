
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef FINITEOMEGA_H
#define FINITEOMEGA_H

#include "MixtureOmegaProcess.h"

class FiniteOmegaProcess : public virtual MixtureOmegaProcess	{

	public: 
	//FiniteProfileProcess() : weight(0), fixncomp(false), empmix(0), Nfixcomp(0), statfix(0), empweight(0), statcons(0), Ncons(0), conscutoff(0)  {
	FiniteOmegaProcess() : omegaweightalpha(1.0), fixnomegacomp(true), empomegamix(0), Nomegafixcomp(0), omegafix(0), empomegaweight(0) {}

	virtual ~FiniteOmegaProcess() {}

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
	void ReadOmegaFix(string name);
	void SetOmegaFix();

	double omegaweightalpha;
	int Nomegafixcomp;
	double* omegafix;
	double* empomegaweight;
	bool fixnomegacomp;
	int empomegamix;
	string omegamixtype;
};

#endif
