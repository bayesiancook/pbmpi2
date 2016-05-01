
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef MULTIGENEPHYLOPROCESS_H
#define MULTIGENEPHYLOPROCESS_H

#include "PhyloProcess.h"

class MultiGeneMPIModule	{

	public:

	MultiGeneMPIModule() : Ngene(0)	{}
	virtual ~MultiGeneMPIModule() {}

	void AllocateAlignments(string datafile, string treefile);

	virtual StateSpace* GetStateSpace() {return statespace;}
	virtual int GetNstate() {return Nstate;}

	protected:

	int Nstate;
	StateSpace* statespace;

	int Ngene;
	int* genealloc;
	int* genesize;
	string* genename;

	int GlobalNsite;
	int* globalnsite;

};

class MultiGeneProfileProcess : public virtual ProfileProcess, public virtual MultiGeneMPIModule	{

};


class MultiGenePhyloProcess : public virtual PhyloProcess, public virtual MultiGeneProfileProcess	{

	public:

	MultiGenePhyloProcess() : genelnL(0) {}
	virtual ~MultiGenePhyloProcess() {}

	virtual void Create();
	virtual void Delete();

	protected:

	virtual SequenceAlignment* GetData() {
		cerr << "error no data in multi gene\n";
		exit(1);
	}

	virtual const SequenceAlignment* GetData() const {
		cerr << "error no data in multi gene\n";
		exit(1);
	}

        virtual void SlaveExecute(MESSAGE);

	void SlaveLikelihood();

	void GlobalSample();
	void SlaveSample();

	void GlobalGeneMove();
	void SlaveGeneMove();

	void SlaveUnfold();
	void SlaveCollapse();
	void SlaveRoot(int n);
	void SlavePropose(int n, double x);
	void SlaveRestore(int n);
	void SlaveReset(int n, bool v);
	void SlaveMultiply(int n, int m, bool v);
	void SlaveSMultiply(int n, bool v);
	void SlaveInitialize(int n, int m, bool v);
	void SlavePropagate(int n, int m, bool v, double t);
	void SlaveGibbsSPRScan(int idown, int iup);

	PhyloProcess** process;
	double* genelnL;
	double* tmpgenelnL;

};

#endif
