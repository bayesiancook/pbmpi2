
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

	// anything here ?
	// make default versions of stupid functions
	virtual double* GetProfile(int site)	{
		cerr << "error: in multi gene getprofile\n";
		exit(1);
		return 0;
	}

	virtual void SampleProfile()	{
		cerr << "error: in multi gene sampleprofile\n";
		exit(1);
	}

	virtual double ProfileSuffStatLogProb()	{
		cerr << "error: in multi gene profilesuffstatlogprob\n";
		exit(1);
		return 0;
	}

	virtual double LogProfilePrior() {
		cerr << "error: in multi gene\n";
		exit(1);
		return 0;
	}
	virtual double LogHyperPrior()	{
		cerr << "error: in multi gene\n";
		exit(1);
		return 0;
	}
	virtual double LogStatPrior()	{
		cerr << "error: in multi gene\n";
		exit(1);
		return 0;
	}
	virtual double LogFrequencyStatPrior(double* prof)	{
		cerr << "error: in multi gene\n";
		exit(1);
		return 0;
	}
	virtual double GetMinStat(double* profile, int site)	{
		cerr << "error: in multi gene\n";
		exit(1);
		return 0;
	}

	virtual void SampleHyper()	{
		cerr << "error: in multi gene\n";
		exit(1);
	}
	virtual void SampleStat()	{
		cerr << "error: in multi gene\n";
		exit(1);
	}
	virtual void SampleFrequencyStat(double* prof)	{
		cerr << "error: in multi gene\n";
		exit(1);
	}
	virtual void UpdateSiteProfileSuffStat()	{
		cerr << "error: in multi gene\n";
		exit(1);
	}
	virtual void GlobalUpdateSiteProfileSuffStat()	{
		cerr << "error: in multi gene\n";
		exit(1);
	}
	virtual void SlaveUpdateSiteProfileSuffStat()	{
		cerr << "error: in multi gene\n";
		exit(1);
	}

};

class MultiGenePhyloProcess : public virtual PhyloProcess, public virtual MultiGeneProfileProcess	{

	public:

	MultiGenePhyloProcess() : genelnL(0) {}
	virtual ~MultiGenePhyloProcess() {}

	virtual void Create();
	virtual void Delete();

	double GetLogLikelihood();

	void GlobalToStream(ostream& os);
	void GlobalFromStream(istream& is);
	void SlaveToStream();
	void SlaveFromStream();

	void WaitLoop();

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

	void GlobalSample();
	void GlobalUnfold();
	void GlobalCollapse();
	void GlobalGeneMove();

	void SlaveSample();
	void SlaveUnfold();
	void SlaveCollapse();
	void SlaveGeneMove();

	void GlobalCollectGeneLikelihoods();
	void SlaveSendGeneLikelihoods();

	void GlobalCollectGeneParameters();
	void SlaveSendGeneParameters();

	PhyloProcess** process;
	double* genelnL;
	double* tmpgenelnL;

};

#endif
