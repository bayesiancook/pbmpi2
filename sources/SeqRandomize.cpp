#include "Random.h"
#include "SequenceAlignment.h"

int main(int argc, char* argv[])	{

	string filename = argv[1];
	string output = argv[2];
	int nrep = atoi(argv[3]);

	SequenceAlignment* from = new FileSequenceAlignment(filename);

	for (int rep=0; rep<nrep; rep++)	{
		cerr << '.';
		from->Randomize();
		ostringstream s;
		s << "rnd" << rep << output << ".ali";
		ofstream os(s.str().c_str());
		from->ToStream(os);	
	}
	cerr << '\n';
}

