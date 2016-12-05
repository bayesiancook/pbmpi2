#include "Random.h"
#include "Parallel.h"
MPI_Datatype Propagate_arg;

using namespace std;
#define IOS_APPEND std::ios_base::app

const string seqheader = "seqheader";

const bool numbering = false;

const int SitePerLine = 50;
const double logoThreshold = 0.001;
const double EntropyEpsilon = 1e-10;

const char AminoAcids[] = {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','-'};

const int LinePerPage = 2;
const double MaxHeight = 3;
const double MaxHeightRelative = 6;
const double cmToPixels = 28.3465;
const double SpaceBetweenLines = -130;
const double MaxKL = 20;
const double KLEpsilon = 1e-10;
const double ratesmoothing = 0.1;

double* rate;
double** Stationary;
int mSiteNumber;

int Nstate = 20;

inline double& Frequency(int i, int j)	{
	return Stationary[i][j];
}

void MakeLogo(ostream& os, int style)	{


	double* entropy = new double[mSiteNumber];
	double MaxEnt = log(Nstate);
	double max = 0;
	double min = 0;


	if(style == 1)	{ // heights proportional to stat information content

		min = MaxEnt;

		for (int i=0; i<mSiteNumber; i++)	{
			double h = 0;
			for (int j=0; j<Nstate; j++)	{
				double temp = Frequency(i,j);
				h += (temp<EntropyEpsilon) ? 0 : -temp * log(temp);
			}
			entropy[i] = MaxEnt -  h;

			if (max < entropy[i])	{
				max = entropy[i];
			}
			if (min > entropy[i])	{
				min= entropy[i];
			}

			// cout << entropy[i] << '\n';
		}
	}

	if (style == 2)	{ // heights all equal

		for (int i=0; i<mSiteNumber; i++)	{
			entropy[i] = 1;
			max = 1;
			min = 0;
		}
	}

	if (style == 3)	{ // heights proportional to rates

		min = 1;
		max = 0;
		
		for (int i=0; i<mSiteNumber; i++)	{

			entropy[i] = - log(rate[i]);
			
			if (max < entropy[i])	{
				max = entropy[i];
			}
			if (min > entropy[i])	{
				min = entropy[i];
			}
		}
	}


	int pageNumber = 1;
	int lineNumber = 0;
	int siteNumber = 0;

	os << "%%Page: " << pageNumber << ' ' << pageNumber << '\n';
	os << "startpage\n";
	os << "startline\n";

	double* array = new double[Nstate];
	int* letters = new int[Nstate];

	for (int i=0; i<mSiteNumber; i++)	{
		entropy[i] = MaxHeight * (entropy[i] - min) /( max - min) + 0.3;

		for (int k=0; k<Nstate; k++)	{
			array[k] = Frequency(i,k) *entropy[i] ;
			letters[k] = k;
		}

		for (int j=0; j<Nstate-1; j++)	{
			for (int k=Nstate-1; k>j; k--)	{
				if (array[k] < array[k-1])	{
					double temp = array[k];
					array[k] = array[k-1];
					array[k-1] = temp;

					int intemp = letters[k];
					letters[k] = letters[k-1];
					letters[k-1] = intemp;
				}
			}
		}

		if (numbering)	{
			os << "numbering {(" << i << ") makenumber} if\n";
		}
		os << "gsave\n";

		for (int k=0; k<Nstate; k++)	{
			double height = (array[k] > 0) ? array[k] : -array[k];
			if (height > logoThreshold)	{
				os << height << " (" << AminoAcids[letters[k]] << ") numchar\n";
			}
		}

		os << "grestore\n";
		os << "shift\n";

		siteNumber++;
		if (siteNumber == SitePerLine)	{
			os << "endline\n";
			os << " 0 " << SpaceBetweenLines << " translate\n";

			siteNumber = 0;
			lineNumber ++;

			if (lineNumber == LinePerPage)	{
				lineNumber = 0;
				pageNumber++;
				os << "endpage\n";
				os << "%%Page: " << pageNumber << " " << pageNumber << "\n";
				os << "startpage\n";
			}
			os << "startline\n";
		}
	}

	os << "endline\n";
	os << "endpage\n";

	os << "%%Trailer\n";
	os << "%%Pages: " << pageNumber << "\n";

	delete[] array;
	delete[] letters;
	delete[] entropy;

}


int main(int argc, char* argv[])	{

	if (argc == 1)	{
		cout << '\n';
		cout << "makelogo <sourcefile> <targetfile> [logo-style]\n";
		cout << '\n';
		cout << "\tsource file : \n";
		cout << "\t\tSiteNumber\n";
		cout << "\tname\trate\t\tstat1\tstat2\t...\tstat20\n";
		cout << '\n';
		cout << "\ttarget file : logo in ps format\n";
		cout << '\n';
		exit(1);
	}
	
	ifstream is(argv[1]);
	string appel = "cp " + seqheader + " " + ((string) argv[2]);
	system(appel.c_str());

	ofstream os(argv[2], IOS_APPEND);

	int style = 1;
	// int style = atoi(argv[3]);

	is >> mSiteNumber;

	// allocating stationary frequency array
	Stationary = new double*[mSiteNumber];
	for (int i=0; i<mSiteNumber; i++)	{
		Stationary[i] = new double[Nstate];
	}

	rate = new double[mSiteNumber];

	// stuffing from infile
	for (int i=0; i<mSiteNumber; i++)	{
		/*
		string temp;
		is >> temp;
		*/
		is >> rate[i];
		for (int j=0; j<Nstate; j++)	{
			is >> Frequency(i,j);
		}
	}

	MakeLogo(os, style);

	for (int i=0; i<mSiteNumber; i++)	{
		delete[] Stationary[i];
	}
	delete[] Stationary;
}
