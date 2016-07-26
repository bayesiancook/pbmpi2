
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <sstream>

using namespace std;

int main(int argc, char* argv[])	{

	string name = "";
	ostringstream ss;
	ostringstream ss2;

	ss << "qsub";

	system("pwd > tmpsub");
	ifstream is("tmpsub");
	string wd;
	is >> wd;

	int mock = 0;

	try	{

		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];

			if (s == "-q")	{
				ss << " -q ";
				i++;
				s = argv[i];
				ss << s;
			}
			else if (s == "-l")	{
				ss << " -l ";
				i++;
				s = argv[i];
				ss << s;
			}
			else if (s == "-n")	{
				i++;
				name = argv[i];
				ss << " -o " << wd << "/" << name << ".out";
				ss << " -e " << wd << "/" << name << ".err";
				ss << " " << wd << "/" << name << ".sub";
			}
			else if (s == "-m")	{
				mock = 1;
			}
			else	{
				while (i < argc)	{
					string s = argv[i];
					ss2 << s << " ";
					i++;
				}
			}
			i++;
		}
	
	}
	catch(...)	{
		cerr << "error in command\n";
		exit(1);
	}

	/*
	cerr << "command\n";
	cerr << ss.str() << '\n';
	cerr << ss2.str() << '\n';
	*/

	ofstream os((name + ".sub").c_str());
	os << "cd " << wd << '\n';
	os << ss2.str() << '\n';
	os.close();

	if (mock)	{
		cerr << '\n';
		cerr << "working dir: " << wd << '\n';
		cerr << '\n';
		cerr << ss.str() << '\n';
		cerr << '\n';
	}
	else	{
		system(ss.str().c_str());
	}
}

