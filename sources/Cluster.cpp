
#include "Random.h"
#include <vector>
#include <map>

const char AAset[] = {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','a','c','d','e','f','g','h','i','k','l','m','n','p','q','r','s','t','v','w','y','-','?','$','.','B','Z','*','X','x'};
const int Nstate = 20;

class Stat	{

	vector<double> v;
	int weight;

	public:

	Stat() : v(Nstate,0), weight(0) {}

	Stat(const Stat& from) : v(from.v), weight(from.weight)	{}

	~Stat()	{}

	int GetWeight()	{
		return weight;
	}

	void Add(const Stat& from)	{
		for (int i=0; i<Nstate; i++)	{
			v[i] += from.v[i];
		}
		weight ++;
	}

	void FromStream(istream& is)	{
		int tmp;
		is >> tmp;
		double tot = 0;
		for (int i=0; i<Nstate; i++)	{
			is >> v[i];
			tot += v[i];
		}
		weight = 1;
		if (fabs(tot -1) > 0.01)	{
			cerr << "error in normalisation : " << tot << "\n";
			for (int i=0; i<Nstate; i++)	{
				cerr << v[i] << '\t';
			}
			cerr << '\n';
			exit(1);
		}
			
	}

	void ToStreamSummary(ostream& os,double c, int type)	{
		os << weight << '\t';
		vector<int> tmp = GetSubset(c,type);
		for (int i=0; i<Nstate; i++)	{
			if (tmp[i])	{
				os << AAset[i];
			}
		}
		os << '\t' << GetMaxMinor(tmp);
		os << '\n';
	}
	
	void ToStream(ostream& os)	{
		os << weight;
		for (int i=0; i<Nstate; i++)	{
			os << '\t' << v[i];
		}
		os << '\n';
	}
	
	int GetNaa(double cutoff, int type)	{
		vector<int> tmp = GetSubset(cutoff,type);
		int naa = 0;
		for (int i=0; i<Nstate; i++)	{
			naa += tmp[i];
		}
		return naa;
	}

    double GetEntropy() {
        double tot = 0;
        for (int i=0; i<Nstate; i++)    {
            if (v[i] > 1e-6)    {
                tot -= v[i]*log(v[i]);
            }
        }
        return tot;
    }

	vector<int> GetSubset(double cutoff, int type)	{

		vector<int> a(Nstate,0);

        // all amino-acids with freq > cutoff
		if (type == 0)	{
			for (int i=0; i<Nstate; i++)	{
				a[i] = ((v[i] > cutoff));
			}
			return a;
		}
	
		else if (type == 1)	{

			double* tmp = new double[Nstate];
			for (int i=0; i<Nstate; i++)	{
				tmp[i] = v[i];
			}
			double tot = 0;
			while (tot < cutoff)	{
				double max = 0;
				int imax = 0;
				for (int i=0; i<Nstate; i++)	{
					if (max < tmp[i])	{
						max = tmp[i];
						imax = i;
					}
				}
				if (max == 0)	{
					cerr << "error: total mass is not 1\n";
					exit(1);
				}
				tot += tmp[imax];
				a[imax] = 1;
				tmp[imax] = 0;
			}
			delete[] tmp;
			return a;
		}
        else if (type == 2) {

            int naa = int(exp(GetEntropy()) + 0.5);
            vector<double> tmp = v;
            for (int k=0; k<naa; k++)   {
                double max = 0;
                int imax = 0;
                for (int i=0; i<Nstate; i++)    {
                    if (max < tmp[i])   {
                        max = tmp[i];
                        imax = i;
                    }
                }
                a[imax] = 1;
                tmp[imax] = 0;
            }
            return a;
        }
		// type == 3

		int naa = 1;
		int found = 0;
        double f = 0.33;
        double g = 3;
		while ((g * (1-cutoff) * naa < f * cutoff * (Nstate-naa)) && (! found))	{
			int n = 0;
            int complies = 1;
			for (int i=0; i<Nstate; i++)	{
				if (v[i] > f*cutoff/naa)	{
					n++;
				}
                else    {
                    if (v[i] > g*(1-cutoff)/(Nstate-naa))    {
                        complies = 0;
                    }
                }
			}
			if ((n == naa) && (complies))	{
				found = 1;
				for (int i=0; i<Nstate; i++)	{
					if (v[i] > f*cutoff/naa)	{
						a[i] = 1;
					}
				}
			}
			else	{
				naa++;
			}
		}
        if (! found)  {
            for (int i=0; i<Nstate; i++)    {
                a[i] = 0;
            }
        }
		return a;
	}

	void Normalise()	{
		for (int i=0; i<Nstate; i++)	{
			v[i] /= weight;
		}
	}

	double GetMaxMinor(const vector<int>& group)	{
		
		double max = 0;
		for (int i=0; i<Nstate; i++)	{
			double tmp = (1 - group[i]) * v[i];
			if (max < tmp)	{
				max = tmp;
			}
		}
		return max;
	}

};

int main(int argc, char* argv[])	{

	ifstream is(argv[1]);
	string name = argv[2];
	double c = atof(argv[3]);
	double weight = atof(argv[4]);
	int type = atoi(argv[5]);

	int N;
	int tmp;
	is >> N >> tmp;
	cerr << N << '\n';

	Stat* profile = new Stat[N];

	double* temp = new double[Nstate];
	for (int i=0; i<N; i++)	{
		profile[i].FromStream(is);
	}

	map<vector<int>,Stat> table;
	for (int i=0; i<N; i++)	{
		table[profile[i].GetSubset(c,type)].Add(profile[i]);
	}
	cerr << '\n';

	size_t Ncluster = table.size();
	cerr << "number of clusters: " << Ncluster << '\n';

	ofstream os((name + ".summary").c_str());
	ofstream sos((name + ".clusters").c_str());
	int ncat = 0;
	int residue = 0;
	for (map<vector<int>,Stat>::iterator i=table.begin(); i!=table.end(); i++)	{
		if (i->second.GetWeight() >= weight)	{
			i->second.Normalise();
			ncat++;
		}
		else	{
			residue += i->second.GetWeight();
		}
	}

	sos << ncat << '\t' << Nstate << '\n';
	os << ncat << '\t' << Nstate << '\n';

	for (int naa = 0; naa<Nstate; naa++)	{
		int n = 0;
		double totweight = 0;
		for (map<vector<int>,Stat>::iterator i=table.begin(); i!=table.end(); i++)	{
			if (i->second.GetWeight() >= weight)	{
				if (i->second.GetNaa(c,type) == naa)	{
					n++;
					totweight += i->second.GetWeight();
					i->second.ToStream(sos);
					i->second.ToStreamSummary(os,c,type);
				}
			}
		}
		if (n)	{
			cerr << naa << '\t' << n << '\t' << totweight << '\n';
		}
	}

	/*
	os << '\n';
	for (map<vector<int>,Stat>::iterator i=table.begin(); i!=table.end(); i++)	{
		if (i->second.GetWeight() < weight)	{
			i->second.Normalise();
			i->second.ToStream(os,c,type);
		}
	}
	*/

	cerr << "after threshold: " << ncat << '\n';
	cerr << "residue : " << residue << '\t' << N << '\n';
}

