
#include "CovMatrix.h"
#include "linalg.h"

int CovMatrix::maxn = 0;
int CovMatrix::ndiag = 0;
double CovMatrix::maxerror = 0;
double CovMatrix::toterror = 0;


void CovMatrix::CreateAllMatrix(){
	value = new double*[GetDim()];
	for (int i=0; i<GetDim(); i++)	{
		value[i] = new double[GetDim()];
		for (int j=0; j<GetDim(); j++)	{
			value[i][j] = 0;
		}
	}

	invvalue = new double*[GetDim()];
	for (int i=0; i<GetDim(); i++)	{
		invvalue[i] = new double[GetDim()];
	}

	u = new double*[GetDim()];
	for (int i=0; i<GetDim(); i++)	{
		u[i] = new double[GetDim()];
	}

	invu = new double*[GetDim()];
	for (int i=0; i<GetDim(); i++)	{
		invu[i] = new double[GetDim()];
	}

	v = new double[GetDim()];
	logv = new double[GetDim()];
}

double CovMatrix::ProposeMove(double tuning)	{
	int j = (int) (GetDim() * (GetDim()+1) * 0.5 * rnd::GetRandom().Uniform());
	int i=0;
	while(i+1<=j){
		i++;
		j -= i;
	}

	if (i == j)	{
		double infbound = 0;
		for (int k=0; k<GetDim(); k++)	{
			if (k != i)	{
				double tmp = value[i][k] * value[i][k] / value[k][k];
				if (infbound < tmp)	{
					infbound = tmp;
				}
			}
		}
		double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
		value[i][i] += m;
		if (value[i][i] < infbound)	{
			value[i][i] = 2 * infbound - value[i][i];
		}
	}
	else	{
		double bound = sqrt(value[i][i] * value[j][j]);
		double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
		value[i][j] += m;
		while ((value[i][j] < -bound) || (value[i][j] > bound))	{
			if (value[i][j] < -bound)	{
				value[i][j] = -2 * bound - value[i][j];
			}
			if (value[i][j] > bound)	{
				value[i][j] = 2 * bound - value[i][j];
			}
		}
		value[j][i] = value[i][j];
	}
	return 0;
}

void CovMatrix::drawVal(double* vec){
	double* randomvalues = new double[GetDim()];
	for (int i=0; i<GetDim(); i++) {
		randomvalues[i] =  rnd::GetRandom().sNormal() * sqrt(GetEigenVal()[i]);
		vec[i]=0;
	}
	for (int i=0; i<GetDim(); i++) {
		for (int j=0; j<GetDim(); j++) {
			vec[j] +=  randomvalues[i] * GetEigenVect()[j][i];
		}
	}
	delete[] randomvalues;
}

void CovMatrix::drawValInv(double* vec){
	double* randomvalues = new double[GetDim()];
	for (int i=0; i<GetDim(); i++) {
		randomvalues[i] =  rnd::GetRandom().sNormal() / sqrt(GetEigenVal()[i]);
		vec[i]=0;
	}
	for (int i=0; i<GetDim(); i++) {
		for (int j=0; j<GetDim(); j++) {
			vec[j] +=  randomvalues[i] * GetEigenVect()[j][i];
		}
	}
	delete[] randomvalues;
}

double CovMatrix::logValProb(double* dval)	{

	double tXSX = 0;
	for (int i=0; i<GetDim() ; i++) {
		tXSX += GetInvMatrix()[i][i] * dval[i] * dval[i];
		for (int j=0; j<i ; j++) {
			tXSX += 2 * GetInvMatrix()[i][j] * dval[j] * dval[i];
		}
	}
	return -0.5 * (GetLogDeterminant() + tXSX);
}

void CovMatrix::Project(int index, double** m)	{

	int k = 0;
	for (int i=0; i<GetDim(); i++)	{
		if (i != index)	{
			int l = 0;
			for (int j=0; j<GetDim(); j++)	{
				if (j != index)	{
					m[k][l] = value[i][j] - value[i][index] * value[j][index] / value[index][index];
					l++;
				}
			}
			k++;
		}
	}
}

int CovMatrix::Invert(){
	double** a = new double*[GetDim()];
	for (int i=0; i<GetDim(); i++)	{
		a[i] = new double[GetDim()];
	}

	// copy value into a :
	for (int i=0; i<GetDim(); i++)	{
		for (int j=0; j<GetDim(); j++)	{
			a[i][j] = value[i][j];
		}
	}

	// invert a into value
	// InvertMatrix(a, GetDim(), w, iw, value);
	double logdet = LinAlg::Gauss(a,GetDim(),value);

	// cerr << "check inverse : " << CheckInverse() << '\n';
	for (int i=0; i<GetDim(); i++)	{
		delete[] a[i];
	}
	delete[] a;
	diagflag = false;
	if (isinf(logdet))	{
		cerr << "error in cov matrix: non invertible\n";
		return 1;
		exit(1);
	}
	return 0;
}

double CovMatrix::CheckInverse()	{
	double max = 0;
	for (int i=0; i<GetDim(); i++)	{
		for (int j=0; j<GetDim(); j++)	{
			double tot = 0;
			for (int k=0; k<GetDim(); k++)	{
				tot += value[i][k] * GetInvMatrix()[k][j];
			}
			if (i == j)	{
				tot --;
			}
			if (max < fabs(tot))	{
				max = fabs(tot);
			}
		}
	}
	return max;
}


int CovMatrix::Diagonalise()	{
	int nmax = 1000;
	double epsilon = 1e-10;

	int n = LinAlg::DiagonalizeSymmetricMatrix(value,dim,nmax,epsilon,v,u);
	if (maxn < n)	{
		maxn = n;
	}
	bool failed = (n == nmax);
	if (failed)	{
		cerr << "diag failed\n";
		cerr << n << '\n';
		for (int i=0; i<dim; i++)	{
			for (int j=0; j<dim; j++)	{
				cerr << value[i][j] << '\t';
			}
			cerr << '\n';
		}
		exit(1);
	}

	// normalise u
	for (int i=0; i<dim; i++)	{
		double total = 0;
		for (int j=0; j<dim; j++)	{
			total += u[j][i] * u[j][i];
		}
	}
	// u-1 = tu
	for (int i=0; i<dim; i++)	{
		for (int j=0; j<dim; j++)	{
			invu[j][i] = u[i][j];
		}
	}

	for (int i=0; i<GetDim(); i++)	{
		logv[i] = log(v[i]);
	}

	LinAlg::Gauss(value,GetDim(),invvalue);

	diagflag = true;
	double tmp = CheckDiag();
	if (maxerror < tmp)	{
		maxerror = tmp;
	}
	toterror += tmp;
	ndiag ++;
	return failed;
}

double CovMatrix::CheckDiag()	{
	double** a = new double*[GetDim()];
	double** b = new double*[GetDim()];
	for (int i=0; i<GetDim(); i++)	{
		a[i] = new double[GetDim()];
		b[i] = new double[GetDim()];
	}

	for (int i=0; i<GetDim(); i++)	{
		for (int j=0; j<GetDim(); j++)	{
			double tot = 0;
			for (int k=0; k<GetDim(); k++)	{
				tot += invu[i][k] * GetMatrix()[k][j];
			}
			a[i][j] = tot;
		}
	}

	double max = 0;

	for (int i=0; i<GetDim(); i++)	{
		for (int j=0; j<GetDim(); j++)	{
			double tot = 0;
			for (int k=0; k<GetDim(); k++)	{
				tot += a[i][k] * u[k][j];
			}
			b[i][j] = tot;
			if (i != j)	{
				if (max < fabs(tot))	{
					max = fabs(tot);
				}
			}
		}
	}

	for (int i=0; i<GetDim(); i++)	{
		delete[] a[i];
		delete[] b[i];
	}
	delete[] a;
	delete[] b;

	return max;
}

void CovMatrix::Scatter(double** inValues, int Nval){
	for (int i=0; i<GetDim(); i++) {
		for (int j=0; j<GetDim(); j++) {
			value[i][j] = 0;
			for (int l=0; l<Nval; l++) {
				value[i][j] += inValues[l][i] * inValues[l][j];
			}
		}
	}
	diagflag = false;
}

void CovMatrix::PrintCorrelationCoefficients(ostream& os)	{
	for (int i=0; i<GetDim(); i++)	{
		for (int j=0; j<GetDim(); j++)	{
			if (i==j) 	{
				os << "\t-";
			}
			else	{
				os << '\t' << GetMatrix()[i][j] / sqrt(GetMatrix()[i][i] * GetMatrix()[j][j]);
			}
		}
		os << '\n';
	}
}

void CovMatrix::PrintEigenVectors(ostream& os)	{
	os << "val";
	for (int j=0; j<GetDim(); j++)	{
		os << '\t' << j;
	}
	os << '\n';
	for (int i=0; i<GetDim(); i++)	{
		os << v[i] << '\t';
		for (int j=0; j<GetDim(); j++)	{
			os << '\t' << u[i][j];
		}
		os << '\n';
	}
	os << '\n';
	os << "inverse eigenvector matrix\n";
	for (int i=0; i<GetDim(); i++)	{
		for (int j=0; j<GetDim(); j++)	{
			os << '\t' << invu[i][j];
		}
		os << '\n';
	}
	os << '\n';

	// proportion of variance explained
	double** prop = new double*[GetDim()];
	for (int i=0; i<GetDim(); i++)	{
		prop[i] = new double[GetDim()];
	}
	for (int i=0; i<GetDim(); i++)	{
		double total = 0;
		for (int j=0; j<GetDim(); j++)	{
			prop[i][j] = invu[j][i] * invu[j][i] * v[j];
			// prop[i][j] = u[i][j] * u[i][j] * v[j];
			// prop[j][i] = invu[j][i] * invu[j][i] * v[j];
			total += prop[i][j];
		}
		for (int j=0; j<GetDim(); j++)	{
			prop[i][j] /= total;
		}
	}
	os << "proportion of variance explained\n";
	for (int i=0; i<GetDim(); i++)	{
		os << i;
		for (int j=0; j<GetDim(); j++)	{
			os << '\t' << prop[i][j];
		}
		os << '\n';
	}
	os << '\n';
	for (int i=0; i<GetDim(); i++)	{
		delete[] prop[i];
	}
	delete[] prop;

}


void InverseWishartMatrix::drawSample(CovMatrix* A)	{

	// algorithm of Odell and Feiveson, 1966
	double* v = new double[GetDim()];
	for (int i=0; i<GetDim(); i++)	{
		v[i] = rnd::GetRandom().Gamma(0.5*(df-i),0.5);
	}
	double** n = new double*[GetDim()];
	double** b = new double*[GetDim()];
	double** a = new double*[GetDim()];
	for (int i=0; i<GetDim(); i++)	{
		n[i] = new double[GetDim()];
		b[i] = new double[GetDim()];
		a[i] = new double[GetDim()];
	}

	for (int i=0; i<GetDim(); i++)	{
		for (int j=i+1; j<GetDim(); j++)	{
			n[i][j] = rnd::GetRandom().sNormal();
		}
	}
	for (int i=0; i<GetDim(); i++)	{
		b[i][i] = v[i];
		for (int k=0; k<i; k++)	{
			b[i][i] += n[k][i] * n[k][i];
		}
	}
	for (int i=0; i<GetDim(); i++)	{
		for (int j=i+1; j<GetDim(); j++)	{
			b[i][j] = n[i][j] * sqrt(v[i]);
			for (int k=0; k<i; k++)	{
				b[i][j] += n[k][i] * n[k][j];
			}
			b[j][i] = b[i][j];
		}
	}

	A->CorruptDiag();
	A->Diagonalise();

	double** p = A->GetEigenVect();
	double** invp = A->GetInvEigenVect();
	double* d = A->GetEigenVal();
	
	
	for (int i=0; i<GetDim(); i++)	{
		for(int j=0; j<GetDim(); j++)	{
			a[i][j] = p[i][j] / sqrt(d[j]);
		}
	}

	for (int i=0; i<GetDim(); i++)	{
		for(int j=0; j<GetDim(); j++)	{
			double tmp = 0;
			for (int k=0; k<GetDim(); k++)	{
				tmp += b[i][k] * a[j][k];
			}
			n[i][j] = tmp;
		}
	}

	for (int i=0; i<GetDim(); i++)	{
		for(int j=0; j<GetDim(); j++)	{
			double tmp = 0;
			for (int k=0; k<GetDim(); k++)	{
				tmp += a[i][k] * n[k][j];
			}
			(*this)[i][j] = tmp;
		}
	}
		
	Invert();
	
	for (int i=0; i<GetDim(); i++)	{
		delete[] n[i];
		delete[] b[i];
		delete[] a[i];
	}
	delete[] n;
	delete[] b;
	delete[] a;
	delete[] v;
}
