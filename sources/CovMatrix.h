
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#ifndef COVMATRIX_H
#define COVMATRIX_H

#include <iostream>
#include <cmath>
using namespace std;

#include "Random.h"

class CovMatrix {

	protected:

	// data members

	int dim;

	// value : the infinitesimal generator matrix
	double ** value;

	bool orthonormal;
	bool diagflag;

	static int maxn;
	static int ndiag;
	static double maxerror;
	static double toterror;

	public:

	static double GetMeanDiagError()	{
		if (ndiag)	{
			return toterror / ndiag;
		}
		return 0;
	}

	static double GetMaxDiagError()	{
		return maxerror;
	}

	static double GetMaxDiagN()	{
		return maxn;
	}

	CovMatrix() : dim(0) , value(0) , orthonormal(false)  {}

	CovMatrix(int indim) : orthonormal(false) {
		dim = indim;
		CreateAllMatrix();
		diagflag = false;
	}

	CovMatrix(double** Q, int indim) : orthonormal(false) {
		dim = indim;
		CreateAllMatrix();
		for (int i=0; i<GetDim(); i++)	{
			value[i][i] = Q[i][i];
			for (int j=0; j<i; j++)	{
				value[i][j] =  Q[i][j];
				value[j][i] = value[i][j];
			}
		}
		diagflag = false;
	}

	CovMatrix(const CovMatrix& from) : orthonormal(false) {
		dim = from.GetDim();
		CreateAllMatrix();
		diagflag = false;
		*this = from;
	}

	void SetOrthoNormal(bool in)	{
		orthonormal = in;
	}

	void CreateAllMatrix();

	~CovMatrix() {

		for (int i=0; i<GetDim(); i++)	{
			delete[] value[i];
			delete[] invvalue[i];
			delete[] u[i];
			delete[] invu[i];
		}
		delete[] value;
		delete[] invvalue;
		delete[] u;
		delete[] invu;

		delete[] v;
		delete[] logv;
	}

	CovMatrix&	operator=(const CovMatrix& from)	{
		if (!dim)	{
			dim = from.GetDim();
			CreateAllMatrix();
			diagflag = false;
			}
		if (dim != from.dim)	{
			cerr << "error : non matching dimenstion for matrix\n";
			cerr << GetDim() << '\t' << from.GetDim() << '\n';
			exit(1);
		}
		for (int i=0; i<GetDim(); i++){
			value[i][i] =from.value[i][i];
			for  (int j=0; j<i; j++){
				value[i][j] =from.value[i][j];
				value[j][i] =from.value[i][j];
			}
		}
		diagflag = false;
		return *this;
	}

	void setVal(double** val, int inDim) {
		if (!dim)	{
			dim = inDim;
			CreateAllMatrix();
			diagflag = false;
			}
		if (dim != inDim)	{
			cerr << "error : non matching dimenstion for matrix\n";
			cerr << GetDim() << '\t' << inDim << '\n';
			exit(1);
		}
		for (int i=0; i<GetDim(); i++){
			value[i][i] =val[i][i];
			for  (int j=0; j<i; j++){
				value[i][j] =val[i][j];
				value[j][i] =val[i][j];
			}
		}
		diagflag = false;

	}

	CovMatrix&	operator+=(const CovMatrix from)	{
				for (int i=0; i<GetDim(); i++){
					value[i][i] += from.value[i][i];
					for  (int j=0; j<i; j++){
						value[i][j] += from.value[i][j];
						value[j][i] += from.value[i][j];
					}
				}
				return *this;
				diagflag = false;
			}

	double*		operator[](int i) {
				double* v = value[i];
				return v;
			}

	double	operator()(int i, int j) const {
		return value[i][j];
	}

	CovMatrix& operator*=(double scal) {
		for (int i=0; i<GetDim(); i++){
			value[i][i] *= scal;
			for  (int j=0; j<i; j++){
				value[i][j] *= scal;
				value[j][i] *= scal;
			}
		}
		diagflag = false;
		return *this;
	}

	CovMatrix operator*(double scal) {

		CovMatrix cm = *this;
		cm*=scal;
		return cm;
	}

	virtual double	ProposeMove(double tuning);

	virtual void drawVal(double* vec);

	void drawValInv(double* vec);

	double logValProb(double* dval);

	private:

	double ** invvalue;

	// v : eigenvalues
	// u : the matrix of eigen vectors
	// invu : the inverse of u

	double ** u;
	double ** invu;
	double * v;
	double * logv;

	public:


	double* GetEigenVal() {
		if (! diagflag) Diagonalise();
		return v;
	}

	double* GetLogEigenVal() {
		if (! diagflag) Diagonalise();
		return logv;
	}

	double** GetEigenVect() {
		if (! diagflag)	{
			Diagonalise();
		}

		return u;
	}

	double GetDeterminant() {
		double ret = 1;
		for (int i=0; i<GetDim(); i++) {
			ret *= GetEigenVal()[i];
		}
		return ret;
	}

	double GetLogDeterminant() {
		double ret = 0;
		for (int i=0; i<GetDim(); i++) {
			ret += GetLogEigenVal()[i];
		}
		/*
		if (isnan(ret))	{
			cerr << "covmatrix det: nan\n";
			exit(1);
		}
		*/
		return ret;
	}

	double** GetInvEigenVect() {
		if (! diagflag) Diagonalise();
		return invu;
	}

	int GetDim()const {
		return dim;
	}

	void Reset()	{
		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				value[i][j] = 0;
			}
		}
	}

	void SetIdentity()	{
		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				value[i][j] = 0;
			}
		}
		for (int i=0; i<GetDim(); i++)	{
			value[i][i] = 1;
		}
	}

	void Project(int index, double** m);

	int ScalarMultiplication(double scal)	{
		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				value[i][j] *= scal;
			}
		}
		diagflag = false;
		return GetDim();
	}

	double** GetMatrix() {
		return value;
	}

	double** GetMatrix() const{
		return value;
	}

	double** GetInvMatrix() {
		if (! diagflag) Diagonalise();
		return invvalue;
	}


	bool isPositive(){
		if (! diagflag) Diagonalise();
		bool r = true;
		for (int i=0; i<GetDim(); i++){
			if(GetEigenVal()[i] <= 1e-6){
				r = false;
			}
		}
		return r;
	}

	void CorruptDiag()	{
		diagflag = false;
	}

	double GetMax()	{
		double max = 0;
		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				double tmp = fabs(value[i][j]);
				if (max < tmp)	{
					max = tmp;
				}
			}
		}
		return max;
	}

	//Set the matrix to it s inverse //loook si diagflag
	int Invert();

	double CheckInverse();

	int Diagonalise();

	double CheckDiag();

	void Scatter(double** inValues, int Nval);

	void PrintCorrelationCoefficients(ostream& os);

	void PrintEigenVectors(ostream& os);

	friend ostream& operator<<(ostream& os, const CovMatrix& r)  {
		for (int i=0; i<r.GetDim(); i++)	{
			for (int j=0; j<r.GetDim(); j++)	{
				os << '\t' << r.GetMatrix()[i][j];
			}
			os << '\n';
		}
		return os;
	}

	friend istream& operator>>(istream& is, CovMatrix& r)  {
		for (int i=0; i<r.GetDim(); i++)	{
			for (int j=0; j<r.GetDim(); j++)	{
				is >> r.value[i][j];
			}
		}
		return is;
	}

};

class InverseWishartMatrix : public virtual CovMatrix	{

	protected :

	double* kappa;
	int df;

	int shapestat;
	CovMatrix scalestat;
	int bkshapestat;
	CovMatrix bkscalestat;

	public :

	int GetDF() {
		return df;
	}

	InverseWishartMatrix(double* inkappa, int indim, int indf) : CovMatrix(indim), scalestat(indim), bkscalestat(indim)	{
		kappa = inkappa;
		df = indf;
		Sample();
	}

	void Sample(){
		int cont = 1;
		double** sample = new double*[df];
		for (int i=0; i<df ; i++) {
			sample[i] = new double[GetDim()];
		}

		while (cont)	{
			for (int i=0; i<df; i++) {
				for (int j=0; j<GetDim(); j++)	{
					sample[i][j] = rnd::GetRandom().sNormal() / kappa[j];
				}
			}
			Scatter(sample,df);
			cont = Invert();
		}


		for (int i=0; i<df ; i++) {
			delete[] sample[i];
		}
		delete[] sample;
	}

	void drawSample(CovMatrix* A);

	void SetDiagonal(){
		for(int i=0; i<GetDim(); i++){
			for(int j=0; j<GetDim(); j++){
				value[i][j] = 0;
			}
		}
		for(int i=0; i<GetDim(); i++){
			value[i][i] = kappa[i];
		}
	}

	double GetLogProb(){
		if(isPositive()){
			double sum = 0;
			for(int i=0; i<GetDim(); i++){
				sum += GetInvMatrix()[i][i] * kappa[i];
			}
			double d = - ((GetLogDeterminant() * (GetDim() + df + 1)) + sum) * 0.5;
			double tmp = 0;
			for(int i=0; i<GetDim(); i++){
				tmp += log(kappa[i]);
			}
			d += tmp * df * 0.5;
			return d;
		}
		else{
			cerr << "singular cov matrix\n";
			exit(1);
			return log(0);
		}
	}

	void ResetSufficientStatistic()	{
		shapestat = 0;
		for( int i=0; i<dim; i++){
			for( int j=0; j<dim; j++){
				scalestat[i][j] = 0;
			}
		}
	}


	void SaveSufficientStatistic()	{
		bkshapestat = shapestat;
		for( int i=0; i<dim; i++){
			bkscalestat[i][i] = scalestat[i][i];
			for( int j=0; j<i; j++){
				bkscalestat[i][j] = scalestat[i][j];
			}
		}
	}

	void RestoreSufficientStatistic()	{
		bkshapestat = shapestat;
		for( int i=0; i<dim; i++){
			scalestat[i][i] = bkscalestat[i][i];
			for( int j=0; j<i; j++){
				scalestat[i][j] = bkscalestat[i][j];
				scalestat[j][i] = bkscalestat[i][j];
			}
		}
	}

	void AddToShape()	{
		shapestat ++;
	}

	void AddToScale(const double* in)	{
		for( int i=0; i<dim; i++){
			scalestat[i][i] += in[i] * in[i];
			for( int j=0; j<i; j++){
				scalestat[i][j] += in[i] * in[j];
				scalestat[j][i] = scalestat[i][j];
			}
		}
	}

	void AddToShape(int N)	{
		shapestat +=N;
	}

	void AddToScale(double** in)	{
		for( int i=0; i<dim; i++){
			for(int j=0; j<dim; j++) {
				 scalestat[i][j] += in[i][j];
			}
		}
	}

	void ComputeSufficientStatistics(double** in, int N)	{

		ResetSufficientStatistic();
		for (int i=0; i<N; i++)	{
			AddToScale(in[i]);
		}
		AddToShape(N);
	}

	double SuffStatLogProb()	{
		double t, trace = 0;
		for( int i=0; i<dim; i++){
			t = 0;
			for( int j=0; j<dim; j++){
				t += scalestat[i][j] * GetInvMatrix()[j][i];
			}
			trace += t;
		}
		return -0.5 * (shapestat * GetLogDeterminant() + trace);
	}

	void GibbsResample(double** in, int N)	{

		ComputeSufficientStatistics(in,N);
		for(int i=0; i<GetDim(); i++){
			scalestat[i][i] += kappa[i];
		}
		df += shapestat;
		scalestat.Diagonalise();
		drawSample(&scalestat);
		for(int i=0; i<GetDim(); i++){
			scalestat[i][i] -= kappa[i];
		}
		df -= shapestat;
		if (!isPositive())	{
			cout << " conjugate Gibbs draw of inverse wishart is singular\n";
			for(int i=0; i<GetDim(); i++){
				cerr << kappa[i] << '\t';
			}
			cerr << '\n';
			cerr << df << '\n';
			for(int i=0; i<GetDim(); i++){
				for(int j=0; j<GetDim(); j++){
					cerr << scalestat[i][j] << '\t';
				}
				cerr << '\n';
			}
			exit(1);
		}

	}
};

#endif
