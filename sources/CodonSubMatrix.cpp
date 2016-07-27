
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#include "CodonSubMatrix.h"

void CodonSubMatrix::ComputeArray(int i)	{

	double total = 0;
	
	for (int j=0; j<GetNstate(); j++)       {
		if (i!=j)       {
			int pos = GetDifferingPosition(i,j);
			if ((pos != -1) && (pos != 3))  {
				int a = GetCodonPosition(pos,i);
				int b = GetCodonPosition(pos,j);
				Q[i][j] = nucrr[GetNucRRIndex(a,b)] * nucstat[b];
				total += Q[i][j];
			}
			else    {
				Q[i][j] = 0;
			}
		}
	}
	Q[i][i] = -total;
}

void CodonSubMatrix::ComputeStationary()	{

	double total = 0;
	for (int i=0; i<GetNstate(); i++)	{
		mStationary[i] = nucstat[GetCodonPosition(0,i)] * nucstat[GetCodonPosition(1,i)] * nucstat[GetCodonPosition(2,i)];
		total += mStationary[i];
	}

	// re-normalize

	for (int i=0; i<GetNstate(); i++)	{
		mStationary[i] /= total;
	}
}


double CodonSubMatrix::RateAwaySyn(int i)	{
	double total = 0;
	int j = 0;		
	while ( (GetCodonNearestNeighbors(i,j) != -1) && (j < (Nnuc-1)*GetCodonStateSpace()->Npos) )	{
		if (Synonymous(i,GetCodonNearestNeighbors(i,j)))	{
			total += Q[i][GetCodonNearestNeighbors(i,j)];
		}
		j++;
	}
	return total;
}

double CodonSubMatrix::RateAwayNonsyn(int i)	{
	double total = 0;
	int j = 0;		
	while ( (GetCodonNearestNeighbors(i,j) != -1) && (j < (Nnuc-1)*GetCodonStateSpace()->Npos) )	{
		if (!Synonymous(i,GetCodonNearestNeighbors(i,j)))	{
			total += Q[i][GetCodonNearestNeighbors(i,j)];
		}
		j++;
	}
	if (total < 0)	{
		cerr << "error: negative rate away non syn\n";
		exit(1);
	}
	return total;
}

double CodonSubMatrix::NonSynRate()	{

	double total = 0;
	for (int i=0; i<GetNstate(); i++)	{
		total += mStationary[i] * RateAwayNonsyn(i);
	}
	return total;
}

void AACodonMutSelProfileSubMatrix::ComputeArray(int i)	{

	double total = 0;
	double deltaF;
	double mutrate = 0;
	for (int j=0; j<GetNstate(); j++)       {
		if (i!=j)       {
			int pos = GetDifferingPosition(i,j);
			if ((pos != -1) && (pos != 3))  {
				int a = GetCodonPosition(pos,i);
				int b = GetCodonPosition(pos,j);
				if (a == b)     {
					cerr << "identical states\n";
					cerr << GetCodonStateSpace()->GetState(i) << '\t' << GetCodonStateSpace()->GetState(j) << '\n';
					cerr << pos << '\n';
					exit(1);
				}
				Q[i][j] = nucrr[GetNucRRIndex(a,b)] * nucstat[b];
				mutrate = Q[i][j];
				if (mutrate < 0)	{
					cerr << "error: negative mutation rate \n";
					cerr << a << '\t' << b << '\t' << GetNucRRIndex(a,b) << '\n';
					cerr << nucrr[GetNucRRIndex(a,b)] << '\t' << nucstat[b] << '\n';
				}
				if (! Synonymous(i,j))  {
					deltaF = log((aaprofile)[GetCodonStateSpace()->Translation(j)] / (aaprofile)[GetCodonStateSpace()->Translation(i)]) +
							log( (codonprofile)[j] / (codonprofile)[i] );
					Q[i][j] *= *omega;

					//cerr << "in ComputeArray, omega is " << *omega << "\n";	
					//cerr << "Q[" << i << "][" << j << "]: " << Q[i][j] << "\n";
					//cerr.flush();
				}
				else	{
					deltaF = log( (codonprofile)[j] / (codonprofile)[i] );

				}

				if (fabs(deltaF) < TOOSMALL)        {
					Q[i][j] /= ( 1.0 - (deltaF / 2) );
				}
				else if (deltaF > TOOLARGE)	{
					Q[i][j] *= deltaF;
				}
				else if (deltaF < TOOLARGENEGATIVE)	{
					Q[i][j] = 0.0;
					//Q[i][j] = 1e-10;
				}
				else    {
					Q[i][j] *=  (deltaF)/(1.0 - exp(-deltaF));
				}	
			}
			else    {
				Q[i][j] = 0;
			}
			total += Q[i][j];



			if (Q[i][j] < 0)        {
				cerr << "negative entry in matrix\n";
				cerr << Q[i][j] << '\n';
				cerr << "mut rate : " << mutrate << '\n';
				cerr << "deltaF: " << deltaF << "\n";
				cerr << "omega : " << *omega << '\n';
				cerr << "codonprofile[" << i << "]: " << codonprofile[i] << "\n";
				cerr << "codonprofile[" << j << "]: " << codonprofile[j] << "\n";
				cerr << "aaprofile[" << GetCodonStateSpace()->Translation(j) << "]: " << (aaprofile)[GetCodonStateSpace()->Translation(j)] << "\n";
				cerr << "aaprofile[" << GetCodonStateSpace()->Translation(i) << "]: " << (aaprofile)[GetCodonStateSpace()->Translation(i)] << "\n";
				exit(1);
			}
			if (isinf(Q[i][j]))	{
				cerr << "inf Q[i][j]\n";
				cerr << "deltaF: " << deltaF << "\n";
				cerr << "codonprofile[" << i << "]: " << codonprofile[i] << "\n";
				cerr << "codonprofile[" << j << "]: " << codonprofile[j] << "\n";
				cerr << "aaprofile[" << GetCodonStateSpace()->Translation(j) << "]: " << (aaprofile)[GetCodonStateSpace()->Translation(j)] << "\n";
				cerr << "aaprofile[" << GetCodonStateSpace()->Translation(i) << "]: " << (aaprofile)[GetCodonStateSpace()->Translation(i)] << "\n";
				exit(1);
			}
			if (isnan(Q[i][j]))	{
				cerr << "nan Q[i][j]\n";
				cerr << "deltaF: " << deltaF << "\n";
				cerr << "codonprofile[" << i << "]: " << codonprofile[i] << "\n";
				cerr << "codonprofile[" << j << "]: " << codonprofile[j] << "\n";
				cerr << "aaprofile[" << GetCodonStateSpace()->Translation(j) << "]: " << (aaprofile)[GetCodonStateSpace()->Translation(j)] << "\n";
				cerr << "aaprofile[" << GetCodonStateSpace()->Translation(i) << "]: " << (aaprofile)[GetCodonStateSpace()->Translation(i)] << "\n";
				exit(1);
			}
			
		}
	}
	Q[i][i] = -total;
	if (total <0)   {
		cerr << "negative rate away\n";
		exit(1);
	}
}

void AACodonMutSelProfileSubMatrix::ComputeStationary()	{

	double total = 0;
	for (int i=0; i<GetNstate(); i++)	{
		mStationary[i] =nucstat[GetCodonPosition(0,i)] * 
				nucstat[GetCodonPosition(1,i)] * 
				nucstat[GetCodonPosition(2,i)] *
				codonprofile[i] *
				aaprofile[GetCodonStateSpace()->Translation(i)];
		//if (mStationary[i] < TOOSMALL)	{
		//	mStationary[i] = 0;
		//}
		total += mStationary[i];
	}

	// re-normalize
	
	//double min=1;
	for (int i=0; i<GetNstate(); i++)	{
		mStationary[i] /= total;
		//if (mStationary[i] < min) min = mStationary[i];
	}
	//cout << "smallest stat: " << min << "\n";
	//cout.flush();	
}

double AACodonMutSelProfileSubMatrix::GetRate()	{
	
	if (! ArrayUpdated())	{
		UpdateStationary();
		for (int k=0; k<Nstate; k++)	{
			ComputeArray(k);
		}
	}
	for (int k=0; k<Nstate; k++)	{
		flagarray[k] = true;
	}
	double norm = 0;
	for (int i=0; i<Nnuc-1; i++)	{
		for (int j=i+1; j<Nnuc; j++)	{
			norm += nucstat[i] * nucstat[j] * nucrr[GetNucRRIndex(i,j)];
		}
	}
	return 2 * (norm * 3);
}

void HBAACodonMutSelProfileSubMatrix::ComputeArray(int i)	{

	double total = 0;
	for (int j=0; j<GetNstate(); j++)       {
		if (i!=j)       {
			int pos = GetDifferingPosition(i,j);
			if ((pos != -1) && (pos != 3))  {
				int a = GetCodonPosition(pos,i);
				int b = GetCodonPosition(pos,j);
				if (a == b)     {
					cerr << "identical states\n";
					cerr << GetCodonStateSpace()->GetState(i) << '\t' << GetCodonStateSpace()->GetState(j) << '\n';
					cerr << pos << '\n';
					exit(1);
				}
				Q[i][j] = nucrr[GetNucRRIndex(a,b)] * nucstat[b];
				double f = nucstat[a] / nucstat[b];
				if (! Synonymous(i,j))  {
					f *= exp((*Ne) * log(aaprofile[GetCodonStateSpace()->Translation(j)] / aaprofile[GetCodonStateSpace()->Translation(i)]));
				}

				if (fabs(f-1) > 1e-8)	{
					Q[i][j] *= log(f) / (1 - 1.0/f);
				}

				if (! Synonymous(i,j))  {
					Q[i][j] *= *omega;
				}
			}
			else    {
				Q[i][j] = 0;
			}
			total += Q[i][j];
		}
	}
	Q[i][i] = -total;
}

void HBAACodonMutSelProfileSubMatrix::ComputeStationary()	{

	for (int i=0; i<GetNstate(); i++)	{
		mStationary[i] = nucstat[GetCodonStateSpace()->GetCodonPosition(0,i)] * nucstat[GetCodonStateSpace()->GetCodonPosition(1,i)] * nucstat[GetCodonStateSpace()->GetCodonPosition(2,i)];
	}
	double aanorm[Naa];
	for (int k=0; k<Naa; k++)	{
		aanorm[k] = 0;
	}
	for (int i=0; i<GetNstate(); i++)	{
		aanorm[GetCodonStateSpace()->Translation(i)] += mStationary[i];
	}
	double total = 0;
	for (int i=0; i<GetNstate(); i++)	{
		int aa = GetCodonStateSpace()->Translation(i);
		mStationary[i] *= exp((*Ne) * log(aaprofile[aa])) / aanorm[aa];
		total += mStationary[i];
	}
	for (int i=0; i<GetNstate(); i++)	{
		mStationary[i] /= total;
	}
}

void CodonMutSelProfileSubMatrix::ComputeArray(int i)	{

	double total = 0;
	for (int j=0; j<GetNstate(); j++)       {
		if (i!=j)       {
			int pos = GetDifferingPosition(i,j);
			if ((pos != -1) && (pos != 3))  {
				int a = GetCodonPosition(pos,i);
				int b = GetCodonPosition(pos,j);
				if (a == b)     {
					cerr << "identical states\n";
					cerr << GetCodonStateSpace()->GetState(i) << '\t' << GetCodonStateSpace()->GetState(j) << '\n';
					cerr << pos << '\n';
					exit(1);
				}
				//Q[i][j] = (*NucMatrix)(a,b);
				Q[i][j] = nucrr[GetNucRRIndex(a,b)] * nucstat[b];
				double deltaF = log((codonprofile)[j] / (codonprofile)[i]);  
				if (fabs(deltaF) < TOOSMALL)        {
					Q[i][j] /= ( 1.0 - (deltaF / 2) );
				}
				else if (deltaF > TOOLARGE)	{
					Q[i][j] *= deltaF;
				}
				else if (deltaF < TOOLARGENEGATIVE)	{
					Q[i][j] = 0;
				}
				else    {
					Q[i][j] *=  (deltaF)/(1.0 - exp(-deltaF));
				}
			}
			else    {
				Q[i][j] = 0;
			}
			total += Q[i][j];

			if (Q[i][j] < 0)        {
				cerr << "negative entry in matrix\n";
				exit(1);
			}
			if (isinf(Q[i][j]))	{
				cerr << "inf Q[i][j]\n";
				exit(1);
			}
			if (isnan(Q[i][j]))	{
				cerr << "nan Q[i][j]\n";
				exit(1);
			}
			
		}
	}
	Q[i][i] = -total;
	if (total <0)   {
		cerr << "negative rate away\n";
		exit(1);
	}
}

void CodonMutSelProfileSubMatrix::ComputeStationary()	{

	double total = 0;
	for (int i=0; i<GetNstate(); i++)	{
		mStationary[i] =nucstat[GetCodonPosition(0,i)] * 
				nucstat[GetCodonPosition(1,i)] * 
				nucstat[GetCodonPosition(2,i)] *
				codonprofile[i];
		total += mStationary[i];
	}

	// re-normalize
	for (int i=0; i<GetNstate(); i++)	{
		mStationary[i] /= total;
	}

}

double CodonMutSelProfileSubMatrix::GetRate()	{
	
	if (! ArrayUpdated())	{
		UpdateStationary();
		ComputeFullArray();
	}
	for (int k=0; k<Nstate; k++)	{
		flagarray[k] = true;
	}
	double norm = 0;
	for (int i=0; i<Nnuc-1; i++)	{
		for (int j=i+1; j<Nnuc; j++)	{
			norm += nucstat[i] * nucstat[j] * nucrr[GetNucRRIndex(i,j)];
		}
	}
	return 2 * (norm * 3);
	/*double mutstatnorm = 0;
	for (int i=0; i<Nstate; i++)	{
		mutstatnorm +=  nucstat[GetCodonPosition(0,i)] *
				nucstat[GetCodonPosition(1,i)] *
				nucstat[GetCodonPosition(2,i)];
	}	

	double norm = 0;
	int a, b;
	for (int i=0; i<Nstate-1; i++)	{
		for (int j=i+1; j<Nstate; j++)	{
			int pos = GetDifferingPosition(i,j);
			if ((pos != -1) && (pos != 3))  {
				a = GetCodonPosition(pos,i);
				b = GetCodonPosition(pos,j);
				norm += ((nucstat[GetCodonPosition(0,i)] * 
					nucstat[GetCodonPosition(1,i)] * 
					nucstat[GetCodonPosition(2,i)]) / mutstatnorm) *
					nucrr[GetNucRRIndex(a,b)] * nucstat[b];
			}
		}
	}
	return 2 * norm;
	*/
}

