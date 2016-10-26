
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "Parallel.h"
#include "PartitionGTRProfileProcess.h"
#include "Random.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* PartitionGTRProfileProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------


void PartitionGTRProfileProcess::Create()	{
	if (! rr)	{
		ProfileProcess::Create();
		Nrr = GetDim() * (GetDim()-1) / 2;
		allocrr = new double[Npart*Nrr];
		rr = new double*[Npart];
		for (int part=0; part<Npart; part++)	{
			rr[part] = allocrr + part*Nrr;
		}
		if (rrtype != "None")	{
			// should be a file with a list
			cerr << "in PartitionGTR: setRR\n";
			exit(1);
			SetRR(rrtype);
		}
		else	{
			SampleRR();
		}
	}
}

void PartitionGTRProfileProcess::Delete()	{
	if (rr)	{
		delete[] rr;
		delete[] allocrr;
		rr = 0;
		ProfileProcess::Delete();
	}
}

double PartitionGTRProfileProcess::LogRRPrior()	{

	double tot = 0;
	for (int part=0; part<Npart; part++)	{
		tot += LogRRPrior(part);
	}
	return tot;
}

double PartitionGTRProfileProcess::LogRRPrior(int part)	{
	double total = 0;
	for (int i=0; i<GetNrr(); i++)	{
		total -= rr[part][i];
	}
	return total;
}

void PartitionGTRProfileProcess::SampleRR()	{
	for (int part=0; part<Npart; part++)	{
		for (int i=0; i<GetNrr(); i++)	{
			rr[part][i] = rnd::GetRandom().sExpo();
		}
	}
}

double PartitionGTRProfileProcess::GlobalParametersMove()	{
	if (! fixrr)	{
		GlobalUpdateParameters();
		GlobalUpdateSiteProfileSuffStat();
		GlobalUpdateModeProfileSuffStat();
		MoveRR();
		GlobalUpdateParameters();
	}
}

double  PartitionGTRProfileProcess::MoveRR(double tuning, int n, int nrep)	{
	double nacc = 0;
	for (int part=0; part<Npart; part++)	{
		nacc += MoveRR(part,tuning,n,nrep);
	}
	return nacc / Npart;
}

double  PartitionGTRProfileProcess::MoveRR(int part, double tuning, int n, int nrep)	{
	GlobalUpdateRRSuffStat();
	int naccepted = 0;

	int* choose = new int[n];
	double* bk = new double[n];

	for (int i=0; i<nrep; i++)	{

		rnd::GetRandom().DrawFromUrn(choose,n,Nrr);
		for (int j=0; j<n; j++)	{
			bk[j] = rr[part][choose[j]];
		}
		
		double deltalogratio = - LogRRPrior(part) - ProfileSuffStatLogProb();

		for (int j=0; j<n; j++)	{
			double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
			double e = exp(m);
			rr[part][choose[j]] *= e;
			deltalogratio += m;
		}

		UpdateMatrices();

		deltalogratio += LogRRPrior(part) + ProfileSuffStatLogProb();
		int accepted = (log(rnd::GetRandom().Uniform()) < deltalogratio);
		if (accepted)	{
			naccepted++;
		}
		else	{
			for (int j=0; j<n; j++)	{
				rr[part][choose[j]] = bk[j];
			}
			UpdateMatrices();
		}
	}

	delete[] bk;
	delete[] choose;

	return ((double) naccepted) / nrep;
}

void PartitionGTRProfileProcess::SetRR(string file)	{
	ifstream is(file.c_str());
	int tmp;
	is >> tmp;
	if (tmp != Npart)	{
		cerr << "error in PartitionGTRProfileProcess::SetRR: non matching number of partitions\n";
		exit(1);
	}
	for (int part=0; part<Npart; part++)	{
		string rrtype;
		is >> rrtype;
		SetRR(rrtype,rr[part]);
	}
}

void PartitionGTRProfileProcess::SetRR(string type, double* rr)	{

	rrtype = type;
	if (type != "None")	{
		fixrr = true;
		if ((type == "WAG") || (type == "wag"))	{
			if (Nrr != 190)	{
				if (! GetMyid())	{
					cerr << "error : wag only applies to amino acid recoded data\n";
					cerr << '\n';
				}
				MPI_Finalize();
				exit(1);
			}
			double total = 0;
			for (int i=0; i<Nrr; i++)	{
				rr[i]= WAG_RR[i];
				total += rr[i];
			}
			for (int i=0; i<Nrr; i++)	{
				rr[i] /= total / Nrr;
			}
		}

		else if ((type == "cg6") || (type == "CG6"))	{
			if (Nrr != 190)	{
				if (! GetMyid())	{
					cerr << "error : wag only applies to amino acid recoded data\n";
					cerr << '\n';
				}
				MPI_Finalize();
				exit(1);
			}
			double total = 0;
			for (int i=0; i<Nrr; i++)	{
				rr[i]= CG6RR[i];
				total += rr[i];
			}
			for (int i=0; i<Nrr; i++)	{
				rr[i] /= total / Nrr;
			}
		}

		else if ((type == "cg10") || (type == "CG10"))	{
			if (Nrr != 190)	{
				if (! GetMyid())	{
					cerr << "error : wag only applies to amino acid recoded data\n";
					cerr << '\n';
				}
				MPI_Finalize();
				exit(1);
			}
			double total = 0;
			for (int i=0; i<Nrr; i++)	{
				rr[i]= CG10RR[i];
				total += rr[i];
			}
			for (int i=0; i<Nrr; i++)	{
				rr[i] /= total / Nrr;
			}
		}

		else if ((type == "cg20") || (type == "CG20"))	{
			if (Nrr != 190)	{
				if (! GetMyid())	{
					cerr << "error : wag only applies to amino acid recoded data\n";
					cerr << '\n';
				}
				MPI_Finalize();
				exit(1);
			}
			double total = 0;
			for (int i=0; i<Nrr; i++)	{
				rr[i]= CG20RR[i];
				total += rr[i];
			}
			for (int i=0; i<Nrr; i++)	{
				rr[i] /= total / Nrr;
			}
		}

		else if ((type == "cg30") || (type == "CG30"))	{
			if (Nrr != 190)	{
				if (! GetMyid())	{
					cerr << "error : wag only applies to amino acid recoded data\n";
					cerr << '\n';
				}
				MPI_Finalize();
				exit(1);
			}
			double total = 0;
			for (int i=0; i<Nrr; i++)	{
				rr[i]= CG30RR[i];
				total += rr[i];
			}
			for (int i=0; i<Nrr; i++)	{
				rr[i] /= total / Nrr;
			}
		}

		else if ((type == "cg40") || (type == "CG40"))	{
			if (Nrr != 190)	{
				if (! GetMyid())	{
					cerr << "error : wag only applies to amino acid recoded data\n";
					cerr << '\n';
				}
				MPI_Finalize();
				exit(1);
			}
			double total = 0;
			for (int i=0; i<Nrr; i++)	{
				rr[i]= CG40RR[i];
				total += rr[i];
			}
			for (int i=0; i<Nrr; i++)	{
				rr[i] /= total / Nrr;
			}
		}

		else if ((type == "cg50") || (type == "CG50"))	{
			if (Nrr != 190)	{
				if (! GetMyid())	{
					cerr << "error : wag only applies to amino acid recoded data\n";
					cerr << '\n';
				}
				MPI_Finalize();
				exit(1);
			}
			double total = 0;
			for (int i=0; i<Nrr; i++)	{
				rr[i]= CG50RR[i];
				total += rr[i];
			}
			for (int i=0; i<Nrr; i++)	{
				rr[i] /= total / Nrr;
			}
		}

		else if ((type == "cg60") || (type == "CG60"))	{
			if (Nrr != 190)	{
				if (! GetMyid())	{
					cerr << "error : wag only applies to amino acid recoded data\n";
					cerr << '\n';
				}
				MPI_Finalize();
				exit(1);
			}
			double total = 0;
			for (int i=0; i<Nrr; i++)	{
				rr[i]= CG60RR[i];
				total += rr[i];
			}
			for (int i=0; i<Nrr; i++)	{
				rr[i] /= total / Nrr;
			}
		}

		else if ((type == "LG") || (type == "lg"))	{
			if (Nrr != 190)	{
				if (! GetMyid())	{
					cerr << "error : lg only applies to amino acid recoded data\n";
					cerr << '\n';
				}
				MPI_Finalize();
				exit(1);
			}
			double total = 0;
			for (int i=0; i<Nrr; i++)	{
				rr[i]= LG_RR[i];
				total += rr[i];
			}
			for (int i=0; i<Nrr; i++)	{
				rr[i] /= total / Nrr;
			}
		}

		else if ((type == "JTT") || (type == "jtt"))	{
			if (Nrr != 190)	{
				if (! GetMyid())	{
					cerr << "error : jtt only applies to amino acid recoded data\n";
					cerr << '\n';
				}
				MPI_Finalize();
				exit(1);
			}
			double total = 0;
			for (int i=0; i<Nrr; i++)	{
				rr[i]= JTT_RR[i];
				total += rr[i];
			}
			for (int i=0; i<Nrr; i++)	{
				rr[i] /= total / Nrr;
			}
		}

		else if ((type == "mtREV") || (type == "mtrev"))	{
			if (Nrr != 190)	{
				if (! GetMyid())	{
					cerr << "error : mtrev only applies to amino acid recoded data\n";
					cerr << '\n';
				}
				MPI_Finalize();
				exit(1);
			}
			double total = 0;
			for (int i=0; i<Nrr; i++)	{
				rr[i]= mtREV_RR[i];
				total += rr[i];
			}
			for (int i=0; i<Nrr; i++)	{
				rr[i] /= total / Nrr;
			}
		}

		else if ((type == "mtZOA") || (type == "mtzoa"))	{
			if (Nrr != 190)	{
				if (! GetMyid())	{
					cerr << "error : mtzoa only applies to amino acid recoded data\n";
					cerr << '\n';
				}
				MPI_Finalize();
				exit(1);
			}
			double total = 0;
			for (int i=0; i<Nrr; i++)	{
				rr[i]= mtZOA_RR[i];
				total += rr[i];
			}
			for (int i=0; i<Nrr; i++)	{
				rr[i] /= total / Nrr;
			}
		}

		else if ((type == "mtART") || (type == "mtart"))	{
			if (Nrr != 190)	{
				if (! GetMyid())	{
					cerr << "error : mtart only applies to amino acid recoded data\n";
					cerr << '\n';
				}
				MPI_Finalize();
				exit(1);
			}
			double total = 0;
			for (int i=0; i<Nrr; i++)	{
				rr[i]= mtART_RR[i];
				total += rr[i];
			}
			for (int i=0; i<Nrr; i++)	{
				rr[i] /= total / Nrr;
			}
		}

		else	{
			ifstream is(type.c_str());
			if (!is)	{
				cerr << "unrecognized file for relative rates : " << type << '\n';
				exit(1);
			}

			int permut[GetDim()];
			for (int k=0; k<GetDim(); k++)	{
				string c;
				is >> c;
				permut[k] = GetStateSpace()->GetState(c);
			}
			for (int k=0; k<GetDim()-1; k++)	{
				for (int l=k+1; l<GetDim(); l++)	{
					double tmp;
					is >> tmp;
					if (tmp < 0)	{
						if (! GetMyid())	{
							cerr << "error when reading exchangeabilities from " << type << '\n';
							cerr << tmp << '\n';
							cerr << '\n';
						}
						MPI_Finalize();
						exit(1);
					}
					rr[rrindex(permut[k],permut[l],GetDim())] = tmp;
					
				}
			}
		}	
	}
}


