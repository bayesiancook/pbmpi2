
#include "Random.h"
#include "Parallel.h"

#include "MPIModule.h"

void MPIModule::CreateMPI(int innsite)	{


	nsite = innsite;
	sitemin = new int[nprocs];
	sitemax = new int[nprocs];
	globalrank = new int[nsite];

	fmin = 0;
	fmax = 1;

	sitemin[0] = 0;
	sitemax[0] = nsite;
	if (nprocs > 1)	{
		maxwidth = 0;
		int width = nsite/(nprocs-1);
		for (int id=1; id<nprocs; id++)	{
			sitemin[id] = (id-1)*width;
			if (id == (nprocs-1)) {
				sitemax[id] = nsite;
			}
			else {
				sitemax[id] = id*width;
			}
			if (maxwidth < (sitemax[id] - sitemin[id]))	{
				maxwidth = sitemax[id] - sitemin[id];
			}
		}
	}
	else	{
		maxwidth = nsite;
	}

	if (nprocs > 1)	{

		if (! myid)	{
			GlobalReshuffleSites(0);
		}
		else	{
			MESSAGE signal;
			MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
			if (signal != RESHUFFLE)	{
				cerr << "error in create mpi\n";
				exit(1);
			}
			SlaveReshuffleSites();
		}
	}
	else	{
		NonMPIReshuffleSites(0);
	}
}

void MPIModule::NonMPIReshuffleSites(int rand)	{

	int* permut = new int[nsite];
	rnd::GetRandom().DrawFromUrn(permut,nsite,nsite);
	for (int i=0; i<nsite; i++)	{
		if (rand)	{
			globalrank[i] = permut[i];
		}
		else	{
			globalrank[i] = i;
		}
	}
	delete[] permut;
}

void MPIModule::GlobalReshuffleSites(int rand)	{

	if (nprocs == 1)	{
		NonMPIReshuffleSites(rand);
	}
	else	{

	MESSAGE signal = RESHUFFLE;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	int* tmprank = new int[nsite];
	int currentsize[nprocs];
	for (int j=1; j<nprocs; j++)	{
		currentsize[j] = 0;
	}
	int i = 0;
	int j = 1;
	while (i < nsite)	{
		if (currentsize[j] < (sitemax[j] - sitemin[j]))	{
			tmprank[sitemin[j] + currentsize[j]] = i;
			i++;
			currentsize[j]++;
		}
		j++;
		if (j == nprocs)	{
			j = 1;
		}
	}

	for (int j=1; j<nprocs; j++)	{
		int localnsite = sitemax[j] - sitemin[j];
		int* permut = new int[localnsite];
		rnd::GetRandom().DrawFromUrn(permut,localnsite,localnsite);
		for (int i=0; i<localnsite; i++)	{
			if (rand)	{
				globalrank[sitemin[j] + permut[i]] = tmprank[sitemin[j] + i];
			}
			else	{
				globalrank[sitemin[j] + i] = tmprank[sitemin[j] + i];
			}
		}
		delete[] permut;
	}

	delete[] tmprank;

	MPI_Bcast(globalrank,nsite,MPI_INT,0,MPI_COMM_WORLD);
	}
}

void MPIModule::SlaveReshuffleSites()	{

	MPI_Bcast(globalrank,nsite,MPI_INT,0,MPI_COMM_WORLD);
}

void MPIModule::GlobalWriteSiteRankToStream(ostream& os)	{

	for (int i=0; i<nsite; i++)	{
		os << globalrank[i] << '\t';
	}
	os << '\n';
}

void MPIModule::GlobalReadSiteRankFromStream(istream& is)	{

	for (int i=0; i<nsite; i++)	{
		is >> globalrank[i];
	}
	if (GetNprocs() > 1)	{
		MESSAGE signal = RESHUFFLE;
		MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(globalrank,nsite,MPI_INT,0,MPI_COMM_WORLD);
	}
}

int MPIModule::GetNactiveSite()	{

	int count = 0;
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		if (ActiveSite(i))	{
			count++;
		}
	}
	return count;
}

void MPIModule::SetPartition(string partitionfile)	{

	if (partalloc)	{
		cerr << "error in MPIModule::SetPartition: partition already allocated\n";
		exit(1);
	}
	ifstream is(partitionfile.c_str());
	int tmpnsite;
	is >> tmpnsite >> Npart;
	if (tmpnsite != nsite)	{
		cerr << "error in MPIModule::SetPartition: number of sites defined by datafile and partition file do not match\n";
		exit(1);
	}
	partalloc = new int[nsite];
	partnsite = new int[Npart];
	for (int part=0; part<Npart; part++)	{
		partnsite[part] = 0;
	}
	for (int i=0; i<nsite; i++)	{
		is >> partalloc[i];
		partnsite[partalloc[i]]++;
	}
}

