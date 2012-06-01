/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  (C) 2001 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include </exports/home/sabs/mirror/mpich2/include/mpi.h>
#include <stdio.h>
#include <math.h>
#include "System.h"

using namespace std;
using namespace MSL;

struct Job {
	int jobid;
	char instrn[10]; // END can be used to stop the slave process - all others may be ignored
};
MPI_Datatype job_MPIDataType;

void constructMessage(Job& _msg, int _id, string _instrn) {
	_msg.jobid = _id;
	_instrn.copy(_msg.instrn,_instrn.size());
	_msg.instrn[_instrn.size() + 1] = '\0';
	//cout << _msg.instrn << endl;

}

void constructJobDatatype () {
	MPI_Datatype array_of_types[2]={MPI_INT, MPI_CHAR};
	MPI_Aint array_of_displacements[2], intExt;
	MPI_Type_extent(MPI_INT, &intExt);
	array_of_displacements[0] = (MPI_Aint) 0; 
	array_of_displacements[1] = intExt;

	int array_of_blocklengths[2] = {1,10};
	MPI_Type_struct(2,array_of_blocklengths,array_of_displacements,array_of_types,&job_MPIDataType);
	MPI_Type_commit(&job_MPIDataType);
}

void master(int myid, int numprocs, int argc, char* argv[]) {
	int num_jobs = 200;
	int sent_jobs = 0;
	Job msg;
	while(sent_jobs < numprocs  -1) {
		// Send one job to each proc
		//MPI_Send(&sent_jobs,1,MPI_INT,sent_jobs + 1,1,MPI_COMM_WORLD);
		constructMessage(msg,sent_jobs,"RUN");
		MPI_Send(&msg,1,job_MPIDataType,sent_jobs + 1,1,MPI_COMM_WORLD);
		sent_jobs++;
	}
	int sent_end = 0;

	int ans = 0;
	MPI_Status status;
	vector<int> jobs(numprocs-1,0);
	while(sent_end < numprocs - 1) {
		// receive from any proc - the tag specifies the job it was received from 
		
		//MPI_Recv(&msg,1,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
		MPI_Recv(&msg,1,job_MPIDataType,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
		jobs[status.MPI_SOURCE -1]++;
		ans += msg.jobid;
		if(sent_jobs == num_jobs) {
			
			//cout << "Sending -1 to " << status.MPI_SOURCE << endl; 
			// if no more jobs send end (here -1)
			//MPI_Send(&msg,1,MPI_INT,status.MPI_SOURCE,1,MPI_COMM_WORLD);
			constructMessage(msg,-1,"END");
			MPI_Send(&msg,1,job_MPIDataType,status.MPI_SOURCE,1,MPI_COMM_WORLD);
			sent_end++;
		} else {
			// send a job to the proc you just received from
			//MPI_Send(&sent_jobs,1,MPI_INT,status.MPI_SOURCE,1,MPI_COMM_WORLD);
			constructMessage(msg,sent_jobs,"RUN");
			MPI_Send(&msg,1,job_MPIDataType,status.MPI_SOURCE,1,MPI_COMM_WORLD);
			sent_jobs++;
		}
	}

	/*
	for(int i = 0; i < jobs.size(); i++) {
		cout << i << " " << jobs[i] << endl;
	}
	*/
	cout << ans << endl;
}

void slave(int myid, int argc, char* argv[]) {
	MPI_Status status;
	Job msg;
	while(1) {
		MPI_Recv(&msg,1,job_MPIDataType,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
		if(msg.jobid == -1) {
			//cout << msg.instrn << endl;
			return;
		}
		// send the msg 
		cout << msg.jobid << endl;
		MPI_Send(&msg,1,job_MPIDataType,0,1,MPI_COMM_WORLD);
	}
}

int main(int argc,char *argv[])
{
    int    myid, numprocs;
    int    namelen;
    char   processor_name[MPI_MAX_PROCESSOR_NAME];

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    MPI_Get_processor_name(processor_name,&namelen);
    constructJobDatatype();

    /*
    fprintf(stdout,"Process %d of %d is on %s\n",
	    myid, numprocs, processor_name);
    fflush(stdout);
    */

    if (myid == 0) {
	    master(myid,numprocs,argc, argv);
    } else {
	    slave(myid,argc,argv);
    }

    MPI_Finalize();
    return 0;
}
