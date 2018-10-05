#include "mpi.h"
#include "../include/constants.h"

extern float *fbuf;

extern void mpisendf(float*, int, int, int);
extern void mpisendrecvf(int, int, int, int);
extern void mpirecvf(float*, int, int, int);

void mpishift(int mi, int& mid, float* fvar) {
  int lsend,lrecv,i;
  MPI_Status istatus;
  MPI_Request ireq1,ireq2;

  lrecv = (myrank+1)%nprocs;
  lsend = (myrank-1+nprocs)%nprocs;

  MPI_Isend(&mi,1,MPI_INT,lrecv,1,MPI_COMM_WORLD,&ireq1);
  MPI_Irecv(&mid,1,MPI_INT,lsend,1,MPI_COMM_WORLD,&ireq2);
  MPI_Wait(&ireq1,&istatus);
  MPI_Wait(&ireq2,&istatus);

  MPI_Isend(fvar,mi,MPI_FLOAT,lrecv,2,MPI_COMM_WORLD,&ireq1);
  MPI_Irecv(fbuf,mid,MPI_FLOAT,lsend,2,MPI_COMM_WORLD,&ireq2);
  MPI_Wait(&ireq1,&istatus);
  MPI_Wait(&ireq2,&istatus);

  for( i=0; i<mid; i++ ) fvar[i] = fbuf[i];
}
