#include "mpi.h"
#include "../include/constants.h"

extern int *nbuf,*nbufd;
extern float *fbuf,*fbufd;
extern double *dbuf,*dbufd;
extern std::complex<float> *cbuf,*cbufd;

// mpi_sendrecv for integer
void mpisendi(int* isend, int ista, int iend, int lsend) {
  int i,ic;

  if( myrank == lsend ) {
    for( i=ista; i<=iend; i++ ) {
      ic = i-ista;
      nbuf[ic] = isend[i];
    }
  }
}

void mpisendrecvi(int nsize, int lsend, int lrecv, int ntag) {
  int ncall,icall,iwork1,iwork2,ista,iend,isize,i;
  MPI_Status istatus;

  nbufd = new int [mpim];

  if( lsend != lrecv ) {
    ncall = nsize/mpim+1;
    for( icall=0; icall<ncall; icall++ ) {
      iwork1 = nsize/ncall;
      iwork2 = nsize%ncall;
      ista = icall*iwork1+std::min(icall,iwork2);
      iend = ista+iwork1-1;
      if( iwork2 > icall ) iend++;
      isize = iend-ista+1;
      for( i=ista; i<=iend; i++ ) nbufd[i-ista] = nbuf[i];
      if( myrank == lsend ) {
        MPI_Send(nbufd,isize,MPI_INT,lrecv,ntag,MPI_COMM_WORLD);
      } else if ( myrank == lrecv ) {
        MPI_Recv(nbufd,isize,MPI_INT,lsend,ntag,MPI_COMM_WORLD,&istatus);
      }
      for( i=ista; i<=iend; i++ ) nbuf[i] = nbufd[i-ista];
    }
  }

  delete[] nbufd;
}

void mpirecvi(int* irecv, int ista, int iend, int lrecv) {
  int i,ic;

  if( myrank == lrecv ) {
    for( i=ista; i<=iend; i++ ) {
      ic = i-ista;
      irecv[i] = nbuf[ic];
    }
  }
}

// mpi_sendrecv for float
void mpisendf(float* fsend, int ista, int iend, int lsend) {
  int i,ic;

  if( myrank == lsend ) {
    for( i=ista; i<=iend; i++ ) {
      ic = i-ista;
      fbuf[ic] = fsend[i];
    }
  }
}

void mpisendrecvf(int nsize, int lsend, int lrecv, int ntag) {
  int ncall,icall,iwork1,iwork2,ista,iend,isize,i;
  MPI_Status istatus;

  fbufd = new float [mpim];

  if( lsend != lrecv ) {
    ncall = nsize/mpim+1;
    for( icall=0; icall<ncall; icall++ ) {
      iwork1 = nsize/ncall;
      iwork2 = nsize%ncall;
      ista = icall*iwork1+std::min(icall,iwork2);
      iend = ista+iwork1-1;
      if( iwork2 > icall ) iend++;
      isize = iend-ista+1;
      for( i=ista; i<=iend; i++ ) fbufd[i-ista] = fbuf[i];
      if( myrank == lsend ) {
        MPI_Send(fbufd,isize,MPI_FLOAT,lrecv,ntag,MPI_COMM_WORLD);
      } else if ( myrank == lrecv ) {
        MPI_Recv(fbufd,isize,MPI_FLOAT,lsend,ntag,MPI_COMM_WORLD,&istatus);
      }
      for( i=ista; i<=iend; i++ ) fbuf[i] = fbufd[i-ista];
    }
  }

  delete[] fbufd;
}

void mpirecvf(float* frecv, int ista, int iend, int lrecv) {
  int i,ic;

  if( myrank == lrecv ) {
    for( i=ista; i<=iend; i++ ) {
      ic = i-ista;
      frecv[i] = fbuf[ic];
    }
  }
}

// mpi_sendrecv for double
void mpisendd(double* dsend, int ista, int iend, int lsend) {
  int i,ic;

  if( myrank == lsend ) {
    for( i=ista; i<=iend; i++ ) {
      ic = i-ista;
      dbuf[ic] = dsend[i];
    }
  }
}

void mpisendrecvd(int nsize, int lsend, int lrecv, int ntag) {
  int ncall,icall,iwork1,iwork2,ista,iend,isize,i;
  MPI_Status istatus;

  dbufd = new double [mpim];

  if( lsend != lrecv ) {
    ncall = nsize/mpim+1;
    for( icall=0; icall<ncall; icall++ ) {
      iwork1 = nsize/ncall;
      iwork2 = nsize%ncall;
      ista = icall*iwork1+std::min(icall,iwork2);
      iend = ista+iwork1-1;
      if( iwork2 > icall ) iend++;
      isize = iend-ista+1;
      for( i=ista; i<=iend; i++ ) dbufd[i-ista] = dbuf[i];
      if( myrank == lsend ) {
        MPI_Send(dbufd,isize,MPI_DOUBLE,lrecv,ntag,MPI_COMM_WORLD);
      } else if ( myrank == lrecv ) {
        MPI_Recv(dbufd,isize,MPI_DOUBLE,lsend,ntag,MPI_COMM_WORLD,&istatus);
      }
      for( i=ista; i<=iend; i++ ) dbuf[i] = dbufd[i-ista];
    }
  }

  delete[] dbufd;
}

void mpirecvd(double* drecv, int ista, int iend, int lrecv) {
  int i,ic;

  if( myrank == lrecv ) {
    for( i=ista; i<=iend; i++ ) {
      ic = i-ista;
      drecv[i] = dbuf[ic];
    }
  }
}

// mpi_sendrecv for complex
void mpisendc(std::complex<float>* csend, int ista, int iend, int lsend) {
  int i,ic;

  if( myrank == lsend ) {
    for( i=ista; i<=iend; i++ ) {
      ic = i-ista;
      cbuf[ic] = csend[i];
    }
  }
}

void mpisendrecvc(int nsize, int lsend, int lrecv, int ntag) {
  int ncall,icall,iwork1,iwork2,ista,iend,isize,i;
  MPI_Status istatus;

  cbufd = new std::complex<float> [mpim];

  if( lsend != lrecv ) {
    ncall = nsize/mpim+1;
    for( icall=0; icall<ncall; icall++ ) {
      iwork1 = nsize/ncall;
      iwork2 = nsize%ncall;
      ista = icall*iwork1+std::min(icall,iwork2);
      iend = ista+iwork1-1;
      if( iwork2 > icall ) iend++;
      isize = iend-ista+1;
      for( i=ista; i<=iend; i++ ) cbufd[i-ista] = cbuf[i];
      if( myrank == lsend ) {
        MPI_Send(cbufd,isize,MPI_COMPLEX,lrecv,ntag,MPI_COMM_WORLD);
      } else if ( myrank == lrecv ) {
        MPI_Recv(cbufd,isize,MPI_COMPLEX,lsend,ntag,MPI_COMM_WORLD,&istatus);
      }
      for( i=ista; i<=iend; i++ ) cbuf[i] = cbufd[i-ista];
    }
  }

  delete[] cbufd;
}

void mpirecvc(std::complex<float>* crecv, int ista, int iend, int lrecv) {
  int i,ic;

  if( myrank == lrecv ) {
    for( i=ista; i<=iend; i++ ) {
      ic = i-ista;
      crecv[i] = cbuf[ic];
    }
  }
}
