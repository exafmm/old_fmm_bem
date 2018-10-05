#include "mpi.h"
#include "../include/constants.h"

extern int *nbuf,*nbufd;
extern float *fbuf,*fbufd;
extern double *dbuf,*dbufd;
extern std::complex<float> *cbuf,*cbufd;

extern void mpisendi(int*, int, int, int);
extern void mpisendrecvi(int, int, int, int);
extern void mpirecvi(int*, int, int, int);
extern void mpisendf(float*, int, int, int);
extern void mpisendrecvf(int, int, int, int);
extern void mpirecvf(float*, int, int, int);
extern void mpisendd(double*, int, int, int);
extern void mpisendrecvd(int, int, int, int);
extern void mpirecvd(double*, int, int, int);
extern void mpisendc(std::complex<float>*, int, int, int);
extern void mpisendrecvc(int, int, int, int);
extern void mpirecvc(std::complex<float>*, int, int, int);

// mpi_alltoallv for integer
void mpialltoallvi(int* nvar, int* nscnt, int* nsdsp, int* nvard, int* nrcnt, int* nrdsp, int nmax) {
  MPI_Barrier(MPI_COMM_WORLD);
  if( nmpi == 1 ) {
    MPI_Alltoallv(nvar,nscnt,nsdsp,MPI_INT,nvard,nrcnt,nrdsp,MPI_INT,MPI_COMM_WORLD);
  } else if( nmpi == 2 ) {
    int loclproc,loclrank,globproc,globrank;
    int mscnt[nprocs],msdsp[nprocs],mrcnt[nprocs],mrdsp[nprocs];
    MPI_Comm MPI_COMM_LOCAL,MPI_COMM_GLOBAL;
    MPI_Comm_split(MPI_COMM_WORLD,myrank/nsplit,myrank%nsplit,&MPI_COMM_LOCAL);
    MPI_Comm_split(MPI_COMM_WORLD,myrank%nsplit,myrank/nsplit,&MPI_COMM_GLOBAL);
    MPI_Comm_size(MPI_COMM_LOCAL,&loclproc);
    MPI_Comm_rank(MPI_COMM_LOCAL,&loclrank);
    MPI_Comm_size(MPI_COMM_GLOBAL,&globproc);
    MPI_Comm_rank(MPI_COMM_GLOBAL,&globrank);
    mrdsp[0] = 0;
    for(int i=0; i<globproc; i++) {
      MPI_Alltoall(&nscnt[i*loclproc],1,MPI_INT,mrcnt,1,MPI_INT,MPI_COMM_LOCAL);
      mscnt[i] = 0;
      for(int j=0; j<loclproc-1; j++) {
        mrdsp[j+1] = mrdsp[j]+mrcnt[j];
        mscnt[i] += mrcnt[j];
      }
      mscnt[i] += mrcnt[loclproc-1];
      MPI_Alltoallv(nvar,&nscnt[i*loclproc],&nsdsp[i*loclproc],MPI_INT,
                    nbuf,mrcnt,mrdsp,MPI_INT,MPI_COMM_LOCAL);
      mrdsp[0] += mscnt[i];
    }
    MPI_Alltoall(mscnt,1,MPI_INT,mrcnt,1,MPI_INT,MPI_COMM_GLOBAL);
    msdsp[0] = 0;
    mrdsp[0] = 0;
    for(int i=0; i<globproc-1; i++) {
      msdsp[i+1] = msdsp[i] + mscnt[i];
      mrdsp[i+1] = mrdsp[i] + mrcnt[i];
    }
    MPI_Alltoallv(nbuf,mscnt,msdsp,MPI_INT,nvard,mrcnt,mrdsp,MPI_INT,MPI_COMM_GLOBAL);
  } else if( nmpi == 3 ) {
    int lsend,lrecv,isize,jsize,ista=0,iend=0,jsta=0,jend=0,irank,i;
    MPI_Status istatus;
    MPI_Request ireq1,ireq2;
  
    for( irank=0; irank<nprocs; irank++ ) {
      lrecv = (myrank+irank)%nprocs;
      lsend = (myrank-irank+nprocs)%nprocs;
      ista = nsdsp[lrecv];
      isize = nscnt[lrecv];
      iend = ista+isize-1;
      jsta = nrdsp[lsend];
      jsize = nrcnt[lsend];
      jend = jsta+jsize-1;
      for( i=ista; i<=iend; i++ ) nbuf[i-ista] = nvar[i];
      if( irank != 0 ) {
        MPI_Isend(nbuf,isize,MPI_INT,lrecv,1,MPI_COMM_WORLD,&ireq1);
        MPI_Irecv(nbufd,jsize,MPI_INT,lsend,1,MPI_COMM_WORLD,&ireq2);
        MPI_Wait(&ireq1,&istatus);
        MPI_Wait(&ireq2,&istatus);
      } else {
        for( i=0; i<isize; i++ ) nbufd[i] = nbuf[i];
      }
      for( i=jsta; i<=jend; i++ ) nvard[i] = nbufd[i-jsta];
    }
  }
}

// mpi_alltoallv for float
void mpialltoallvf(float* fvar, int* nscnt, int* nsdsp, float* fvard, int* nrcnt, int* nrdsp, int nmax) {
  MPI_Barrier(MPI_COMM_WORLD);
  if( nmpi == 1 ) {
    MPI_Alltoallv(fvar,nscnt,nsdsp,MPI_FLOAT,fvard,nrcnt,nrdsp,MPI_FLOAT,MPI_COMM_WORLD);
  } else if( nmpi == 2 ) {
    int loclproc,loclrank,globproc,globrank;
    int mscnt[nprocs],msdsp[nprocs],mrcnt[nprocs],mrdsp[nprocs];
    MPI_Comm MPI_COMM_LOCAL,MPI_COMM_GLOBAL;
    MPI_Comm_split(MPI_COMM_WORLD,myrank/nsplit,myrank%nsplit,&MPI_COMM_LOCAL);
    MPI_Comm_split(MPI_COMM_WORLD,myrank%nsplit,myrank/nsplit,&MPI_COMM_GLOBAL);
    MPI_Comm_size(MPI_COMM_LOCAL,&loclproc);
    MPI_Comm_rank(MPI_COMM_LOCAL,&loclrank);
    MPI_Comm_size(MPI_COMM_GLOBAL,&globproc);
    MPI_Comm_rank(MPI_COMM_GLOBAL,&globrank);
    mrdsp[0] = 0;
    for(int i=0; i<globproc; i++) {
      MPI_Alltoall(&nscnt[i*loclproc],1,MPI_INT,mrcnt,1,MPI_INT,MPI_COMM_LOCAL);
      mscnt[i] = 0;
      for(int j=0; j<loclproc-1; j++) {
        mrdsp[j+1] = mrdsp[j]+mrcnt[j];
        mscnt[i] += mrcnt[j];
      }
      mscnt[i] += mrcnt[loclproc-1];
      MPI_Alltoallv(fvar,&nscnt[i*loclproc],&nsdsp[i*loclproc],MPI_FLOAT,
                    fbuf,mrcnt,mrdsp,MPI_FLOAT,MPI_COMM_LOCAL);
      mrdsp[0] += mscnt[i];
    }
    MPI_Alltoall(mscnt,1,MPI_INT,mrcnt,1,MPI_INT,MPI_COMM_GLOBAL);
    msdsp[0] = 0;
    mrdsp[0] = 0;
    for(int i=0; i<globproc-1; i++) {
      msdsp[i+1] = msdsp[i]+mscnt[i];
      mrdsp[i+1] = mrdsp[i]+mrcnt[i];
    }
    MPI_Alltoallv(fbuf,mscnt,msdsp,MPI_FLOAT,fvard,mrcnt,mrdsp,MPI_FLOAT,MPI_COMM_GLOBAL);
  } else if( nmpi == 3 ) {
    int lsend,lrecv,isize,jsize,ista=0,iend=0,jsta=0,jend=0,irank,i;
    MPI_Status istatus;
    MPI_Request ireq1,ireq2;

    for( irank=0; irank<nprocs; irank++ ) {
      lrecv = (myrank+irank)%nprocs;
      lsend = (myrank-irank+nprocs)%nprocs;
      ista = nsdsp[lrecv];
      isize = nscnt[lrecv];
      iend = ista+isize-1;
      jsta = nrdsp[lsend];
      jsize = nrcnt[lsend];
      jend = jsta+jsize-1;
      for( i=ista; i<=iend; i++ ) fbuf[i-ista] = fvar[i];
      if( irank != 0 ) {
        MPI_Isend(fbuf,isize,MPI_FLOAT,lrecv,1,MPI_COMM_WORLD,&ireq1);
        MPI_Irecv(fbufd,jsize,MPI_FLOAT,lsend,1,MPI_COMM_WORLD,&ireq2);
        MPI_Wait(&ireq1,&istatus);
        MPI_Wait(&ireq2,&istatus);
      } else {
        for( i=0; i<isize; i++ ) fbufd[i] = fbuf[i];
      }
      for( i=jsta; i<=jend; i++ ) fvard[i] = fbufd[i-jsta];
    } 
  }
}

// mpi_alltoallv for double
void mpialltoallvd(double* dvar, int* nscnt, int* nsdsp, double* dvard, int* nrcnt, int* nrdsp, int nmax) {
  MPI_Barrier(MPI_COMM_WORLD);
  if( nmpi == 1 ) {
    MPI_Alltoallv(dvar,nscnt,nsdsp,MPI_DOUBLE,dvard,nrcnt,nrdsp,MPI_DOUBLE,MPI_COMM_WORLD);
  } else if( nmpi == 2 ) {
    int loclproc,loclrank,globproc,globrank;
    int mscnt[nprocs],msdsp[nprocs],mrcnt[nprocs],mrdsp[nprocs];
    MPI_Comm MPI_COMM_LOCAL,MPI_COMM_GLOBAL;
    MPI_Comm_split(MPI_COMM_WORLD,myrank/nsplit,myrank%nsplit,&MPI_COMM_LOCAL);
    MPI_Comm_split(MPI_COMM_WORLD,myrank%nsplit,myrank/nsplit,&MPI_COMM_GLOBAL);
    MPI_Comm_size(MPI_COMM_LOCAL,&loclproc);
    MPI_Comm_rank(MPI_COMM_LOCAL,&loclrank);
    MPI_Comm_size(MPI_COMM_GLOBAL,&globproc);
    MPI_Comm_rank(MPI_COMM_GLOBAL,&globrank);
    mrdsp[0] = 0;
    for(int i=0; i<globproc; i++) {
      MPI_Alltoall(&nscnt[i*loclproc],1,MPI_INT,mrcnt,1,MPI_INT,MPI_COMM_LOCAL);
      mscnt[i] = 0;
      for(int j=0; j<loclproc-1; j++) {
        mrdsp[j+1] = mrdsp[j]+mrcnt[j];
        mscnt[i] += mrcnt[j];
      }
      mscnt[i] += mrcnt[loclproc-1];
      MPI_Alltoallv(dvar,&nscnt[i*loclproc],&nsdsp[i*loclproc],MPI_DOUBLE,
                    dbuf,mrcnt,mrdsp,MPI_DOUBLE,MPI_COMM_LOCAL);
      mrdsp[0] += mscnt[i];
    }
    MPI_Alltoall(mscnt,1,MPI_INT,mrcnt,1,MPI_INT,MPI_COMM_GLOBAL);
    msdsp[0] = 0;
    mrdsp[0] = 0;
    for(int i=0; i<globproc-1; i++) {
      msdsp[i+1] = msdsp[i] + mscnt[i];
      mrdsp[i+1] = mrdsp[i] + mrcnt[i];
    }
    MPI_Alltoallv(dbuf,mscnt,msdsp,MPI_DOUBLE,dvard,mrcnt,mrdsp,MPI_DOUBLE,MPI_COMM_GLOBAL);
  } else if( nmpi == 3 ) {
    int lsend,lrecv,isize,jsize,ista=0,iend=0,jsta=0,jend=0,irank,i;
    MPI_Status istatus;
    MPI_Request ireq1,ireq2;

    for( irank=0; irank<nprocs; irank++ ) {
      lrecv = (myrank+irank)%nprocs;
      lsend = (myrank-irank+nprocs)%nprocs;
      ista = nsdsp[lrecv];
      isize = nscnt[lrecv];
      iend = ista+isize-1;
      jsta = nrdsp[lsend];
      jsize = nrcnt[lsend];
      jend = jsta+jsize-1;
      for( i=ista; i<=iend; i++ ) dbuf[i-ista] = dvar[i];
      if( irank != 0 ) {
        MPI_Isend(dbuf,isize,MPI_DOUBLE,lrecv,1,MPI_COMM_WORLD,&ireq1);
        MPI_Irecv(dbufd,jsize,MPI_DOUBLE,lsend,1,MPI_COMM_WORLD,&ireq2);
        MPI_Wait(&ireq1,&istatus);
        MPI_Wait(&ireq2,&istatus);
      } else {
        for( i=0; i<isize; i++ ) dbufd[i] = dbuf[i];
      }
      for( i=jsta; i<=jend; i++ ) dvard[i] = dbufd[i-jsta];
    }
  }
}

// mpi_alltoallv for complex
void mpialltoallvc(std::complex<float>* cvar, int* nscnt, int* nsdsp, std::complex<float>* cvard, int* nrcnt, int* nrdsp, int nmax) {
  MPI_Barrier(MPI_COMM_WORLD);
  if( nmpi == 1 ) {
    MPI_Alltoallv(cvar,nscnt,nsdsp,MPI_COMPLEX,cvard,nrcnt,nrdsp,MPI_COMPLEX,MPI_COMM_WORLD);
  } else if( nmpi == 2 ) {
    int loclproc,loclrank,globproc,globrank;
    int mscnt[nprocs],msdsp[nprocs],mrcnt[nprocs],mrdsp[nprocs];
    MPI_Comm MPI_COMM_LOCAL,MPI_COMM_GLOBAL;
    MPI_Comm_split(MPI_COMM_WORLD,myrank/nsplit,myrank%nsplit,&MPI_COMM_LOCAL);
    MPI_Comm_split(MPI_COMM_WORLD,myrank%nsplit,myrank/nsplit,&MPI_COMM_GLOBAL);
    MPI_Comm_size(MPI_COMM_LOCAL,&loclproc);
    MPI_Comm_rank(MPI_COMM_LOCAL,&loclrank);
    MPI_Comm_size(MPI_COMM_GLOBAL,&globproc);
    MPI_Comm_rank(MPI_COMM_GLOBAL,&globrank);
    mrdsp[0] = 0;
    for(int i=0; i<globproc; i++) {
      MPI_Alltoall(&nscnt[i*loclproc],1,MPI_INT,mrcnt,1,MPI_INT,MPI_COMM_LOCAL);
      mscnt[i] = 0;
      for(int j=0; j<loclproc-1; j++) {
        mrdsp[j+1] = mrdsp[j]+mrcnt[j];
        mscnt[i] += mrcnt[j];
      }
      mscnt[i] += mrcnt[loclproc-1];
      MPI_Alltoallv(cvar,&nscnt[i*loclproc],&nsdsp[i*loclproc],MPI_COMPLEX,
                    cbuf,mrcnt,mrdsp,MPI_COMPLEX,MPI_COMM_LOCAL);
      mrdsp[0] += mscnt[i];
    }
    MPI_Alltoall(mscnt,1,MPI_INT,mrcnt,1,MPI_INT,MPI_COMM_GLOBAL);
    msdsp[0] = 0; 
    mrdsp[0] = 0; 
    for(int i=0; i<globproc-1; i++) {
      msdsp[i+1] = msdsp[i] + mscnt[i];
      mrdsp[i+1] = mrdsp[i] + mrcnt[i];
    } 
    MPI_Alltoallv(cbuf,mscnt,msdsp,MPI_COMPLEX,cvard,mrcnt,mrdsp,MPI_COMPLEX,MPI_COMM_GLOBAL);
  } else if( nmpi == 3 ) {
    int lsend,lrecv,isize,jsize,ista=0,iend=0,jsta=0,jend=0,irank,i;
    MPI_Status istatus;
    MPI_Request ireq1,ireq2;

    for( irank=0; irank<nprocs; irank++ ) {
      lrecv = (myrank+irank)%nprocs;
      lsend = (myrank-irank+nprocs)%nprocs;
      ista = nsdsp[lrecv];
      isize = nscnt[lrecv];
      iend = ista+isize-1;
      jsta = nrdsp[lsend];
      jsize = nrcnt[lsend];
      jend = jsta+jsize-1;
      for( i=ista; i<=iend; i++ ) cbuf[i-ista] = cvar[i];
      if( irank != 0 ) {
        MPI_Isend(cbuf,isize,MPI_COMPLEX,lrecv,1,MPI_COMM_WORLD,&ireq1);
        MPI_Irecv(cbufd,jsize,MPI_COMPLEX,lsend,1,MPI_COMM_WORLD,&ireq2);
        MPI_Wait(&ireq1,&istatus);
        MPI_Wait(&ireq2,&istatus);
      } else {
        for( i=0; i<isize; i++ ) cbufd[i] = cbuf[i];
      }
      for( i=jsta; i<=jend; i++ ) cvard[i] = cbufd[i-jsta];
    }
  }
}
