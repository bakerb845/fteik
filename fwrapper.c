#include <stdio.h>
#include <stdlib.h>
#include "fteik_fortran.h"

//----------------------------------------------------------------------------//
//                                 Model                                      //
//----------------------------------------------------------------------------//
int fteik_model_grid2index(const int iz, const int ix, const int iy,
                           const int nz, const int nzx)
{
    int indx;
    indx = fteik_model_grid2indexF(iz+1, ix+1, iy+1, nz, nzx) - 1;
    return indx;
}
//----------------------------------------------------------------------------//
//                                Receiver                                    //
//----------------------------------------------------------------------------//
#ifdef FTEIK_USE_MPI
int fteik_receiver_broadcast(const int root, const MPI_Comm comm)
{
    void fteik_receiver_broadcastF(const int root, const int comm, int *mpierr);
    int mpierr;
    MPI_Fint fComm = MPI_Comm_c2f(comm);
    fteik_receiver_broadcastF(fComm, comm, &mpierr);
    return mpierr;
}
#endif

//----------------------------------------------------------------------------//
//                                Locate                                      //
//----------------------------------------------------------------------------//
void locate_setObservation64f(const int evnmbr, const int nttimes,
                              const bool lhaveOT, const double t0In,
                              const int *__restrict__ obs2tf,
                              const double *__restrict__ pickTimes,
                              const double *__restrict__ wts, int *ierr)
{
    int *obs2tfF = (int *) calloc((size_t) nttimes, sizeof(int));
    int i;
    for (i=0; i<nttimes; i++){obs2tfF[i] = obs2tf[i] + 1;} 
    locate_setObservation64fF(evnmbr+1, nttimes, lhaveOT, t0In, 
                              obs2tfF, pickTimes, wts, ierr); 
    free(obs2tfF);
    return; 
}

