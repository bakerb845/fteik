#ifndef FTEIK_FORTRAN_H__
#define FTEIK_FORTRAN_H__ 1
#include "fteik/fteik_config.h"
#ifdef FTEIK_USE_MPI
#include <mpi.h>
#endif
#include <stdbool.h>

#ifdef __cplusplus
extern "C"
{
#endif

//----------------------------------------------------------------------------//
//                              Model Module                                  //
//----------------------------------------------------------------------------//
void fteik_model_free(void);
void fteik_model_getVelocityModel64f(const int nwork, int *nv, 
                                     double *__restrict__ vel, int *ierr);
void fteik_model_intializeGeometry(
    const int nz, const int nx, const int ny,
    const double dz, const double dx, const double dy,
    const double z0, const double x0, const double y0,
    int *ierr);
void fteik_model_setGridSize(const int nz, const int nx, const int ny,
                             int *ierr);
void fteik_model_setOrigin(const double z0, const double x0, const double y0);
void fteik_model_setGridSpacing( 
    const double dz, const double dx, const double dy, int *ierr);
void fteik_model_setVelocityModel64f(
    const int nv, const double *__restrict__ vel, int *ierr);
void fteik_model_setVelocityModel32f(
    const int nv, const float *__restrict__ vel4, int *ierr);
void fteik_model_isModel3D(bool *lis3d, int *ierr);
#pragma omp declare simd uniform(nz, nzx)
int fteik_model_grid2indexF(const int iz, const int ix, const int iy, 
                            const int nz, const int nzx);
#pragma omp declare simd uniform(nz, nzx)
int fteik_model_grid2index(const int iz, const int ix, const int iy,
                           const int nz, const int nzx);
#pragma omp declare simd uniform (ierr)
void fteik_model_index2gridF(const int igrd, int *iz, int *ix, int *iy,
                             int *ierr);
#pragma omp declare simd uniform (nzm1, nzm1_nxm1)
int fteik_model_velGrid2indexF(const int iz, const int ix, const int iy,
                               const int nzm1, const int nzm1_nxm1);
//----------------------------------------------------------------------------//
//                            Receiver Module                                 //
//----------------------------------------------------------------------------//
void fteik_receiver_initialize64f(const int nrec,
                                  const double *__restrict__ z,
                                  const double *__restrict__ x,
                                  const double *__restrict__ y,
                                  const int verbose, int *ierr);
void fteik_receiver_getTravelTimes64f(const int nrec, const int ngrd,
                                      const double ttimes[], double ttr[],
                                      int *ierr);
void fteik_receiver_free(void);

#ifdef FTEIK_USE_MPI
int fteik_receiver_broadcast(const int root, const MPI_Comm comm);
#endif
//----------------------------------------------------------------------------//
//                            Source Module                                   //
//----------------------------------------------------------------------------//
void fteik_source_free(void);
void fteik_source_initialize64f(const int nsrcIn,
                                const double *__restrict__ zsrcIn,
                                const double *__restrict__ xsrcIn,
                                const double *__restrict__ ysrcIn,
                                const int verbose, int *ierr);
void fteik_source_getSolverInfo64fF(const int isrc, 
                                    int *zsi, int *xsi, int *ysi,
                                    double *zsa, double *xsa, double *ysa,
                                    double *szero, double *szero2,
                                    int *ierr);
/*
void fteik_source_getSourceIndices64fF(const int isrc,
                                       double *zsa, double *xsa, double *ysa,
                                       int *ierr);
void fteik_source_getSourceIndices32iF(const int isrc,
                                       int *zsi, int *xsi, int *ysi, int *ierr);
void fteik_source_getSzero64fF(const int isrc,
                               double *szero, double *szero2, int *ierr);
*/
//----------------------------------------------------------------------------//
//                             Solver Module                                  //
//----------------------------------------------------------------------------//
void fteik_solver3d_initialize64f(
    const int nzIn, const int nxIn, const int nyIn,
    const double z0In, const double x0In, const double y0In,
    const double dzIn, const double dxIn, const double dyIn,
    const int nsweepIn, const double epsIn,
    const int verbose, int *ierr);
void fteik_solver3d_free(void);
void fteik_solver3d_setVelocityModel64f(const int ncell, 
                                        const double vel[],
                                        int *ierr);
void fteik_solver3d_setSources64f(const int nsrc,
                                  const double zsrc[],
                                  const double xsrc[],
                                  const double ysrc[],
                                  int *ierr);
void fteik_solver3d_getTravelTimes64f(const int nrec, double ttr[], int *ierr);
void fteik_solver3d_getNumberOfSources(int *nsrc, int *ierr);
void fteik_solver3d_getNumberOfReceivers(int *nrec, int *ierr);
void fteik_solver3d_setReceivers64f(const int nrec,
                                    const double zrec[],
                                    const double xrec[],
                                    const double yrec[],
                                    int *ierr);
void fteik_solver3d_solveSourceLSM(const int isrc, int *ierr);
//----------------------------------------------------------------------------//
//                           2D Solver Module                                 //
//----------------------------------------------------------------------------//
void fteik_solver2d_initialize64f(const int nz, const int nx,
                                  const double z0, const double x0,
                                  const double dz, const double dx,
                                  const int nsweep, const double eps,
                                  const int verbose, int *ierr);
void fteik_solver2d_setVelocityModel64f(const int ncell,
                                        const double vel[],
                                        int *ierr);
void fteik_solver2d_setReceivers64f(const int nrec,
                                    const double zrec[],
                                    const double xrec[],
                                    int *ierr);
void fteik_solver2d_solveSourceLSM(const int isrc, int *ierr);
void fteik_solver2d_solveSourceFSM(const int isrc, int *ierr);
void fteik_solver2d_getTravelTimeField64f(const int ngin,
                                          double ttout[], int *ierr);
void fteik_solver2d_setSources64f(const int nsrc,
                                  const double zsrc[],
                                  const double xsrc[],
                                  int *ierr);
void fteik_solver2d_getTravelTimes64f(const int nrec, double ttr[], int *ierr);
void fteik_solver2d_getNumberOfReceivers(int *nrec, int *ierr);
void fteik_solver2d_free(void);
//----------------------------------------------------------------------------//
//                             Locator Module                                 //
//----------------------------------------------------------------------------//
void locate_initializeF(const int nEventsIn, const int ntfIn, const int ngrdIn,
                        int *ierr);
void locate_setTravelTimeField64fF(const int ngrdIn, const int itf,
                                   const double *__restrict__ ttIn,
                                   int *ierr);
void locate_setTravelTimeField32fF(const int ngrdIn, const int itf,
                                   const double *__restrict__ ttIn,
                                   int *ierr);
void locate_finalizeF(void);
void locate_setObservation64fF(const int evnmbr, const int nttimes,
                               const bool lhaveOT, const double t0In,
                               const int *__restrict__ obs2tf,
                               const double *__restrict__  pickTimes,
                               const double *__restrict__  wts,
                               int *ierr);
void locate_locateEventF(const int evnmbr, int *optindx,
                         double *t0Opt, double *objOpt);
void locate_setObservation64f(const int evnmbr, const int nttimes,
                              const bool lhaveOT, const double t0In,
                              const int *__restrict__ obs2tf,
                              const double *__restrict__ pickTimes,
                              const double *__restrict__ wts, int *ierr);
#ifdef FTEIK_USE_MPI
void locate_initializeMPI(const int root, const MPI_Comm comm,
                          const int nEventsIn, const int ntfIn,
                          const int ngrdIn, int *ierr);
void locate_initializeMPIF(const int root, const int comm,
                           const int nEventsIn, const int ntfIn,
                           const int ngrdIn, int *ierr);
void locate_setTravelTimeField64fMPIF(const int root, const int ngrd,
                                      const int itf,
                                      const double *__restrict__ tt,
                                      int *ierr);
void locate_locateEventF(const int evnmbr, int *optindx,
                         double *t0Opt, double *objOpt);
void locate_locateEventMPIF(const int evnmbr, int *optIndx,
                            double *t0Opt, double *objOpt);

#endif
//----------------------------------------------------------------------------//
//                            3D Graph Reordering                             //
//----------------------------------------------------------------------------//
/* Initializes the 3D graph structure. */
void fteik_graph3d_initialize(const int nzIn, const int nxIn, const int nyIn,
                              const int verboseIn, int *ierr);
/* Makes the level structure */
void fteik_graph3d_makeLevelStructures(int *ierr);
/* Frees memory on the 3D graph. */
void fteik_graph3d_free(void);
void fteik_graph3d_getLevelPointerF(const int nwork, int levelPtr[], int *ierr);
int fteik_graph3d_getMaxLevelSize(int *ierr);
int fteik_graph3d_getNumberOfLevels(int *ierr);
void fteik_graph3d_getIJKVF(const int nwork, const int sweep,
                            int ijkv[], int *ierr);

//----------------------------------------------------------------------------//

// These are the functions worth knowing
void fteik_initializeF(const int *nz, const int *nx, const int *ny,
                       const double *z0, const double *x0, const double *y0,
                       const double *dz, const double *dx, const double *dy,
                       const int *nsweepIn, const double *epsIn, int *ierr);
void fteik_finalizeF(void);
void fteik_setSourceLocationF(const double *zs, const double *xs,
                              const double *ys, int *ierr);
void fteik_setReceivers64fF(const int nrec,
                            const double *__restrict__ z,
                            const double *__restrict__ x,
                            const double *__restrict__ y,
                            int *ierr);
void fteik_setVelocityModel32fF(const int *nv, const float *__restrict__ vel4,
                                int *ierr);
void fteik_setVelocityModel64fF(const int *nv, const double *__restrict__ vel,
                                int *ierr);
void fteik_solveEikonalLSMF(int *ierr);
void fteik_solveEikonalFSMF(int *ierr);
void fteik_solveEikonalDebugF(int *ierr);
void fteik_getTravelTimes64fF(const int ng, double *__restrict__ tt, int *ierr);
void fteik_getTravelTimesAtReceivers64fF(const int nrecIn,
                                         const double *__restrict__ trec, int *ierr); 
void fteik_getGridSizeF(int *nzOut, int *nxOut, int *nyOut,
                        int *ngrdOut, int *ncellOut, int *ierr);
void fteik_getNumberOfLevelsF(int *nLevelsOut);
void fteik_getGraphPointersF(const int sweep, const int nLevelsIn, const int ngrd4,
                             int *__restrict__ levelPtrOut, int *__restrict__ ijkvOut,
                             int *ierr);
// You can have more fine grained control
void fteik_setNumberOfSweepsF(const int *nsweepIn, int *ierr);
void fteik_setGridSpacingF(const double *dzIn, const double *dxIn,
                           const double *dyIn, int *ierr);
void fteik_setModelOriginF(const double *z0In, const double *x0In,
                           const double *y0In);
void fteik_setSphericalToCartesianEpsilonF(const double *epsIn, int *ierr);
void fteik_getSlownessPermF(const int *sweep, int *__restrict__ perm,
                            int *ierr);
void fteik_getTravelTimePermF(const int *sweep, int *__restrict__ perm,
                              int *ierr);
void fteik_getSweepSigns64fF(const int *sweep,
                             int *sgntz, int *sgntx, int *sgnty,
                             int *sgnvz, int *sgnvx, int *sgnvy,
                             double *sgnrz, double *sgnrx, double *sgnry); 
void fteik_getSweepLimitsF(const int *sweep, const bool *linitk,
                           const int *nz, const int *nx, const int *ny,
                           const int *zsi, const int *xsi, const int *ysi,
                           int *z1, int *z2, int *x1, int *x2, int *y1, int *y2,
                           int *ierr);
// I wouldn't recommend calling these functions unless you've taken particular
// care to initialize initialize the fteik Fortran module and are sure it has
// not gone out of scope.
void fteik_setGridSizeF(const int *nzIn, const int *nxIn, const int *nyIn,
                        int *ierr); 
void fteik_initializeTravelTimesNearSourceF(double *__restrict__ ttimes,
                                            int *ierr);
void fteik_computeGraphF(int *ierr);
int fteik_grid2indexF(const int i, const int j, const int k,
                      const int nz, const int nzx);
void fteik_index2gridF(const int igrd, int *i, int *j, int *k, int *ierr);

#ifdef __cplusplus
}
#endif
#endif
