#ifndef FTEIK_FORTRAN_H__
#define FTEIK_FORTRAN_H__ 1
#include <stdbool.h>
#ifdef __cplusplus
extern "C"
{
#endif

// These are the functions worth knowing
void fteik_initializeF(const int *nz, const int *nx, const int *ny,
                       const double *z0, const double *x0, const double *y0,
                       const double *dz, const double *dx, const double *dy,
                       const int *nsweepIn, const double *epsIn, int *ierr);
void fteik_finalizeF(void);
void fteik_setSourceLocationF(const double *zs, const double *xs,
                              const double *ys, int *ierr);
void fteik_setVelocityModel32fF(const int *nv, const float *__restrict__ vel4,
                                int *ierr);
void fteik_setVelocityModel64fF(const int *nv, const double *__restrict__ vel,
                                int *ierr);
void fteik_solveEikonalLSMF(int *ierr);
void fteik_solveEikonalFSMF(int *ierr);
void fteik_solveEikonalDebugF(int *ierr);
void fteik_getTravelTimes64fF(const int *ng, double *__restrict__ tt, int *ierr);
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
