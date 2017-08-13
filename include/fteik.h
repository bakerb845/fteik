#ifndef FTEIK_H__
#define FTEIK_H__ 1
#include "fteik_config.h"
#include "fteik_struct.h"
#include "fteik_graph.h"
#include "fteik_fortran.h"
#include <stdbool.h>

#ifdef __cplusplus
extern "C"
{
#endif

int fteik_free(struct fteikSolver_struct *solver);
int fteik_legacySolver(struct fteikSolver_struct *solver);
int fteik_initializeGraph(struct fteikSolver_struct *solver);
int fteik_solveEikonalEquation(struct fteikSolver_struct *solver);
int fteik_levelSetSolver(struct fteikSolver_struct *solver);
int fteik_setNumberOfSweeps(const int nsweep, struct fteikSolver_struct *solver);
int fteik_setVelocityModelPointer64f(const int ncell, double *__restrict__ vel, 
                                     struct fteikSolver_struct *solver);
int fteik_dereferenceVelocityModelPointer(struct fteikSolver_struct *solver);
int fteik_setTravelTimesPointer64f(const int nxyz, double *__restrict__ tt,
                                   struct fteikSolver_struct *solver);
int fteik_dereferenceTravelTimePointer(struct fteikSolver_struct *solver);
int fteik_setGridSize(const int nz, const int nx, const int ny,
                      struct fteikSolver_struct *solver);
int fteik_setGridSpacing(const double dz, const double dx, const double dy,
                         struct fteikSolver_struct *solver);
int fteik_setSlownessModel64f(const int ncell, const double *__restrict__ vel,
                              struct fteikSolver_struct *solver);
int fteik_setSourceIndex(const double zs, const double xs, const double ys,
                         struct fteikSolver_struct *solver);
int fteik_setModelOrigin(const double z0, const double x0, const double y0,
                         struct fteikSolver_struct *solver);
int fteik_setSphericalToCartesianEpsilon(const double eps,
                                         struct fteikSolver_struct *solver);
int fteik_performSweeps64f(const int kiter,
                           struct fteikSolver_struct *solver);
int fteik_setInitialUpdateMask(struct fteikSolver_struct *solver);
int fteik_performSweep64f(const int kiter, const int sweep,
                          struct fteikSolver_struct *solver);
int fteik_prefetchTravelTimes(const int sweep,
                              const int i, const int j, const int k,
                              const int sgntz, const int sgntx, const int sgnty,
                              const int nz, const int nx, const int ny,
                              const double *__restrict__ ttimes,
                              const int *__restrict__ ttPerm,
                              double *__restrict__ ttLoc);
int fteik_prefetchSlowness(const int sweep,
                           const int ic, const int jc, const int kc,
                           const int sgnvz, const int sgnvx, const int sgnvy,
                           const int nz, const int nx, const int ny,
                           const double *__restrict__ slow,
                           const int *__restrict__ velPerm,
                           double *__restrict__ slowStencils);

int fteik_applyEvaluateLevelSetSweeps(
    const int n1, const int n2, const int n3);

bool fteik_haveValidSourceIndex(const struct fteikSolver_struct solver);
bool fteik_haveValidGridSpacing(const struct fteikSolver_struct solver);
bool fteik_haveValidGridSize(const struct fteikSolver_struct solver);
bool fteik_haveValidEpsilon(const struct fteikSolver_struct solver);

#ifdef __cplusplus
}
#endif
#endif
