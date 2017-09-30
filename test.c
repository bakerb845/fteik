#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <omp.h>
#include "fteik.h"
#include "fteik_analytic.h"
#include "fteik_h5io.h"
#include "fteik_graph.h"
#include "fteik_fortran.h"


void fteikTestGraph(const int nz, const int nx, const int ny);

void fteik_(const double *__restrict__ vel,
            double *__restrict__ tt,
            const int *nz, const int *nx, const int *ny,
            const double *zsin, const double *xsin, const double *ysin,
            const double *dzin, const double *dxin, const double *dyin,
            const int *nsweep, const double *epsin);
int graph_testGrd2ijk(const int n1, const int n2, const int n3);
int graph_createGraph3D(const int n1, const int n2, const int n3,
                        int **xadj, int **adjncy);
int graph_createLevelStructure(const int n, const int *__restrict__ xadj,
                               const int *__restrict__ adjncy,
                               int *nLevels, int **levelPtrOut,
                               int **ijkvOut, int **n2lOut);
int fteik_io_readVelocityModel(const char *fileName,
                               int *ncellz, int *ncellx, int *ncelly,
                               double *dz, double *dx, double *dy,
                               double *z0, double *x0, double *y0,
                               double *__restrict__ vel[]);
int fteik_xdmf_writeVelocityModel(const int job, 
                                  const char *h5flName,
                                  const char *velModel,
                                  const bool lflip,
                                  const int nz, const int nx, const int ny, 
                                  const double dz, const double dx, 
                                  const double dy, 
                                  const double z0, const double x0, 
                                  const double y0);
/*!
 */
int main()
{
    const char *projnm = "testing";
    struct levelSets_struct levelSets;
    struct xdmf_struct xdmf;
    struct fteikSolver_struct solver;
    double *tt, *ttRef, *vel, xrms;
    clock_t t0, t1;
    int *xadj, *adjncy, *levelPtr, *ijkv, *n2l, *n2l8;
    int ierr, ncell, nLevels, ngrd, nrec;
    const double eps = 2.0;
    const int nsweep = 5;
    int nz = 13, nx = 10, ny = 15;
nz = 23; nx = 20; ny = 25;
//nz = 35; nx = 76; ny = 76;
//nz = 55; nx = 76; ny = 76;
    const double x0 = 0.0, y0 = 0.0, z0 = 0.0;
    const double dx = 1.1, dy = 1.2, dz = 1.3;
    const double xs = 0.5*(double) (nx - 1)*dx - 0.1;
    const double ys = 0.5*(double) (ny - 1)*dy + 0.1;
    const double zs = 0.5*(double) (nz - 1)*dz + 0.2;
    const double vconst = 6.e3;
    hid_t fileID;

int ncellx, ncelly, ncellz;
double dxt, dyt, dzt, x0t, y0t, z0t;
double *vwork = NULL;
 fteik_io_readVelocityModel("/home/bakerb25/C/SEG_EAGE_3DOverhrustModel/Disk1/3D-Velocity-Grid/segOverthrust.h5",
                            &ncellz, &ncellx, &ncelly,
                            &dxt, &dyt, &dzt, &z0t, &x0t, &y0t,
                            &vwork);
int nzFact = (ncellz - ncellz%(nz - 2))/(nz - 2); //floor((double) ncellz/(double) nz + 0.5);
int nxFact = (ncellx - ncellx%(nx - 2))/(nx - 2); //floor((double) ncellx/(double) nx + 0.5);
int nyFact = (ncelly - ncelly%(ny - 2))/(ny - 2); //floor((double) ncelly/(double) ny + 0.5);
int ix, iy, iz;
int jz = 0;
int jx = 0;
int jy = 0;
    ngrd = nx*ny*nz;
    ncell = (nx - 1)*(ny - 1)*(nz - 1); 
    vel = (double *) aligned_alloc(64, (size_t) (ncell)*sizeof(double));
    ttRef = (double *) aligned_alloc(64, (size_t) (nx*ny*nz)*sizeof(double));
    tt = (double *) aligned_alloc(64, (size_t) (nx*ny*nz)*sizeof(double));
for (iy=0; iy<ny-1; iy++)
{
  for (ix=0; ix<nx-1; ix++)
  {
     for (iz=0; iz<nz-1; iz++)
     {
        jz = MAX(0, MIN(ncellz-1, iz*nzFact));
        jx = MAX(0, MIN(ncellx-1, ix*nxFact));
        jy = MAX(0, MIN(ncelly-1, iy*nyFact));
        vel[iy*(nz-1)*(nx-1)     + ix*(nz-1)   + iz] 
     =(double) vwork[jy*(ncellz)*(ncellx) + jx*(ncellz) + jz];
//int loc = iy*(nz-1)*(nx-1)         + ix*(nz-1)     + iz;
//printf("%d %d %d %d %d %d %d %f\n", iz, ix, iy, jz, jx, jy, loc, vel[loc]);
     }
  }
}
free(vwork);
//printf("Press (ENTER)\n");
//getchar();
//omp_set_num_threads(1);
fteik_xdmf_writeVelocityModel(1, "debug.h5", "debugModel", false,
                              nz, nx, ny, dz, dx, dy, z0, x0, y0);
fteik_xdmf_writeVelocityModel(3, "debug.h5", "debugTimes", false,
                              nz, nx, ny, dz, dx, dy, z0, x0, y0);
fteik_xdmf_writeVelocityModel(4, "debug.h5", "debugModel", false,
                              nz, nx, ny, dz, dx, dy, z0, x0, y0);
fteik_xdmf_writeVelocityModel(5, "debug.h5", "debugModel", false,
                              nz, nx, ny, dz, dx, dy, z0, x0, y0);

    memset(&xdmf, 0, sizeof(struct xdmf_struct));
    memset(&solver, 0, sizeof(struct fteikSolver_struct));
omp_set_num_threads(1);
//graph_testGrd2ijk(34, 10, 21);
    // Compute an analytic solution in a constant velocity model
/*
    ngrd = nx*ny*nz;
    ncell = (nx - 1)*(ny - 1)*(nz - 1); 
    vel = (double *) aligned_alloc(64, (size_t) (ncell)*sizeof(double));
    ttRef = (double *) aligned_alloc(64, (size_t) (nx*ny*nz)*sizeof(double));
    tt = (double *) aligned_alloc(64, (size_t) (nx*ny*nz)*sizeof(double));
*/
/*
    analyticSolution_wholeSpace(nz, nx ,ny, dz, dx, dy, vconst, vel);
    analyticSolution_wholeSpaceAnalyticSolution(nz, nx, ny, 
                                                dz, dx, dy, 
                                                zs, xs, ys, 
                                                vconst, ttRef);
*/
    // Set some receiver locations at the free surface
    nrec = MIN3(nx, ny, nz)/4;
    printf("Number of receivers: %d\n", nrec);
    double *xrec = (double *) calloc((size_t) nrec, sizeof(double));
    double *yrec = (double *) calloc((size_t) nrec, sizeof(double));
    double *zrec = (double *) calloc((size_t) nrec, sizeof(double));
    double *trec = (double *) calloc((size_t) nrec, sizeof(double));
    for (int ir=0; ir<nrec; ir++)
    {
        xrec[ir] = x0 + 3.0/4.0*ir;
        yrec[ir] = y0 + 3.0/4.0*ir;
        zrec[ir] = z0;
    }
    // Initialize the solver
    printf("Initializing solver...\n");
//printf("SOMETHING IS AMISS WITH THE LOCAL SOLVERS!!!!!!!!!!!\n");
//printf("INSPECTOR SAYS slow(?), i.e., prefetchSlowness ARE UNINITIALIZED - CHECK THERE FIRST\n");
//getchar();
    fteik_solver_initialize64fF(nz, nx, ny,
                                z0, x0, y0,
                                dz, dx, dy,
                                nsweep, eps, &ierr);
    // Set the sources
    printf("Setting sources...\n");
    fteik_solver_setSources64fF(1, &zs, &xs, &ys, &ierr);
    // Set the velocity model
    printf("Setting velocity model...\n");
    fteik_solver_setVelocityModel64fF(ncell, vel, &ierr);
    // Solve the Eikonal equation
    printf("Solving the Eikonal equation...\n");
    fteik_solver_solveSourceLSMF(1, &ierr);
    printf("Finalizing the solver\n");
//    fteik_solver_finalizeF();
//return 0;

    // Intialize the solver
    printf("Initializing solver...\n");
    fteik_initializeF(&nz, &nx, &ny, 
                      &z0, &x0, &y0,
                      &dz, &dx, &dy,
                      &nsweep, &eps, &ierr); 
    printf("Setting receiver locations...\n");
    fteik_setReceivers64fF(nrec, zrec, xrec, yrec, &ierr);
    printf("Setting velocity model...\n");
    fteik_setVelocityModel64fF(&ncell, vel, &ierr);
    printf("Setting source...\n");
    fteik_setSourceLocationF(&zs, &xs, &ys, &ierr);
    printf("Using debugging solver...\n");
    //fteik_solveEikonalDebugF(&ierr);
    t0 = clock();
    fteik_solveEikonalDebugF(&ierr);
    t1 = clock();
    printf("Debug solver time %f (s)\n", ((float)(t1-t0))/CLOCKS_PER_SEC);
    t0 = clock();
    fteik_solveEikonalFSMF(&ierr);
    t1 = clock();
    printf("FSM solver time %f (s)\n", ((float)(t1-t0))/CLOCKS_PER_SEC);
    t0 = clock();
    fteik_solveEikonalLSMF(&ierr);
    t1 = clock();
    printf("LSM solver time %f (s)\n", ((float)(t1-t0))/CLOCKS_PER_SEC);
    t0 = clock();
    fteik_solveEikonalLSMF(&ierr);
    t1 = clock();
    printf("LSM solver time %f (s)\n", ((float)(t1-t0))/CLOCKS_PER_SEC);
 
    printf("Copying travel-times...\n");
    fteik_getTravelTimes64fF(ngrd, tt, &ierr);
    printf("Getting traveltimes at stations...\n");
    fteik_getTravelTimesAtReceivers64fF(nrec, trec, &ierr); 

    printf("trying some output\n");
    fteik_h5io_initializeF("./debug.h5");
    fteik_h5io_writeLevelSchedulesF();
    fteik_h5io_writeVelocityModelF("debugModel");
    fteik_h5io_writeTravelTimesF("debugTimes");

    printf("Finalizing solver\n");
    fteik_finalizeF();

    printf("Finalize the solver\n");
    fteik_solver_finalizeF();


float zs4 = zs;
float xs4 = xs;
float ys4 = ys;
float dz4 = dz;
float dx4 = dx;
float dy4 = dy;
float eps4 = eps;
double *tori = (double *) calloc((size_t) (nx*ny*nz), sizeof(double));
t0=clock();
int ldo = 0;
omp_set_num_threads(1);
if (ldo == 0)
{
fteik_(vel, tori,
       &nz, &nx, &ny,
       &zs, &xs, &ys,
       &dz, &dx, &dy, &nsweep, &eps);
//return 0;
}
t1=clock() - t0;
printf("solver time %f (s)\n", ((float)t1)/CLOCKS_PER_SEC);
xrms = 0.0;
for (int i=0; i<ngrd; i++)
{
 xrms = xrms + pow(tt[i] - tori[i], 2);
}
printf("xrms: %e\n", sqrt(xrms)/(double) (ngrd));


free(tori);
free(tt);
free(xrec);
free(yrec);
free(zrec);
free(trec);
free(vel);
free(ttRef);
return 0;
}

int analyticSolution_wholeSpace(
    const int nz, const int nx, const int ny,
    const double dz, const double dx, const double dy,
    const double vconst, double *__restrict__ vel)
{
    int i, ielem, j, k;
    for (k=0; k<ny-1; k++)
    {
        for (j=0; j<nx-1; j++)
        {
            for (i=0; i<nz-1; i++)
            {
                ielem = k*(nx - 1)*(nz - 1) + j*(nz - 1) + i;
                vel[ielem] = vconst; //(double) (ielem+1); //vconst;
            }
        }
    }
    return 0;
}

int analyticSolution_wholeSpaceAnalyticSolution(
    const int nz, const int nx, const int ny,
    const double dz, const double dx, const double dy,
    const double zs, const double xs, const double ys,
    const double vconst, double *__restrict__ ttimes)
{
    double dist, delx, dely, delz, slowness, x, y, z;
    int i, indx, j, k;
    slowness = 1.0/vconst;
    for (k=0; k<ny; k++)
    {
        for (j=0; j<nx; j++)
        {
            for (i=0; i<nz; i++)
            {
                indx = k*nx*nz + j*nz + i;
                z = (double) i*dz; 
                x = (double) j*dx;
                y = (double) k*dy;
                delz = z - zs;
                delx = x - xs;
                dely = y - ys;
                dist = sqrt(delx*delx + dely*dely + delz*delz);
                ttimes[indx] = dist*slowness;
            }
        }
    }
    return 0;
}
