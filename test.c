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
            const float *zsin, const float *xsin, const float *ysin,
            const float *dzin, const float *dxin, const float *dyin,
            const int *nsweep, const float *epsin);
int graph_testGrd2ijk(const int n1, const int n2, const int n3);
int graph_createGraph3D(const int n1, const int n2, const int n3,
                        int **xadj, int **adjncy);
int graph_createLevelStructure(const int n, const int *__restrict__ xadj,
                               const int *__restrict__ adjncy,
                               int *nLevels, int **levelPtrOut,
                               int **ijkvOut, int **n2lOut);
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
    int ierr, ncell, nLevels, ngrd;
    const double eps = 2.0;
    const int nsweep = 5;
    int nz = 13, nx = 10, ny = 15;
nz = 35; nx = 36; ny = 37;
    const double x0 = 0.0, y0 = 0.0, z0 = 0.0;
    const double dx = 1.1, dy = 1.2, dz = 1.3;
    const double xs = 0.5*(double) (nx - 1)*dx - 0.1;
    const double ys = 0.5*(double) (ny - 1)*dy + 0.1;
    const double zs = 0.5*(double) (nz - 1)*dz + 0.2;
    const double vconst = 6.e3;
    hid_t fileID;

    memset(&xdmf, 0, sizeof(struct xdmf_struct));
    memset(&solver, 0, sizeof(struct fteikSolver_struct));
omp_set_num_threads(2);
//graph_testGrd2ijk(34, 10, 21);
    // Compute an analytic solution in a constant velocity model
    ngrd = nx*ny*nz;
    ncell = (nx - 1)*(ny - 1)*(nz - 1); 
    vel = (double *) aligned_alloc(64, (size_t) (ncell)*sizeof(double));
    ttRef = (double *) aligned_alloc(64, (size_t) (nx*ny*nz)*sizeof(double));
    tt = (double *) aligned_alloc(64, (size_t) (nx*ny*nz)*sizeof(double));
    analyticSolution_wholeSpace(nz, nx ,ny, dz, dx, dy, vconst, vel);
    analyticSolution_wholeSpaceAnalyticSolution(nz, nx, ny, 
                                                dz, dx, dy, 
                                                zs, xs, ys, 
                                                vconst, ttRef);
    // Intialize the solver
    printf("Initializing solver...\n");
    fteik_initializeF(&nz, &nx, &ny, 
                      &z0, &x0, &y0,
                      &dz, &dx, &dy,
                      &nsweep, &eps, &ierr); 
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

    fteik_getTravelTimes64fF(&ngrd, tt, &ierr);
    printf("Finalizing solver\n");
    fteik_finalizeF();
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
       &zs4, &xs4, &ys4,
       &dz4, &dx4, &dy4, &nsweep, &eps4);
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
return 0;
getchar();
    // Create the level set structure for the 8 different sweeps 
    ierr = graph_createLevelStructs(nz, nx, ny, &levelSets);
    if (ierr != 0)
    {
        printf("Error creating level set structures\n");
        return EXIT_FAILURE;
    }
    // Create the output H5 file
    fileID = H5Fcreate("levels.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    // Dump the model geometry
    fteik_h5io_writeGeometry(fileID, 
                             nz, nx, ny,
                             dz, dx, dy,
                             z0, x0, y0);
    // Initialize the XDMF file for plotting and dump the levels
    fteik_xdmf_initialize("./\0", "levels.h5", projnm, nz, nx, ny, &xdmf);
    for (int k=0; k<8; k++)
    {
        char name[128];
        memset(name, 0, 128*sizeof(char));
        sprintf(name, "levelSet_%d", k+1);
        fteik_h5io_writeLevelSet32i(fileID, name,
                                    nz, nx, ny,
                                    levelSets.levelSet[k].n2l);
        fteik_xdmf_addLevelSet(xdmf, k+1, name);
    }
/*
    // Compute an analytic solution in a constant velocity model
    ncell = (nx - 1)*(ny - 1)*(nz - 1);
    vel = (double *) aligned_alloc(64, (size_t) (ncell)*sizeof(double));
    ttRef = (double *) aligned_alloc(64, (size_t) (nx*ny*nz)*sizeof(double));
    analyticSolution_wholeSpace(nz, nx ,ny, dz, dx, dy, vconst, vel);
    analyticSolution_wholeSpaceAnalyticSolution(nz, nx, ny,
                                                dz, dx, dy,
                                                zs, xs, ys,
                                                vconst, ttRef);
*/
    // Write the analytic solutions to the H5 archive
    fteik_h5io_writeVelocityModel64f(fileID, "HomogeneousModel\0",
                                     nz, nx, ny, vel);
    fteik_h5io_writeTravelTimes64f(fileID, "AnalyticSolution\0",
                                   nz, nx ,ny, ttRef);
    fteik_xdmf_addVelocityModel(xdmf, "HomogeneousModel\0");
    fteik_xdmf_addTravelTimes(xdmf, "AnalyticSolution\0");


    // Set the constants for fteik
    fteik_setGridSize(nz, nx, ny, &solver);
    fteik_setGridSpacing(dz, dx, dy, &solver);
    fteik_setSphericalToCartesianEpsilon(eps, &solver);
    fteik_setNumberOfSweeps(nsweep, &solver);
    // Set the level-set elimination tree
    fteik_initializeGraph(&solver); 
    // Set the velocity model
    fteik_setSlownessModel64f(ncell, vel, &solver);
    // Set the source
    fteik_setSourceIndex(zs, xs, ys, &solver);
    // solve
printf("fix here\n"); solver.nsweeps = 0;
/*
t0=clock();
fteik_levelSetSolver(&solver);
t1=clock() - t0;
*/
printf("levelset %f (s)\n", ((float)t1/CLOCKS_PER_SEC));
fteik_h5io_writeTravelTimes64f(fileID, "LevelSetSolution\0",
                                   nz, nx ,ny, solver.tt);
double *temp = (double *) calloc((size_t) (nx*ny*nz), sizeof(double));
memcpy(temp, solver.tt, (size_t) (nx*ny*nz)*sizeof(double));
fteik_xdmf_addTravelTimes(xdmf, "LevelSetSolution\0");
//getchar();

printf("solving\n");
t0=clock();
ldo = 0;
if (ldo == 0)
{
fteik_(vel, tori, 
       &nz, &nx, &ny,
       &zs4, &xs4, &ys4,
       &dz4, &dx4, &dy4, &nsweep, &eps4);
//return 0;
}
t1=clock() - t0;
printf("solver time %f (s)\n", ((float)t1)/CLOCKS_PER_SEC);
    fteik_h5io_writeTravelTimes64f(fileID, "FortranSolution\0",
                                   nz, nx ,ny, tori);
    fteik_xdmf_addTravelTimes(xdmf, "FortranSolution\0");
printf("solver 2\n");
t0 = clock();
    //fteik_solveEikonalEquation(&solver);
 fteik_legacySolver(&solver);
t1 = clock() - t0;
printf("solver time %f (s)\n", ((float)t1)/CLOCKS_PER_SEC);

fteik_h5io_writeTravelTimes64f(fileID, "cSolution\0",
                               nz, nx ,ny, solver.tt);
fteik_xdmf_addTravelTimes(xdmf, "cSolution\0");

xrms = 0.0;
double xrmsNew = 0.0;
for (int i=0; i<solver.ngrd; i++)
{
 xrms = xrms + pow(solver.tt[i] - tori[i], 2);
 xrmsNew = xrmsNew + pow(solver.tt[i] - temp[i], 2);
}
printf("%e\n", sqrt(xrms)/(double) solver.ngrd);
printf("debug rms %e\n", sqrt(xrmsNew)/(double) solver.ngrd);

    // Close the archive file and finalize the XDMF file
    H5Fclose(fileID);
    fteik_xdmf_finalize(&xdmf);
    // Free memory
    graph_freeLevelsStruct(&levelSets);
    free(ttRef);
    free(vel);
    //fteik_applyEvaluateLevelSetSweeps(3, 5, 5);
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
                vel[ielem] = (double) (ielem+1); //vconst;
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
