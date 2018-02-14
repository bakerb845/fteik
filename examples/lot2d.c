#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fteik/fteik.h"

void solver2d_c(const double *slow,
                double *tt,
                const int nz, const int nx,
                const double zsrc, const double xsrc,
                const double dz, const double dx,
                const int nsrc);

int create2DLayeredVelocityModel(const int nlayer, const int nz, const int nx, 
                                 double dz, 
                                 const double interface[], const double v1d[],
                                 double **vmod);
/*!
 * @brief 2D layer over halfspace problem.
 */
int main( )
{
    int nlayer = 3;                     // Number of layers in model.
    double dx = 10.0;                   // Grid spacing (m) in x.
    double dz = 10.0;                   // Grid spacing (m) in z.
    double x0 = 0.0;                    // x model origin (m).
    double x1 = 2000.0;                 // x model extent (m).
    double z0 = 0.0;                    // z model origin (m); e.g., free surface
    double z1 = 800.0;                  // z model extent (m).
    double v1d[3] = {3500.0, 4200.0, 5500.0};   // Velocity (m/s) in layers; z+ down
    double interface[2] = {100.0, 250.0};      // Interface depth (m); z+ down
    // Source info
    int nsrc = 1;                       // Number of sources
    double zsrc[1] = {z0+(z1-z0)*0.5 + 3.0}; // z source depth (m)
    double xsrc[1] = {x0+(x1-x0)*0.5 + 3.0};  // x source offset (m).
    // Solver parameters
    int nsweep = 2;                     // Number of Gauss-Seidel iterations.
    double eps = 5.0;                   // Number of grid points around source to
                                        // use spherical approximation.
    double convTol = 1.e-10;            // Convergence tolerance for Gauss-Seidel
    // Other parameters
    const char *archiveFile = "loh.h5";
    const char *modelName = "vpLayerOverHalfSpaceModel";
    const char *travelTimeName = "vpLayerOverHalfSpaceTravelTimes";
    double *vel, *ttimes, xwidth, zwidth;
    int ierr, ncell, nx, nz;
   
    xwidth = fabs(x1 - x0);
    zwidth = fabs(z1 - z0);
    nx = (int) (round(xwidth/dx)) + 1;
    nz = (int) (round(zwidth/dz)) + 1;
    // Some checks 
    if (dx <= 0.0 || dz <= 0.0)
    {
        fprintf(stderr, "%s: Invalid grid spacing\n", __func__);
        return EXIT_FAILURE;
    }
printf("%d %d\n", nx, nz);
    if (nx < 3 || nz < 3)
    {
        fprintf(stderr, "%s: Insufficient number of grid points: %d %d\n",
                __func__, nx, nz);
        return EXIT_FAILURE;
    }
    // Create a 1D velocity model in the solver coordinate system
    ierr = create2DLayeredVelocityModel(nlayer, nz, nx, dz,
                                        interface, v1d, &vel);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Failed to make velocity model\n", __func__);
        return EXIT_FAILURE;
    }
    // Initialize the solver
    fprintf(stdout, "%s: Initializing solver...\n", __func__);
    fteik_solver2d_initialize64f(nz, nx, 
                                 z0, x0,  
                                 dz, dx, 
                                 nsweep, eps,
                                 convTol, 0,
                                 &ierr);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error initializing solver\n", __func__);
        return EXIT_FAILURE;
    }
    // Set the velocity model
    fprintf(stdout, "%s: Setting the velocity model...\n", __func__);
    ncell = (nz - 1)*(nx - 1);
    fteik_solver2d_setVelocityModel64f(ncell, vel, &ierr);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error setting velocity model\n", __func__);
        return EXIT_FAILURE;
    } 
    fprintf(stdout, "%s: Setting the source locations...\n", __func__);
    fteik_solver2d_setSources64f(nsrc, zsrc, xsrc, &ierr); 
ttimes = (double *) calloc((size_t) (nz*nx), sizeof(double));
double *slow = (double *) calloc((size_t) ((nz-1)*(nx-1)), sizeof(double));
for (int k=0; k<(nz-1)*(nx-1); k++){slow[k] = 1.0/vel[k];}
printf("calling ref\n");
solver2d_c(slow, ttimes, nz, nx, zsrc[0], xsrc[0], dz, dx, nsweep); 
    fprintf(stdout, "%s: Calling debug solver...\n", __func__);
fteik_solver2d_solveSourceFSM(1, &ierr);
printf("calling lsm\n"); 
    fprintf(stdout, "%s: Calling solver...\n", __func__);
    fteik_solver2d_solveSourceLSM(1, &ierr);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error calling solver\n", __func__);
        return EXIT_FAILURE;
    }
    // Initialize the output
    fprintf(stdout, "%s: Writing results...\n", __func__);
    fteik_h5io_initializeF(archiveFile);
    fteik_h5io_writeVelocityModelF(modelName);
    fteik_h5io_writeTravelTimesF(travelTimeName); 
    fteik_h5io_finalizeF();

    fprintf(stdout, "%s: Finalizing solver...\n", __func__);
    fteik_solver2d_free();
    // Free resources
    free(vel);
    return EXIT_SUCCESS;
}
//============================================================================//

int create2DLayeredVelocityModel(const int nlayer, const int nz, const int nx,
                                 double dz,
                                 const double interface[], const double v1d[],
                                 double **vmod)
{
    double *vel = NULL;
    double v, z;
    int ix, iz, k, ncellx, ncellz;
    *vmod = NULL;
    if (nlayer < 1)
    {
        fprintf(stderr, "%s: Error no layers\n", __func__);
        return EXIT_FAILURE;
    }
    ncellz = nz - 1;
    ncellx = nx - 1;
    vel = (double *) calloc((size_t) (ncellz*ncellx), sizeof(double));
    for (iz=0; iz<ncellz; iz++)
    {
        z = (double) iz*dz;
        for (k=0; k<nlayer; k++)
        {
            ix*ncellz + iz;
        }
        v = v1d[0];
        for (k=0; k<nlayer-1; k++)
        {
            if (z >= interface[k] - 1.e-10){v = v1d[k+1];}
        }
printf("%d %f\n", iz, v);
        //printf("%f %f\n", z, v);
        // Copy velocity for all x locations
        #pragma omp simd
        for (ix=0; ix<ncellx; ix++)
        {
            vel[ix*ncellz + iz] = v; //ix*ncellz + iz + 2; //v;
        }
    }
    // Ensure that I got all points 
    for (k=0; k<ncellx*ncellz; k++)
    {
        if (vel[k] <= 0.0)
        {
            fprintf(stderr, "%s: Index %d has invalid velocity\n", __func__, k);
            free(vel);
            return EXIT_FAILURE;
        }
    }
    *vmod = vel;
    return EXIT_SUCCESS;
}
