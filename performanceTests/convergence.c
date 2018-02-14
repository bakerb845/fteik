#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "fteik/fteik_fortran.h"

int converge2D(const int job);
static void computeNorms(const int n,
                         const double x[],
                         const double y[],
                         double *l1, double *l2, double *l8);

int main( )
{
    converge2D(1);
    converge2D(2);
    return EXIT_SUCCESS;
}

int converge2D(const int job)
{
    const double x0 = 0.0;
    const double y0 = 0.0;
    const double z0 = 0.0;
    const double x1 = 1.e3;
    const double z1 = 1.e3;
    double dx0 = 100.0;
    double dy0 = 1.0;
    double dz0 = 100.0;
    const double vconst = 1.e3;
    const double vmin = 1.e3;
    const double vmax = 2.e3;
    double xsrc[1], ysrc[1], zsrc[1]; 
    double l1, l2, l8;
    double *ttimes, *ttref, *vel, vint;
    const double convTol = 0.0;
    const int verbose = 0;
    const int nsweep = 2;
    const double eps = 3.0;
    const int isrc = 1;
    int i, ierr, indx, ix, iz, ncell, ngrd, nx, nz;
    const int nsrc = 1;
    clock_t start, diff;
    int n = 8;
    for (i=0; i<n; i++)
    {
        nx = (int) ((x1 - x0)/dx0 + 0.5) + 1;
        nz = (int) ((z1 - z0)/dz0 + 0.5) + 1;
        // Set the source
        xsrc[0] = (x0 + (double) (nx - 1)*dx0)/2.0;
        ysrc[0] = 0.0;
        zsrc[0] = (z0 + (double) (nz - 1)*dz0)/2.0;
        // Allocate space
        ngrd = nx*nz;
        ttimes = (double *) calloc((size_t) ngrd, sizeof(double));
        ttref  = (double *) calloc((size_t) ngrd, sizeof(double));
        // Build a velocity model
        ncell = (nx - 1)*(nz - 1);
        vel = (double *) calloc((size_t) ncell, sizeof(double));
        if (job == 1)
        {
            for (iz=0; iz<nz-1; iz++)
            {
                indx = iz*(nx - 1);
                for (ix=0; ix<nx-1; ix++){vel[indx+ix] = vconst;}
            }
        }
        // Linear gradient
        else
        {
            for (iz=0; iz<nz-1; iz++)
            {
                vint = vmin + (vmax - vmin)/(z1 - z0)*(iz*dz0 + z0);
                indx = iz*(nx - 1);
                for (ix=0; ix<nx-1; ix++){vel[indx+ix] = vint;}
            }
        }
        // Run the eikonal solver
        fteik_solver2d_initialize64f(nz, nx, z0, x0, dz0, dx0, 
                                     nsweep, eps, convTol, verbose, &ierr); 
        if (ierr != 0){return EXIT_FAILURE;}
        fteik_solver2d_setSources64f(nsrc, zsrc, xsrc, &ierr);
        if (ierr != 0){return EXIT_FAILURE;}
        fteik_solver2d_setCellVelocityModel64f(ncell, FTEIK_XZ_ORDERING,
                                               vel, &ierr);
        if (ierr != 0){return EXIT_FAILURE;}
        start = clock();
        fteik_solver2d_solveSourceLSM(isrc, &ierr);
        diff = clock() - start;
        if (ierr != 0){return EXIT_FAILURE;}
        fteik_solver2d_getTravelTimeField64f(ngrd, FTEIK_ZX_ORDERING,
                                             ttimes, &ierr);
        if (ierr != 0){return EXIT_FAILURE;}
        fteik_solver2d_free();
        // Run the analytic solver
        fteik_analytic_initialize64f(nz, nx, 1,
                                     z0, x0, y0,
                                     dz0, dx0, dy0,
                                     verbose, &ierr);
        if (ierr != 0){return EXIT_FAILURE;} 
        fteik_analytic_setSources64f(nsrc, zsrc, xsrc, ysrc, &ierr);
        if (ierr != 0){return EXIT_FAILURE;}
        fteik_analytic_setConstantVelocity64f(vconst, &ierr);
        fteik_analytic_setLinearVelocityGradient64f(vmin, vmax, &ierr);
        if (job == 1)
        {
            fteik_analytic_solveSourceConstantVelocity64f(isrc, &ierr);
        }
        else
        {
            fteik_analytic_solveSourceLinearVelocityGradient64f(isrc, &ierr);
        }
        if (ierr != 0){return EXIT_FAILURE;}
        fteik_analytic_getTravelTimeField64f(ngrd, FTEIK_ZYX_ORDERING,
                                             ttref, &ierr);
        if (ierr != 0){return EXIT_FAILURE;}
        fteik_analytic_free();
        // Compute the norms
        computeNorms(ngrd, ttimes, ttref, &l1, &l2, &l8);
        l1 = l1/(double) ngrd; // Compare apples and apples
        l2 = l2/(double) ngrd; // Compare apples and apples
        l8 = l8;               // This doesn't need normalization
printf("%8d %e %e %e %e %e %e\n", ngrd, dz0, dx0, (double) diff/CLOCKS_PER_SEC, l1, l2, l8);
        // Clean up the workspace
        free(ttimes);
        free(ttref);
        free(vel);
        // Update my grid spacings
        dx0 = dx0/2.0;
        dz0 = dz0/2.0;
    }
    return EXIT_SUCCESS;
}

static void computeNorms(const int n,
                         const double x[],
                         const double y[],
                         double *l1, double *l2, double *l8)
{
    double l1norm, l2norm, l8norm, res;
    int i;
    l1norm = 0.0;
    l2norm = 0.0;
    l8norm = 0.0;
//    #pragma omp simd
    for (i=0; i<n; i++)
    {
        res = x[i] - y[i];
        l1norm = l1norm + fabs(res);
        l2norm = l2norm + res*res;
        l8norm = fmax(l8norm, res);
    }
    *l1 = l1norm;
    *l2 = sqrt(l2norm);
    *l8 = l8norm;
    return; 
}
