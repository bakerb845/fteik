#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "fteik_fortran.h"

/*!
 * @brief Tests the locator. 
 */
int main()
{
    clock_t t0, t1;
    const double slowConst = 1.0/5.e3;
/*
    const double z0 = 0.0;
    const double y0 = 0.0;
    const double x0 = 0.0;
    const double dz = 1.0;
    const double dx = 1.0;
    const double dy = 1.0;
*/
    const int nz = 101;
    const int nx = 95;
    const int ny = 98;
    const int nsrc = 1;
    const int nrec = 10;
    const int ngrd = nz*nx*ny;
    const double tori[nsrc] = {1.0}; // Origin time
    int zsi[nsrc] = {nz/2};
    int xsi[nsrc] = {nx/2};
    int ysi[nsrc] = {ny/2};
    // Create receiver positions and observed times
    double *ttimes = (double *) calloc((size_t) ngrd, sizeof(double));
    double *tobs = (double *) calloc((size_t) (nrec*nsrc), sizeof(double));
    double *wts = (double *) calloc((size_t) (nrec*nsrc), sizeof(double));
    int *zri = (int *) calloc((size_t) nrec, sizeof(int));
    int *xri = (int *) calloc((size_t) nrec, sizeof(int));
    int *yri = (int *) calloc((size_t) nrec, sizeof(int));
    double diffz, diffx, diffy, t;
    int i, ierr, indx, ix, iy, iz;
    for (i=0; i<nrec; i++)
    {
        zri[i] = 0;
        xri[i] = nrec/4 + i;
        yri[i] = nrec/4 + i;
    }
    int isrc;
    // Initialize the locator
    printf("Initializing locator...\n");
    locate_initializeF(nsrc, nrec, ngrd, &ierr);
    if (ierr != 0)
    {
        printf("Error initializing locator\n");
        return EXIT_FAILURE;
    }
    // Loop on the sources
    printf("Tabulating observed times...\n");
    for (isrc=0; isrc<nsrc; isrc++)
    {
        // Compute the theorietical travel times
        for (i=0; i<nrec; i++)
        {
            diffz = (double) (zri[i] - zsi[isrc]);
            diffx = (double) (xri[i] - xsi[isrc]);
            diffy = (double) (yri[i] - ysi[isrc]);
            t = sqrt(diffz*diffz + diffx*diffx + diffy*diffy)*slowConst;
            tobs[isrc*nrec+i] = t + tori[isrc];
            wts[isrc*nrec+i] = 1.0;
        }
    }
    // Loop on receivers and tabulate reciprocal travel-times
    printf("Tabulating travel times from receivers to all points...\n");
    int irec;
    t0 = clock();
    for (irec=0; irec<nrec; irec++)
    {
        // Compute travel-times in homogeneous velocity model
        for (iy=0; iy<ny; iy++)
        {
            for (ix=0; ix<nx; ix++)
            {
                for (iz=0; iz<nz; iz++) 
                {
                    indx = fteik_model_grid2index(iz, ix, iy, nz, nz*nx);
                    diffz = (double) (iz - zri[irec]);
                    diffx = (double) (ix - xri[irec]);
                    diffy = (double) (iy - yri[irec]);
                    t = sqrt(diffz*diffz + diffx*diffx + diffy*diffy)*slowConst;
                    ttimes[indx] = t; //ttimes[irec*ngrd+indx] = t;
                }
            }
        }
        // Set the travel-times
        locate_setTravelTimeField64fF(ngrd, irec+1, ttimes, &ierr);
        if (ierr != 0)
        {
            printf("Failed to set travel time field\n");
            return EXIT_FAILURE;
        }
    } // Loop on sources
    t1 = clock();
    printf("Travel-time computation time %lf (s)\n", ((double)(t1-t0))/CLOCKS_PER_SEC);
    // Set the observed travel times
    printf("Setting observed travel times...\n");
    t0 = clock();

    t1 = clock();
    printf("Time to set observations %lf (s)\n", ((double)(t1-t0))/CLOCKS_PER_SEC);
    // Free memory
    printf("Finalizing solver...\n");
    locate_finalizeF();
    printf("Freeing memory...\n");
    free(ttimes);
    free(wts);
    free(tobs);
    free(zri);
    free(xri);
    free(yri);
    return EXIT_SUCCESS;
}
