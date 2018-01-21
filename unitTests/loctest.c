#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "fteik/fteik_fortran.h"

/*!
 * @brief Tests the locator. 
 */
int main(int argc, char *argv[])
{
    const int master = 0;
    int myid = 0;
#ifdef FTEIK_USE_MPI
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
#endif
    clock_t t0, t1;
    const double velConst = 5.e3;
    const double slowConst = 1.0/velConst;
    const double slowConstS = 1.0/(velConst/sqrt(3.0));
/*
    const double z0 = 0.0;
    const double y0 = 0.0;
    const double x0 = 0.0;
*/
    const double dz = 100.0;
    const double dx = 100.0;
    const double dy = 100.0;

    const int nz = 101;
    const int nx = 95;
    const int ny = 98;
    const int nsrc = 1;
    const int nrec = 10;
    const int ngrd = nz*nx*ny;
    const double tori[1] = {+1.5}; // Origin time
    int zsi[1] = {nz/2};
    int xsi[1] = {nx/2};
    int ysi[1] = {ny/2};
    // Create receiver positions and observed times
    double *ttimes = (double *) calloc((size_t) ngrd, sizeof(double));
    double *tobs = (double *) calloc((size_t) (nrec*nsrc), sizeof(double));
    double *wts = (double *) calloc((size_t) (nrec*nsrc), sizeof(double));
    int *obs2tf = (int *) calloc((size_t) (nrec*nsrc), sizeof(int));
    int *zri = (int *) calloc((size_t) nrec, sizeof(int));
    int *xri = (int *) calloc((size_t) nrec, sizeof(int));
    int *yri = (int *) calloc((size_t) nrec, sizeof(int));
    double diffz, diffx, diffy, t, slowUse;
    int i, ierr, indx, ix, iy, iz;
    for (i=0; i<nrec-1; i++)
    {
        zri[i] = 0;
        xri[i] = nrec/4 + i;
        yri[i] = nrec/4 + i;
    }
    // Set the final receiver location twice
    zri[nrec-1] = 0;
    zri[nrec-1] = nrec/4 + nrec - 2;
    zri[nrec-1] = nrec/4 + nrec - 2;
    int isrc;
    // Initialize the locator
    if (myid == master){printf("Initializing locator...\n");}
#ifdef FTEIK_USE_MPI
    locate_initializeMPI(master, MPI_COMM_WORLD, nsrc, nrec, ngrd, &ierr); 
#else
    locate_initializeF(nsrc, nrec, ngrd, &ierr);
#endif
    if (ierr != 0)
    {
        printf("Error initializing locator\n");
        return EXIT_FAILURE;
    }
    // Loop on the sources
    printf("Tabulating observed times...\n");
    for (isrc=0; isrc<nsrc; isrc++)
    {
        // Compute the theoretical travel times
        for (i=0; i<nrec; i++)
        {
            diffz = (double) (zri[i] - zsi[isrc])*dz;
            diffx = (double) (xri[i] - xsi[isrc])*dx;
            diffy = (double) (yri[i] - ysi[isrc])*dy;
            slowUse = slowConst;
            if (i == nrec - 1){slowUse = slowConstS;}
            t = sqrt(diffz*diffz + diffx*diffx + diffy*diffy)*slowUse;
            tobs[isrc*nrec+i] = t + tori[isrc];
            obs2tf[isrc*nrec+i] = i;
            wts[isrc*nrec+i] = 1.0;
        }
    }
    // Loop on receivers and tabulate reciprocal travel times
    if (myid == master)
    {
        printf("Tabulating travel times from receivers to all points...\n");
    }
    int irec;
    t0 = clock();
    for (irec=0; irec<nrec; irec++)
    {
        slowUse = slowConst;
        if (irec == nrec - 1){slowUse = slowConstS;}
        // Compute travel-times in homogeneous velocity model
        for (iy=0; iy<ny; iy++)
        {
            for (ix=0; ix<nx; ix++)
            {
                for (iz=0; iz<nz; iz++) 
                {
                    indx = fteik_model_grid2index(iz, ix, iy, nz, nz*nx);
                    diffz = (double) (iz - zri[irec])*dz;
                    diffx = (double) (ix - xri[irec])*dx;
                    diffy = (double) (iy - yri[irec])*dy;
                    t = sqrt(diffz*diffz + diffx*diffx + diffy*diffy)*slowUse;
                    ttimes[indx] = t; //ttimes[irec*ngrd+indx] = t;
                }
            }
        }
        // Set the travel-times
#ifdef FTEIK_USE_MPI
        locate_setTravelTimeField64fMPIF(master, ngrd, irec+1, ttimes, &ierr);
#else 
        locate_setTravelTimeField64fF(ngrd, irec+1, ttimes, &ierr);
#endif
        if (ierr != 0)
        {
            printf("Failed to set travel time field\n");
            return EXIT_FAILURE;
        }
    } // Loop on sources
    free(ttimes);
    t1 = clock();
    if (myid == master)
    {
        printf("Travel-time computation time %lf (s)\n",
               ((double)(t1-t0))/CLOCKS_PER_SEC);
    }
    // Set the observed travel times
    if (myid == master){printf("Setting observed travel times...\n");}
    t0 = clock();
    locate_setObservation64f(0, nrec,  false, 0.0,
                             obs2tf, tobs, wts, &ierr); 
    if (ierr != 0)
    {
        printf("Failed to set observations");
        return EXIT_FAILURE;
    }
    t1 = clock();
    if (myid == master)
    {
        printf("Time to set observations %lf (s)\n", ((double)(t1-t0))/CLOCKS_PER_SEC);
    }
    // Locate events
double t0Opt, objOpt;
int optIndx; 
    t0 = clock();
#ifdef FTEIK_USE_MPI
    locate_locateEventMPIF(1, &optIndx, &t0Opt, &objOpt);
#else
    locate_locateEventF(1, &optIndx, &t0Opt, &objOpt);
#endif
    t1 = clock();
if (myid == master)
{
printf("optindex, t0, obj: %d %e %e\n", optIndx, t0Opt, objOpt);
printf("True index: %d\n", fteik_model_grid2index(zsi[0], xsi[0], ysi[0], nz, nz*nx) + 1);
}
    if (myid == master)
    {
        printf("Time to locate event %lf (s)\n", ((double)(t1-t0))/CLOCKS_PER_SEC);
    }
    // Free memory
    if (myid == master){printf("Finalizing solver...\n");}
    locate_finalizeF();
printf("back\n");
    if (myid == master){printf("Freeing memory...\n");}
    free(obs2tf);
    free(wts);
    free(tobs);
    free(zri);
    free(xri);
    free(yri);
#ifdef FTEIK_USE_MPI
    MPI_Finalize();
#endif
    return EXIT_SUCCESS;
}
