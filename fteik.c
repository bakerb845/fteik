#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <getopt.h>
#include <limits.h>
#include "fteik.h"
#include <iniparser.h>

static void printUsage(void);
static int parseArguments(int argc, char *argv[],
                          char iniFile[PATH_MAX]);
void fteik_(const double *__restrict__ vel,
            double *__restrict__ tt, 
            const int *nz, const int *nx, const int *ny,
            const double *zsin, const double *xsin, const double *ysin,
            const double *dzin, const double *dxin, const double *dyin,
            const int *nsweep, const double *epsin);;
static int decimateVelocityModel(
    const int ncellz, const int ncellx, const int ncelly,
    const int nzGrid, const int nxGrid, const int nyGrid,
    const double dz, const double dx, const double dy, 
    const double *__restrict__ vel,
    double *dzGrid, double *dxGrid, double *dyGrid,
    double **vdecim);

struct fteikParms_struct
{
    char **modelFiles; /*!< Name of model files to read.  This is an array
                            of dimension [nmodels]. */
    double *xr;        /*!< x receiver locations in meters relative to origin.
                            This is an array of dimension [nrec]. */
    double *yr;        /*!< y receiver locations in meters relative to origin.
                            This is an array of dimension [nrec]. */
    double *zr;        /*!< z receiver locations in meters relative to origin.
                            This is an array of dimension [nrec]. */
    double *xs;        /*!< x source locations in meters relative to origin.
                            This is an array of dimension [nsrc]. */
    double *ys;        /*!< y source locations in meters relative to origin. 
                            This is an array of dimension [nsrc]. */
    double *zs;        /*!< z source locations in meters relative to origin.
                            This is an array of dimension [nsrc]. */
    double eps;        /*!< Radius, in grid points, where the finite 
                            differencing switches from spherical to Cartesian
                            operators. */
    int nsweep;        /*!< Number of Gauss-Seidel sweeping iterations. */
    int nmodels;       /*!< Number of velocity models. */
    int nrec;          /*!< Number of receivers. */
    int nsrc;          /*!< Number of sources. */
    int nzGrid;        /*!< Number of grid points in z. */
    int nxGrid;        /*!< Number of grid points in x. */
    int nyGrid;        /*!< Number of grid points in y. */
    bool ldecim;       /*!< If true then decimate the grid from the input model
                            to (nzGrid, nxGrid, nyGrid). */ 
    bool lreciprocal ; /*!< If false then performing source to receiver
                            travel time modeling. \n
                            If true then performing receiver to source
                            travel time modeling. */
};

int main(int argc, char *argv[])
{
    char iniFile[PATH_MAX];
    int ierr, isrc, m, ncell, ncellx, ncelly, ncellz,
        nrec, nsrc, nx, nx0, ny, ny0, nz, nz0;
    double *vel, *velDecim, *ttr, dx, dy, dz, x0, y0, z0;
    bool linit;
    const int myid = 0;
    const int master = 0;
    size_t nbytes;
    struct fteikSolver_struct solver;
    struct fteikParms_struct parms;
    ttr = NULL;
    memset(&parms,  0, sizeof(struct fteikParms_struct));
    memset(&solver, 0, sizeof(struct fteikSolver_struct)); 
    ierr = parseArguments(argc, argv, iniFile);
    if (ierr != 0)
    {
        if (ierr ==-2){return EXIT_SUCCESS;}
        fprintf(stderr, "%s: Error starting program\n", __func__);
        return EXIT_FAILURE;
    }
    // Read the ini file
    ierr = fteik_readIniFile(iniFile, &parms); 
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error reading ini file\n", __func__);
        return EXIT_FAILURE;
    }
    linit = false;
    // Loop on the models
    for (m=0; m<parms.nmodels; m++)
    {
        vel = NULL;
        ncell = 0;
        nsrc = 0;
        // Load the model
        if (myid == master)
        {
            fprintf(stdout, "%s: Loading model: %s\n",
                    __func__, parms.modelFiles[m]);
            // Read the model where the velocities are considered cellular
            ierr = fteik_io_readVelocityModel(parms.modelFiles[m],
                                              &ncellz, &ncellx, &ncelly,
                                              &dz, &dx, &dy,
                                              &z0, &x0, &y0,
                                              &vel);
            if (ierr != 0)
            {
                fprintf(stderr, "%s: Error reading file\n", __func__);
                return EXIT_FAILURE;
            }
            // Turns cellular data into grid points for travel time field
            nz = ncellz + 1;
            nx = ncellx + 1;
            ny = ncelly + 1;
int nzGrid = nz; //75;
int nxGrid = nx; //210;
int nyGrid = ny; //150;
            // Potentially decimate grid
            if (nzGrid != nz || nxGrid != nx || nyGrid != ny)
            {
                fprintf(stdout, "%s: Decimating grid...\n", __func__);
                nz0 = ncellz + 1;
                nx0 = ncellx + 1;
                ny0 = ncelly + 1;
                velDecim = NULL;
                decimateVelocityModel(ncellz, ncellx, ncelly,
                                      nzGrid, nxGrid, nyGrid,
                                      dz, dx, dy,
                                      vel,
                                      &dz, &dx, &dy,
                                      &velDecim);
                nz = nzGrid;
                nx = nxGrid;
                ny = nyGrid;
                ncellz = nz - 1;
                ncellx = nx - 1;
                ncelly = ny - 1;
                nbytes = (size_t) (ncellz*ncellx*ncelly)*sizeof(double);
                free(vel);
                vel = (double *) aligned_alloc(64, nbytes);
                memset(vel, 0, nbytes);
                memcpy(vel, velDecim, nbytes);
                free(velDecim);
                fprintf(stdout,
                        "%s: Grid decimated from (%d,%d,%d) to (%d,%d,%d)\n",
                         __func__, nz, nx, ny, nzGrid, nxGrid, nyGrid);
            }
            ncell = ncellz*ncellx*ncelly;
        }
        //MPI_Barrier(MPI_COMM_WORLD);
        // Initialize the solver for the first model
        if (m == 0)
        {
            fprintf(stdout, "%s: Initializing solver...\n", __func__);
            fteik_solver_initialize64fF(nz, nx, ny,
                                        z0, x0, y0,
                                        dz, dx, dy,
                                        parms.nsweep, parms.eps, &ierr);
            if (ierr != 0)
            {
                fprintf(stderr, "%s: Error initializing solver\n", __func__);
                return EXIT_FAILURE;
            }
            // Set the sources and receivers
            if (!parms.lreciprocal)
            {
                fprintf(stdout, "%s: Setting the sources...\n", __func__);
                fteik_solver_setSources64fF(parms.nsrc,
                                            parms.zs, parms.xs, parms.ys,
                                            &ierr);
                if (ierr != 0)
                {
                    fprintf(stderr, "%s: Error setting sources\n", __func__);
                    return EXIT_FAILURE;
                }
                fprintf(stdout, "%s: Setting the receivers...\n", __func__);
                fteik_solver_setReceivers64fF(parms.nrec,
                                              parms.zr, parms.xr, parms.yr,
                                              &ierr);
                if (ierr != 0)
                {
                    fprintf(stderr, "%s: Error setting receivers\n", __func__);
                    return EXIT_FAILURE;
                }
            }
            // Reciprocal modeling - set receivers as sources and vice versa
            else
            {
                fprintf(stdout, "%s: Setting receivers as sources...\n",
                        __func__);
                fteik_solver_setSources64fF(parms.nrec,
                                            parms.zr, parms.xr, parms.yr,
                                            &ierr); 
                if (ierr != 0)
                {
                    fprintf(stderr, "%s: Error setting sources\n", __func__);
                    return EXIT_FAILURE;
                }
                fprintf(stdout, "%s: Setting sources as receiver...\n",
                            __func__); 
                fteik_solver_setReceivers64fF(parms.nsrc,
                                              parms.zs, parms.xs, parms.ys,
                                              &ierr);
                if (ierr != 0)
                {
                    fprintf(stderr, "%s: Error setting receiver locations\n", 
                            __func__);
                    return EXIT_FAILURE;
                }
            }
            // Get the number of sources 
            fteik_solver_getNumberOfSources(&nsrc, &ierr);
            fteik_solver_getNumberOfReceivers(&nrec, &ierr);
            if (nrec > 0)
            {
                nbytes = nrec*sizeof(double);
                ttr = aligned_alloc(64, nbytes); 
                memset(ttr, 0, nbytes);
            }
            
            linit = true;
        }
        // Set the model
        fteik_solver_setVelocityModel64fF(ncell, vel, &ierr);
        if (vel != NULL){free(vel);}
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error setting velocity model\n", __func__);
            return EXIT_FAILURE;
        }
        // Solve eikonal equation for each source
        for (isrc=0; isrc<nsrc; isrc++)
        {
/*
            double *tt = (double *) calloc((size_t) (nz*nx*ny), sizeof(double));
            double zsin, xsin, ysin;
            if (!parms.lreciprocal)
            {
                 zsin = parms.zs[isrc] - z0;
                 xsin = parms.xs[isrc] - x0;
                 ysin = parms.ys[isrc] - y0;
            }
            else
            {
                 zsin = parms.zr[isrc] - z0; 
                 xsin = parms.xr[isrc] - x0; 
                 ysin = parms.yr[isrc] - y0;
            }
            printf("solving...\n");
printf("%e %e %e\n", zsin, xsin, ysin);
            fteik_(vel, tt, 
                   &nz, &nx, &ny, 
                   &zsin, &xsin, &ysin,
                   &dz, &dx, &dy,
                   &parms.nsweep, &parms.eps);
             int ngrd = nx*ny*nz;
 double *ttRef = calloc((size_t) nrec, sizeof(double));
             fteik_receiver_getTravelTimes64fF(nrec, ngrd, tt, ttRef, &ierr);
            free(tt);
*/

            fprintf(stdout, "%s: Solving source %d with LSMF...\n", __func__, isrc+1); 
            fteik_solver_solveSourceLSMF(isrc+1, &ierr);
            if (ierr != 0)
            {
                fprintf(stderr, "%s: Error solving eikonal equation\n", __func__);
                return EXIT_FAILURE;
            }
            // Get the travel times from the source to the receiver
            if (nrec > 0)
            {
                fteik_solver_getTravelTimes64fF(nrec, ttr, &ierr);
for (int i=0; i<nrec; i++)
{
//printf("%10.5f %10.5f\n", ttr[i], ttRef[i]);
printf("%10.5f\n", ttr[i]);
}
/*
    fteik_h5io_initializeF("./debug.h5");
    fteik_h5io_writeLevelSchedulesF();
    fteik_h5io_writeVelocityModelF("debugModel");
    fteik_h5io_writeTravelTimesF("debugTimes");

fteik_xdmf_writeVelocityModel(1, "debug.h5", "debugModel", false,
                              nz, nx, ny, dz, dx, dy, z0, x0, y0);
fteik_xdmf_writeVelocityModel(3, "debug.h5", "debugTimes", false,
                              nz, nx, ny, dz, dx, dy, z0, x0, y0);
fteik_xdmf_writeVelocityModel(4, "debug.h5", "debugModel", false,
                              nz, nx, ny, dz, dx, dy, z0, x0, y0);
fteik_xdmf_writeVelocityModel(5, "debug.h5", "debugModel", false,
                              nz, nx, ny, dz, dx, dy, z0, x0, y0);
*/
            }
        }  // Loop on sources
break;
    }
    // Release the memory on the eikonal solver
    //MPI_Barrier(MPI_COMM_WORLD);
    if (ttr != NULL){free(ttr);} 
    if (linit)
    {
        fprintf(stdout, "%s: Finalizing solver...\n", __func__);
        fteik_finalizeF();
    }
    return EXIT_SUCCESS;
}
//============================================================================//
/*!
 * @brief Reads the initialization file.
 */
int fteik_readIniFile(const char *iniFile,
                       struct fteikParms_struct *parms)
{
    const char *s;
    char varname[128], swork[128];
    dictionary *ini;
    int i, ierr, k;
    ierr = 1;
    memset(parms, 0, sizeof(struct fteikParms_struct));
    ini = iniparser_load(iniFile);
    parms->nmodels = iniparser_getint(ini, "fteik:nmodels\0", 0);
    if (parms->nmodels < 1)
    {
        fprintf(stderr, "%s: No models\n", __func__);
        goto ERROR;
    }
    parms->nsweep = iniparser_getint(ini, "fteik:nsweep\0", 2);
    if (parms->nsweep < 1)
    {
        fprintf(stdout, "%s: No Gauss-Seidel sweeps\n", __func__);
        goto ERROR;
    }
    parms->eps = iniparser_getdouble(ini, "fteik:eps\0", 2.0);
    if (parms->eps < 0.0)
    {
        fprintf(stdout, "%s: spherical to cartesian eps cannot be negative\n",
                __func__);
        goto ERROR;
    }
    parms->modelFiles
        = (char **) calloc((size_t) (parms->nmodels), sizeof(char *));
    for (i=0; i<parms->nmodels; i++)
    {
        parms->modelFiles[i] = (char *) calloc(PATH_MAX, sizeof(char));
        memset(varname, 0, 128*sizeof(char));
        sprintf(varname, "fteik:model_%d", i+1);
        s = iniparser_getstring(ini, varname, NULL);
        if (s != NULL)
        {
            strcpy(parms->modelFiles[i], s);
            if (!fteik_os_isfile(parms->modelFiles[i]))
            {
                fprintf(stderr, "%s: model file %s doesn't exist\n",
                        __func__, parms->modelFiles[i]);
                goto ERROR;
            }
        }
        else
        {
            fprintf(stderr, "%s: Error h5 file is NULL\n", __func__);
            goto ERROR;
        }
    }
    // Now read the receiver locations
    parms->lreciprocal
        = iniparser_getboolean(ini, "fteik:lreciprocal\0", false);
    parms->nsrc = iniparser_getint(ini, "fteik:nsrc\0", 0);
    parms->nrec = iniparser_getint(ini, "fteik:nrec\0", 0);
    if (parms->lreciprocal)
    {
        if (parms->nrec < 1)
        {
            fprintf(stderr, "%s: No receivers for reciprocal modeling\n",
                     __func__);
            goto ERROR;
        }
    }
    else
    {
        if (parms->nsrc < 1)
        {
            fprintf(stderr, "%s: No sources\n", __func__);
            goto ERROR;
        }
    }
    // Read the receivers
    if (parms->nrec > 0)
    {
        parms->xr = (double *) calloc((size_t) parms->nrec, sizeof(double));
        parms->yr = (double *) calloc((size_t) parms->nrec, sizeof(double));
        parms->zr = (double *) calloc((size_t) parms->nrec, sizeof(double));
        for (i=0; i<parms->nrec; i++)
        {
            memset(varname, 0, 128*sizeof(char));
            sprintf(varname, "fteik:xr_%d", i+1);
            s = iniparser_getstring(ini, varname, NULL);
            if (s == NULL)
            {
                fprintf(stderr, "%s: Receiver %d is not defined\n",
                        __func__, i + 1);
                goto ERROR;
            }
            memset(swork, 0, 128*sizeof(char));
            strcpy(swork, s);
            k = 0;
            char *token = strtok(swork, " ;,"); 
            while (token)
            {
                if (k == 0){parms->xr[i] = atof(token);}
                if (k == 1){parms->yr[i] = atof(token);}
                if (k == 2){parms->zr[i] = atof(token);}
                k = k + 1;
                token = strtok(NULL, " ;,");
            }
        }
    }
    // Read the sources
    if (parms->nsrc > 0)
    {
        parms->xs = (double *) calloc((size_t) parms->nsrc, sizeof(double));
        parms->ys = (double *) calloc((size_t) parms->nsrc, sizeof(double));
        parms->zs = (double *) calloc((size_t) parms->nsrc, sizeof(double));
        for (i=0; i<parms->nsrc; i++)
        {
            memset(varname, 0, 128*sizeof(char));
            sprintf(varname, "fteik:xs_%d", i+1);
            s = iniparser_getstring(ini, varname, NULL);
            if (s == NULL)
            {
                fprintf(stderr, "%s: Source %d is not defined\n",
                        __func__, i + 1); 
                goto ERROR;
            }
            memset(swork, 0, 128*sizeof(char));
            strcpy(swork, s); 
            k = 0;
            char *token = strtok(swork, " ;,"); 
            while (token)
            {
                if (k == 0){parms->xs[i] = atof(token);}
                if (k == 1){parms->ys[i] = atof(token);}
                if (k == 2){parms->zs[i] = atof(token);}
                k = k + 1;
                token = strtok(NULL, " ;,");
            }
        }
    }
    ierr = 0;
ERROR:;
    iniparser_freedict(ini);
    return ierr;
}
//============================================================================//
static int decimateVelocityModel(
    const int ncellz, const int ncellx, const int ncelly,
    const int nzGrid, const int nxGrid, const int nyGrid,
    const double dz, const double dx, const double dy,
    const double *__restrict__ vel,
    double *dzGrid, double *dxGrid, double *dyGrid,
    double **vdecim)
{
    double *vwork;
    int ix, iy, iz, jx, jy, jz, ncell, nxFact, nyFact, nzFact;
    nzFact = (ncellz - ncellz%(nzGrid - 2))/(nzGrid - 2);
    nxFact = (ncellx - ncellx%(nxGrid - 2))/(nxGrid - 2);
    nyFact = (ncelly - ncelly%(nyGrid - 2))/(nyGrid - 2);
    *dzGrid = dz*(double) nzFact;
    *dxGrid = dx*(double) nxFact;
    *dyGrid = dy*(double) nyFact;
    ncell = (nzGrid - 1)*(nxGrid - 1)*(nyGrid - 1);
    vwork = (double *) aligned_alloc(64, (size_t) ncell*sizeof(double));
    memset(vwork, 0, (size_t) ncell*sizeof(double));
    for (iy=0; iy<nyGrid-1; iy++)
    {
        for (ix=0; ix<nxGrid-1; ix++)
        {
            for (iz=0; iz<nzGrid-1; iz++)
            {
                jz = MAX(0, MIN(ncellz-1, iz*nzFact));
                jx = MAX(0, MIN(ncellx-1, ix*nxFact));
                jy = MAX(0, MIN(ncelly-1, iy*nyFact));
                vwork[iy*(nzGrid-1)*(nxGrid-1) + ix*(nzGrid-1) + iz]
                   = vel[jy*(ncellz)*(ncellx) + jx*(ncellz) + jz];
            }
        }
    }
    *vdecim = vwork;
    return 0;
}
//============================================================================//
/*!
 * @brief Extracts the ini file from the input arguments.
 *
 */
static int parseArguments(int argc, char *argv[],
                          char iniFile[PATH_MAX])
{
    bool liniFile;
    memset(iniFile, 0, PATH_MAX*sizeof(char));
    liniFile = false;
    while (true)
    {
        static struct option longOptions[] =
        {
            {"help", no_argument, 0, '?'},
            {"help", no_argument, 0, 'h'},
            {"iniFile", required_argument, 0, 'i'},
            {0, 0, 0, 0}
        };
        int c, optionIndex;
        c = getopt_long(argc, argv, "?hi:",
                        longOptions, &optionIndex);
        if (c ==-1){break;}
        if (c == 'i')
        {
            strcpy(iniFile, (const char *) optarg);
            liniFile = true;
        }
        else if (c == 'h' || c == '?')
        {
            printUsage();
            return -2;
        }
        else
        {
            fprintf(stderr, "%s: Unknown options: %s\n",
                    __func__, argv[optionIndex]);
        }
    }
    if (liniFile)
    {
        if (!fteik_os_isfile(iniFile))
        {
            fprintf(stderr, "%s: ini file %s does not exist\n",
                    __func__, iniFile);
            return -1;
        }
    }
    else
    {
        fprintf(stderr, "%s: Invalid arguments\n", __func__);
        printUsage(); 
        return -1;
    }
    return 0;
}

static void printUsage(void)
{
    fprintf(stdout, "Usage:\n    fteik -i iniFile\n\n");
    fprintf(stdout, "Require arguments:\n");
    fprintf(stdout, "   -i iniFile is the initialization file\n\n");
    fprintf(stdout, "Optional arguments:\n");
    fprintf(stdout, "   -h displays this message\n");
    return; 
}
