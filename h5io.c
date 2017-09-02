#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "fteik_h5io.h"
#include "h5io.h"
#include "fteik_os.h"

#define VELMODEL_GROUP "/VelocityModels\0"
#define LEVEL_SCHEDULER_GROUP "/LevelSchedules\0"
#define TRAVELTIMES_GROUP "/TravelTimeFields\0"

int fteik_h5io_initialize(const char *fileName)
{
    hid_t h5fl, groupID;
    char dirName[PATH_MAX];
    // Check the file name makes sense
    if (fileName == NULL)
    {
        fprintf(stderr, "%s: Error file name is NULL\n", __func__);
        return -1;
    }
    if (strlen(fileName) == 0)
    {
        fprintf(stderr, "%s: Error file name is blank\n", __func__);
        return -1;
    }
    // Does the output directory exist?
    fteik_os_dirname_work(fileName, dirName);
    if (!fteik_os_path_isdir(dirName))
    {
        // Try to make it
        if (fteik_os_makedirs(dirName) != 0)
        {
            fprintf(stderr, "%s: Failed to make dirctory: %s\n",
                    __func__, dirName);
            return -1;
        }
    }
    // Warn the user I'm about to zap their file
    if (fteik_os_isfile(fileName))
    {
        fprintf(stdout, "%s: Overwriting %s\n", __func__, fileName);
    }
    h5fl = H5Fcreate(fileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    // Directory for level sets
    groupID = H5Gcreate2(h5fl, LEVEL_SCHEDULER_GROUP,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(groupID);
    // Directory for velocity models
    groupID = H5Gcreate2(h5fl, VELMODEL_GROUP,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(groupID);
    // Directory for the travel times
    groupID = H5Gcreate2(h5fl, TRAVELTIMES_GROUP,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
 
    H5Fclose(h5fl);
    return 0;    

}

int fteik_h5io_writeGeometry(const hid_t fileID,
                             const int nz, const int nx, const int ny,
                             const double dz, const double dx, const double dy,
                             const double z0, const double x0, const double y0)
{
    const char *fcnm = "fteik_h5io_writeGeometry\0";
    double *xlocs, *ylocs, *zlocs;
    hid_t groupID;
    int *connect, i, ielem, ierr, ix, iy, iz, nelem, nnodes;
    // Create and open the geometry group for writing
    groupID = H5Gcreate2(fileID, "/Geometry\0",
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); 
    // Set space
    nnodes = nx*ny*nz;
    nelem = (nx - 1)*(ny - 1)*(nz - 1);
    xlocs = (double *) calloc((size_t) nnodes, sizeof(double));
    ylocs = (double *) calloc((size_t) nnodes, sizeof(double));
    zlocs = (double *) calloc((size_t) nnodes, sizeof(double));
    connect = (int *) calloc((size_t) (8*nelem), sizeof(int));
    // Define the hexahedral connectivity
    for (iy=0; iy<ny-1; iy++)
    {
        for (ix=0; ix<nx-1; ix++)
        {
            for (iz=0; iz<nz-1; iz++)
            {
                ielem = iy*(nx - 1)*(nz - 1) + ix*(nz - 1) + iz;
                connect[8*ielem+0] = (iy + 0)*nx*nz + (ix + 0)*nz + iz;
                connect[8*ielem+1] = (iy + 0)*nx*nz + (ix + 1)*nz + iz;
                connect[8*ielem+2] = (iy + 1)*nx*nz + (ix + 1)*nz + iz;
                connect[8*ielem+3] = (iy + 1)*nx*nz + (ix + 0)*nz + iz;
                connect[8*ielem+4] = (iy + 0)*nx*nz + (ix + 0)*nz + iz + 1;
                connect[8*ielem+5] = (iy + 0)*nx*nz + (ix + 1)*nz + iz + 1;
                connect[8*ielem+6] = (iy + 1)*nx*nz + (ix + 1)*nz + iz + 1;
                connect[8*ielem+7] = (iy + 1)*nx*nz + (ix + 0)*nz + iz + 1;
            }
        }
    }
    // Create the physical nodal locations
    for (iy=0; iy<ny; iy++)
    {
        for (ix=0; ix<nx; ix++)
        {
            for (iz=0; iz<nz; iz++)
            {
                i = iy*nx*nz + ix*nz + iz;
                xlocs[i] = x0 + dx*(double) ix;
                ylocs[i] = y0 + dy*(double) iy;
                zlocs[i] = z0 + dz*(double) iz;
            }
        }
    }
    // Write it
    ierr = 0;
    ierr += h5io_writeArray64f(groupID, "xLocations", nnodes, xlocs);
    if (ierr != 0)
    {
        printf("%s: Failed to write xLocations\n", fcnm);
        goto ERROR;
    }
    ierr += h5io_writeArray64f(groupID, "yLocations", nnodes, ylocs);
    if (ierr != 0)
    {
        printf("%s: Failed to write yLocations\n", fcnm);
        goto ERROR;
    }
    ierr += h5io_writeArray64f(groupID, "zLocations", nnodes, zlocs);
    if (ierr != 0)
    {
        printf("%s: Failed to write zLocations\n", fcnm);
        goto ERROR;
    }
    ierr += h5io_writeArray32i(groupID, "Connectivity", 8*nelem, connect);
    if (ierr != 0)
    {
        printf("%s: Failed to write Connectivity\n", fcnm);
        goto ERROR;
    }
ERROR:; 
    // Release workspace
    free(connect);
    free(xlocs);
    free(ylocs);
    free(zlocs); 
    H5Gclose(groupID);
    return ierr;
}
//============================================================================//
int fteik_h5io_writeLevelScheduleF(const int64_t h5fl, const int sweep,
                                   const int nz, const int nx, const int ny,
                                   const int16_t *__restrict__ levelSchedule)
{
    char dataName[64];
    hid_t fileID, groupID;
    herr_t status;
    int ierr, ngrd;
    // Require the group exists and if not make it 
    fileID = (hid_t) h5fl;
    if (H5Lexists(fileID, LEVEL_SCHEDULER_GROUP, H5P_DEFAULT) > 0)
    {
        groupID = H5Gopen2(fileID, LEVEL_SCHEDULER_GROUP, H5P_DEFAULT);
    }
    else
    {
        groupID = H5Gcreate2(fileID, LEVEL_SCHEDULER_GROUP,
                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }
    memset(dataName, 0, 64*sizeof(char));
    sprintf(dataName, "Sweep_%d", sweep);
    ngrd = nz*nx*ny;
    ierr = h5io_writeArray16i(groupID, dataName, ngrd, levelSchedule);
    if (ierr != 0)
    {   
        fprintf(stderr, "%s: Error writing %s\n", __func__, dataName);
    }   
    status = H5Gclose(groupID);
    if (status < 0)
    {   
        fprintf(stderr, "%s: Error closing group\n", __func__);
        ierr = 1;
    }
    return ierr;
}
//============================================================================//
int fteik_h5io_writeTravelTimes32fF(const int64_t h5fl, const char *ttName,
                                    const int nz, const int nx, const int ny,
                                    const float *__restrict__ tt)
{
    hid_t fileID, groupID;
    herr_t status;
    int ngrd;
    int ierr;
    // Require the group exists and if not make it 
    fileID = (hid_t) h5fl;
    if (H5Lexists(fileID, TRAVELTIMES_GROUP, H5P_DEFAULT) > 0)
    {
        groupID = H5Gopen2(fileID, TRAVELTIMES_GROUP, H5P_DEFAULT);
    }
    else
    {
        groupID = H5Gcreate2(fileID, TRAVELTIMES_GROUP,
                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }
    // Check this velocity model doesn't exist
    if (H5Lexists(groupID, ttName, H5P_DEFAULT) > 0)
    {
        fprintf(stderr, "%s: Error %s already exists\n", __func__, ttName);
        H5Gclose(groupID);
        return -1;
    }
    ngrd = nx*ny*nz;
    // Write it
    ierr = h5io_writeArray32f(groupID, ttName, ngrd, tt);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error writing %s\n", __func__, ttName);
    }
    status = H5Gclose(groupID);
    if (status < 0)
    {
        fprintf(stderr, "%s: Error closing group\n", __func__);
        ierr = 1;
    }
    return ierr;
}
//============================================================================//
int fteik_h5io_writeVelocityModel16iF(const int64_t h5fl, const char *velName,
                                      const int nz, const int nx, const int ny,
                                      const int16_t *__restrict__ vel)
{
    hid_t fileID, groupID;
    herr_t status;
    int nelem;
    int ierr;
    // Require the group exists and if not make it 
    fileID = (hid_t) h5fl;
    if (H5Lexists(fileID, VELMODEL_GROUP, H5P_DEFAULT) > 0)
    {   
        groupID = H5Gopen2(fileID, VELMODEL_GROUP, H5P_DEFAULT);
    }   
    else
    {
        groupID = H5Gcreate2(fileID, VELMODEL_GROUP,
                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }
    // Check this velocity model doesn't exist
    if (H5Lexists(groupID, velName, H5P_DEFAULT) > 0)
    {
        fprintf(stderr, "%s: Error %s already exists\n", __func__, velName);
        H5Gclose(groupID);
        return -1;
    }
    nelem = (nx - 1)*(ny - 1)*(nz - 1);
    // Write it
    ierr = h5io_writeArray16i(groupID, velName, nelem, vel);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error writing %s\n", __func__, velName);
    } 
    status = H5Gclose(groupID);
    if (status < 0)
    {
        fprintf(stderr, "%s: Error closing group\n", __func__);
        ierr = 1;
    }
    return ierr;
}
//============================================================================//
int fteik_h5io_writeVelocityModel64f(const hid_t fileID, const char *velName,
                                     const int nz, const int nx, const int ny,
                                     const double *__restrict__ vel)
{
    const char *fcnm = "fteik_h5io_writeVelocityModel64f\0";
    int ierr, nelem;
    hid_t groupID;
    if (H5Lexists(fileID, "/VelocityModel", H5P_DEFAULT) > 0)
    {
        groupID = H5Gopen2(fileID, "/VelocityModel\0", H5P_DEFAULT); 
    }
    else
    {
        groupID = H5Gcreate2(fileID, "/VelocityModel\0",
                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }
    nelem = (nx - 1)*(ny - 1)*(nz - 1); 
    ierr = h5io_writeArray64f(groupID, velName, nelem, vel); 
    if (ierr != 0)
    {
        printf("%s: Failed to write velocity model %s\n", fcnm, velName);
    }
    H5Gclose(groupID);
    return ierr;
}
//============================================================================//
/*!
 * @brief Writes the travel times for the given model.
 */
int fteik_h5io_writeTravelTimes64f(const hid_t fileID, const char *ttName,
                                   const int nz, const int nx, const int ny, 
                                   const double *__restrict__ tt)
{
    const char *fcnm = "fteik_h5io_writeTravelTimes64f\0";
    int ierr, nnodes;
    hid_t groupID;
    if (H5Lexists(fileID, "/TravelTimes", H5P_DEFAULT) > 0)
    {   
        groupID = H5Gopen2(fileID, "/TravelTimes\0", H5P_DEFAULT); 
    }   
    else
    {   
        groupID = H5Gcreate2(fileID, "/TravelTimes\0",
                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }   
    nnodes = nx*ny*nz;
    ierr = h5io_writeArray64f(groupID, ttName, nnodes, tt); 
    if (ierr != 0)
    {   
        printf("%s: Failed to write travel time model %s\n", fcnm, ttName);
    }   
    H5Gclose(groupID);
    return ierr;
}
//============================================================================//
int fteik_h5io_writeLevelSet32i(const hid_t fileID, const char *levelSetName,
                                const int nz, const int nx, const int ny,
                                const int *__restrict__ level)
{
    const char *fcnm = "fteik_h5io_writeLevelSet32i\0";
    int ierr, nnodes;
    hid_t groupID;
    if (H5Lexists(fileID, "/LevelSet", H5P_DEFAULT) > 0)
    {   
        groupID = H5Gopen2(fileID, "/LevelSet\0", H5P_DEFAULT); 
    }   
    else
    {   
        groupID = H5Gcreate2(fileID, "/LevelSet\0",
                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }
    nnodes = nx*ny*nz;
    ierr = h5io_writeArray32i(groupID, levelSetName, nnodes, level);
    if (ierr != 0)
    {
        printf("%s: Failed to write levelset %s\n", fcnm, levelSetName);
    }
    H5Gclose(groupID);
    return ierr;
}
//============================================================================//

