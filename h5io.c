#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fteik_h5io.h"
#include "h5io.h"

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

