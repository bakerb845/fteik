#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <float.h>
#include <libgen.h>
#include <hdf5.h>
#include <stdbool.h>
#include "fteik_os.h"

static void unpackInt16ModelToDouble(const bool lzDown,
                                     const int nz, const int nx, const int ny, 
                                     const int16_t *__restrict__ vi2,
                                     double *__restrict__ vin);
static void unpackInt32ModelToDouble(const bool lzDown,
                                     const int nz, const int nx, const int ny,
                                     const int32_t *__restrict__ vi4,
                                     double *__restrict__ vin);
static void unpackFloatModelToDouble(const bool lzDown,
                                     const int nz, const int nx, const int ny, 
                                     const float *__restrict__ v4,
                                     double *__restrict__ vin);
static void unpackDoubleModelToDouble(const bool lzDown,
                                      const int nz, const int nx, const int ny, 
                                      const double *__restrict__ v8, 
                                      double *__restrict__ vin);

 
/*
int fteik_io_writeVelocityModel(const char *fileName,
                            const int ncellz, const int ncellx, const int nelly,
                            const double dz, const double dx, const double dy,
                            const double z0, const double x0, const double y0,
                            const double *__restrict__ vel)
{

}
*/
//============================================================================//
/*!
 * @brief Reads the velocity model from an HDF5 file.  This is a finite
 *        difference model approriate for viewing with XDMF and Paraview.
 *
 * @param[in] fileName     Name of file to read.
 * @param[out] ncellz      Number of cells in z.
 * @param[out] ncellx      Number of cells in x.
 * @param[out] ncelly      Number of cells in y.
 * @param[out] dz          Cell spacing in z (meters).
 * @param[out] dx          Cell spacing in x (meters).
 * @param[out] dy          Cell spacing in y (meters).
 * @param[out] z0          z origin (meters).
 * @param[out] x0          x origin (meters).
 * @param[out] y0          y origin (meters).
 * @param[out] vel         Velocity model (m/s).  This is an array of
 *                         dimension [ncellz x ncellx x ncelly] where
 *                         the (iz,ix,iy)'th cell is accessed by:
 *                         iy*ncellz*ncellx + ix*ncellz + iz.  It has
 *                         convention +x west-to-east, +y north-to-south,
 *                         and +z down.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker
 *
 * @copyright MIT
 *
 */
int fteik_io_readVelocityModel(const char *fileName,
                               int *ncellz, int *ncellx, int *ncelly,
                               double *dz, double *dx, double *dy,
                               double *z0, double *x0, double *y0,
                               double *__restrict__ vel[])
{
    hid_t attrID, attrSpace, dataSet, dataSpace, dataType, fileID, memSpace;
    int zDown, rightHanded;
    hsize_t dims[3];
    herr_t status;
    int i, ierr, nCell, nx, ny, nz, rank;
    double *vin, *v8, vmin, vmax;
    float *v4;
    int32_t *vi4;
    int16_t *vi2;
    bool lzDown;
    status = 0;
    ierr = 0;
    vi2 = NULL;
    vi4 = NULL;
    v4 = NULL;
    v8 = NULL; 
    *ncellz = 0;
    *ncellx = 0;
    *ncelly = 0;
    if (!fteik_os_isfile(fileName))
    {
        fprintf(stderr, "%s: File %s doesn't exist\n", __func__, fileName);
        return -1;
    }
    printf("%s: Loading %s...\n", __func__, fileName); //basename(fileName));
    // Open file and check the p velocity model exists
    fileID = H5Fopen(fileName, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (H5Lexists(fileID, "PVelocity\0", H5P_DEFAULT) <= 0)
    {
        fprintf(stderr, "%s: PVelocity model doesn't exist\n", __func__);
        H5Fclose(fileID);
        return -1;
    }
    // Open the dataset
    dataSet = H5Dopen(fileID, "PVelocity\0", H5P_DEFAULT);
    // Pick off its attributes; dx
    attrID = H5Aopen(dataSet, "dx\0", H5P_DEFAULT);
    attrSpace = H5Aget_space(attrID);
    H5Aread(attrID, H5T_NATIVE_DOUBLE, dx);
    H5Sclose(attrSpace);
    H5Aclose(attrID);
    // dy
    attrID = H5Aopen(dataSet, "dy\0", H5P_DEFAULT);
    attrSpace = H5Aget_space(attrID);
    H5Aread(attrID, H5T_NATIVE_DOUBLE, dy);
    H5Sclose(attrSpace);
    H5Aclose(attrID);
    // dz
    attrID = H5Aopen(dataSet, "dz\0", H5P_DEFAULT);
    attrSpace = H5Aget_space(attrID);
    H5Aread(attrID, H5T_NATIVE_DOUBLE, dz);
    H5Sclose(attrSpace);
    H5Aclose(attrID);
    // and the origin
    attrID = H5Aopen(dataSet, "x0\0", H5P_DEFAULT);
    attrSpace = H5Aget_space(attrID);
    H5Aread(attrID, H5T_NATIVE_DOUBLE, x0);
    H5Sclose(attrSpace);
    H5Aclose(attrID);
    // dy
    attrID = H5Aopen(dataSet, "y0\0", H5P_DEFAULT);
    attrSpace = H5Aget_space(attrID);
    H5Aread(attrID, H5T_NATIVE_DOUBLE, y0);
    H5Sclose(attrSpace);
    H5Aclose(attrID);
    // dz
    attrID = H5Aopen(dataSet, "z0\0", H5P_DEFAULT);
    attrSpace = H5Aget_space(attrID);
    H5Aread(attrID, H5T_NATIVE_DOUBLE, z0);
    H5Sclose(attrSpace);
    H5Aclose(attrID);
    // z up or down
    attrID = H5Aopen(dataSet, "zPositiveDown\0", H5P_DEFAULT);
    attrSpace = H5Aget_space(attrID);
    H5Aread(attrID, H5T_NATIVE_INT, &zDown);
    H5Sclose(attrSpace);
    H5Aclose(attrID);
    // left-handed coordinate system
    attrID = H5Aopen(dataSet, "rightHandedCoordinateSystem\0", H5P_DEFAULT);
    attrSpace = H5Aget_space(attrID);
    H5Aread(attrID, H5T_NATIVE_INT, &rightHanded);
    H5Sclose(attrSpace);
    H5Aclose(attrID);
    // Get the data sizes
    dataSpace = H5Dget_space(dataSet);
    rank = H5Sget_simple_extent_ndims(dataSpace);
    if (rank != 3)
    {
        printf("Error model rank=%d must be 3D\n", rank);
        ierr = 1;
        goto ERROR;
    }
    H5Sget_simple_extent_dims(dataSpace, dims, NULL);
    *ncellz = (int) dims[0];
    *ncelly = (int) dims[1];
    *ncellx = (int) dims[2];
    lzDown = false;
    if (zDown == 1){lzDown = true;}
    nCell = (*ncellx)*(*ncelly)*(*ncellz);
    if (nCell < 1)
    {
        printf("%s: No data to read!\n", __func__);
        ierr = 1;
        goto ERROR;
    } 
    // Let's read it
    dataType = H5Dget_type(dataSet);
    memSpace = H5Screate_simple(rank, dims, NULL);
    vin = (double *) calloc((size_t) nCell, sizeof(double));
    if (H5Tequal(dataType, H5T_NATIVE_INT16))
    {
        printf("Reading int16_t data\n");
        vi2 = (int16_t *) calloc((size_t) (nCell), sizeof(int16_t));
        status = H5Dread(dataSet, H5T_NATIVE_INT16, memSpace,
                         dataSpace, H5P_DEFAULT, vi2);
        unpackInt16ModelToDouble(lzDown, *ncellz, *ncellx, *ncelly, vi2, vin);
    }
    else if (H5Tequal(dataType, H5T_NATIVE_INT32))
    {
        printf("Reading int32_t data\n");
        vi4 = (int32_t *) calloc((size_t) (nCell), sizeof(int32_t));
        status = H5Dread(dataSet, H5T_NATIVE_INT32, memSpace, 
                         dataSpace, H5P_DEFAULT, vi4);
        unpackInt32ModelToDouble(lzDown, *ncellz, *ncellx, *ncelly, vi4, vin);
    }
    else if (H5Tequal(dataType, H5T_NATIVE_FLOAT))
    {
        printf("Reading float data\n");
        v4 = (float *) calloc((size_t) (nCell), sizeof(float));
        status = H5Dread(dataSet, H5T_NATIVE_FLOAT, memSpace, 
                         dataSpace, H5P_DEFAULT, v4);
        unpackFloatModelToDouble(lzDown, *ncellz, *ncellx, *ncelly, v4, vin);
    }
    else if (H5Tequal(dataType, H5T_NATIVE_DOUBLE))
    {
        printf("Double not yet done\n");
        v8 = (double *) calloc((size_t) (nCell), sizeof(double));
        status = H5Dread(dataSet, H5T_NATIVE_DOUBLE, memSpace, 
                         dataSpace, H5P_DEFAULT, v8);
        unpackDoubleModelToDouble(lzDown, *ncellz, *ncellx, *ncelly, v8, vin);
    }
    else
    {
        printf("Unknown datatype\n");
        status =-1;
    }
    H5Sclose(memSpace);
    H5Tclose(dataType);
    if (status < 0)
    {
        printf("Error reading velocity model\n");
        ierr = 1;
        goto ERROR;
    }
    printf("Number of cells (nz,nx,ny)=(%d,%d,%d)\n",
           *ncellz, *ncellx, *ncelly);
    printf("Model origin (z,x,y)=(%f,%f,%f) meters\n", *z0, *x0, *y0);
    printf("Cell spacing in (z,x,y)=(%f,%f,%f) meters\n", *dz ,*dx, *dy);
    lzDown = true;
    if (!zDown){lzDown = false;}
    if ((rightHanded && !lzDown) || (!(rightHanded) && lzDown))
    {
        printf("I don't know how to unpack this model\n");
        ierr = 1;
        goto ERROR;
    }
    // Unpack the model
    vmin = DBL_MAX;
    vmax =-DBL_MAX;
    #pragma omp simd reduction(min: vmin) reduction(max: vmax)
    for (i=0; i<nCell; i++)
    {
        vmin = fmin(vmin, vin[i]);
        vmax = fmax(vmax, vin[i]);
    }
    if (vmin <= 0.0)
    {
        printf("Error minimum velocity %f is not positive\n", vmin);
        ierr = 1;
    }
    else
    {
         printf("Maximum velocity is %f m/s\n", vmax);
         printf("Minimum velocity is %f m/s\n", vmin);
    }
ERROR:;
    if (vi2 != NULL){free(vi2);}
    if (vi4 != NULL){free(vi4);}
    if (v4  != NULL){free(v4);}
    if (v8  != NULL){free(v8);}
    H5Sclose(dataSpace);
    H5Dclose(dataSet); 
    H5Fclose(fileID); 
    *vel = vin;
    return ierr; 
}

//============================================================================//

static void unpackInt16ModelToDouble(const bool lzDown,
                                     const int nz, const int nx, const int ny,
                                     const int16_t *__restrict__ vi2,
                                     double *__restrict__ vin)
{
    int indx, ix, iy, iz, jndx;
    if (!lzDown)
    {
        for (iy=0; iy<ny; iy++)
        {
            for (ix=0; ix<nx; ix++)
            {
                for (iz=0; iz<nz; iz++)
                {
                    indx = iz*nx*ny + iy*nx + ix;
                    jndx = iy*nz*nx + ix*nz + iz;
                    vin[jndx] = (double) vi2[indx];
                }
            }
        }
    }
    else
    {
        printf("Flipping model...\n");
        for (iy=0; iy<ny; iy++)
        {
            for (ix=0; ix<nx; ix++)
            {
                for (iz=0; iz<nz; iz++)
                {
                    indx = (nz-1-iz)*nx*ny + iy*nx + ix;
                    jndx = iy*nz*nx + ix*nz + iz;
                    vin[jndx] = (double) vi2[indx];
                }
            }
        }
    }
}

static void unpackInt32ModelToDouble(const bool lzDown,
                                     const int nz, const int nx, const int ny, 
                                     const int32_t *__restrict__ vi4,
                                     double *__restrict__ vin)
{
    int indx, ix, iy, iz, jndx;
    if (lzDown)
    {   
        for (iy=0; iy<ny; iy++)
        {
            for (ix=0; ix<nx; ix++)
            {
                for (iz=0; iz<nz; iz++)
                {
                    indx = iz*nx*ny + iy*nx + ix; 
                    jndx = iy*nz*nx + ix*nz + iz; 
                    vin[jndx] = (double) vi4[indx];
                }
            }
        }
    }   
    else
    {   
        printf("Flipping model...\n");
        for (iy=0; iy<ny; iy++)
        {
            for (ix=0; ix<nx; ix++)
            {
                for (iz=0; iz<nz; iz++)
                {
                    indx = (nz-1-iz)*nx*ny + iy*nx + ix; 
                    jndx = iy*nz*nx + ix*nz + iz; 
                    vin[jndx] = (double) vi4[indx];
                }
            }
        }
    }   
}

static void unpackFloatModelToDouble(const bool lzDown,
                                     const int nz, const int nx, const int ny, 
                                     const float *__restrict__ v4,
                                     double *__restrict__ vin)
{
    int indx, ix, iy, iz, jndx;
    if (lzDown)
    {   
        for (iy=0; iy<ny; iy++)
        {
            for (ix=0; ix<nx; ix++)
            {
                for (iz=0; iz<nz; iz++)
                {
                    indx = iz*nx*ny + iy*nx + ix; 
                    jndx = iy*nz*nx + ix*nz + iz; 
                    vin[jndx] = (double) v4[indx];
                }
            }
        }
    }   
    else
    {   
        printf("Flipping model...\n");
        for (iy=0; iy<ny; iy++)
        {
            for (ix=0; ix<nx; ix++)
            {
                for (iz=0; iz<nz; iz++)
                {
                    indx = (nz-1-iz)*nx*ny + iy*nx + ix; 
                    jndx = iy*nz*nx + ix*nz + iz; 
                    vin[jndx] = (double) v4[indx];
                }
            }
        }
    }   
}

static void unpackDoubleModelToDouble(const bool lzDown,
                                      const int nz, const int nx, const int ny,
                                      const double *__restrict__ v8, 
                                      double *__restrict__ vin)
{
    int indx, ix, iy, iz, jndx;
    if (lzDown)
    {   
        for (iy=0; iy<ny; iy++)
        {
            for (ix=0; ix<nx; ix++)
            {
                for (iz=0; iz<nz; iz++)
                {
                    indx = iz*nx*ny + iy*nx + ix; 
                    jndx = iy*nz*nx + ix*nz + iz; 
                    vin[jndx] = v8[indx];
                }
            }
        }
    }   
    else
    {   
        printf("Flipping model...\n");
        for (iy=0; iy<ny; iy++)
        {
            for (ix=0; ix<nx; ix++)
            {
                for (iz=0; iz<nz; iz++)
                {
                    indx = (nz-1-iz)*nx*ny + iy*nx + ix; 
                    jndx = iy*nz*nx + ix*nz + iz; 
                    vin[jndx] = v8[indx];
                }
            }
        }
    }   
}
