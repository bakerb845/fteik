#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>
#include "fteik/h5io.h"


int h5io_closeFileF(const int64_t fileID)
{
    hid_t fid;
    herr_t status;
    int ierr;
    fid = (hid_t) fileID;
    status = H5Fclose(fid);
    ierr = (int) status;
    return ierr;
}
//============================================================================//
int64_t h5io_openFileReadOnlyF(const char *fname)
{
    hid_t fid;
    int64_t fileID;
    fid = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);    
    fileID = (int) fid; 
    return fid;
}
//============================================================================//
int64_t h5io_openFileReadWriteF(const char *fname)
{
    hid_t fid;
    int64_t fileID;
    fid = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);    
    fileID = (int64_t) fid; 
    return fileID;
}
//============================================================================//
int64_t h5io_createFileF(const char *fname)
{
    hid_t fid;
    int64_t fileID;
    fid = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    fileID = (int64_t) fid;
    return fileID;
}
//============================================================================//
int64_t h5io_createGroup_finter(const int64_t *fid, const char *groupName)
{
    hid_t fileID, groupID;
    int64_t group;
    fileID = (hid_t) *fid;
    if (H5Lexists(fileID, groupName, H5P_DEFAULT) > 0)
    {
        fprintf(stderr, "%s: Group %s already exists\n", __func__, groupName);
        groupID = H5Gopen2(fileID, groupName, H5P_DEFAULT);
    }
    else
    {
        groupID = H5Gcreate2(fileID, groupName,
                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }
    group = (int64_t) groupID; 
    return group;
}
//============================================================================//
int64_t h5io_openGroup_finter(const int64_t *fid, const char *groupName)
{
    hid_t fileID, groupID;
    int64_t group;
    fileID = (hid_t) *fid;
    groupID = H5Gopen2(fileID, groupName, H5P_DEFAULT);
    group = (int64_t) groupID;
    return group;
}
//============================================================================//
int64_t h5io_openDataset_finter(const int64_t *fid, const char *dataSetName)
{
    hid_t fileID, datasetID;
    int64_t dataset; 
    fileID = (hid_t) *fid;
    datasetID = H5Dopen(fileID, dataSetName, H5P_DEFAULT);
    dataset = (int64_t) datasetID;
    return dataset;
}
//============================================================================//
int h5io_closeDataset_finter(const int64_t *dataset)
{
    hid_t datasetID;
    herr_t status;
    int ierr;
    datasetID = (hid_t) *dataset;
    status = H5Dclose(datasetID);
    ierr = (int) status;
    return status;
}
//============================================================================//
int h5io_closeGroup_finter(const int64_t *group)
{
    hid_t groupID;
    herr_t status;
    int ierr;
    groupID = (hid_t) *group;
    status = H5Gclose(groupID);
    ierr = (int) status;
    return status;
}
//============================================================================//
int h5io_writeAttribute32i_finter(const int64_t *dataSetID,
                                  const char *attributeName,
                                  const int *n, const int *__restrict__ x)
{
    hid_t dataSet, fileID;
    int ierr;
    dataSet = (hid_t) *dataSetID;
    ierr = h5io_writeAttribute32i(dataSet, attributeName, *n, x); 
    return ierr;
}
//============================================================================//
int h5io_writeAttribute64f_finter(const int64_t *dataSetID,
                                  const char *attributeName,
                                  const int *n, const double *__restrict__ x)
{
    hid_t dataSet, fileID;
    int ierr;
    dataSet = (hid_t) *dataSetID;
    ierr = h5io_writeAttribute64f(dataSet, attributeName, *n, x);
    return ierr;
}
//============================================================================//
int h5io_writeArray32iF(const int64_t fid, const char *dataName,
                        const int n, const int *__restrict__ x)
{
    hid_t fileID;
    int ierr;
    fileID = (hid_t) fid;
    ierr = h5io_writeArray32i(fileID, dataName, n, x);
    return ierr;
}
//============================================================================//
int h5io_writeArray64fF(const int64_t fid, const char *dataName,
                        const int n, const double *__restrict__ x)
{
    hid_t fileID;
    int ierr;
    fileID = (hid_t) fid;
    ierr = h5io_writeArray64f(fileID, dataName, n, x);
    return ierr;
}
//============================================================================//
int h5io_readArray32i_finter(const int64_t *fid, const char *dataName,
                             const int *nwork, int *n, int *__restrict__ x)
{
    hid_t fileID;
    int ierr;
    fileID = (hid_t) *fid;
    ierr = h5io_readArray32i(fileID, dataName, *nwork, n, x);
    return ierr; 
}
//============================================================================//
int h5io_readArray64f_finter(const int64_t *fid, const char *dataName,
                             const int *nwork, int *n, double *__restrict__ x)
{
    hid_t fileID;
    int ierr;
    fileID = (hid_t) *fid;
    ierr = h5io_readArray64f(fileID, dataName, *nwork, n, x); 
    return ierr; 
}
//============================================================================//
int h5io_writeArray16iF(const int64_t fileID, const char *dataName,
                        const int n, const int16_t *__restrict__ x)
{
    hid_t fid = (hid_t) fileID;
    return h5io_writeArray16iF(fid, dataName, n, x);
}
int h5io_writeArray16i(const hid_t fileID, const char *dataName,
                       const int n, const int16_t *__restrict__ x)
{
    char *citem;
    int ierr;
    hsize_t dims[1];
    hid_t dataSet, dataSpace;
    herr_t status;
    //------------------------------------------------------------------------//
    ierr = 0;
    if (n < 1 || x == NULL)
    {   
        if (n < 1){fprintf(stderr, "%s: Error no points to write\n", __func__);}
        if (x == NULL){fprintf(stderr, "%s: Error x is NULL\n", __func__);}
        return -1; 
    }   
    // Copy the data name
    citem = (char *) calloc(strlen(dataName)+1, sizeof(char));
    if (H5Lexists(fileID, dataName, H5P_DEFAULT) > 0)
    {
        fprintf(stderr, "%s: item %s already exists\n", __func__, citem);
        free(citem);
        return -1;
    }
    strcpy(citem, dataName);
    dims[0] = (hsize_t) n;
    // Create the dataspace
    dataSpace = H5Screate_simple(1, dims, NULL);
    // Create the dataset
    dataSet = H5Dcreate2(fileID, citem, H5T_NATIVE_SHORT,
                         dataSpace,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    free(citem);
    // Write it
    status = H5Dwrite(dataSet, H5T_NATIVE_SHORT, H5S_ALL,
                      H5S_ALL, H5P_DEFAULT, x); 
    if (status != 0)
    {   
        fprintf(stderr, "%s: Failed to write data\n", __func__);
        ierr =-1;
    }   
    // Free H5 resources
    status  = H5Sclose(dataSpace);
    status += H5Dclose(dataSet);
    if (status != 0)
    {   
        fprintf(stderr, "%s: Failed to close dataset\n", __func__);
        ierr =-1;
    }   
    return ierr;
}
//============================================================================//
int h5io_writeArray32i(const hid_t fileID, const char *dataName,
                       const int n, const int *__restrict__ x)
{
    char *citem;
    int ierr;
    hsize_t dims[1];
    hid_t dataSet, dataSpace;
    herr_t status;
    //------------------------------------------------------------------------//
    ierr = 0;
    if (n < 1 || x == NULL)
    {
        if (n < 1){fprintf(stderr, "%s: Error no points to write\n", __func__);}
        if (x == NULL){fprintf(stderr, "%s: Error x is NULL\n", __func__);}
        return -1;
    }
    // Copy the data name
    citem = (char *) calloc(strlen(dataName)+1, sizeof(char));
    strcpy(citem, dataName);
    dims[0] = (hsize_t) n;
    // Create the dataspace
    dataSpace = H5Screate_simple(1, dims, NULL);
    // Create the dataset
    dataSet = H5Dcreate2(fileID, citem, H5T_NATIVE_INT,
                         dataSpace,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    free(citem);
    // Write it
    status = H5Dwrite(dataSet, H5T_NATIVE_INT, H5S_ALL,
                      H5S_ALL, H5P_DEFAULT, x);
    if (status != 0)
    {   
        fprintf(stderr, "%s: Failed to write data\n", __func__);
        ierr =-1;
    }
    // Free H5 resources
    status  = H5Sclose(dataSpace);
    status += H5Dclose(dataSet);
    if (status != 0)
    {
        fprintf(stderr, "%s: Failed to close dataset\n", __func__);
        ierr =-1;
    }
    return ierr;
}
//============================================================================//
int h5io_writeArray32f(const hid_t fileID, const char *dataName,
                       const int n, const float *__restrict__ x)
{
    char *citem;
    int ierr;
    hsize_t dims[1];
    hid_t dataSet, dataSpace;
    herr_t status;
    //------------------------------------------------------------------------//
    ierr = 0;
    if (n < 1 || x == NULL)
    {
        if (n < 1){fprintf(stderr, "%s: Error no points to write\n", __func__);}
        if (x == NULL){fprintf(stderr, "%s: Error x is NULL\n", __func__);}
        return -1;
    }
    // Copy the data name
    citem = (char *) calloc(strlen(dataName)+1, sizeof(char));
    strcpy(citem, dataName);
    dims[0] = (hsize_t) n;
    // Create the dataspace
    dataSpace = H5Screate_simple(1, dims, NULL);
    // Create the dataset
    dataSet = H5Dcreate2(fileID, citem, H5T_NATIVE_FLOAT,
                         dataSpace,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    free(citem);
    // Write it
    status = H5Dwrite(dataSet, H5T_NATIVE_FLOAT, H5S_ALL,
                      H5S_ALL, H5P_DEFAULT, x);
    if (status != 0)
    {
        fprintf(stderr, "%s: Failed to write data\n", __func__);
        ierr =-1;
    }
    // Free H5 resources
    status  = H5Sclose(dataSpace);
    status += H5Dclose(dataSet);
    if (status != 0)
    {
        fprintf(stderr, "%s: Failed to close dataset\n", __func__);
        ierr =-1;
    }
    return ierr;
}
//============================================================================//
int h5io_writeAttribute64f(const hid_t dataSet, const char *attributeName,
                           const int n, const double *__restrict__ x)
{
    char *citem;
    hid_t attribute, dataSpace;
    herr_t status;
    hsize_t dims[1];
    int ierr;

    ierr = 0; 
    citem = (char *) calloc(strlen(attributeName)+1, sizeof(char));
    strcpy(citem, attributeName); 
    dims[0] = n;
    dataSpace = H5Screate_simple(1, dims, NULL);
    attribute = H5Acreate2(dataSet, citem, H5T_NATIVE_DOUBLE,
                           dataSpace, H5P_DEFAULT, H5P_DEFAULT); 
    status = H5Awrite(attribute, H5T_NATIVE_DOUBLE, x);
    status = H5Aclose(attribute);
    status = H5Sclose(dataSpace);
    free(citem);
    return ierr;
}
//============================================================================//
int h5io_readAttribute32i(const hid_t dataSet, const char *attributeName,
                          const int nwork, int *n, int *__restrict__ x)
{
    char *citem;
    hid_t aspace, attr;
    herr_t status;
    hsize_t dims[1];
    int ierr, rank;
    // open the attribute space and get the size 
    *n = 0;
    citem = (char *) calloc(strlen(attributeName)+1, sizeof(char));
    strcpy(citem, attributeName);
    attr = H5Aopen(dataSet, citem, H5P_DEFAULT);
    aspace = H5Aget_space(attr);
    rank = H5Sget_simple_extent_ndims(aspace);
    if (rank != 1)
    {
        printf("%s: Only rank 1 attributes done\n", __func__);
        return -1;
    }
    H5Sget_simple_extent_dims(aspace, dims, NULL);
    *n = (int) dims[0];
    if (nwork < 0)
    {
        H5Sclose(aspace);
        H5Aclose(attr);
        return 0;
    }
    if (*n > nwork)
    {
        printf("%s: Insufficient space in x\n", __func__);
        H5Sclose(aspace);
        H5Aclose(attr);
        return -1;
    }
    // read it
    status = H5Aread(attr, H5T_NATIVE_INT, x);
    if (status != 0)
    {
        fprintf(stderr, "%s: Failed to read attribute\n", __func__);
        ierr = 1;
    } 
    // close it up
    status  = H5Sclose(aspace);
    status += H5Aclose(attr);
    free(citem);
    ierr = (int) status;
    return ierr;
}
//============================================================================//
int h5io_writeAttribute32i(const hid_t dataSet, const char *attributeName,
                           const int n, const int *__restrict__ x)
{
    char *citem;
    hid_t attribute, dataSpace;
    herr_t status;
    hsize_t dims[1];
    int ierr;

    ierr = 0;
    citem = (char *) calloc(strlen(attributeName)+1, sizeof(char));
    strcpy(citem, attributeName);
    dims[0] = n;
    dataSpace = H5Screate_simple(1, dims, NULL);
    attribute = H5Acreate2(dataSet, citem, H5T_NATIVE_INT,
                           dataSpace, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(attribute, H5T_NATIVE_INT, x);
    if (status < 0)
    {
        fprintf(stderr, "%s: Failed to write attribute\n", __func__);
    }
    status = H5Aclose(attribute);
    status = H5Sclose(dataSpace);
    free(citem);
    return ierr;
}
//============================================================================//
int h5io_writeArray64f(const hid_t fileID, const char *dataName,
                       const int n, const double *__restrict__ x)
{
    char *citem;
    int ierr;
    hsize_t dims[1];
    hid_t dataSet, dataSpace;
    herr_t status;
    //------------------------------------------------------------------------//
    ierr = 0;
    if (n < 1 || x == NULL)
    {
        if (n < 1){printf("%s: Error no points to write\n", __func__);}
        if (x == NULL){printf("%s: Error x is nULL\n", __func__);}
        return -1;
    }
    // Copy the data name
    citem = (char *) calloc(strlen(dataName)+1, sizeof(char));
    strcpy(citem, dataName);
    dims[0] = (hsize_t) n;
    // Create the dataspace
    dataSpace = H5Screate_simple(1, dims, NULL);
    // Create the dataset
    dataSet = H5Dcreate2(fileID, citem, H5T_NATIVE_DOUBLE,
                         dataSpace,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    free(citem);
    // Write it
    status = H5Dwrite(dataSet, H5T_NATIVE_DOUBLE, H5S_ALL,
                      H5S_ALL, H5P_DEFAULT, x);
    if (status != 0)
    {
        printf("%s: Failed to write data\n", __func__);
        ierr =-1;
    }
    // Free H5 resources
    status  = H5Sclose(dataSpace);
    status += H5Dclose(dataSet);
    if (status != 0)
    {
        printf("%s: Failed to close dataset\n", __func__);
        ierr =-1;
    }
    return ierr;
}
//============================================================================//
int h5io_readArray64f(const hid_t fileID, const char *dataName, 
                      const int nwork, int *n, double *__restrict__ x)
{
    char *citem;
    int i, ierr, nw, rank;
    hid_t dataSet, dataSpace, memSpace;
    herr_t status;
    hsize_t *dims;
    ierr = 0;
    if (H5Lexists(fileID, dataName, H5P_DEFAULT) < 1)
    {
        printf("%s: Dataset %s does not exist\n", __func__, dataName);
        return -1; 
    }
    // Open the dataset
    citem = (char *) calloc(strlen(dataName)+1, sizeof(char));
    strcpy(citem, dataName);
    dataSet = H5Dopen(fileID, citem, H5P_DEFAULT);
    free(citem);
    // Create the dataspace
    dataSpace = H5Dget_space(dataSet);
    rank = H5Sget_simple_extent_ndims(dataSpace);
    dims = (hsize_t *) calloc(rank, sizeof(hsize_t));
    H5Sget_simple_extent_dims(dataSpace, dims, NULL);
    // Get sizes
    nw = 1;
    for (i=0; i<rank; i++)
    {   
        nw = nw*(int) dims[i];
    }   
    *n = nw; 
    // Memory query
    if (nwork ==-1){goto EXIT1;}
    // Avoid a segfault
    if (nw > nwork || x == NULL)
    {
        if (nw > nwork)
        {
            fprintf(stderr, "%s: Insufficient workspace %d - require %d\n",
                     __func__, nwork, nw);
        }
        if (x == NULL){fprintf(stderr, "%s: Error x is NULL\n", __func__);}
        ierr = 1;
        goto EXIT1;
    }
    // Create the memory space and load the data
    memSpace = H5Screate_simple(rank, dims, NULL);
    status = H5Dread(dataSet, H5T_NATIVE_DOUBLE, memSpace,
                     dataSpace, H5P_DEFAULT, x);
    if (status != 0)
    {
        fprintf(stderr, "%s: Failed to read dataset %s\n", __func__, dataName);
        ierr = 1;
    }
    H5Sclose(memSpace);
EXIT1:;
    free(dims);
    H5Sclose(dataSpace);
    H5Dclose(dataSet);
    return ierr;
}
//============================================================================//
int h5io_readArray32i(const hid_t fileID, const char *dataName, 
                      const int nwork, int *n, int *__restrict__ x)
{
    char *citem;
    int i, ierr, nw, rank;
    hid_t dataSet, dataSpace, memSpace;
    herr_t status;
    hsize_t *dims;
    ierr = 0;
    if (H5Lexists(fileID, dataName, H5P_DEFAULT) < 1)
    {
        printf("%s: Dataset %s does not exist\n", __func__, dataName);
        return -1;
    }
    // Open the dataset
    citem = (char *) calloc(strlen(dataName)+1, sizeof(char));
    strcpy(citem, dataName);
    dataSet = H5Dopen(fileID, citem, H5P_DEFAULT);
    free(citem);
    // Create the dataspace
    dataSpace = H5Dget_space(dataSet);
    rank = H5Sget_simple_extent_ndims(dataSpace);
    dims = (hsize_t *) calloc(rank, sizeof(hsize_t));
    H5Sget_simple_extent_dims(dataSpace, dims, NULL);
    // Get sizes
    nw = 1;
    for (i=0; i<rank; i++)
    {
        nw = nw*(int) dims[i];
    }
    *n = nw;
    // Memory query
    if (nwork ==-1){goto EXIT1;}
    // Avoid a segfault
    if (nw > nwork || x == NULL)
    {
        if (nw > nwork)
        {
            fprintf(stderr, "%s: Insufficient workspace %d - require %d\n",
                    __func__, nwork, nw);
        }
        if (x == NULL){fprintf(stderr, "%s: Error x is NULL\n", __func__);}
        ierr = 1;
        goto EXIT1;
    }
    // Create the memory space and load the data
    memSpace = H5Screate_simple(rank, dims, NULL);
    status = H5Dread(dataSet, H5T_NATIVE_INT, memSpace,
                     dataSpace, H5P_DEFAULT, x);
    if (status != 0)
    {
        fprintf(stderr, "%s: Failed to read dataset %s\n", __func__, dataName);
        ierr = 1;
    }
    H5Sclose(memSpace);
EXIT1:;
    free(dims);
    H5Sclose(dataSpace);
    H5Dclose(dataSet);
    return ierr;
}
