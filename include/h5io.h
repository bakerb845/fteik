#ifndef H5IO_H__
#define H5IO_H__ 1
#include <stdint.h>
#include <hdf5.h>

#ifdef __cplusplus
extern "C"
{
#endif
int64_t h5io_openDataset_finter(const int64_t *fid, const char *dataSetName);
int h5io_writeArray32i_finter(const int64_t *fid, const char *dataName,
                              const int *n, const int *__restrict__ x);
int h5io_writeArray32i(const hid_t fileID, const char *dataName,
                       const int n, const int *__restrict__ x);
int h5io_writeArray64f_finter(const int64_t *fid, const char *dataName,
                              const int *n, const double *__restrict__ x);
int h5io_writeArray64f(const hid_t fileID, const char *dataName,
                       const int n, const double *__restrict__ x);
int h5io_writeArray32f(const hid_t fileID, const char *dataName,
                       const int n, const float *__restrict__ x);
int h5io_writeAttribute64f_finter(const int64_t *dataSetID,
                                  const char *attributeName,
                                  const int *n, const double *__restrict__ x);
int h5io_writeAttribute64f(const hid_t dataSet, const char *attributeName,
                           const int n, const double *__restrict__ x);
int h5io_writeAttribute32i_finter(const int64_t *dataSetID,
                                  const char *attributeName,
                                  const int *n, const int *__restrict__ x);
int h5io_readAttribute32i(const hid_t dataSet, const char *attributeName,
                          const int nwork, int *n, int *__restrict__ x);
int h5io_writeAttribute32i(const hid_t dataSet, const char *attributeName,
                           const int n, const int *__restrict__ x);
int h5io_readArray32i_finter(const int64_t *fid, const char *dataName,
                             const int *nwork, int *n, int *__restrict__ x);
int h5io_readArray32i(const hid_t fileID, const char *dataName,
                      const int nwork, int *n, int *__restrict__ x);
int h5io_readArray64f_finter(const int64_t *fid, const char *dataName,
                             const int *nwork, int *n, double *__restrict__ x);
int h5io_readArray64f(const hid_t fileID, const char *dataName,
                      const int nwork, int *n, double *__restrict__ x);

#ifdef __cplusplus
}
#endif
#endif
