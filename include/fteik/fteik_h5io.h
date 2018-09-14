#ifndef FTEIK_H5IO_H__
#define FTEIK_H5IO_H__ 1
#include <stdbool.h>
#include <limits.h>
#include "fteik/fteik_config.h"
#include <hdf5.h>

struct xdmf_struct
{
    FILE *xdmfl;
    char h5fl[PATH_MAX];
    char projnm[PATH_MAX];
    int nx; 
    int ny; 
    int nz; 
    int nelem;
    int nnodes;
    int ipart;
    bool linit;
    bool lfileOpen;
};

#ifdef __cplusplus
extern "C"
{
#endif

int fteik_h5io_initializeF(const char *fileName);
int fteik_h5io_finalizeF(void);
int fteik_h5io_initialize(const char *fileName,
                          const int nz, const int nx, const int ny,
                          const double dz, const double dx, const double dy,
                          const double z0, const double x0, const double y0);
int fteik_h5io_writeVelocityModelF(const char *dataName);
int fteik_h5io_writeTravelTimesF(const char *dataName);
int fteik_h5io_writeLevelSchedulesF(void);


int fteik_h5io_writeTravelTimes64f(const hid_t fileID, const char *ttName,
                                   const int nz, const int nx, const int ny, 
                                   const double *__restrict__ tt);
int fteik_h5io_writeVelocityModel64f(const hid_t fileID, const char *velName,
                                     const int nz, const int nx, const int ny,
                                     const double *__restrict__ vel);
int fteik_h5io_writeGeometry(const hid_t fileID,
                             const int nz, const int nx, const int ny,
                             const double dz, const double dx, const double dy,
                             const double z0, const double x0, const double y0);
int fteik_h5io_writeLevelSet32i(const hid_t fileID, const char *levelSetName,
                                const int nz, const int nx, const int ny,
                                const int *__restrict__ level);


int fteik_xdmf_initialize(const char *h5dir,
                          const char *h5flName, const char *projnm,
                          const int nz, const int nx, const int ny, 
                          struct xdmf_struct *xdmf);
int fteik_xdmf_addLevelSet(const struct xdmf_struct xdmf,
                           const int sweep, const char *levelName);
int fteik_xdmf_addVelocityModel(const struct xdmf_struct xdmf,
                                const char *vmodelName);
int fteik_xdmf_addTravelTimes(const struct xdmf_struct xdmf,
                              const char *ttName);
int fteik_xdmf_finalize(struct xdmf_struct *xdmf);

#ifdef __cplusplus
}
#endif
#endif
