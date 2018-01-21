#ifndef FTEIK_IO_H__
#define FTEIK_IO_H__ 1

#ifdef __cplusplus
extern "C"
{
#endif

int fteik_io_readVelocityModel(const char *fileName,
                               int *ncellz, int *ncellx, int *ncelly,
                               double *dz, double *dx, double *dy,
                               double *z0, double *x0, double *y0,
                               double *__restrict__ vel[]);
int fteik_io_writeVelocityModel(const char *fileName,
                            const bool lisVP,
                            const int ncellz, const int ncellx, const int ncelly,
                            const double dz, const double dx, const double dy, 
                            const double z0, const double x0, const double y0, 
                            const double *__restrict__ vel);
#ifdef __cplusplus
}
#endif
#endif
