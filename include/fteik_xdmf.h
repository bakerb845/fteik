#ifndef FTEIK_XDMF_H__
#define FTEIK_XDMF_H__ 1

#ifdef __cplusplus
extern "C"
{
#endif

int fteik_xdmf_writeVelocityModel(const int job, 
                                  const char *h5flName,
                                  const char *velModel,
                                  const bool lflip,
                                  const int nz, const int nx, const int ny, 
                                  const double dz, const double dx, 
                                  const double dy, 
                                  const double z0, const double x0, 
                                  const double y0);

#ifdef __cplusplus
}
#endif

#endif
