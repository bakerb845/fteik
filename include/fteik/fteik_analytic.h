#ifndef FTEIK_ANALYTIC_H__
#define FTEIK_ANALYTIC_H__ 1

#ifdef __cplusplus
extern "C"
{
#endif

int analyticSolution_wholeSpace(
    const int nz, const int nx, const int ny, 
    const double dz, const double dx, const double dy, 
    const double vconst, double *__restrict__ vel);
int analyticSolution_wholeSpaceAnalyticSolution(
    const int nz, const int nx, const int ny, 
    const double dz, const double dx, const double dy, 
    const double zs, const double xs, const double ys, 
    const double vconst, double *__restrict__ ttimes);

#ifdef __cplusplus
}
#endif
#endif
