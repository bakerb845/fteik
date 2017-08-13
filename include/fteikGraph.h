#ifndef FTEIK_CPP_GRAPHINTER
#define FTEIK_CPP_GRAPHINTER 1
//#include "fteikGraph.hpp"
#ifdef __cplusplus
extern "C"
{
#endif

void *fteik_graph_initializeF(const int *nz, const int *nx,
                                    const int *ny, int *ierr);
void *fteik_graph_intialize(const int nz, const int nx, const int ny, 
                            int *ierr);
void fteik_graph_finalizeF(fteikGraph *graphIn);
void fteik_graph_finalize(void *graphIn);
int fteik_graph_getNumberOfLevels(void *graphIn);
int fteik_graph_getNumberOfGridPoints(void *graphIn);
int fteik_graph_getIJKVF(void *graphIn, const int *nwork,
                         const int *sweep, int *__restrict__ ijkv);
int fteik_graph_getIJKV(void *graphIn, const int nwork,
                        const int sweep, int *__restrict__ ijkv);
extern "C" int fteik_graph_getLevelPointer(void *graphIn, const int nwork,
                                           const int sweep,
                                           int *__restrict__ levelPtr);


#ifdef __cplusplus
}
#endif
#endif
