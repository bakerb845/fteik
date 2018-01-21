#ifndef FTEIK_GRAPH_H__
#define FTEIK_GRAPH_H__ 1
#include "fteik/fteik_config.h"
#include "fteik/fteik_struct.h"

#ifdef __cplusplus
extern "C"
{
#endif

void graph_sweep2sgnx64f(const int sweep,
                         int *sgntz, int *sgntx, int *sgnty,
                         int *sgnvz, int *sgnvx, int *sgnvy,
                         double *sgnrz, double *sgnrx, double *sgnry);
int graph_createTravelTimeAccessPointers(struct levelSet_struct *levelSet);
int graph_createUpdateMask(const int sweep,
                           struct levelSet_struct *levelSet);
int graph_createUpdateInitMask(const int sweep,
                               const int zsi, const int xsi, const int ysi,
                               struct levelSet_struct *levelSet);
int graph_getSweepLoopLimits(const int sweep, const bool linitk,
                             const int nz, const int nx, const int ny,
                             const int zsi, const int xsi, const int ysi,
                             int *z1, int *z2,
                             int *x1, int *x2,
                             int *y1, int *y2);
int graph_createLevelStructs(const int nz, const int nx, const int ny, 
                             struct levelSets_struct *levelSets);
int graph_createGraph3D(const int nz, const int nx, const int ny,
                        int **xadjOut, int **adjncyOut);
int graph_getVelocityPermutation(const int sweep,
                                 int *__restrict__ perm);
int graph_getTravelTimePermutation(const int sweep,
                                   int *__restrict__ perm);
int graph_freeLevelsStruct(struct levelSets_struct *levelSets);
int graph_freeLevelStruct(struct levelSet_struct *levelSet);
int graph_createLevelStructure(const int n, const int *__restrict__ xadj,
                               const int *__restrict__ adjncy,
                               int *nLevels, int **levelPtrOut,
                               int **ijkvOut, int **n2lOut);
int graph_permuteLevelStructure(const int nz, const int nx, const int ny, 
                                const int sweep, const int nLevels,
                                const int *__restrict__ n2l,
                                int **ijkvOut,
                                int **levelPtrOut,
                                int **n2lPermOut);


#ifdef __cplusplus
}
#endif
#endif
