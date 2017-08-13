#ifndef FTEIK_STRUCT_H__
#define FTEIK_STRUCT_H__ 1
#include <stdbool.h>
#include <stdint.h>
#define DEFAULT_FTEIK_FD_CHUNK_SIZE 16 /*< Default chunk size. */

typedef double *__attribute__((aligned(64))) double64f;

struct levelSet_struct
{
    int *ijkv;      /*!< Maps from the i'th node in the elimination
                         tree to the global node index [nnodes]. */
    int *levelPtr;  /*!< Maps from the level'th level to the start
                         index of ijkv [nLevels+1]. */
    int *n2ijk;     /*!< Maps from the i'th node to the (i,j,k) indices
                         in the grid where i corresponds to z,
                         j corresponds to x, and k corresponds to y
                         [3*nnodes]. */
    int *n2l;       /*!< Maps from the i'th node to the level
                         number [nnodes]. */
    int *ttPerm;    /*!< Maps from travel-times to a travel time workspace array
                         that for better cache-use in travel time accesses
                         [8*nnodes]. */
    int16_t *ttLocPerm; /*!< Maps from travel time workspace array to a local
                         travel time array appropriate for finite differencing.
                         [8*nnodes]. */
    bool *lupdate;  /*!< If true then this node is able to be updated [nnodes]. */
    bool *lupdateInit; /*!< If true then this node is unable to be update
d                           during the initialization phase because it is
                            upwind for the source [nnodes]. */
    int nLevels;    /*!< Number of levels. */
    int nnodes;     /*!< Number of nodes in travel time grid (=nx*ny*nz). */
    int ny;         /*!< Number of y nodes in travel time grid. */
    int nx;         /*!< Number of x nodes in travel time grid. */
    int nz;         /*!< Number of z nodes in travel time grid. */
    int sweep;      /*!< Sweep number of level set structure. */
    int chunkSize;  /*!< Chunk size for vectorizing finite difference
                         operator. */
};

struct levelSets_struct
{
    struct levelSet_struct
      *levelSet;           /*!< Level set structure for each sweep [nsweeps]. */
    int nsweeps;           /*!< Number of sweeps (should be 8). */
};

struct fteikSolver_struct
{
    struct levelSets_struct levelSets; /*!< Level set structures. */
    double *tt;  /*!< Travel times at each grid points (seconds).  This
                      has dimensions [nz x nx x ny] and is packed in 
                      column major format. */
    double *slow; /*!< Slowness model (meters/seconds).  This has 
                       dimensions [nz-1 x nx-1 x ny-1] and is packed in
                       column major format. */
    double dz; /*!< Grid spacing (meters) in z. */
    double dx; /*!< Grid spacing (meters) in x. */
    double dy; /*!< Grid spacing (meters) in y. */
    double z0; /*!< Z origin (meters). Default is 0. */
    double x0; /*!< X origin (meters). Default is 0. */
    double y0; /*!< Y origin (meters). Default is 0. */
    double eps; /*!< Radius (meters) where the solver transisitions
                     from a spherical approximation to a cartesian
                     finite difference. */
    double _zsa;
    double _xsa;
    double _ysa;
    double _dzi;
    double _dxi;
    double _dyi;
    double _dz2i;
    double _dx2i;
    double _dy2i;
    double _dz2dx2;
    double _dz2dy2;
    double _dx2dy2;
    double _dsum;
    double _dsumi;
    double _szero;
    int zsi;   /*!< Source z index in grid. */
    int xsi;   /*!< Source x index in grid. */
    int ysi;   /*!< Source y index in grid. */
    int nz;    /*!< Number of z grid points. */
    int nx;    /*!< Number of x grid points. */
    int ny;    /*!< Number of y grid points. */
    int ncell; /*!< Number of cells in mesh (= (nz-1)*(nx-1)*(ny-1)). */
    int ngrd;  /*!< Number of grid points (=nz*nx*ny). */
    int nsweeps; /*!< Number of sweeps to perform over model. */
    bool lhaveGraph; /*!< If true then the graph analysis has been 
                          performed. */
    bool luseVelPtr; /*!< If true then the velocity model is a pointer.
                          Otherwise it has been copied to the structure and 
                          exists in a private memory address space. */
    bool luseTTPtr;  /*!< If true then the travel times is a pointer.
                          Otherwise it has been copied to the structure
                          and exists in a private memory address space. */
};

#endif
