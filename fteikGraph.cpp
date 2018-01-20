#include <boost/config.hpp>
#include <iostream>
#include <new>
#include <vector>
#include <utility>
#include <stdio.h>
#include <stdlib.h>
#include <boost/graph/properties.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/edge_coloring.hpp>
#include "fteikGraph.hpp"
#include "fteikGraph.h"
#include "sorting.h"

using namespace std;
using namespace boost;

extern "C"
void fteikTestGraph(const int nz, const int nx, const int ny)
{
     size_t ngrd;
     int *levels, nLevels, sweep;

     fteikGraph G = fteikGraph(nz, nx, ny);
     ngrd = G.getNumberOfGridPoints(); 
     printf("Computing node levels...\n");
     sweep = 1;
     G.computeNodeLevels(sweep);
     nLevels = G.getNumberOfLevels();
     printf("Number of levels: %d\n", nLevels);
     levels = G.getLevelPointer(sweep);
     return;
}
//============================================================================//
/*!
 * @brief Initializes the graph class.
 *
 * @param[in] nzIn    Number of z grid points.  Must be at least 3.
 * @param[in] nxIn    Number of x grid points.  Must be at least 3.
 * @param[in] nyIn    Number of y grid points.  Must be at least 3.
 *
 */
fteikGraph::fteikGraph(const int nzIn, const int nxIn, const int nyIn,
                       const int verboseIn)
{
    const char *fcnm = "fteikGraph\0";
    int ierr;
    linit = false;
    levels1 = NULL;
    levelPtr1 = NULL;
/*
    levels2 = NULL; levels3 = NULL; levels4 = NULL;
    levels5 = NULL; levels6 = NULL; levels7 = NULL; levels8 = NULL;
    levelPtr2 = NULL; levelPtr3 = NULL; levelPtr4 = NULL;
    levelPtr5 = NULL; levelPtr6 = NULL; levelPtr7 = NULL; levelPtr8 = NULL;
*/
    nx = 0;
    ny = 0; 
    nz = 0;
    ngrd = 0;
    nnz = 0;
    nLevels = 0;
    verbose = verboseIn;
    if (nzIn < 3)
    {
        printf("%s: Error nz=%d must be at least 3\n", fcnm, nzIn);
        return;
    }    
    if (nxIn < 3)
    {
        printf("%s: Error nx=%d must be at least 3\n", fcnm, nxIn);
        return; 
    }
    if (nyIn < 3)
    {
        printf("%s: Error ny=%d must be at least 3\n", fcnm, nyIn);
        return; 
    }
    nz = nzIn;
    nx = nxIn;
    ny = nyIn;
    ngrd = static_cast<size_t> (nz)
          *static_cast<size_t> (nx)
          *static_cast<size_t> (ny);
    levels1 = (int *) aligned_alloc(64, ngrd*sizeof(int));
    memset(levels1, -1, ngrd*sizeof(int));
/*
    levels2 = (int *) aligned_alloc(64, ngrd*sizeof(int));
    levels3 = (int *) aligned_alloc(64, ngrd*sizeof(int));
    levels4 = (int *) aligned_alloc(64, ngrd*sizeof(int));
    levels5 = (int *) aligned_alloc(64, ngrd*sizeof(int));
    levels6 = (int *) aligned_alloc(64, ngrd*sizeof(int));
    levels7 = (int *) aligned_alloc(64, ngrd*sizeof(int));
    levels8 = (int *) aligned_alloc(64, ngrd*sizeof(int));
    memset(levels2, -1, ngrd*sizeof(int));
    memset(levels3, -1, ngrd*sizeof(int));
    memset(levels4, -1, ngrd*sizeof(int));
    memset(levels5, -1, ngrd*sizeof(int));
    memset(levels6, -1, ngrd*sizeof(int));
    memset(levels7, -1, ngrd*sizeof(int));
    memset(levels8, -1, ngrd*sizeof(int));
*/
    ijkv1 = (int *) aligned_alloc(64, 4*ngrd*sizeof(int));
    ijkv2 = (int *) aligned_alloc(64, 4*ngrd*sizeof(int));
    ijkv3 = (int *) aligned_alloc(64, 4*ngrd*sizeof(int));
    ijkv4 = (int *) aligned_alloc(64, 4*ngrd*sizeof(int));
    ijkv5 = (int *) aligned_alloc(64, 4*ngrd*sizeof(int));
    ijkv6 = (int *) aligned_alloc(64, 4*ngrd*sizeof(int));
    ijkv7 = (int *) aligned_alloc(64, 4*ngrd*sizeof(int));
    ijkv8 = (int *) aligned_alloc(64, 4*ngrd*sizeof(int));
    memset(ijkv1, -1, 4*ngrd*sizeof(int));
    memset(ijkv2, -1, 4*ngrd*sizeof(int));
    memset(ijkv3, -1, 4*ngrd*sizeof(int));
    memset(ijkv4, -1, 4*ngrd*sizeof(int));
    memset(ijkv5, -1, 4*ngrd*sizeof(int));
    memset(ijkv6, -1, 4*ngrd*sizeof(int));
    memset(ijkv7, -1, 4*ngrd*sizeof(int));
    memset(ijkv8, -1, 4*ngrd*sizeof(int));
    // Some initialization
    generateSweepOrderings();
    setSweepSigns();
    // Compute the levels; if the first sweep fails the others will too
    linit = true;
    ierr = computeNodeLevels(1);
    if (ierr != 0)
    {
        printf("%s: Failed to compute first level\n", fcnm);
        finalize();
        return;
    }
/*
    for (int sweep=2; sweep<=8; sweep++)
    {
        ierr = computeNodeLevels(sweep);
        if (ierr != 0)
        {
            printf("%s: Failed to create level structure for sweep: %d\n",
                   fcnm, sweep);
            finalize();
            return;
        }
    }
*/
    return;
}
//============================================================================//
/*!
 * @brief Class destructor.
 */
fteikGraph::~fteikGraph(void)
{
    finalize();
    return;
}
//============================================================================//
/*!
 * @brief Resets internal variables and releases memory.
 */
void fteikGraph::finalize(void)
{
    linit = false;
    nz = 0;
    nx = 0;
    ny = 0;
    ngrd = 0;
    nnz = 0;
    nLevels = 0;
    verbose = 0;
    if (levels1 != NULL){free(levels1);}
/*
    if (levels2 != NULL){free(levels2);}
    if (levels3 != NULL){free(levels3);}
    if (levels4 != NULL){free(levels4);}
    if (levels5 != NULL){free(levels5);}
    if (levels6 != NULL){free(levels6);}
    if (levels7 != NULL){free(levels7);}
    if (levels8 != NULL){free(levels8);}
*/
    if (levelPtr1 != NULL){free(levelPtr1);}
/*
    if (levelPtr2 != NULL){free(levelPtr2);}
    if (levelPtr3 != NULL){free(levelPtr3);}
    if (levelPtr4 != NULL){free(levelPtr4);}
    if (levelPtr5 != NULL){free(levelPtr5);}
    if (levelPtr6 != NULL){free(levelPtr6);}
    if (levelPtr7 != NULL){free(levelPtr7);}
    if (levelPtr8 != NULL){free(levelPtr8);}
*/
    if (ijkv1 != NULL){free(ijkv1);}
    if (ijkv2 != NULL){free(ijkv2);}
    if (ijkv3 != NULL){free(ijkv3);}
    if (ijkv4 != NULL){free(ijkv4);}
    if (ijkv5 != NULL){free(ijkv5);}
    if (ijkv6 != NULL){free(ijkv6);}
    if (ijkv7 != NULL){free(ijkv7);}
    if (ijkv8 != NULL){free(ijkv8);}
    levels1 = NULL;
/*
    levels2 = NULL; levels3 = NULL; levels4 = NULL;
    levels5 = NULL; levels6 = NULL; levels7 = NULL; levels8 = NULL;
*/
    levelPtr1 = NULL;
/*
    levelPtr2 = NULL; levelPtr3 = NULL; levelPtr4 = NULL;
    levelPtr5 = NULL; levelPtr6 = NULL; levelPtr7 = NULL; levelPtr8 = NULL;
*/
    ijkv1 = NULL; ijkv2 = NULL; ijkv3 = NULL; ijkv4 = NULL;
    ijkv5 = NULL; ijkv6 = NULL; ijkv7 = NULL; ijkv8 = NULL;
    return;
}
//============================================================================//
/*!
 * @brief Gets the number of grid points (nodes) in the graph.
 */
size_t fteikGraph::getNumberOfGridPoints(void)
{
    return ngrd;
}
//============================================================================//
bool fteikGraph::isInitialized(void)
{
    return linit;
}
//============================================================================//
/*!
 * @brief Gets the number of levels in the graph.
 */
int fteikGraph::getNumberOfLevels(void)
{
    return nLevels;
}
//============================================================================//
void fteikGraph::generateSweepOrderings(void)
{
    // Sweep 1
    y1Sweep1 = 1;  x1Sweep1 = 1;  z1Sweep1 = 1; 
    y2Sweep1 = ny; x2Sweep1 = nx; z2Sweep1 = nz;
    // Sweep 2
    y1Sweep2 = 1;  x1Sweep2 = nx-2; z1Sweep2 = 1; 
    y2Sweep2 = ny; x2Sweep2 =-1;    z2Sweep2 = nz;
    // Sweep 3
    y1Sweep3 = ny-2; x1Sweep3 = 1;  z1Sweep3 = 1;
    y2Sweep3 =-1;    x2Sweep3 = nx; z2Sweep3 = nz;
    // Sweep 4
    y1Sweep4 = ny-2; x1Sweep4 = nx-2; z1Sweep4 = 1;
    y2Sweep4 =-1;    x2Sweep4 =-1;    z2Sweep4 = nz;
    // Sweep 5
    y1Sweep5 = 1;  x1Sweep5 = 1;  z1Sweep5 = nz-2;
    y2Sweep5 = ny; x2Sweep5 = nx; z2Sweep5 =-1;
    // Sweep 6
    y1Sweep6 = 1;  x1Sweep6 = nx-2; z1Sweep6 = nz-2;
    y2Sweep6 = ny; x2Sweep6 =-1;    z2Sweep6 =-1;
    // Sweep 7
    y1Sweep7 = ny-2; x1Sweep7 = 1;  z1Sweep7 = nz-2;
    y2Sweep7 =-1;    x2Sweep7 = nx; z2Sweep7 =-1;
    // Sweep 8
    y1Sweep8 = ny-2; x1Sweep8 = nx-2; z1Sweep8 = nz-2;
    y2Sweep8 =-1;    x2Sweep8 =-1;    z2Sweep8 =-1;
    return;
}

/* TODO delete */
void fteikGraph::setSweepSigns(void)
{
    sweep2sgnx64f(1,
                  &sgntzSweep1, &sgntxSweep1, &sgntySweep1,
                  &sgnvzSweep1, &sgnvxSweep1, &sgnvySweep1,
                  &sgnrzSweep1, &sgnrxSweep1, &sgnrySweep1);
    sweep2sgnx64f(2,
                  &sgntzSweep2, &sgntxSweep2, &sgntySweep2,
                  &sgnvzSweep2, &sgnvxSweep2, &sgnvySweep2,
                  &sgnrzSweep2, &sgnrxSweep2, &sgnrySweep2);
    sweep2sgnx64f(3,
                  &sgntzSweep3, &sgntxSweep3, &sgntySweep3,
                  &sgnvzSweep3, &sgnvxSweep3, &sgnvySweep3,
                  &sgnrzSweep3, &sgnrxSweep3, &sgnrySweep3);
    sweep2sgnx64f(4,
                  &sgntzSweep4, &sgntxSweep4, &sgntySweep4,
                  &sgnvzSweep4, &sgnvxSweep4, &sgnvySweep4,
                  &sgnrzSweep4, &sgnrxSweep4, &sgnrySweep4);
    sweep2sgnx64f(5,
                  &sgntzSweep5, &sgntxSweep5, &sgntySweep5,
                  &sgnvzSweep5, &sgnvxSweep5, &sgnvySweep5,
                  &sgnrzSweep5, &sgnrxSweep5, &sgnrySweep5);
    sweep2sgnx64f(6,
                  &sgntzSweep6, &sgntxSweep6, &sgntySweep6,
                  &sgnvzSweep6, &sgnvxSweep6, &sgnvySweep6,
                  &sgnrzSweep6, &sgnrxSweep6, &sgnrySweep6);
    sweep2sgnx64f(7,
                  &sgntzSweep7, &sgntxSweep7, &sgntySweep7,
                  &sgnvzSweep7, &sgnvxSweep7, &sgnvySweep7,
                  &sgnrzSweep7, &sgnrxSweep7, &sgnrySweep7);
    sweep2sgnx64f(8,
                  &sgntzSweep8, &sgntxSweep8, &sgntySweep8,
                  &sgnvzSweep8, &sgnvxSweep8, &sgnvySweep8,
                  &sgnrzSweep8, &sgnrxSweep8, &sgnrySweep8);
}
//============================================================================//
void fteikGraph::sweep2sgnx64f(const int sweep,
                               int *sgntz, int *sgntx, int *sgnty,
                               int *sgnvz, int *sgnvx, int *sgnvy,
                               double *sgnrz, double *sgnrx, double *sgnry)
{
    const int sgntzv[8] = {1, 1, 1, 1,-1,-1,-1,-1};
    const int sgntxv[8] = {1,-1, 1,-1, 1,-1, 1,-1};
    const int sgntyv[8] = {1, 1,-1,-1, 1, 1,-1,-1};
    const int sgnvzv[8] = {1, 1, 1, 1, 0, 0, 0, 0};
    const int sgnvxv[8] = {1, 0, 1, 0, 1, 0, 1, 0};
    const int sgnvyv[8] = {1, 1, 0, 0, 1, 1, 0, 0};
    int k;
    k = sweep - 1;
    *sgntz = sgntzv[k];
    *sgntx = sgntxv[k];
    *sgnty = sgntyv[k];

    *sgnvz = sgnvzv[k];
    *sgnvx = sgnvxv[k];
    *sgnvy = sgnvyv[k];

    *sgnrz = static_cast<double> (*sgntz);
    *sgnrx = static_cast<double> (*sgntx);
    *sgnry = static_cast<double> (*sgnty);
    return;
}

/*
int fteikGraph::generateLevelPointers(const int nLevels,
                                      const int *__restrict__ levels,
                                      int *__restrict__ levelPtr,
                                      int *__restrict__ ijkv)
{
    int *work;
    work = (int *) aligned_alloc(64, ngrd*sizeof(int)); 
    return 0;
}
*/
//============================================================================//
/*!
 * @brief Returns the node to level map for the given sweep.
 *
 * @param[in] sweep    Sweep number.  This must be in the interval [1,8].
 *
 * @result Pointer that given the inode'th node returns the corresponding
 *         level.  The dimension of this array is [nz*nx*ny].
 * 
 */
int *fteikGraph::getLevelPointer(const int sweep)
{
    const char *fcnm = "getLevelPointer\0";
    int *levels = NULL;
    levels = levels1;
    return levels;
/*
    if (sweep == 1)
    {
        levels = levels1;
    }
    else if (sweep == 2)
    {
        levels = levels2;
    }
    else if (sweep == 3)
    {
        levels = levels3; 
    }
    else if (sweep == 4)
    {
        levels = levels4;
    }
    else if (sweep == 5)
    {
        levels = levels5;
    }
    else if (sweep == 6)
    {   
        levels = levels6;
    }
    else if (sweep == 7)
    {
        levels = levels7; 
    }   
    else if (sweep == 8)
    {   
        levels = levels8;
    }
    else
    {
        printf("%s: Error invalid sweep number: %d\n", fcnm, sweep);
    } 
*/
    return levels;
}
//============================================================================//
/*!
 * @brief Returns a pointer to the IJKV pointer.
 *
 */
int *fteikGraph::getIJKVPointer(const int sweep)
{
    const char *fcnm = "getIJKVPointer\0";
    int *ijkv = NULL;
    if (sweep == 1)
    {   
        ijkv = ijkv1;
    }   
    else if (sweep == 2)
    {   
        ijkv = ijkv2;
    }   
    else if (sweep == 3)
    {   
        ijkv = ijkv3; 
    }   
    else if (sweep == 4)
    {   
        ijkv = ijkv4;
    }   
    else if (sweep == 5)
    {   
        ijkv = ijkv5;
    }   
    else if (sweep == 6)
    {   
        ijkv = ijkv6;
    }   
    else if (sweep == 7)
    {   
        ijkv = ijkv7; 
    }
    else if (sweep == 8)
    {
        ijkv = ijkv8;
    }
    else
    {   
        printf("%s: Error invalid sweep number: %d\n", fcnm, sweep);
    }   
    return ijkv;
}

int fteikGraph::copyLevelPtr(const int nwork, const int sweep,
                             int *__restrict__ levelPtrOut)
{
    const char *fcnm = "copyLevelPtr\0";
    int *levelPtr;
    if (nwork < nLevels || levelPtrOut == NULL)
    {
        if (nwork < nLevels + 1)
        {
            printf("%s: Insufficient space for output\n", fcnm);
        }
        if (levelPtrOut == NULL){printf("%s: levelPtrOut is NULL\n", fcnm);}
        return -1;
    }
    memset(levelPtrOut, -1, static_cast<size_t> (nwork)*sizeof(int)); 
    levelPtr = getLevelPtrPointer(sweep);
    if (levelPtr == NULL)
    {
        for (int i=0; i<nwork; i++){levelPtrOut[i] =-1;}
        printf("%s: Error getting pointer\n", fcnm);
        return -1;
    }
    #pragma omp simd
    for (int i=0; i<nLevels + 1; i++){levelPtrOut[i] = levelPtr[i];}
    levelPtr = NULL;
    return 0; 
}
//============================================================================//
/*!
 * @brief Copies the IJKV vector to ijkvOut.
 *
 * @param[in] nwork     Max space for ijkvOut.  This must be at least
 *                      [4*nx*ny*nz].
 * @param[in] sweep     Sweep number.  This must be in the range [1,8].
 *
 * @param[out] ijkvOut  This is an array that maps from the n'th node
 *                      in the level structure fo the (iz, ix, iy, meshNode).
 *                      This is an array of dimension [nwork].
 *
 * @result 0 indicates success.
 *
 */
int fteikGraph::copyIJKV(const int nwork, const int sweep,
                         int *__restrict__ ijkvOut)
{
    const char *fcnm = "copyIJKV\0";
    int *ijkv;
    if (nwork < 4*ngrd || ijkvOut == NULL)
    {
        if (nwork < 4*ngrd)
        {
            printf("%s: Insufficient space for output\n", fcnm);
        }
        if (ijkvOut == NULL){printf("%s: Error ijkvOut is NULL\n", fcnm);}
        return -1;
    }
    memset(ijkvOut, -1, static_cast<size_t> (nwork)*sizeof(int));
    ijkv = getIJKVPointer(sweep);
    if (ijkv == NULL)
    {
        for (int i=0; i<4*nwork; i++){ijkvOut[i] =-1;}
        printf("%s: Error getting pointer\n", fcnm);
        return -1;
    }
    #pragma omp simd
    for (int i=0; i<4*ngrd; i++){ijkvOut[i] = ijkv[i];}
    ijkv = NULL;
    return 0;
}
//============================================================================//
/*!
 * @brief Utility routine for allocating the level pointers.
 */
void fteikGraph::allocateLevelPtrs(void)
{
    size_t nbytes = static_cast<size_t> (nLevels+1)*sizeof(int);
    if (levelPtr1 == NULL){levelPtr1 = (int *) aligned_alloc(64, nbytes);}
/*
    if (levelPtr2 == NULL){levelPtr2 = (int *) aligned_alloc(64, nbytes);}
    if (levelPtr3 == NULL){levelPtr3 = (int *) aligned_alloc(64, nbytes);}
    if (levelPtr4 == NULL){levelPtr4 = (int *) aligned_alloc(64, nbytes);}
    if (levelPtr5 == NULL){levelPtr5 = (int *) aligned_alloc(64, nbytes);}
    if (levelPtr6 == NULL){levelPtr6 = (int *) aligned_alloc(64, nbytes);}
    if (levelPtr7 == NULL){levelPtr7 = (int *) aligned_alloc(64, nbytes);}
    if (levelPtr8 == NULL){levelPtr8 = (int *) aligned_alloc(64, nbytes);}
*/
    return;
}
//============================================================================//
/*!
 * @brief Gets a pointer to the levelPtr array for the given sweep.
 *
 * @param[in] sweep     Sweep number.  This must be in the range [1,8].
 *
 * @result Pointer to levelPtr in given sweep.  If NULL then the sweep number 
 *         was invalid.
 *
 */
int *fteikGraph::getLevelPtrPointer(const int sweep)
{
    const char *fcnm = "getLevelPtrPointer\0";
    int *levelPtr = NULL;
    levelPtr = levelPtr1;
    return levelPtr;
/*
    if (sweep == 1)
    {   
        levelPtr = levelPtr1;
    }
    else if (sweep == 2)
    {   
        levelPtr = levelPtr2;
    }
    else if (sweep == 3)
    {
        levelPtr = levelPtr3;
    }
    else if (sweep == 4)
    {   
        levelPtr = levelPtr4;
    }
    else if (sweep == 5)
    {
        levelPtr = levelPtr5;
    }
    else if (sweep == 6)
    {
        levelPtr = levelPtr6;
    }
    else if (sweep == 7)
    {
        levelPtr = levelPtr7;
    }    
    else if (sweep == 8)
    {
        levelPtr = levelPtr8;
    }
    else
    {
        printf("%s: Error invalid sweep number: %d\n", fcnm, sweep);
    } 
    return levelPtr;
*/
}
//============================================================================//
/*!
 * @brief Computes the level that each node belongs to for the given sweep.
 *
 * @param[in] sweep     Sweep number.  This is must be in the interval [1,8].
 *
 * @result 0 indicates success.
 *
 */
int fteikGraph::computeNodeLevels(const int sweep)
{
    const char *fcnm = "computeNodeLevels\0";
    double sgnrx, sgnry, sgnrz;
    int *ijkv, *levels, *levelPtr, *nbinsAll, *work;
    int iloc[8], ierr, grd1, grd2, ix, iy, iz, maxWork,
        nLevels0, sgntx, sgnty, sgntz, sgnvx, sgnvy, sgnvz;
    typedef adjacency_list<vecS, vecS, directedS,
                           property<vertex_color_t, default_color_type> > Graph;
    typedef graph_traits<Graph>::vertex_descriptor Vertex;
/*
    typedef adjacency_list<vecS, vecS, directedS, 
                           property<vertex_color_t, default_color_type,
                           property<vertex_degree_t,int> > > Graph;
*/
    typedef std::pair<int, int> Pair;
    nLevels0 = nLevels;
    nLevels = 0;
    if (!linit)
    {
        printf("%s: Error fteikGraph not initialized\n", fcnm);
        return -1;
    }
    levels = getLevelPointer(sweep);
    if (levels == NULL)
    {
        printf("%s: Error getting level pointer\n", fcnm);
        return -1;
    }
    memset(levels, -1, ngrd*sizeof(int));
    // Figure out sweep ordering
    sweep2sgnx64f(sweep,
                  &sgntz, &sgntx, &sgnty,
                  &sgnvz, &sgnvx, &sgnvy,
                  &sgnrz, &sgnrx, &sgnry);
    size_t maxEdges = ngrd*7;
    Pair *edges = (Pair *) calloc(maxEdges, sizeof(Pair)); //Pair edges[maxEdges];
    size_t edge = 0;
    vector<int> xadj(ngrd+1, 0);
    vector<int> adjncy(maxEdges, 0);
    int vtx = 0;
    xadj[vtx] = 0;
    //for (int iy = y1Sweep1; iy < y2Sweep1; iy++)
    for (iy=0; iy<ny; iy++)
    {
        //for (int ix = x1Sweep1; ix < x2Sweep1; ix++)
        for (ix=0; ix<nx; ix++)
        {
            //for (int iz = z1Sweep1; iz < z2Sweep1; iz++)
            for (iz=0; iz<nz; iz++)
            {
                bool laddz = false;
                bool laddx = false;
                bool laddy = false;
                if (iz-sgntz >= 0 && iz-sgntz < nz){laddz = true;}
                if (ix-sgntx >= 0 && ix-sgntx < nx){laddx = true;}
                if (iy-sgnty >= 0 && iy-sgnty < ny){laddy = true;}
                // Get connections to grd1
                grd1 = grid2index(iz, ix, iy);
                if (laddz)
                {
                    grd2 = grid2index(iz-sgntz, ix, iy);
                    edges[edge] = Pair(grd1, grd2);
                    adjncy[edge] = grd2;
                    edge = edge + 1;
                }
                if (laddx)
                {
                    grd2 = grid2index(iz, ix-sgntx, iy);
                    edges[edge] = Pair(grd1, grd2);
                    adjncy[edge] = grd2;
                    edge = edge + 1;
                }
                if (laddy)
                {
                    grd2 = grid2index(iz, ix, iy-sgnty);
                    edges[edge] = Pair(grd1, grd2); 
                    adjncy[edge] = grd2;
                    edge = edge + 1;
                }
                if (laddz && laddx)
                {
                    grd2 = grid2index(iz-sgntz, ix-sgntx, iy);
                    edges[edge] = Pair(grd1, grd2); 
                    adjncy[edge] = grd2;
                    edge = edge + 1;
                }
                if (laddx && laddy)
                {
                    grd2 = grid2index(iz, ix-sgntx, iy-sgnty);
                    edges[edge] = Pair(grd1, grd2); 
                    adjncy[edge] = grd2;
                    edge = edge + 1;
                }
                if (laddz && laddy)
                {
                    grd2 = grid2index(iz-sgntz, ix, iy-sgnty);
                    edges[edge] = Pair(grd1, grd2); 
                    adjncy[edge] = grd2;
                    edge = edge + 1;
                }
                if (laddz && laddx && laddy)
                {
                    grd2 = grid2index(iz-sgntz, ix-sgntx, iy-sgnty);
                    edges[edge] = Pair(grd1, grd2);
                    adjncy[edge] = grd2;
                    edge = edge + 1;
                }
                vtx = vtx + 1;
                xadj[vtx] = edge;
                int nsort = xadj[vtx] - xadj[vtx-1];
                if (nsort > 1)
                {
                    sorting_sort32i_work(nsort, &adjncy[xadj[vtx-1]],
                                         SORT_ASCENDING);
                }
            }
            if (edge > maxEdges)
            {
                printf("%s: Invalid space estimate\n", fcnm);
                return -1;
            }
        }
    }
    Graph G(ngrd); 
    for (int i=0; i<edge; i++)
    {
        //cout << edges[i].first << " " <<  edges[i].second << endl;
        add_edge(edges[i].first, edges[i].second, G); 
    }
    typedef std::vector<Vertex> container;
    container c;
    topological_sort(G, std::back_inserter(c));
    nnz = ngrd + edge;
/*
cout << "write result" << endl;
    for (container::reverse_iterator ii=c.rbegin(); ii != c.rend(); ++ii)
    {
        cout << c[*ii] << " " ;
    }
cout << endl;
cout << "nedges: " << edge << endl;
cout << "ngrd: " << ngrd << " " << nx*ny*nz << endl;
*/
    free(edges);
    nLevels = 1;
    // Now take the topological sort and build up a dependency list
    container::iterator firstNode = c.begin();
    levels[*firstNode] = nLevels;
    for (container::iterator kk=c.begin(); kk != c.end(); ++kk)
    {
        if (levels[*kk] !=-1){continue;}
        // Find highest level in this row 
        int maxLevel = 1;
        for (int j=xadj[*kk]; j<xadj[*kk+1]; j++)
        {
            int jcol = adjncy[j];
            maxLevel = MAX(maxLevel, levels[jcol]);
        }
        levels[*kk] = maxLevel + 1; 
        if (maxLevel + 1 > nLevels){nLevels = maxLevel + 1;}
    }
    // Ensure the number of levels is consistent with other sweeps
    if (nLevels0 == 0){nLevels0 = nLevels;}
    if (nLevels != nLevels0)
    {   
        printf("%s: Error - inconsistent number of levels in sweeps\n", fcnm);
        return -1;
    }
    // Ensure all levels have been visited
    ierr = 0;
    nbinsAll = (int *) calloc(static_cast<size_t> (nLevels), sizeof(int));
    for (int k=0; k<ngrd; k++)
    {
        if (levels[k] < 1 || levels[k] > nLevels)
        {
            printf("%s: Failed to initialize node: %d\n", fcnm, k);
            ierr = ierr + 1;
        }
        nbinsAll[levels[k]-1] = nbinsAll[levels[k]-1] + 1;
    }
    // Build the level sets pointer
    allocateLevelPtrs();
    levelPtr = getLevelPtrPointer(sweep);
    levelPtr[0] = 0;
    maxWork = 0;
    for (int k=0; k<nLevels; k++)
    {
        levelPtr[k+1] = levelPtr[k] + nbinsAll[k];
        maxWork = MAX(maxWork, levelPtr[k+1] - levelPtr[k]);
    }
    if (verbose > 1 || (verbose == 1 && sweep == 1))
    {
        for (int k=0; k<nLevels; k++)
        {
            printf("%s: There are %d nodes in level %d for sweep %d\n",
                   fcnm, nbinsAll[k], k+1, sweep);
        }
    }
    maxLevelSize = maxWork;
    free(nbinsAll);
    // Build the (i, j, k, node) arrays
    generateIJKV(1);
    generateIJKV(2);
    generateIJKV(3);
    generateIJKV(4);
    generateIJKV(5);
    generateIJKV(6);
    generateIJKV(7);
    generateIJKV(8);
/*
    // Build the (i,j,k,node) 
    ierr = 0;
    size_t gotEmAll = 0;
    work = (int *) calloc(static_cast<size_t> (maxWork), sizeof(int));
    ijkv = getIJKVPointer(sweep);
    for (int il=1; il <= nLevels; il++)
    {
        // Get the nodes in this level
        int j = 0;
        for (int i=0; i<ngrd; i++)
        {
            if (levels[i] == il)
            {
                work[j] = i;
                j = j + 1;
            }
        }
        // Sort the work array 
        int nsort = j;
        if (levelPtr[il] - levelPtr[il-1] != nsort)
        {
            printf("%s: levelPtr not set correctly %d %d\n",
                   fcnm, levelPtr[il] - levelPtr[il-1], nsort);
        }
        if (nsort > 1)
        {
            sorting_sort32i_work(nsort, work, SORT_ASCENDING);
        }
        // Put into the array
        for (int j=0; j<nsort; j++)
        {
            int indx = levelPtr[il-1] + j;
            ierr += index2grid(work[j], &iz, &ix, &iy);
            ijkv[4*indx+0] = iz;
            ijkv[4*indx+1] = ix; 
            ijkv[4*indx+2] = iy; 
            ijkv[4*indx+3] = work[j];
            gotEmAll = gotEmAll + 1;
        }
    }
    if (ierr != 0)
    {
        printf("%s: Failed to convert index to grid\n", fcnm);
        return -1;
    }
    if (gotEmAll != ngrd)
    {
        printf("%s: Failed to initialize all nodes\n", fcnm);
        return -1;
    }
    //free(edges);
    free(work);
printf("i think this is the spot that's broken (4,1,3) depnds on (4,2,3)\n");
*/
//getchar();
/*
    // Compute the level pointer
    for (container::iterator ii=c.begin(); ii != c.end(); ii++)
    {
       cout << c[*ii] << " " << levels[*ii] << " " << endl;
    }
    cout << "Number of levels " << nLevels << endl;
*/
    return 0;
}
//============================================================================//
int fteikGraph::generateIJKV(const int sweep)
{
    int *work, *ijkv, i, ierr, indx, il, ix, iy, iz, j, jx, jy, jz, nsort;
    size_t gotEmAll;
    // Build the (i,j,k,node) 
    ierr = 0;
    gotEmAll = 0;
    work = (int *) calloc(static_cast<size_t> (maxLevelSize), sizeof(int));
    ijkv = getIJKVPointer(sweep);
    for (il=1; il <= nLevels; il++)
    {
        // Get the nodes in this level
        j = 0;
        for (i=0; i<ngrd; i++)
        {
            if (levels1[i] == il) 
            {
                // Permute the grid point from the first sweep
                ierr = index2grid(i, &iz, &ix, &iy);
                if (sweep == 1)
                {
                    jz = iz;
                    jx = ix;
                    jy = iy;
                }
                else if (sweep == 2)
                {
                    jz = iz;
                    jx = nx - 1 - ix;
                    jy = iy;
                }
                else if (sweep == 3) 
                {
                    jz = iz;
                    jx = ix;
                    jy = ny - 1 - iy;
                }
                else if (sweep == 4)
                {
                    jz = iz;
                    jx = nx - 1 - ix;
                    jy = ny - 1 - iy;
                }
                else if (sweep == 5)
                {    
                    jz = nz - 1 - iz;
                    jx = ix;
                    jy = iy;
                }    
                else if (sweep == 6)
                {    
                    jz = nz - 1 - iz;
                    jx = nx - 1 - ix;
                    jy = iy;
                }    
                else if (sweep == 7)
                {
                    jz = nz - 1 - iz;
                    jx = ix;
                    jy = ny - 1 - iy;
                }
                else //if (sweep == 8)
                {
                    jz = nz - 1 - iz;
                    jx = nx - 1 - ix;
                    jy = ny - 1 - iy;
                }
                work[j] = grid2index(jz, jx, jy);
                j = j + 1;
            }   
        }
        // Sort the work array 
        nsort = j;
        if (levelPtr1[il] - levelPtr1[il-1] != nsort)
        {
            printf("%s: levelPtr not set correctly %d %d\n",
                   __func__, levelPtr1[il] - levelPtr1[il-1], nsort);
        }
        if (nsort > 1)
        {
            sorting_sort32i_work(nsort, work, SORT_ASCENDING);
        }
        // Put into the array
        for (j=0; j<nsort; j++)
        {
            indx = levelPtr1[il-1] + j;
            ierr += index2grid(work[j], &iz, &ix, &iy);
            ijkv[4*indx+0] = iz;
            ijkv[4*indx+1] = ix;
            ijkv[4*indx+2] = iy;
            ijkv[4*indx+3] = work[j];
            gotEmAll = gotEmAll + 1;
        }
    }
    free(work);
    if (ierr != 0)
    {
        printf("%s: Failed to convert index to grid\n", __func__);
        return -1;
    }
    if (gotEmAll != ngrd)
    {
        printf("%s: Failed to initialize all nodes\n", __func__);
        return -1;
    }
    return 0;
}
//============================================================================//
extern "C" void *fteik_graph_initializeF(const int *nz, const int *nx,
                                         const int *ny, int *ierr)
{
    void *graph;
    graph = fteik_graph_initialize(*nz, *nx, *ny, ierr);
    return graph;
}

extern "C" void *fteik_graph_initialize(const int nz, const int nx, const int ny,
                                        int *ierr)
{
    const char *fcnm = "fteik_graphIntialize\0";
    fteikGraph *graph = NULL; 
    *ierr = 0;
    if (nz < 3 || nx < 3 || ny < 3)
    {
        if (nz < 3){printf("%s: Error nz=%d must be at least 3\n", fcnm, nz);}
        if (nx < 3){printf("%s: Error nx=%d must be at least 3\n", fcnm, nx);}
        if (ny < 3){printf("%s: Error ny=%d must be at least 3\n", fcnm, ny);}
        *ierr = 1;
        return graph;
    }
    graph = new fteikGraph(nz, nx, ny);
    return static_cast<void *> (graph);
}

extern "C" void fteik_graph_finalize(void *graphIn)
{
    fteikGraph *graph = NULL;
    graph = static_cast<fteikGraph *> (graphIn); 
    graph->finalize();
    delete graph;
    return;
}

extern "C" int fteik_graph_getNumberOfLevels(void *graphIn)
{
    int nLevels;
    fteikGraph *graph = NULL;
    graph = static_cast<fteikGraph *> (graphIn);
    nLevels = graph->getNumberOfLevels();
    graph = NULL;
    return nLevels;
}

extern "C" int fteik_graph_getNumberOfGridPoints(void *graphIn)
{
    int ngrd;
    fteikGraph *graph = static_cast<fteikGraph *> (graphIn);
    ngrd = static_cast<int> (graph->getNumberOfGridPoints());
    graph = NULL;
    return ngrd;
}

extern "C" int fteik_graph_getLevelPointerF(void *graphIn, const int *nwork,
                                            const int *sweep,
                                            int *__restrict__ levelPtr)
{
    int i, ierr;
    ierr = fteik_graph_getLevelPointer(graphIn, *nwork, *sweep, levelPtr);
    if (ierr == 0)
    {
        #pragma omp simd
        for (i=0; i<*nwork; i++){levelPtr[i] = levelPtr[i] + 1;} // C 2 Fortran
    }
    else
    {   
        printf("%s: Error copying levelPtr\n", __func__);
    }   
    return ierr;       
}

extern "C" int fteik_graph_getLevelPointer(void *graphIn, const int nwork,
                                           const int sweep,
                                           int *__restrict__ levelPtr)
{
    int ierr;
    fteikGraph *graph = static_cast<fteikGraph *> (graphIn);
    ierr = graph->copyLevelPtr(nwork, sweep, levelPtr);
    if (ierr != 0){printf("%s: Error getting levelPtr\n", __func__);}
    graph = NULL;
    return ierr;
}

extern "C" int fteik_graph_getIJKVF(void *graphIn, const int *nwork,
                                    const int *sweep, int *__restrict__ ijkv)
{
    const char *fcnm = "fteik_graph_getIJKVF\0";
    int i, ierr;
    ierr = fteik_graph_getIJKV(graphIn, *nwork, *sweep, ijkv);
    if (ierr == 0)
    {
        #pragma omp simd
        for (i=0; i<*nwork; i++){ijkv[i] = ijkv[i] + 1;} // Fortran numbering
    }
    else
    {
        printf("%s: Error copying ijkv\n", fcnm);
    }
    return ierr;
}

extern "C" int fteik_graph_getIJKV(void *graphIn, const int nwork,
                                   const int sweep, int *__restrict__ ijkv)
{
    const char *fcnm = "fteik_graph_getIJKV\0";
    int ierr;
    fteikGraph *graph = NULL;
    graph = static_cast<fteikGraph *> (graphIn);
    ierr = graph->copyIJKV(nwork, sweep, ijkv);
    if (ierr != 0){printf("%s: Error getting ijkv\n", fcnm);}
    graph = NULL; 
    return ierr;
}

extern "C" int fteik_graph_getMaxLevelSize(void *graphIn)
{
    int maxLevelSize;
    fteikGraph *graph = static_cast<fteikGraph *> (graphIn);
    maxLevelSize = graph->getMaxLevelSize();
    graph = NULL;
    return maxLevelSize;
}
