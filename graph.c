#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include "sorting.h"
#include "fteik_graph.h"

static int grd2indx(const int i, const int j, const int k,
                    const int n1, const int n2, const int n3);
static int grd2ijk(const int igrd,
                   const int n1, const int n2, const int n3, 
                   int *i, int *j, int *k);

/*!
 * @brief This is the driver function that creates the level set
 *        structures.  This significance of these structures is that
 *        the update at each point in a level is independent of all
 *        other points in that level and is therefore safe to paralellize
 *        as there are no data races.  The algorithm: \n
 *         (1) Defines the graph for the first forward sweep
 *             (increases in z, x, and y).  This graph corresponds
 *             to a lower triangular of a sparse matrix. \n
 *         (2) From the lower triangle sparse graph generates a 
 *             level set structure (e.g., Saad Iterative Methods for 
 *             Sparse Matrices Ed. 2. pg. 389).
 *         (3) Finally, it permutes the level set structure for the
 *             first sweep to all the other sweeps in a fashion that
 *             is consistent with the Noble et. al. (2014).
 * 
 * @param[in] nz          Number of z points in travel time grid (> 2).
 * @param[in] nx          Number of x points in travel time grid (> 2).
 * @param[in] ny          Number of y points in travel time grid (> 2).
 *
 * @param[out] levelSets  Structure containing the level set structure
 *                        for all sweeps.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker
 *
 * @copyright Ben Baker distributed under the MIT license.
 *
 */
int graph_createLevelStructs(const int nz, const int nx, const int ny,
                             struct levelSets_struct *levelSets)
{
    const char *fcnm = "graph_createLevelStructs\0";
    int *xadj, *adjncy;
    int i, ierr, k, n, sweep;
    if (nz < 3 || nx < 3 || ny < 3)
    {
        printf("%s: Error (nx,ny,nz) must all be at least 2\n", fcnm);
        return -1;
    }
    n = nz*nx*ny;
    xadj = NULL;
    adjncy = NULL;
    memset(levelSets, 0, sizeof(struct levelSets_struct));
    levelSets->levelSet = (struct levelSet_struct *)
                          calloc(8, sizeof(struct levelSet_struct));
    levelSets->nsweeps = 8;
    // Set the chunk size for the finite difference
    for (k=0; k<8; k++)
    {
        levelSets->levelSet[k].chunkSize = DEFAULT_FTEIK_FD_CHUNK_SIZE;
    }
    // Define the graph structure
    ierr = graph_createGraph3D(nz, nx, ny, &xadj, &adjncy);
    if (ierr != 0)
    {
        printf("%s: Error generating graph\n", fcnm);
        return -1;
    }
    // Define the original level structure
    levelSets->levelSet[0].nnodes = n;
    levelSets->levelSet[0].nx = nx;
    levelSets->levelSet[0].ny = ny;
    levelSets->levelSet[0].nz = nz;
    ierr = graph_createLevelStructure(n, xadj, adjncy,
                                      &levelSets->levelSet[0].nLevels,
                                      &levelSets->levelSet[0].levelPtr,
                                      &levelSets->levelSet[0].ijkv, 
                                      &levelSets->levelSet[0].n2l);
    if (ierr != 0)
    {
        printf("%s: Error defining level set structure!\n", fcnm);
        return -1;
    }
    free(xadj);
    free(adjncy);
    // Now permute the level set structures to conform to the
    // various sweep orderings
    for (k=1; k<8; k++)
    {
        levelSets->levelSet[k].nLevels = levelSets->levelSet[0].nLevels;
        levelSets->levelSet[k].nnodes = levelSets->levelSet[0].nnodes;
        levelSets->levelSet[k].nx = levelSets->levelSet[0].nx;
        levelSets->levelSet[k].ny = levelSets->levelSet[0].ny;
        levelSets->levelSet[k].nz = levelSets->levelSet[0].nz;
        sweep = k + 1;
        ierr = graph_permuteLevelStructure(nz, nx, ny,
                                          sweep, levelSets->levelSet[k].nLevels,
                                          levelSets->levelSet[0].n2l,
                                          &levelSets->levelSet[k].ijkv,
                                          &levelSets->levelSet[k].levelPtr,
                                          &levelSets->levelSet[k].n2l);
/*
printf("%d %d %d %d\n", levelSets->levelSet[k].levelPtr[1],
                        levelSets->levelSet[k].levelPtr[2],
                        levelSets->levelSet[k].levelPtr[3],
                        levelSets->levelSet[k].levelPtr[4]);
*/
        if (ierr != 0)
        {
            printf("%s: Error permuting sweep %d\n", fcnm, sweep);
            return -1;
        }
    }
    // Convert the nodes to (i,j,k) indices
    for (k=0; k<8; k++)
    {
        levelSets->levelSet[k].n2ijk
             = (int *) aligned_alloc(64, (size_t) (3*n)*sizeof(int)); 
        for (int level=0; level<levelSets->levelSet[k].nLevels; level++)
        {
            int nodeStart = levelSets->levelSet[k].levelPtr[level];
            int nodeEnd   = levelSets->levelSet[k].levelPtr[level+1];
            for (i=nodeStart; i<nodeEnd; i++)
            {
                grd2ijk(levelSets->levelSet[k].ijkv[i],
                        nz, nx, ny,
                        &levelSets->levelSet[k].n2ijk[3*i+0],
                        &levelSets->levelSet[k].n2ijk[3*i+1],
                        &levelSets->levelSet[k].n2ijk[3*i+2]);
/*
            if (k == 2){printf("%d %d %d %d %d\n", 
                               level+1,
                               levelSets->levelSet[k].n2ijk[3*i+0],
                               levelSets->levelSet[k].n2ijk[3*i+1],
                               levelSets->levelSet[k].n2ijk[3*i+2],
                               levelSets->levelSet[k].ijkv[i]);}
*/
            }
        }
//printf("\n");
/*
        for (i=0; i<levelSets->levelSet[k].nnodes; i++)
        {
            grd2ijk(levelSets->levelSet[k].ijkv[i],
                    nz, nx, ny,
                    &levelSets->levelSet[k].n2ijk[3*i+0],
                    &levelSets->levelSet[k].n2ijk[3*i+1],
                    &levelSets->levelSet[k].n2ijk[3*i+2]);
        }
*/
    }
    // Create the nodes that cannot be updated and sweep signs
    for (k=0; k<8; k++)
    {
        sweep = k + 1;
        levelSets->levelSet[k].sweep = sweep;
        ierr = graph_createUpdateMask(sweep, &levelSets->levelSet[k]);
        if (ierr != 0)
        {
            printf("%s: Error creating update mask\n", fcnm);
            return -1;
        }
    }
    // Create the masks for the nodes that cannot be udpated
    // in the first sweep because of the source location 
    for (k=0; k<8; k++)
    {
         graph_createTravelTimeAccessPointers(&levelSets->levelSet[k]);
    }
    // May use this later
    for (k=0; k<8; k++)
    {
        levelSets->levelSet[k].lupdateInit
           = (bool *) aligned_alloc(64, (size_t) n*sizeof(bool));
    }
    return 0;
}
//============================================================================//
/*!
 */
int graph_createTravelTimeAccessPointers(struct levelSet_struct *levelSet )
{
    const char *fcnm = "graph_createTravelTimeAccessPointers\0";
    double sgnrx, sgnry, sgnrz;
    size_t nwork;
    int *idest, *iloc, *jdest, *perm, *permLocal,
        chunk, i, inode, j, jnode, k, knode, level, maxNodesInLevel,
        node, nnode, nnodesSort, nodeEnd, nodeStart, nsort, nx, ny, nz,
        sgntx, sgnty, sgntz, sgnvx, sgnvy, sgnvz, sweep;
    sweep = levelSet->sweep;
    chunk = levelSet->chunkSize;
    nz = levelSet->nz;
    nx = levelSet->nx;
    ny = levelSet->ny;
    graph_sweep2sgnx64f(sweep,
                        &sgntz, &sgntx, &sgnty,
                        &sgnvz, &sgnvx, &sgnvy,
                        &sgnrz, &sgnrx, &sgnry);
    maxNodesInLevel = 0;
    for (level=0; level<levelSet->nLevels; level++)
    {
        maxNodesInLevel = MAX(levelSet->levelPtr[level+1]
                             -levelSet->levelPtr[level], maxNodesInLevel);
    }
    nwork = (size_t) (8*levelSet->nnodes);
    levelSet->ttPerm = (int *) calloc(nwork, sizeof(int));
    levelSet->ttLocPerm = (int16_t *) calloc(nwork, sizeof(int16_t));
    for (i=0; i<(int) nwork; i++)
    {
        levelSet->ttPerm[i] =-1;
        levelSet->ttLocPerm[i] =-1;
    }
    // Set the workspace
    iloc = (int *) calloc((size_t) (8*chunk), sizeof(int)); 
    perm = (int *) calloc((size_t) (8*chunk), sizeof(int));
    permLocal = (int *) calloc((size_t) (8*chunk), sizeof(int));
    idest = (int *) calloc((size_t) (8*chunk), sizeof(int));
    jdest = (int *) calloc((size_t) (8*chunk), sizeof(int));
    for (level=0; level<levelSet->nLevels; level++)
    {
        nodeStart = levelSet->levelPtr[level];
        nodeEnd   = levelSet->levelPtr[level+1];
        for (node=nodeStart; node<nodeEnd; node=node+chunk)
        {
            // Pre-compute the travel time global memory accesses
            nnode = MIN(nodeEnd, node+chunk - node);
            knode = 0;
            for (jnode=node; jnode<MIN(nodeEnd, node+chunk); jnode++)
            {
                if (levelSet->lupdate[jnode])
                {
                    i = levelSet->n2ijk[3*jnode+0];
                    j = levelSet->n2ijk[3*jnode+1];
                    k = levelSet->n2ijk[3*jnode+2];
                    iloc[8*knode+0] = grd2indx(i-sgntz, j,       k,
                                               nz, nx, ny);
                    iloc[8*knode+1] = grd2indx(i,       j-sgntx, k,
                                               nz, nx, ny);
                    iloc[8*knode+2] = grd2indx(i,       j,       k-sgnty,
                                               nz, nx, ny);
                    iloc[8*knode+3] = grd2indx(i-sgntz, j-sgntx, k,
                                               nz, nx, ny);
                    iloc[8*knode+4] = grd2indx(i,       j-sgntx, k-sgnty,
                                               nz, nx, ny);
                    iloc[8*knode+5] = grd2indx(i-sgntz, j,       k-sgnty,
                                               nz, nx, ny);
                    iloc[8*knode+6] = grd2indx(i-sgntz, j-sgntx, k-sgnty,
                                               nz, nx, ny);
                    iloc[8*knode+7] = grd2indx(i,       j,       k,
                                               nz, nx, ny);
                    for (i=0; i<8; i++){idest[8*knode+i] = 8*knode + i;}
                    knode = knode + 1;
                }
            }
            // Sort the travel-time array memory accesses in this chunk
            nsort = knode;
            nnodesSort = knode;
            if (nnodesSort > 0)
            {
                //printf("nsort: %d\n", nsort);
                sorting_argsort32i_work(8*nnodesSort, iloc,
                                        SORT_ASCENDING, perm);
                for (i=0; i<8*nnodesSort; i++)
                {
                    jdest[i] = idest[perm[i]];
                    //printf("%d %d\n", iloc[perm[i]], idest[perm[i]]);
                    if (iloc[perm[i]] < 0)
                    {
                        printf("%s: This node can't be valid\n", fcnm);
                        return -1;
                    }
                }
                sorting_argsort32i_work(8*nnodesSort, jdest,
                                        SORT_ASCENDING, permLocal);
                // Save nodes in manner to mimic accesses
                knode = 0;
                for (jnode=node; jnode<MIN(nodeEnd, node+chunk); jnode++)
                {
                    if (levelSet->lupdate[jnode])
                    {
                        for (i=0; i<8; i++)
                        {
                            levelSet->ttPerm[8*jnode+i] = iloc[perm[8*knode+i]];
                            levelSet->ttLocPerm[8*jnode+i] = permLocal[8*knode+i];
/*
                            printf("%d %d %d\n", iloc[perm[8*knode+i]], 
                                                 idest[perm[i]],
                                                 permLocal[8*knode+i]);
*/
                        }
                        knode = knode + 1;
                    }
                }
                if (knode != nnodesSort)
                {
                    printf("%s: Lost count somewhere\n", fcnm);
                    return -1;
                }
            }
        }
        //getchar();
    }
    free(iloc);
    free(idest);
    free(jdest);
    free(perm);
    free(permLocal);
    return 0;
}
//============================================================================//
/*!
 * @brief Frees the memory on the levelSet structure.
 *
 * @param[out] levelSet   On exit all memory has been freed and parameters
 *                        set to NULL or 0.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker
 *
 * @copyright Ben Baker distributed under the MIT license.
 *
 */
int graph_freeLevelStruct(struct levelSet_struct *levelSet)
{
    if (levelSet->ijkv != NULL){free(levelSet->ijkv);}
    if (levelSet->levelPtr != NULL){free(levelSet->levelPtr);}
    if (levelSet->n2ijk != NULL){free(levelSet->n2ijk);}
    if (levelSet->n2l != NULL){free(levelSet->n2l);}
    if (levelSet->lupdate != NULL){free(levelSet->lupdate);}
    if (levelSet->lupdateInit != NULL){free(levelSet->lupdateInit);}
    memset(levelSet, 0, sizeof(struct levelSet_struct));
    return 0;
}
//============================================================================//
/*!
 * @brief Frees the memory on the levelSets structure.
 *
 * @param[out] levelSets  On exit all memory has been freed and parameters
 *                        set to NULL or 0.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker
 *
 * @copyright Ben Baker distributed under the MIT license.
 *
 */
int graph_freeLevelsStruct(struct levelSets_struct *levelSets)
{
    int k;
    if (levelSets->nsweeps > 0 && levelSets->levelSet != NULL)
    {
        for (k=0; k<levelSets->nsweeps; k++)
        {
            graph_freeLevelStruct(&levelSets->levelSet[k]);
        }
        free(levelSets->levelSet);
    }
    memset(levelSets, 0, sizeof(struct levelSets_struct));
    return 0;
}
//============================================================================//
/*!
 * @brief Creates the graph for the sweep 1 which will be 
 *        increasing ny, increasing, nx, and increasing nz
 *        so that the sparse matrix is lower triangular.
 *
 * @param[in] nz          Number of z nodal points in travel time grid.
 * @param[in] nx          Number of x nodal points in travel time grid.
 * @param[in] ny          Number of y nodal points in travel time grid.
 *
 * @param[out] xadjOut    Maps from the i'th node to the start index
 *                        of adjncy [nz*nx*ny + 1].
 * @param[out] adjncyOut  Maps to the lower triangle non-zero columns
 *                        in a given row [xadjOut[nz*nx*ny]].
 *
 * @result 0 indicates success. 
 *
 * @author Ben Baker
 *
 * @copyright Ben Baker distributed under MIT license.
 *
 */
int graph_createGraph3D(const int nz, const int nx, const int ny,
                        int **xadjOut, int **adjncyOut)
{
    const char *fcnm = "graph_createGraph3D\0";
    int *xadj, *adjncy, *adjTrim;
    int i, i1, i2, id, ierr, iWall, j, j1, j2, jd, jndx, jWall,
        k, k1, k2, kd, kWall, mnodes, nadj, nnodes, nsort, nwork;
    // Set space
    printf("%s: Generating graph...\n", fcnm);
    mnodes = nz*nx*ny;
    nwork = 7*mnodes; // stencil is a box (but don't count update point) 
    xadj = (int *) calloc((size_t) (mnodes+1), sizeof(int));
    adjncy = (int *) calloc((size_t) nwork, sizeof(int)); 
    nadj = 0;
    // Set properties for sweep 1:
    //   top to bottom, west to east, and south to north
    kWall = 0;    // y plane not updated in sweep
    jWall = 0;    // x plane not updated in sweep
    iWall = 0;    // z plane not updated in sweep
    kd =-1;       // y sweeps north to south
    jd =-1;       // x sweeps west to east
    id =-1;       // z top to bottom 
    k1 = 0;       // y start
    k2 = ny;      // y end
    j1 = 0;       // x start
    j2 = nx;      // x end
    i1 = 0;       // z start
    i2 = nz;      // z end
    nnodes = 0;
    for (k=0; k<ny; k++)  //k=k1; k != k2; k=k-kd)
    {
        for (j=0; j<nx; j++) //j=j1; j != j2; j=j-jd)
        {
            for (i=0; i<nz; i++) //i=i1; i != i2; i=i-id)
            {
                // The finite difference stencil
                nnodes = nnodes + 1;
                if (i > iWall && j > jWall && k > kWall)
                {
                    adjncy[nadj+0] = grd2indx(i+id, j  ,  k,    nz, nx, ny);
                    adjncy[nadj+1] = grd2indx(i+id, j+jd, k,    nz, nx, ny);
                    adjncy[nadj+2] = grd2indx(i   , j+jd, k,    nz, nx, ny);
                    adjncy[nadj+3] = grd2indx(i+id, j   , k+kd, nz, nx, ny);
                    adjncy[nadj+4] = grd2indx(i   , j   , k+kd, nz, nx, ny);
                    adjncy[nadj+5] = grd2indx(i+id, j+jd, k+kd, nz, nx, ny);
                    adjncy[nadj+6] = grd2indx(i   , j+jd, k+kd, nz, nx, ny);
                    nadj = nadj + 7;
                }
                // Boundary condition
                else
                {
                    if (i != iWall)
                    {
                        adjncy[nadj] = grd2indx(i+id, j,    k,    nz, nx, ny); 
                        nadj = nadj + 1;
                    }
                    if (i != iWall && j != jWall)
                    {
                        adjncy[nadj] = grd2indx(i+id, j+jd, k,    nz, nx, ny);
                        nadj = nadj + 1;
                    }
                    if (j != jWall)
                    {
                        adjncy[nadj] = grd2indx(i,    j+jd, k,    nz, nx, ny);
                        nadj = nadj + 1;
                    }
                    if (i != iWall && k != kWall)
                    {
                        adjncy[nadj] = grd2indx(i+id, j,    k+kd, nz, nx, ny);
                        nadj = nadj + 1;
                    }
                    if (k != kWall)
                    {
                        adjncy[nadj] = grd2indx(i,    j,    k+kd, nz, nx, ny);
                        nadj = nadj + 1;
                    }
                    if (i != iWall && j != jWall && k != kWall)
                    {
                        adjncy[nadj] = grd2indx(i+id, j+jd, k+kd, nz, nx, ny);
                        nadj = nadj + 1;
                    } 
                    if (j != jWall && k != kWall)
                    {
                        adjncy[nadj] = grd2indx(i,    j+jd, k+kd, nz, nx, ny);
                        nadj = nadj + 1;
                    }
                } // End check on boundary conditions
                xadj[nnodes] = nadj;
                nsort = xadj[nnodes] - xadj[nnodes-1];
                jndx = xadj[nnodes-1];
                if (nsort > 1)
                {
                    ierr = sorting_sort32i_work(nsort, &adjncy[jndx],
                                                SORT_ASCENDING);
                    if (ierr != 0)
                    {
                        printf("%s: Error sorting array\n", fcnm);
                        return -1;
                    }
                }
            }
        }
    }
    if (nnodes != nz*nx*ny)
    {
        printf("%s: Error - missed some nodes\n", fcnm);
        return -1;
    }
    if (nadj > nwork)
    {
        printf("%s: Error - exceeded workspace\n", fcnm);
        return -1;
    }
    //printf("%d %d %d %d\n", nnodes, nz*nx*ny, nadj, 7*nnodes);
    printf("%s: Number of nodes and non-zeros %d %d\n",
            fcnm, nnodes, nadj+nnodes);
    // Create a tighter array
    adjTrim = (int *) calloc((size_t) nadj, sizeof(int));
    //*adjncyOut = (int *) calloc((size_t) nadj, sizeof(int));
    for (i=0; i<nadj; i++)
    {
        adjTrim[i] = adjncy[i];
    }
    *xadjOut = xadj;
    *adjncyOut = adjTrim;
/*
FILE *graph = fopen("graph.txt", "w");
for (i=0; i<mnodes; i++)
{
  for (j=xadj[i]; j<xadj[i+1]; j++)
  {
   fprintf(graph, "%d %d\n", mnodes-i, adjncy[j]);
  }
  fprintf(graph, "\n");
}
fclose(graph);
*/
    free(adjncy);
    return 0;
}
//============================================================================//
/*!
 * @brief Creates the loop indices for the legacy sweeping so that one can
 *        do something like:
 *         for (k = y1; k != k2; k=k+kdir)
 *             for (j = x1; j != x2; j=j+jdir)
 *                 for (i= z1; i != z2; i=i+idir)
 *
 * @param[in] sweep     Sweep number (=1,2,...,8)
 * @param[in] linitk    If true then points downwind of the source will
 *                      not be in the loop grid.
 *                      Otherwise, appropriate points upwind and downwind
 *                      of source will be in loop grid.
 * @param[in] nz        Number of z grid points.
 * @param[in] nx        Number of x grid points.
 * @param[in] ny        Number of y grid points.
 * @param[in] zsi       z source index (C numbered).
 * @param[in] xsi       x source index (C numbered).
 * @param[in] ysi       y source index (C numbered).
 * @param[out] z1       z start index.
 * @param[out] z2       z stop index (not inclusive).
 * @param[out] x1       x start index.
 * @param[out] x2       x stop index (not inclusive).
 * @param[out] y1       y start index.
 * @param[out] y2       y stop index (not inclusive).
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker
 *
 * @copyright MIT
 *
 */
int graph_getSweepLoopLimits(const int sweep, const bool linitk,
                             const int nz, const int nx, const int ny,
                             const int zsi, const int xsi, const int ysi,
                             int *z1, int *z2,
                             int *x1, int *x2,
                             int *y1, int *y2)
{
    const char *fcnm = "graph_getSweepLoopLimits\0";
    if (sweep == 1)
    {
        *y1 = 1;  *x1 = 1;  *z1 = 1;
        *y2 = ny; *x2 = nx; *z2 = nz;
        if (linitk)
        {
            *y1 = MAX(1, ysi); *x1=MAX(1, xsi); *z1=MAX(1,zsi);
            *y2 = ny;          *x2=nx;          *z2=nz;
        }
    }
    else if (sweep == 2)
    {
        *y1 = 1;  *x1 = nx-2; *z1 = 1;
        *y2 = ny; *x2 =-1;    *z2 = nz;
        if (linitk)
        {
            *y1 = MAX(1, ysi); *x1=xsi+1;       *z1=MAX(1,zsi);
            *y2 = ny;          *x2=-1;          *z2=nz;
        }
    }
    else if (sweep == 3)
    {
        *y1 = ny-2; *x1 = 1;  *z1 = 1;
        *y2 =-1;    *x2 = nx; *z2 = nz;
        if (linitk)
        {
            *y1 = ysi+1;       *x1=MAX(1, xsi); *z1=MAX(1,zsi);
            *y2 =-1;           *x2=nx;          *z2=nz;
        }
    }
    else if (sweep == 4)
    {
        *y1 = ny-2; *x1 = nx-2; *z1 = 1;
        *y2 =-1;    *x2 =-1;    *z2 = nz;
        if (linitk)
        {
            *y1 = ysi+1;       *x1=xsi+1;       *z1=MAX(1,zsi);
            *y2 =-1;           *x2=-1;          *z2=nz;
        }
    }
    else if (sweep == 5)
    {
        *y1 = 1;  *x1 = 1;  *z1 = nz-2;
        *y2 = ny; *x2 = nx; *z2 =-1;
        if (linitk)
        {
            *y1 = MAX(1, ysi); *x1=MAX(1, xsi); *z1=zsi+1;
            *y2 = ny;          *x2=nx;          *z2=-1;
        }
    }
    else if (sweep == 6)
    {
        *y1 = 1;  *x1 = nx-2; *z1 = nz-2;
        *y2 = ny; *x2 =-1;    *z2 =-1;
        if (linitk)
        {
            *y1 = MAX(1, ysi); *x1=xsi+1;       *z1=zsi+1;
            *y2 = ny;          *x2=-1;          *z2=-1;
        }
    }
    else if (sweep == 7)
    {
        *y1 = ny-2; *x1 = 1;  *z1 = nz-2;
        *y2 =-1;    *x2 = nx; *z2 =-1;
        if (linitk)
        {
            *y1 = ysi+1;       *x1=MAX(1, xsi); *z1=zsi+1;
            *y2 =-1;           *x2=nx;          *z2=-1;
        }
    }
    else if (sweep == 8)
    {
        *y1 = ny-2; *x1 = nx-2; *z1 = nz-2;
        *y2 =-1;    *x2 =-1;    *z2 =-1;
        if (linitk)
        {
            *y1 = ysi+1;       *x1=xsi+1;       *z1=zsi+1;
            *y2 =-1;           *x2=-1;          *z2=-1;
        }
    }
    else
    {
        printf("%s: Error; invalid sweep %d\n", fcnm, sweep);
        return -1;
    }
    return 0;
}
//============================================================================//
int graph_createUpdateInitMask(const int sweep,
                               const int zsi, const int xsi, const int ysi,
                               struct levelSet_struct *levelSet)
{
    const char *fcnm = "graph_createUpdateInitMask\0";
    int ix, iy, iz, level, node, nodeEnd, nodeStart, nx, ny, nz,
        xh, xl, yh, yl, zh, zl;
    const bool linitk = true;
    //------------------------------------------------------------------------//
    nz = levelSet->nz;
    nx = levelSet->nx;
    ny = levelSet->ny;
    if (sweep == 1)
    {
        zl = MAX(1, zsi);
        zh = nz - 1;
        xl = MAX(1, xsi);
        xh = nx - 1;
        yl = MAX(1, ysi);
        yh = ny - 1;
    }
    else if (sweep == 2)
    {
        zl = MAX(1, zsi);
        zh = nz - 1;
        xl = 0;
        xh = MIN(nx - 1, xsi + 1);
        yl = MAX(1, ysi);
        yh = ny - 1;
    }
    else if (sweep == 3)
    {
        zl = MAX(1, zsi);
        zh = nz - 1;
        xl = MAX(1, xsi);
        xh = nx - 1;
        yl = 0; 
        yh = MIN(ny - 1, ysi + 1);
    }
    else if (sweep == 4)
    {
        zl = MAX(1, zsi);
        zh = nz - 1;
        xl = 0;
        xh = MIN(nx - 1, xsi + 1);
        yl = 0;
        yh = MIN(ny - 1, ysi + 1);
    }
    else if (sweep == 5)
    {
        zl = 0;
        zh = MIN(nz - 1, zsi + 1);
        xl = MAX(1, xsi);
        xh = nx - 1;
        yl = MAX(1, ysi);
        yh = ny - 1;
    }
    else if (sweep == 6)
    {
        zl = 0;
        zh = MIN(nz - 1, zsi + 1);
        xl = 0;
        xh = MIN(nx - 1, xsi + 1);
        yl = MAX(1, ysi);
        yh = ny - 1;
    }
    else if (sweep == 7)
    {
        zl = 0;
        zh = MIN(nz - 1, zsi + 1);
        xl = MAX(1, xsi);
        xh = nx - 1;
        yl = 0;
        yh = MIN(ny - 1, ysi + 1);
    }
    else if (sweep == 8)
    {   
        zl = 0;
        zh = MIN(nz - 1, zsi + 1);
        xl = 0;
        xh = MIN(nx - 1, xsi + 1); 
        yl = 0;
        yh = MIN(ny - 1, ysi + 1);
    }
    else
    {
        printf("%s: Invalid sweep: %d\n", fcnm, sweep);
        return -1;
    }
printf("%d %d %d %d %d %d\n", zl, zh, xl, xh, yl, yh);
    for (level=0; level<levelSet->nLevels; level++)
    {   
        nodeStart = levelSet->levelPtr[level];
        nodeEnd   = levelSet->levelPtr[level+1];
        for (node=nodeStart; node<nodeEnd; node++)
        {
            iz = levelSet->n2ijk[3*node+0]; // z
            ix = levelSet->n2ijk[3*node+1]; // x
            iy = levelSet->n2ijk[3*node+2]; // y
            levelSet->lupdateInit[node] = levelSet->lupdate[node];
            if (iz < zl || iz > zh){levelSet->lupdateInit[node] = false;}
            if (ix < xl || ix > xh){levelSet->lupdateInit[node] = false;}
            if (iy < yl || iy > yh){levelSet->lupdateInit[node] = false;}
        }
    }
    return 0;
}
//============================================================================//
/*!
 * @brief Sets a mask on the updates that emerges because we do not want
 *        to selectively used the 1D, 2D, and 3D stencils based on the
 *        location in the mesh.
 *
 * @param[in] sweep         Sweep number (=1,2,...,8).
 *
 * @param[in,out] levelSet  On input contains the node to (i,j,k) pointer
 *                          and level set structure.
 *                          On output lupdate nodes marked as false will
 *                          be updated in the sweep.
 *
 * @author Ben Baker
 *
 * @copyright MIT
 *
 */
int graph_createUpdateMask(const int sweep,
                           struct levelSet_struct *levelSet)
{
    const char *fcnm = "graph_createUpdateMask\0";
    int i, indx, j, k, level, node, nodeStart, nodeEnd, wallx, wally, wallz;
    size_t nbytes;
    if (sweep == 1)
    {
        wallz = 0;
        wallx = 0;
        wally = 0;
    }
    else if (sweep == 2)
    {
        wallz = 0;
        wallx = levelSet->nx - 1;
        wally = 0;
    }
    else if (sweep == 3)
    {
        wallz = 0;
        wallx = 0;
        wally = levelSet->ny - 1;
    }
    else if (sweep == 4)
    {
        wallz = 0;
        wallx = levelSet->nx - 1;
        wally = levelSet->ny - 1;
    }
    else if (sweep == 5)
    {   
        wallz = levelSet->nz - 1;
        wallx = 0;
        wally = 0;
    }   
    else if (sweep == 6)
    {   
        wallz = levelSet->nz - 1;
        wallx = levelSet->nx - 1;
        wally = 0;
    }   
    else if (sweep == 7)
    {   
        wallz = levelSet->nz - 1;
        wallx = 0;
        wally = levelSet->ny - 1;
    }
    else if (sweep == 8)
    {
        wallz = levelSet->nz - 1;
        wallx = levelSet->nx - 1;
        wally = levelSet->ny - 1;
    }
    else
    {
        printf("%s: Invalid sweep number %d\n", fcnm, sweep);
        return -1;
    }
    nbytes = (size_t) (levelSet->nx*levelSet->ny*levelSet->nz)*sizeof(bool);
    levelSet->lupdate = (bool *) aligned_alloc(64, nbytes);
    for (level=0; level<levelSet->nLevels; level++)
    {
        nodeStart = levelSet->levelPtr[level];
        nodeEnd   = levelSet->levelPtr[level+1];
        for (node=nodeStart; node<nodeEnd; node++)
        {
            i = levelSet->n2ijk[3*node+0];
            j = levelSet->n2ijk[3*node+1];
            k = levelSet->n2ijk[3*node+2];
            indx = node; //indx = levelSet->ijkv[node];
            levelSet->lupdate[indx] = true;
            if (k == wally){levelSet->lupdate[indx] = false;}
            if (j == wallx){levelSet->lupdate[indx] = false;}
            if (i == wallz){levelSet->lupdate[indx] = false;}
        }
    } 
/*
    for (k=0; k<levelSet->ny; k++)
    {
        for (j=0; j<levelSet->nx; j++)
        {
            for (i=0; i<levelSet->nz; i++)
            {
                indx = grd2indx(i, j, k,
                                levelSet->nz, levelSet->nx, levelSet->ny);
                levelSet->lupdate[indx] = true;
                if (k == wally){levelSet->lupdate[indx] = false;}
                if (j == wallx){levelSet->lupdate[indx] = false;}
                if (i == wallz){levelSet->lupdate[indx] = false;}
            }
        }
    }
*/
    return 0; 
}
//============================================================================//
/*!
 * @brief Permutes the level set structure for sweep 1 to the corresponding
 *        level set structure for the given sweep.
 *
 * @param[in] nz            Number of z grid points in travel time grid.
 * @param[in] nx            Number of x grid points in travel time grid.
 * @param[in] ny            Number of y grid points in travel time grid.
 * @param[in] sweep         Sweep number (=2,3,4,...8).
 * @param[in] nLevels       Number of levels.
 *
 * @param[out] ijkvOut      Maps from the node in the sweep to the global
 *                          nodal number [nz*nx*ny].
 * @param[out] levelPtrOut  Maps from the level'th level to the start index of
 *                          ijkvOut [nLevels+1].
 * @param[out] n2lPermOut   The level to which the i'th node belongs to
 *                          [nz*nx*ny].
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker
 *
 * @copyright Ben Baker distributed under the MIT license.
 *
 */
int graph_permuteLevelStructure(const int nz, const int nx, const int ny,
                                const int sweep, const int nLevels,
                                const int *__restrict__ n2l,
                                int **ijkvOut,
                                int **levelPtrOut,
                                int **n2lPermOut)
{
    const char *fcnm = "graph_permuteLevelStructure\0";
    bool lrevx, lrevy, lrevz;
    int *ijkv, *levelPtr, *n2lPerm, *perm, i, ierr, ii, indx, j, jj, jndx,
        k, kk, n, levelCtr, nnodes, nnodesAll;
    n = nx*ny*nz;
    ijkv     = (int *) aligned_alloc(64, (size_t) n*sizeof(int));
    levelPtr = (int *) aligned_alloc(64, (size_t) (nLevels+1)*sizeof(int));
    n2lPerm  = (int *) aligned_alloc(64, (size_t) n*sizeof(int));
    lrevx = false;
    lrevy = false;
    lrevz = false;
    if (sweep == 1)
    {
        lrevz = false; // top to bottom
        lrevx = false; // west to east
        lrevy = false; // south to north
    }
    else if (sweep == 2)
    {
        lrevz = false; // top to bottom
        lrevx = true;  // east to west 
        lrevy = false; // south to north
    }
    else if (sweep == 3)
    {
        lrevz = false; // top to bottom
        lrevx = false; // west to east 
        lrevy = true;  // north to south 
    }
    else if (sweep == 4)
    {
        lrevz = false; // top to bottom
        lrevx = true;  // east to west 
        lrevy = true;  // north to south
    }
    else if (sweep == 5)
    {
        lrevz = true;  // bottom to top 
        lrevx = false; // west to east
        lrevy = false; // south to north
    }
    else if (sweep == 6)
    {
        lrevz = true;  // bottom to top
        lrevx = true;  // east to west 
        lrevy = false; // south to north
    }
    else if (sweep == 7)
    {
        lrevz = true;  // bottom to top
        lrevx = false; // west to east 
        lrevy = true;  // north to south 
    }
    else if (sweep == 8)
    {
        lrevz = true;  // bottom to top
        lrevx = true;  // east to west 
        lrevy = true;  // north to south
    }
    for (k=0; k<ny; k++)
    {
        for (j=0; j<nx; j++)
        {
            for (i=0; i<nz; i++)
            {
                kk = k;
                jj = j;
                ii = i;
                if (lrevz){ii = nz - 1 - i;}
                if (lrevx){jj = nx - 1 - j;}
                if (lrevy){kk = ny - 1 - k;}
                indx = grd2indx(i,  j,  k,  nz, nx, ny);
                jndx = grd2indx(ii, jj, kk, nz, nx, ny);
                n2lPerm[jndx] = n2l[indx];
            }
        }
    }
    // Now I need to re-sort to make an appropriate levelPtr
    perm = (int *) calloc((size_t) (n+1), sizeof(int));
    sorting_argsort32i_work(n, n2lPerm, SORT_ASCENDING, perm);
    // Gather all the nodes in each level
    nnodesAll = 0;
    levelPtr[0] = 0;
    for (levelCtr=1; levelCtr<=nLevels; levelCtr++)
    {
        // Copy all nodes in this level
        nnodes = 0;
        for (i=nnodesAll; i<n; i++)
        {
            if (n2lPerm[perm[i]] == levelCtr)
            {
                ijkv[nnodesAll+nnodes] = perm[i];
                nnodes = nnodes + 1;
            }
            if (n2lPerm[perm[i+1]] != levelCtr){break;}
        }
        // Put the (i,j,k) nodes in order to raise chances of cache reuse
        ierr = sorting_sort32i_work(nnodes, &ijkv[nnodesAll], SORT_ASCENDING);
        if (ierr != 0)
        {
            printf("%s: Failed to sort nodes\n", fcnm);
            return -1; 
        }
        nnodesAll = nnodesAll + nnodes;
        levelPtr[levelCtr] = nnodesAll;
    }
    if (nnodesAll != n)
    {
        printf("%s: Failed to put all nodes in a level\n", fcnm);
    }
    *ijkvOut = ijkv;
    *levelPtrOut = levelPtr;
    *n2lPermOut = n2lPerm;
    return 0;
}
/*
int graph_node2Level(const int nLevels, const int *__restrict__ levelPtr,
                     int *__restrict__ node2level)
{
    int i;
    for (level=0; level<nLevels; level++)
    {

    }
}
*/
 
/*!
 * @brief Computes the nodes at each level
 *
 * @param[in] n                Number of nodes in mesh.
 * @param[in] xadj             Points from i'th node to start of adjncy [n+1].
 * @param[in] adjncy           Contains the non-zero columns of the
 *                             lower-trianglular sparse matrix [xadj[n]].
 * @param[out] nLevels         The number of levels (probably MAX(nx, ny, nz)).
 * @param[in,out] levelPtrOut  On intput this is NULL.
 *                             On output maps from from the level'th level
 *                             to the start index of ijkvOut [nLevels+1].
 * @param[in,out] ijkvOut      Given the start and indices from levelPtrOut
 *                              l1 = levelPtrOut[level] and
 *                              l2 = levelPtrOut[level+1]
 *                             ijkvOut[l1:l2] (l2 not included) is the sorted
 *                             nodes to update at this level.
 *                    
 * @author Ben Baker
 *
 * @copyright Ben Baker distributed under the MIT license.
 *
 */ 
int graph_createLevelStructure(const int n, const int *__restrict__ xadj,
                               const int *__restrict__ adjncy,
                               int *nLevels, int **levelPtrOut,
                               int **ijkvOut, int **n2lOut)
{
    const char *fcnm = "graph_createLevelStructure\0";
    int *level, *perm, *levelPtr, *ijkv, *n2l,
        i, ierr, j, jcol, k, levelCtr, lfirst,
        nnodes, nnodesAll, noff;
    bool lallDone;
    // Initialize the first level which has no dependencies (i.e. diagonal) 
    // and set all other nodes to -1 
    printf("%s: Generating levels...\n", fcnm);
    *nLevels = 0;
    level = (int *) calloc((size_t) (n + 1), sizeof(int));
    lfirst = 0;
    levelCtr = 1;
    for (i=0; i<n; i++)
    {
        noff = xadj[i+1] - xadj[i];
        level[i] =-1;
        if (noff == 0)
        {
            level[i] = levelCtr;
            lfirst = lfirst + 1;
        }
    }
    // Not sure how to deal with multiple root nodes
    if (lfirst != 1)
    {
        if (lfirst < 1)
        {
            printf("%s: Can't create level 1 - DAG is cyclic\n", fcnm);
        }
        else
        {
            printf("%s: Warning - multiple root nodes %d\n", fcnm, lfirst);
        }
        return -1;
    }
    // Worse case is each node is each node depends on the previous node 
    for (k=0; k<n; k++)
    {
        // Loop on the nodes in the graph
        lallDone = true;
        for (i=0; i<n; i++)
        {
            if (level[i] !=-1){continue;} // Node is done
            lallDone = false;
            for (j=xadj[i]; j<xadj[i+1]; j++)
            {
                jcol = adjncy[j];
                if (level[jcol] == levelCtr)
                {
                    //printf("%d %d %d %d %d\n", k, i, j, jcol, levelCtr);
                    level[i] = levelCtr + 1;
                    break;
                }
            }
        }
        if (lallDone){break;} 
        levelCtr = levelCtr + 1;
    }
    *nLevels = levelCtr;
    // Put into a level structure
    printf("%s: Number of levels: %d\n", fcnm, *nLevels);
    //ijkv = (int *) calloc((size_t) n, sizeof(int));
    //n2l = (int *) calloc((size_t) n, sizeof(int));
    //levelPtr = (int *) calloc((size_t) (*nLevels+1), sizeof(int));
    ijkv = (int *) aligned_alloc(64, ((size_t) n*sizeof(int)));
    n2l = (int *) aligned_alloc(64, ((size_t) n*sizeof(int)));
    levelPtr = (int *) aligned_alloc(64, ((size_t) (*nLevels+1)*sizeof(int)));
    perm = (int *) calloc((size_t) (n + 1), sizeof(int));
    perm[n] = n;
    sorting_argsort32i_work(n, level, SORT_ASCENDING, perm);
    nnodesAll = 0;
    levelPtr[0] = 0;
    for (levelCtr=1; levelCtr<=*nLevels; levelCtr++)
    {
        // Copy all nodes in this level
        nnodes = 0;
        for (i=nnodesAll; i<n; i++)
        {
            if (level[perm[i]] == levelCtr)
            {
                ijkv[nnodesAll+nnodes] = perm[i];
                n2l[perm[i]] = levelCtr;
                nnodes = nnodes + 1;
            }
            if (level[perm[i+1]] != levelCtr){break;}
        }
        // Put the (i,j,k) nodes in order to raise chances of cache reuse
        ierr = sorting_sort32i_work(nnodes, &ijkv[nnodesAll], SORT_ASCENDING);
        if (ierr != 0)
        {
            printf("%s: Failed to sort nodes\n", fcnm);
            return -1;
        }
        nnodesAll = nnodesAll + nnodes;
        levelPtr[levelCtr] = nnodesAll;
    }
    if (nnodesAll != n)
    {
        printf("%s: Failed to put all nodes in a level\n", fcnm);
    }
/*
    // As an example: loop on the levels
    for (levelCtr=0; levelCtr<*nLevels; levelCtr++)
    {
        // Set an output file name
        int j1 = levelPtr[levelCtr];
        int j2 = levelPtr[levelCtr+1];
        char fname[128];
        memset(fname, 0, 128*sizeof(char));
        sprintf(fname, "level_%d.txt", levelCtr+1); 
        FILE *ofl = fopen(fname, "w");
        // Loop on the nodes in the level
        for (j=j1; j<j2; j++)
        {
            i = ijkv[j]; // Extract the node number
            // Print out the connectivity
            for (int k=xadj[i]; k<xadj[i+1]; k++)
            {
                fprintf(ofl, "%d %d\n", n - i, adjncy[k]+1);
            }
            fprintf(ofl, "%d %d\n", n - i, i + 1);
        }
        fclose(ofl);
    }
*/
    *levelPtrOut = levelPtr;
    *ijkvOut = ijkv;
    *n2lOut = n2l;
    free(level);
    free(perm);
    return 0;
}
//============================================================================//
int graph_getTravelTimePermutation(const int sweep,
                                   int *__restrict__ perm)
{
    int i;
    const int perm1[8] = {6, 4, 5, 2, 3, 1,  0,  7};
    const int perm2[8] = {5, 2, 6, 4, 0, 7,  3,  1};
    const int perm3[8] = {3, 1, 0, 7, 6, 4,  5,  2};
    const int perm4[8] = {0, 7, 3, 1, 5, 2,  6,  4};
    const int perm5[8] = {4, 6, 2, 5, 1, 3,  7,  0};
    const int perm6[8] = {2, 5, 4, 6, 7, 0,  1,  3};
    const int perm7[8] = {1, 3, 7, 0, 4, 6,  2,  5};
    const int perm8[8] = {7, 0, 1, 3, 2, 5,  4,  6};
    if (sweep == 1)
    {
        for (i=0; i<8; i++){perm[i] = perm1[i];}
    } 
    else if (sweep == 2)
    {
        for (i=0; i<8; i++){perm[i] = perm2[i];}
    }
    else if (sweep == 3)
    {
        for (i=0; i<8; i++){perm[i] = perm3[i];}
    }
    else if (sweep == 4)
    {
        for (i=0; i<8; i++){perm[i] = perm4[i];}
    }
    else if (sweep == 5)
    {
        for (i=0; i<8; i++){perm[i] = perm5[i];}
    }
    else if (sweep == 6)
    {
        for (i=0; i<8; i++){perm[i] = perm6[i];}
    }
    else if (sweep == 7)
    {
        for (i=0; i<8; i++){perm[i] = perm7[i];}
    }
    else if (sweep == 8)
    {
        for (i=0; i<8; i++){perm[i] = perm8[i];}
    }
    else
    {
         return -1;
    }
    return 0;
}
//============================================================================//
int graph_getVelocityPermutation(const int sweep,
                                 int *__restrict__ perm)
{
    int i;
    const int *permToLoc;
    const int *locToLoc;
    const int permToLoc1[19] = {0, 4, 8, 12, 14, 16, 18, 5, 10, 17, 2, 9, 15, 11, 1, 6, 13, 7, 3};
    const int locToLoc1[19] = {0, 14, 10, 18, 1, 7, 15, 17, 2, 11, 8, 13, 3, 16, 4, 12, 5, 9, 6};
    const int permToLoc2[19] = {0, 8, 14, 10, 2, 4, 9, 12, 15, 16, 18, 5, 11, 17, 1, 3, 6, 13, 7};
    const int locToLoc2[19] = {0, 14, 4, 15, 5, 11, 16, 18, 1, 6, 3, 12, 7, 17, 2, 8, 9, 13, 10};
    const int permToLoc3[19] = {0, 4, 12, 5, 2, 1, 6, 8, 13, 14, 16, 18, 7, 10, 17, 3, 9, 15, 11};
    const int locToLoc3[19] = {0, 5, 4, 15, 1, 3, 6, 12, 7, 16, 13, 18, 2, 8, 9, 17, 10, 14, 11};
    const int permToLoc4[19] = {0, 2, 4, 12, 5, 1, 8, 14, 10, 3, 6, 9, 13, 15, 16, 18, 7, 11, 17};
    const int locToLoc4[19] = {0, 5, 1, 9, 2, 4, 10, 16, 6, 11, 8, 17, 3, 12, 7, 13, 14, 18, 15};
    const int permToLoc5[19] = {4, 8, 16, 0, 5, 10, 12, 14, 17, 18, 9, 2, 11, 15, 6, 1, 7, 13, 3};
    const int locToLoc5[19] = {3, 15, 11, 18, 0, 4, 14, 16, 1, 10, 5, 12, 6, 17, 7, 13, 2, 8, 9};
    const int permToLoc6[19] = {8, 0, 10, 14, 4, 9, 16, 2, 5, 11, 12, 15, 17, 18, 1, 6, 3, 7, 13};
    const int locToLoc6[19] = {1, 14, 7, 16, 4, 8, 15, 17, 0, 5, 2, 9, 10, 18, 3, 11, 6, 12, 13};
    const int permToLoc7[19] = {4, 0, 5, 12, 2, 6, 8, 16, 1, 7, 10, 13, 14, 17, 18, 9, 3, 11, 15};
    const int locToLoc7[19] = {1, 8, 4, 16, 0, 2, 5, 9, 6, 15, 10, 17, 3, 11, 12, 18, 7, 13, 14};
    const int permToLoc8[19] = {0, 4, 2, 5, 12, 8, 1, 10, 14, 6, 9, 16, 3, 7, 11, 13, 15, 17, 18};
    const int locToLoc8[19] = {0, 6, 2, 12, 1, 3, 9, 13, 5, 10, 7, 14, 4, 15, 8, 16, 11, 17, 18};
    if (sweep == 1)
    {
        permToLoc = permToLoc1;
        locToLoc = locToLoc1;
    }
    else if (sweep == 2)
    {    
        permToLoc = permToLoc2;
        locToLoc = locToLoc2;
    }
    else if (sweep == 3)
    {    
        permToLoc = permToLoc3;
        locToLoc = locToLoc3;
    }
    else if (sweep == 4)
    {    
        permToLoc = permToLoc4;
        locToLoc = locToLoc4;
    }
    else if (sweep == 5)
    {    
        permToLoc = permToLoc5;
        locToLoc = locToLoc5;
    }
    else if (sweep == 6)
    {    
        permToLoc = permToLoc6;
        locToLoc = locToLoc6;
    }
    else if (sweep == 7)
    {    
        permToLoc = permToLoc7;
        locToLoc = locToLoc7;
    }
    else if (sweep == 8)
    {    
        permToLoc = permToLoc8;
        locToLoc = locToLoc8;
    }
    else
    {
        return -1;
    }
    for (i=0; i<19; i++)
    {
        perm[i] = permToLoc[i];
    }
    return 0;
}
//============================================================================//
/*!
 */
void graph_sweep2sgnx64f(const int sweep,
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

    *sgnrz = (double) *sgntz;
    *sgnrx = (double) *sgntx;
    *sgnry = (double) *sgnty;
/*
printf("%d %d %d %d %d %d\n", *sgntz,*sgntx,*sgnty, *sgnvz,*sgnvx,*sgnvy);
    *sgntz = 1; *sgntx = 1; *sgnty = 1;
    *sgnvz = 1; *sgnvx = 1; *sgnvy = 1;
    if (sweep == 1)
    {
        *sgntz = 1; *sgntx = 1; *sgnty = 1;
        *sgnvz = 1; *sgnvx = 1; *sgnvy = 1;
    }
    else if (sweep == 2)
    {
        *sgntz = 1; *sgntx =-1; *sgnty = 1;
        *sgnvz = 1; *sgnvx = 0; *sgnvy = 1;
    }
    else if (sweep == 3)
    {
        *sgntz = 1; *sgntx = 1; *sgnty =-1;
        *sgnvz = 1; *sgnvx = 1; *sgnvy = 0;
    }
    else if (sweep == 4)
    {
        *sgntz = 1; *sgntx =-1; *sgnty =-1;
        *sgnvz = 1; *sgnvx = 0; *sgnvy = 0;
    }
    else if (sweep == 5)
    {
        *sgntz =-1; *sgntx = 1; *sgnty = 1;
        *sgnvz = 0; *sgnvx = 1; *sgnvy = 1;
    }
    else if (sweep == 6)
    {
        *sgntz =-1; *sgntx =-1; *sgnty = 1;
        *sgnvz = 0; *sgnvx = 0; *sgnvy = 1;
    }
    else if (sweep == 7)
    {
        *sgntz =-1; *sgntx = 1; *sgnty =-1;
        *sgnvz = 0; *sgnvx = 1; *sgnvy = 0;
    }
    else if (sweep == 8)
    {
        *sgntz =-1; *sgntx =-1; *sgnty =-1;
        *sgnvz = 0; *sgnvx = 0; *sgnvy = 0;
    }
    *sgnrz = (double) *sgntz;
    *sgnrx = (double) *sgntx;
    *sgnry = (double) *sgnty;
 printf("%d %d %d %d %d %d\n", *sgntz,*sgntx,*sgnty, *sgnvz,*sgnvx,*sgnvy);
*/
    return;
}

int graph_testGrd2ijk(const int n1, const int n2, const int n3)
{
    int i, i1, ierr, ierrAll, igrd, j, j1, k, k1;
    ierrAll = 0;
    for (k=0; k<n3; k++)
    {
        for (j=0; j<n2; j++)
        {
            for (i=0; i<n1; i++)
            {
                igrd = grd2indx(i, j, k, n1, n2, n3);
                ierr = grd2ijk(igrd, n1, n2, n3, &i1, &j1, &k1);
                if (ierr != 0)
                {
                    ierrAll = ierrAll + 1;
                    printf("Algorithmic failure\n");
                }
                if (i1 != i || j1 != j || k1 != k)
                {
                    ierrAll = ierrAll + 1;
                    printf("%d %d %d %d %d %d\n", i, i1, j, j1, k, k1);
                }
            }
        }
    }
    return ierrAll;
}

static int grd2indx(const int i, const int j, const int k,
                    const int n1, const int n2, const int n3)
{
    int indx;
    indx = k*n1*n2 + j*n1 + i;
    return indx;
}

static int grd2ijk(const int igrd,
                   const int n1, const int n2, const int n3,
                   int *i, int *j, int *k)
{
    int ierr, n12;
    ierr = 0;
    n12 = n1*n2;
    *k = (igrd)/n12;
    *j = (igrd - *k*n12)/n1;
    *i =  igrd - *k*n12 - *j*n1;
    if (*i < 0 || *i > n1 - 1){ierr = ierr + 1;}
    if (*j < 0 || *j > n2 - 1){ierr = ierr + 1;}
    if (*k < 0 || *k > n3 - 1){ierr = ierr + 1;}
    return ierr;
}
