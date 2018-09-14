#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fteik/fteik_fortran.h>
#include "./fteik_graph.h"
#include "./fteikGraph.h"

//int test_slowPerm(void);
//int test_ttPerm(void);
//int test_sweepSigns(void);
int test_gridConversions(void);
int test_graph3d();

int main()
{
    int ierr;
    ierr = test_gridConversions();
    if (ierr != EXIT_SUCCESS)
    {
        printf("Failed gridConversions test\n");
        return EXIT_FAILURE;
    }
/*
    ierr = test_slowPerm();
    if (ierr != EXIT_SUCCESS)
    {
        printf("Failed slowPerm test\n"); 
        return EXIT_FAILURE;
    } 
    ierr = test_ttPerm();
    if (ierr != EXIT_SUCCESS)
    {
        printf("Failed ttPerm test\n");
        return EXIT_FAILURE;
    }
    ierr = test_sweepSigns();
    if (ierr != EXIT_SUCCESS)
    {
        printf("Failed sweepSigns test\n");
        return EXIT_FAILURE;
    }
*/
    ierr = test_graph3d();
    if (ierr != EXIT_SUCCESS)
    {
        printf("Failed graph test\n");
        return EXIT_FAILURE;
    }
    printf("Passed tests\n");
    return EXIT_SUCCESS;
}

int test_graph3d(void)
{
#ifdef FTEIK_USE_BOOST
    void *graph;
    int nz = 35;
    int nx = 68;
    int ny = 23;
    int ngrd = nx*ny*nz;
    int *ijkv1, *ijkv2, *ijkv3, *ijkv4, *ijkv5, *ijkv6, *ijkv7, *ijkv8;
    int *ijkv1R, *ijkv2R, *ijkv3R, *ijkv4R, *ijkv5R, *ijkv6R, *ijkv7R, *ijkv8R;
    int *levelPtr, *levelPtrRef;
    int i, ierr, maxLevelSize, nLevels, maxLevelSizeRef, nLevelsRef, nwork, sweep;
    fprintf(stdout, "%s: Generating graph...\n", __func__);
    fteik_graph3d_initialize(nz, nx, ny, 1, &ierr); 
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error initilializing graph\n", __func__);
        return EXIT_FAILURE;
    }
    fteik_graph3d_makeLevelStructures(&ierr);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error making the level structures\n", __func__);
        return EXIT_FAILURE;
    }
    // Pick up results
    maxLevelSize = fteik_graph3d_getMaxLevelSize(&ierr);
    nLevels = fteik_graph3d_getNumberOfLevels(&ierr);
    levelPtr = (int *) calloc((size_t) nLevels+1, sizeof(int));
    fteik_graph3d_getLevelPointerF(nLevels+1, levelPtr, &ierr);
    ijkv1 = (int *) calloc((size_t) (4*ngrd), sizeof(int));
    ijkv2 = (int *) calloc((size_t) (4*ngrd), sizeof(int));
    ijkv3 = (int *) calloc((size_t) (4*ngrd), sizeof(int));
    ijkv4 = (int *) calloc((size_t) (4*ngrd), sizeof(int));
    ijkv5 = (int *) calloc((size_t) (4*ngrd), sizeof(int));
    ijkv6 = (int *) calloc((size_t) (4*ngrd), sizeof(int));
    ijkv7 = (int *) calloc((size_t) (4*ngrd), sizeof(int));
    ijkv8 = (int *) calloc((size_t) (4*ngrd), sizeof(int));
    fteik_graph3d_getIJKVF(4*ngrd, 1, ijkv1, &ierr);
    fteik_graph3d_getIJKVF(4*ngrd, 2, ijkv2, &ierr);
    fteik_graph3d_getIJKVF(4*ngrd, 3, ijkv3, &ierr);
    fteik_graph3d_getIJKVF(4*ngrd, 4, ijkv4, &ierr);
    fteik_graph3d_getIJKVF(4*ngrd, 5, ijkv5, &ierr);
    fteik_graph3d_getIJKVF(4*ngrd, 6, ijkv6, &ierr);
    fteik_graph3d_getIJKVF(4*ngrd, 7, ijkv7, &ierr);
    fteik_graph3d_getIJKVF(4*ngrd, 8, ijkv8, &ierr);
    fteik_graph3d_free();

    // Do the reference problem
    graph = fteik_graph_initialize(nz, nx, ny, &ierr);
    nLevelsRef = fteik_graph_getNumberOfLevels(graph);
    maxLevelSizeRef = fteik_graph_getMaxLevelSize(graph);
    levelPtrRef = (int *) calloc((size_t) nLevelsRef+1, sizeof(int));
    ijkv1R = (int *) calloc((size_t) (4*ngrd), sizeof(int)); 
    ijkv2R = (int *) calloc((size_t) (4*ngrd), sizeof(int));
    ijkv3R = (int *) calloc((size_t) (4*ngrd), sizeof(int));
    ijkv4R = (int *) calloc((size_t) (4*ngrd), sizeof(int));
    ijkv5R = (int *) calloc((size_t) (4*ngrd), sizeof(int));
    ijkv6R = (int *) calloc((size_t) (4*ngrd), sizeof(int));
    ijkv7R = (int *) calloc((size_t) (4*ngrd), sizeof(int));
    ijkv8R = (int *) calloc((size_t) (4*ngrd), sizeof(int));
    nwork = nLevelsRef + 1;
    sweep = 1;
    fteik_graph_getLevelPointerF(graph, &nwork, &sweep, levelPtrRef); 
    nwork = 4*ngrd;
    sweep = 1; 
    fteik_graph_getIJKVF(graph, &nwork, &sweep, ijkv1R);
    sweep = 2;  
    fteik_graph_getIJKVF(graph, &nwork, &sweep, ijkv2R);
    sweep = 3;  
    fteik_graph_getIJKVF(graph, &nwork, &sweep, ijkv3R);
    sweep = 4;  
    fteik_graph_getIJKVF(graph, &nwork, &sweep, ijkv4R);
    sweep = 5;  
    fteik_graph_getIJKVF(graph, &nwork, &sweep, ijkv5R);
    sweep = 6;
    fteik_graph_getIJKVF(graph, &nwork, &sweep, ijkv6R);
    sweep = 7;
    fteik_graph_getIJKVF(graph, &nwork, &sweep, ijkv7R);
    sweep = 8;
    fteik_graph_getIJKVF(graph, &nwork, &sweep, ijkv8R);
    fteik_graph_finalize(graph);
    // Test it
    if (nLevelsRef != nLevels)
    {
        fprintf(stderr, "%s: nLevelsRef=%d != nLevels=%d\n", __func__,
                nLevelsRef, nLevels);
        return EXIT_FAILURE;
    }
    if (maxLevelSizeRef != maxLevelSize)
    {
        fprintf(stderr, "%s: maxLevelSizeRef=%d != maxLevelSize=%d\n", __func__,
                maxLevelSizeRef, maxLevelSize);
        return EXIT_FAILURE;
    }
    for (i=0; i<nLevels+1; i++)
    {
        if (levelPtrRef[i] - levelPtr[i] != 0)
        {
            fprintf(stderr, "%s: levelPtrRef[%d]=%d != levelPtr[%d]=%d\n",
                    __func__, i, levelPtrRef[i], i, levelPtr[i]);
            return EXIT_FAILURE;
        }
    }
    for (i=0; i<4*ngrd; i++)
    {
        if (ijkv1[i] - ijkv1R[i] != 0)
        {
            fprintf(stderr, "%s: ijkv1[%d]=%d != ijkv1R[%d]=%d\n",
                    __func__, i, ijkv1[i], i, ijkv1R[i]);
            return EXIT_FAILURE;
        } 
        if (ijkv2[i] - ijkv2R[i] != 0)
        {
            fprintf(stderr, "%s: ijkv2[%d]=%d != ijkv2R[%d]=%d\n",
                    __func__, i, ijkv2[i], i, ijkv2R[i]);
            return EXIT_FAILURE;
        }
        if (ijkv3[i] - ijkv3R[i] != 0)
        {
            fprintf(stderr, "%s: ijkv3[%d]=%d != ijkv3R[%d]=%d\n",
                    __func__, i, ijkv3[i], i, ijkv3R[i]);
            return EXIT_FAILURE;
        }
        if (ijkv4[i] - ijkv4R[i] != 0)
        {
            fprintf(stderr, "%s: ijkv4[%d]=%d != ijkv4R[%d]=%d\n",
                    __func__, i, ijkv4[i], i, ijkv4R[i]);
            return EXIT_FAILURE;
        }
        if (ijkv5[i] - ijkv5R[i] != 0)
        {
            fprintf(stderr, "%s: ijkv5[%d]=%d != ijkv5R[%d]=%d\n",
                    __func__, i, ijkv5[i], i, ijkv5R[i]);
            return EXIT_FAILURE;
        }
        if (ijkv6[i] - ijkv6R[i] != 0)
        {
            fprintf(stderr, "%s: ijkv6[%d]=%d != ijkv6R[%d]=%d\n",
                    __func__, i, ijkv6[i], i, ijkv6R[i]);
            return EXIT_FAILURE;
        }
        if (ijkv7[i] - ijkv7R[i] != 0)
        {
            fprintf(stderr, "%s: ijkv7[%d]=%d != ijkv7R[%d]=%d\n",
                    __func__, i, ijkv7[i], i, ijkv7R[i]);
            return EXIT_FAILURE;
        }
        if (ijkv8[i] - ijkv8R[i] != 0)
        {
            fprintf(stderr, "%s: ijkv1[%d]=%d != ijkv8R[%d]=%d\n",
                    __func__, i, ijkv8[i], i, ijkv8R[i]);
            return EXIT_FAILURE;
        }
    }
    free(levelPtr);
    free(levelPtrRef);
    free(ijkv1); free(ijkv2); free(ijkv3); free(ijkv4);
    free(ijkv5); free(ijkv6); free(ijkv7); free(ijkv8);
    free(ijkv1R); free(ijkv2R); free(ijkv3R); free(ijkv4R); 
    free(ijkv5R); free(ijkv6R); free(ijkv7R); free(ijkv8R);
    return EXIT_SUCCESS;
#else
    return EXIT_SUCCESS;
#endif
}

/*
int test_sweepSigns(void)
{
    const char *fcnm = "test_sweepSigns\0";
    int sweep; 
    double sgnrzC, sgnrxC, sgnryC, sgnrzF, sgnrxF, sgnryF;
    int sgntzC, sgntxC, sgntyC, sgntzF, sgntxF, sgntyF;
    int sgnvzC, sgnvxC, sgnvyC, sgnvzF, sgnvxF, sgnvyF;

    for (sweep=1; sweep<=8; sweep++)
    {
        fteik_getSweepSigns64fF(&sweep,
                                &sgntzF, &sgntxF, &sgntyF,
                                &sgnvzF, &sgnvxF, &sgnvyF,
                                &sgnrzF, &sgnrxF, &sgnryF);
        graph_sweep2sgnx64f(sweep,
                            &sgntzC, &sgntxC, &sgntyC,
                            &sgnvzC, &sgnvxC, &sgnvyC,
                            &sgnrzC, &sgnrxC, &sgnryC);
        if (sgntzF != sgntzC || sgntxF != sgntxC || sgntyF != sgntyC)
        {
            printf("%s: Incorrect t sign in sweep %d\n", fcnm, sweep);
            return EXIT_FAILURE;
        }
        if (sgnvzF != sgnvzC || sgnvxF != sgnvxC || sgnvyF != sgnvyC)
        {
            printf("%s: Incorrect v sign in sweep %d\n", fcnm, sweep);
            return EXIT_FAILURE;
        }
        if (fabs(sgnrzF - sgnrzC) > 1.e-14 ||
            fabs(sgnrxF - sgnrxC) > 1.e-14 ||
            fabs(sgnryF - sgnryC) > 1.e-14)
        {
            printf("%s: Incorrect double cast in sweep %d\n", fcnm, sweep);
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}
*/

int test_gridConversions(void)
{
    const char *fcnm = "test_gridConversions\0";
    int ierr, ierrAll, nx, ny, nz, nzx;
    int indx, i, i1, j, j1, k, k1;
    nz = 34;
    nx = 10;
    ny = 21;
    fteik_setGridSizeF(&nz, &nx, &ny, &ierr);
    if (ierr != 0)
    {
        printf("%s: Failed to set grid size\n", fcnm);
        return EXIT_FAILURE;
    }
    nzx = nz*nx;
    ierrAll = 0;
    for (k=1; k<=ny; k++)
    {
        for (j=1; j<=nx; j++)
        {
            for (i=1; i<=nz; i++)
            {
                indx = fteik_model_grid2indexF(i, j, k, nz, nzx);
                fteik_index2gridF(indx, &i1, &j1, &k1, &ierr);
                if (ierr != 0)
                {
                    printf("%s: algorithmic failure\n", fcnm);
                    ierrAll = ierrAll + 1;
                }
                if (i1 != i || j1 != j || k1 != k)
                {
                    printf("%s: failed inverse %d %d %d %d %d %d \n", fcnm, i1, i, j1, j, k1, k); 
                    ierrAll = ierrAll + 1;
                }
            }
        }
    }
    fteik_finalizeF();
    if (ierrAll != 0){return EXIT_FAILURE;}
    return EXIT_SUCCESS;
}

/*
int test_slowPerm(void)
{
    const char *fcnm = "test_slowPerm\0";
    int permF[19], permC[19];
    int i, ierr, k;
    for (k=1; k<=8; k++)
    {
        fteik_getSlownessPermF(&k, permF, &ierr);
        graph_getVelocityPermutation(k, permC);
        for (i=0; i<19; i++)
        {
            if (permF[i] - 1 != permC[i])
            {
                printf("%s: Failed on sweep %d (%d /= %d)\n", fcnm, k,
                       permF[i], permC[i]);
                return EXIT_FAILURE;
            }
        }
    }
    return EXIT_SUCCESS;
}

int test_ttPerm(void)
{
    const char *fcnm = "test_ttPerm\0";
    int permF[8], permC[8];
    int i, ierr, k;
    for (k=1; k<=8; k++)
    {
        fteik_getTravelTimePermF(&k, permF, &ierr);
        graph_getTravelTimePermutation(k, permC);
        for (i=0; i<8; i++)
        {
            if (permF[i] - 1 != permC[i])
            {
                printf("%s: Failed on sweep %d (%d /= %d)\n", fcnm, k,
                       permF[i], permC[i]);
                return EXIT_FAILURE;
            }
        }
    }
    return EXIT_SUCCESS;
}
*/
