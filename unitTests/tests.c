#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fteik_fortran.h>
#include <fteik_graph.h>

int test_slowPerm(void);
int test_ttPerm(void);
int test_sweepSigns(void);
int test_gridConversions(void);

int main()
{
    int ierr;
    ierr = test_gridConversions();
    if (ierr != EXIT_SUCCESS)
    {
        printf("Failed gridConversions test\n");
        return EXIT_FAILURE;
    }
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
    return EXIT_SUCCESS;
}

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
                indx = fteik_grid2indexF(i, j, k, nz, nzx);
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
