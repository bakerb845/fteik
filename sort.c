#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include "fteik_sorting.h"
#ifdef FTEIK_USE_INTEL
#include <ipps.h>
#endif

struct double2d_struct
{
    double val;
    int indx;
    char pad4[4];
};

struct int2d_struct
{
    int val;
    int indx;
};

static int cmp_int_ascending(const void *x, const void *y);
static int cmp_int_descending(const void *x, const void *y);


static int cmp_double_array(const void *x, const void *y)
{
    struct double2d_struct xx = *(struct double2d_struct *) x;
    struct double2d_struct yy = *(struct double2d_struct *) y;
    if (xx.val < yy.val) return -1;
    if (xx.val > yy.val) return  1;
    return 0;
}

static int cmp_int_array(const void *x, const void *y)
{
    struct int2d_struct xx = *(struct int2d_struct *) x;
    struct int2d_struct yy = *(struct int2d_struct *) y;
    if (xx.val < yy.val) return -1;
    if (xx.val > yy.val) return  1;
    return 0;
}

static int cmp_int_ascending(const void *x, const void *y)
{
    int xx = *(int *) x;
    int yy = *(int *) y;
    if (xx < yy) return -1;
    if (xx > yy) return  1;
    return 0;
}

static int cmp_int_descending(const void *x, const void *y)
{
    int xx = *(int *) x;
    int yy = *(int *) y;
    if (xx < yy) return  1;
    if (xx > yy) return -1;
    return 0;
}

#ifndef FTEIK_USE_INTEL
static int cmp_double(const void *x, const void *y)
{
    double xx = *(double *) x;
    double yy = *(double *) y;
    if (xx < yy) return -1;
    if (xx > yy) return  1;
    return 0;
}

static int cmp_int(const void *x, const void *y) 
{
    int xx = *(int *) x;
    int yy = *(int *) y;
    if (xx < yy) return -1; 
    if (xx > yy) return  1;  
    return 0;
}
#endif

int32_t sorting_sort32i_finter(const int32_t n, int32_t *__restrict__ a,
                               const int32_t orderIn)
{
    int32_t ierr;
    enum sortOrder_enum order;
    order = SORT_ASCENDING;
    if (orderIn == 1){order = SORT_DESCENDING;} 
    ierr = sorting_sort32i_work(n, a, order);
    return ierr; 
}

int32_t sorting_sort64f_finter(const int32_t n, double *__restrict__ a,
                               const int32_t orderIn)
{
    int32_t ierr;
    enum sortOrder_enum order;
    order = SORT_ASCENDING;
    if (orderIn == 1){order = SORT_DESCENDING;}
    ierr = sorting_sort64f_work(n, a, order);
    return ierr; 
}

int32_t sorting_bsearch32i_finter(const int32_t *key,
                                  const int32_t *__restrict__ values,
                                  const int32_t *n,
                                  const int32_t *hintIn)
{
    int32_t indx;
    enum sortOrder_enum hint;
    hint = (enum sortOrder_enum) *hintIn;
    indx = sorting_bsearch32i(*key, values, *n, hint);
    if (indx >-1){indx = indx + 1;} // C -> F numbering
    return indx;
}

int32_t sorting_argsort64f_finter(const int32_t *n,
                                  const double *__restrict__ a,
                                  const int32_t *orderIn,
                                  int *__restrict__ perm)
{
    int i;
    int32_t ierr;
    enum sortOrder_enum order;
    order = (enum sortOrder_enum) *orderIn;
    ierr = sorting_argsort64f_work(*n, a, order, perm);
    if (ierr != 0)
    {
        for (i=0; i<*n; i++)
        {
            perm[i] = perm[i] + 1;
        }
    }
    return ierr;
}

int32_t sorting_sort2(const int32_t *n,
                      double *__restrict__ x, double *__restrict__ y)
{
    double *work;
    int *perm;
    int i, ierr;
    perm = (int *) calloc((size_t) *n, sizeof(int));
    work = (double *) calloc((size_t) *n, sizeof(double));
    ierr = sorting_argsort64f_work(*n, x, SORT_ASCENDING, perm);
    if (ierr != 0)
    {
        printf("sorting_sort2: Sort failure\n");
        return -1; 
    }
    for (i=0; i<*n; i++){work[i] = x[i];}
    for (i=0; i<*n; i++)
    {
        x[i] = work[perm[i]];
    }
    for (i=0; i<*n; i++){work[i] = y[i];}
    for (i=0; i<*n; i++)
    {
        y[i] = work[perm[i]];
    }

    free(perm);
    free(work);
    return 0;  
}

int32_t sorting_bsearch64f_finter(const int32_t *n, const int32_t *job,
                                  const double *tol,
                                  const double *x,
                                  const double *__restrict__ a)
{
    int i, j, k;
    i = 0;
    j = *n - 1;
    while (i <= j)
    {
        k = (i + j) / 2;
        if (fabs(a[k] - *x) < *tol) //== x) {
        {
            i = k;
            if (*job == 1)
            {
                for (j=k; j>=1; j--)
                {
                    if (fabs(a[j] - a[j-1]) > *tol)
                    {
                        return j + 1;
                    }
                }
               return 1;
            }
            return i;
        }
        else if (a[k] < *x)
        {
            i = k + 1;
        }
        else
        {
            j = k - 1;
        }
    }
    return -1;
}
/*! 
 * @brief Perform an argument sort with GSL for a double precision array.
 * 
 * @param[in] n        length of array a
 * @param[in] a        array to be argsorted [n]
 *
 * @param[out] perm    map so that a[perm[:]] is in ascending order [n]
 *
 * @result 0 indicates success 
 *
 * @author Ben Baker, ISTI
 *
 */
int sorting_argsort64f_work(const int n,
                            const double *__restrict__ a,
                            const enum sortOrder_enum order,
                            int *__restrict__ perm)
{
    struct double2d_struct *vals;
    int i;
    //------------------------------------------------------------------------//
    //  
    // Error check
    if (n < 1)
    {   
        printf("%s: Warning: No data to sort!\n", __func__);
        return -1;
    }
    if (a == NULL || perm == NULL)
    {   
        if (a == NULL){printf("%s: Error a is NULL\n", __func__);}
        if (perm == NULL){printf("%s: Error perm is NULL\n", __func__);}
        return -1;
    }
    // Special case
    if (n == 1)
    {   
        perm[0] = 0;
        return 0;
    }
    // Set workspace
    vals = (struct double2d_struct *)
           calloc((size_t) n, sizeof(struct double2d_struct));
    for (i=0; i<n; i++)
    {
        vals[i].val = a[i];
        vals[i].indx = i;
    }
    qsort((void *) vals, (size_t) n,
          sizeof(struct double2d_struct), cmp_double_array);
    // Copy answer back
    if (order == SORT_ASCENDING)
    {
        for (i=0; i<n; i++)
        {
            perm[i] = vals[i].indx;
        }
    }
    else
    {
        for (i=0; i<n; i++)
        {
            perm[i] = vals[n-1-i].indx;
        }
    }
    free(vals);
    return 0;
}
//============================================================================//
/*! 
 * @brief Perform an argument sort with GSL for an integer array.
 * 
 * @param[in] n        length of array a
 * @param[in] ia       array to be argsorted [n]
 *
 * @param[out] perm    map so that a[perm[:]] is in ascending order [n]
 *
 * @result 0 indicates success
 *
 * @author Ben Baker, ISTI
 *
 */
int sorting_argsort32i_work(const int n,
                            const int *__restrict__ ia,
                            const enum sortOrder_enum order,
                            int *__restrict__ perm)
{
    struct int2d_struct *vals;
    int i;
    //------------------------------------------------------------------------//
    //  
    // Error check
    if (n < 1)
    {
        printf("%s: Warning: No data to sort!\n", __func__);
        return -1;
    }
    if (ia == NULL || perm == NULL)
    {
        if (ia == NULL){printf("%s: Error ia is NULL\n", __func__);}
        if (perm == NULL){printf("%s: Error perm is NULL\n", __func__);}
        return -1;
    }
    // Special case
    if (n == 1)
    {
        perm[0] = 0;
        return 0;
    }
    // Set workspace
    vals = (struct int2d_struct *)
           calloc((size_t) n, sizeof(struct int2d_struct));
    for (i=0; i<n; i++)
    {
        vals[i].val = ia[i];
        vals[i].indx = i;
    }
    qsort((void *) vals, (size_t) n,
          sizeof(struct int2d_struct), cmp_int_array);
    // Copy answer back
    if (order == SORT_ASCENDING)
    {
        for (i=0; i<n; i++)
        {
            perm[i] = vals[i].indx;
        }
    }
    else
    {
        for (i=0; i<n; i++)
        {
            perm[i] = vals[n-1-i].indx;
        }
    }
    free(vals);
    return 0;
}
//============================================================================//
/*!
 * @brief Sorts a double array in place
 *
 * @param[in] n       length of array a
 *
 * @param[inout] a    On input the array to sort.
 *                    On successful the sorted array.
 *
 * @param[in] order   ASCENDING -> sort array in ascending order (default)
 *                    DESCENDING -> sort array in descending order
 *
 * @author Ben Baker, ISTI
 *
 */
#ifdef FTEIK_USE_INTEL
int sorting_sort64f_work(const int n,
                         double *__restrict__ a,
                         const enum sortOrder_enum order)
{
    IppStatus status;
    // Error checking
    if (n < 1)
    {
        printf("%s: Error no data to sort\n", __func__);
        return -1;
    }
    if (order == SORT_ASCENDING)
    {
        status = ippsSortAscend_64f_I(a, n);
    }
    else
    {
        status = ippsSortDescend_64f_I(a, n);
    }
    if (ippStsNoErr != status)
    {
        printf("%s: Sort failed!\n", __func__);
        return -1;
    }
    return 0;
}
#else
static int array_reverse64f_work(const int n, double *a, double *b);
int sorting_sort64f_work(const int n, double *__restrict__ a,
                         const enum sortOrder_enum order)
{
    int ierr;
    size_t ns = (size_t) n;
    //------------------------------------------------------------------------//
    //  
    // Error checking
    if (n < 1)
    {
        printf("%s: Error no data to sort\n", __func__);
        return -1;
    }
    // Early return
    if (n == 1){return 0;}
    qsort((void *)a, (size_t) ns, sizeof(double), cmp_double);
    if (order == SORT_DESCENDING)
    {
        ierr = array_reverse64f_work(n, a, a);
        if (ierr != 0)
        {
            printf("%s: Error reversing array!\n", __func__);
            return -1;
        }
    }
    return 0;
}
/*!
 * @brief Reverse a double array
 *
 * @param[in] n     size of arrays
 * @param[in] a     array to be reversed
 *
 * @param[out] b    reversed form of a (can be called with same argument as a)
 *
 * @result 0 indicates success
 *
 * @author Ben Baker, ISTI
 *
 */
static int array_reverse64f_work(const int n, double *a, double *b)
{
#ifdef FTEIK_USE_INTEL
    IppStatus status;
#else
    double bj;
    int i, j;
#endif
    if (n < 1 || a == NULL || b == NULL)
    {
        if (n < 1){printf("%s: Invalid length: %d\n", __func__, n);}
        if (a == NULL){printf("%s: Error a is NULL\n", __func__);}
        if (b == NULL){printf("%s: Error b is NULL\n", __func__);}
        return -1;
    }
    // In place
    if (a == b)
    {
#ifdef FTEIK_USE_INTEL
        status = ippsFlip_64f_I(a, n);
        if (status != ippStsNoErr)
        {
            if (status == ippStsNullPtrErr)
            {
                if (a == NULL)
                {
                    printf("%s: Array is NULL\n", __func__);
                }
            }
            if (status == ippStsSizeErr)
            {
                printf("%s: Invalid size\n", __func__);
            }
            return -1; 
        }
#else
        memcpy(b, a, (size_t)n*sizeof(double));
        j = n - 1;
        for (i=0; i<n/2; i++)
        {
            bj = b[j];
            b[j] = b[i];
            b[i] = bj;
            j = j - 1;
        }
#endif
    }
    // Different memory addresses
    else
    {
#ifdef FTEIK_USE_INTEL
        status = ippsFlip_64f(a, b, n); 
        if (status != ippStsNoErr)
        {
            if (status == ippStsNullPtrErr)
            {
                if (a == NULL || b == NULL)
                {
                    printf("%s: Array a or b is NULL\n", __func__);
                }
            }
            if (status == ippStsSizeErr)
            {
                printf("%s: Invalid size\n", __func__);
            }
            return -1; 
        }
#else
        #pragma omp simd
        for (i=0; i<n; i++)
        {
            j = n - 1 - i;
            b[j] = a[i];
        }
#endif
    }
    return 0;
}
#endif
/*!
 * @brief Sorts an integer array in place
 *
 * @param[in] n       length of array a
 *
 * @param[inout] a    On input the array to sort.
 *                    On successful the sorted array.
 *
 * @param[in] order   ASCENDING -> sort array in ascending order (default)
 *                    DESCENDING -> sort array in descending order
 *
 * @author Ben Baker, ISTI
 *
 */
#ifdef FTEIK_USE_INTEL
int sorting_sort32i_work(const int n, int *__restrict__ a,
                         const enum sortOrder_enum order)
{
    IppStatus status;
    // Error checking
    if (n < 1)
    {
        printf("%s: Error no data to sort\n", __func__);
        return -1; 
    }
    if (order == SORT_ASCENDING)
    {
        status = ippsSortAscend_32s_I(a, n);
    }
    else
    {
        status = ippsSortDescend_32s_I(a, n);
    }
    if (ippStsNoErr != status)
    {
        printf("%s: Sort failed!\n", __func__);
        return -1;
    }
    return 0;
}
#else
static int array_reverse32i_work(const int n, int *a, int *b);
int sorting_sort32i_work(const int n, int *__restrict__ a,
                         const enum sortOrder_enum order)
{
    int ierr;
    size_t ns = (size_t) n;
    //------------------------------------------------------------------------//
    //  
    // Error checking
    if (n < 1)
    {
        printf("%s: Error no data to sort\n", __func__);
        return -1;
    }
    // Early return
    if (n == 1){return 0;}
    qsort((void *)a, (size_t) ns, sizeof(int), cmp_int);
    if (order == SORT_DESCENDING)
    {
        ierr = array_reverse32i_work(n, a, a);
        if (ierr != 0)
        {
            printf("%s: Error reversing array!\n", __func__);
            return -1;
        }
    }
    return 0;
}
/*!
 * @brief Reverse an integer array
 *
 * @param[in] n     size of arrays
 * @param[in] a     array to be reversed
 *
 * @param[out] b    reversed form of a (can be called with same argument as a)
 *
 * @result 0 indicates success
 *
 * @author Ben Baker, ISTI
 *
 */
static int array_reverse32i_work(const int n, int *a, int *b)
{
    int bj;
    int i, j;
    if (n < 1 || a == NULL || b == NULL)
    {   
        if (n < 1){printf("%s: Invalid length: %d\n", __func__, n);}
        if (a == NULL){printf("%s: Error a is NULL\n", __func__);}
        if (b == NULL){printf("%s: Error b is NULL\n", __func__);}
        return -1;
    }       
    // In place
    if (a == b)
    {   
        memcpy(b, a, (size_t) n*sizeof(int));
        j = n - 1;
        for (i=0; i<n/2; i++)
        {
            bj = b[j];
            b[j] = b[i];
            b[i] = bj; 
            j = j - 1; 
        }       
    }
    // Different memory addresses
    else
    {
        #pragma omp simd
        for (i=0; i<n; i++)
        {
            j = n - 1 - i;
            b[j] = a[i];
        }
    }
    return 0;
}
#endif
//============================================================================//
/*!
 * @brief Searches a sorted array, values, for the key.
 *
 * @param[in] key     item to sort for in values
 * @param[in] values  array sorted in ascending or descending order [n]
 * @param[in] n       number of elements in values
 * @param[in] hint    describes values as sorted in ascending or descending
 *                    order.
 *
 * @result if -1 then the key was not found in values.
 *         otherwise it is the index of an element in value which matches
 *         the key.
 *
 * @author Ben Baker
 *
 * @copyright Apache 2
 *
 */
int sorting_bsearch32i(const int key, const int *__restrict__ values,
                       const int n, const enum sortOrder_enum hint)
{
    int *item, indx;
    size_t np;
    bool lascending;
    // check the inputs
    if (n < 1 || values == NULL)
    {
        if (n < 1){printf("%s: Error no points\n", __func__);}
        if (values == NULL){printf("%s: Error values is NULL\n", __func__);}
        return -1;
    }
    // handle the edge case
    if (n == 1)
    {
        if (key == values[0])
        {
            return 0;
        }
        else
        {
            return -1;
        }
    }
    // figure out the order
    np = (size_t) n;
    lascending = false;
    if (hint == SORT_ASCENDING)
    {
        lascending = true;
    }
    else if (hint == SORT_DESCENDING)
    {
        lascending = false;
    }
    else
    {
        if (values[0] < values[n-1])
        {
            lascending = true;
        }
        else
        {
            lascending = false;
        }
    }
    // look for it
    item = NULL;
    if (lascending)
    {
        // out of bounds and edge case
        if (key < values[0]){return -1;}
        if (key > values[n-1]){return -1;}
        // look for it
        item = (int *) bsearch((const void *) &key, (const void *) values,
                               np, sizeof(int),
                               cmp_int_ascending);
    }
    else
    {
        // out of bounds and edge case
        if (key < values[n-1]){return -1;}
        if (key > values[0]){return -1;}
        // look for it
        item = (int *) bsearch((const void *) &key, (const void *) values, np,
                               sizeof(int),
                               cmp_int_descending);
    }
    if (item == NULL)
    {
        indx =-1;
    }
    else
    {
        indx = item - values;
    }
    item = NULL;
    return indx;
}
