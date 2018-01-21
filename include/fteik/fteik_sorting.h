#ifndef FTEIK_SORTING_H__
#define FTEIK_SORTING_H__ 1
#include "fteik/fteik_config.h"
enum sortOrder_enum 
{
    SORT_ASCENDING = 0,  /*!< Sort in ascending order */
    SORT_DESCENDING = 1, /*!< Sort in descending order */
    SORT_EITHER = 2      /*!< When checking if sorted - array can be ascending
                              or descending order */
};

#ifdef __cplusplus
extern "C"
{
#endif
int sorting_sort32i_work(const int n, int *__restrict__ a,
                         const enum sortOrder_enum order);
int sorting_sort64f_work(const int n,
                         double *__restrict__ a,
                         const enum sortOrder_enum order);
int sorting_bsearch32i(const int key, const int *__restrict__ values,
                       const int n, const enum sortOrder_enum hint);
int sorting_argsort64f_work(const int n,
                            const double *__restrict__ a,
                            const enum sortOrder_enum order,
                            int *__restrict__ perm);
int sorting_argsort32i_work(const int n,
                            const int *__restrict__ ia,
                            const enum sortOrder_enum order,
                            int *__restrict__ perm);
#ifdef __cplusplus
}
#endif
#endif
