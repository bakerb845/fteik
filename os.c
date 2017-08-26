#if defined WINNT || defined WIN32 || defined WIN64
#include <windows.h>
#endif
#include <libgen.h> // Might not be a windows thing
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <limits.h>
#include <sys/stat.h>
#include <errno.h>
#include "fteik_os.h"

#ifndef MIN
#define MIN(x,y) (((x) < (y)) ? (x) : (y))
#endif
/*! 
 * @brief Tests if filenm is a file.
 * 
 * @param[in] filenm    Name of file to test.
 * 
 * @result  If true then filenm is an existing file. \n
 *          If false then filenm is not a file.
 *
 * @author Ben Baker, ISTI
 *
 * @copyright ISTI distributed under Apache 2
 *
 */
bool fteik_os_isfile(const char *filenm)
{
    struct stat info;
    if (filenm == NULL){return false;}
    if (strlen(filenm) == 0){return false;}
    // Doesn't exist
    if (stat(filenm, &info) ==-1)
    {   
        return false;
    }   
    // Exists -> check it is a file
    else
    {   
        if (S_ISREG(info.st_mode))
        {
            return true;
        }
        else
        {
            return false;
        }
    }   
}
//============================================================================//
/*! 
 * @brief Tests if dirnm is a directory.
 * 
 * @param[in] dirnm    Name of directory to test.
 *
 * @result If true true then dirnm is an existing directory. \n
 *         If false then dirnm is not a directory.
 * 
 * @author Ben Baker
 *
 * @copyright ISTI distributed under Apache 2.
 *
 */
bool fteik_os_path_isdir(const char *dirnm)
{
    struct stat s;
    int err;
    if (dirnm == NULL){return false;}
    if (strlen(dirnm) == 0){return false;}
    err = stat(dirnm, &s);
    // Doesn't exist
    if (err == -1)
    {
        if (ENOENT == errno)
        {
            return false;
        }
        // Exists
        else
        {
            return true;
        }
    }
    // Exists
    else
    {
        // Test it is a directory
        if (S_ISDIR(s.st_mode))
        {
            return true;
        }
        else
        {
            return false;
        }
    }
}
//============================================================================//
/*!
 * @brief Breaks a null-terminated pathname string into a directory.
 *        For example, /usr/include/lapacke.h will return 
 *        /usr/include
 *
 * @param[in] path      Name of file path to be parsed.
 *
 * @param[out] dirName  Name of directory containing path.
 * 
 * @result 0 indicates success.
 *
 * @author Ben Baker
 *
 * @copyright ISTI distributed under Apache 2.
 *
 */
int fteik_os_dirname_work(const char *path, char dirName[PATH_MAX])
{
    char *dname, *temp;
    size_t lenos;
    // Initialize output
    memset(dirName, 0, PATH_MAX*sizeof(char));
    // Input string is NULL - this may yield some weirdness
    if (path == NULL)
    {
        strcpy(dirName, ".\0");
        return -1;
    }
    // Save copy of string - glibc can destroy it
    dname = NULL;
    lenos = strlen(path);
    temp = (char *) calloc(lenos+1, sizeof(char));
    strcpy(temp, path);
    dname = dirname(temp);
    if (dname != NULL)
    {
        lenos = strlen(dname);
        strncpy(dirName, dname, MIN(lenos, PATH_MAX));
    }
    free(temp);
    return 0;
}
//============================================================================//
/*! 
 * Makes a directory named dirnm with full permissions.
 *
 * @param[in] dirnm    Name of directory to make.
 *
 * @result 0 if success.
 *
 * @author Ben Baker
 *
 * @copyright ISTI distributed under Apache 2
 *
 */
int fteik_os_mkdir(const char *dirnm)
{
    int ierr;
    if (dirnm == NULL)
    {   
        fprintf(stderr, "Directory name is NULL");
        return -1;
    }   
    if (strlen(dirnm) == 0)
    {   
        fprintf(stderr, "Directory name is empty");
        return -1;
    }   
#if defined WINNT || defined WIN32 || defined WIN64
    ierr = mkdir(dirnm);
#else
    ierr = mkdir(dirnm, 0777);
#endif
    if (ierr != 0)
    {   
        fprintf(stderr, "Error making directory %s", dirnm);
        return -1;
    }   
    return 0;
}
//============================================================================//
/*!
 * @brief Recursive directory creation function.
 *
 * @param[in] path    Directory tree to make.
 *
 * @result ISCL error code where 0 indicates success.
 *
 * @author Ben Baker
 *
 * @copyright ISTI distributed under Apache 2.
 *
 */
int fteik_os_makedirs(const char *path)
{
    char *where, *dname, *work, directory[PATH_MAX];
#if defined WINNT || defined WIN32 || defined WIN64
    const char find[] = "\\/";
#else
    const char find[] = "/\0";
#endif
    int indx, lenos;
    int ierr; //enum isclError_enum ierr;
    //------------------------------------------------------------------------//    //
    // Errors 
    ierr = 0; //ISCL_SUCCESS;
    if (path == NULL)
    {   
        fprintf(stderr, "Directory name is NULL");
        return -1; //ISCL_NULL_PATH;
    }
    lenos = strlen(path);
    if (lenos == 0)
    {   
        fprintf(stderr, "Directory name is empty");
        return -1; //ISCL_EMPTY_PATH;
    }
    if (lenos > PATH_MAX - 2)
    {   
        fprintf(stderr, "Directory %s is too long", path);
        return -1; //ISCL_INVALID_INPUT;
    }
    // Already exists
    if (fteik_os_path_isdir(path)){return 0;}
    // Initialize
    work = (char *) calloc(lenos+2, sizeof(char));
    strcpy(work, path); 
    memset(directory, 0, PATH_MAX*sizeof(char));
    dname = work;
#if defined WINNT || defined WIN32 || defined WIN64
    dname[lenos] = '\0'; // No idea what to do
#else
    dname[lenos] = '/'; // Try to catch the final case
#endif
    where = strpbrk(dname, find);
    while ((where != NULL))
    {
        indx = where - dname;
        lenos = strlen(directory);
        strncat(directory, dname, indx);
        // Windows does nonsense like: C:
#if defined WINNT || defined WIN32 || defined WIN64
        if (directory[strlen(directory) - 1] != ':')
#endif
        {
            // If directory doesn't exist then make it
            if (!fteik_os_path_isdir(directory))
            {
                ierr = fteik_os_mkdir(directory);
                if (ierr != 0) //ISCL_SUCCESS)
                {
                    fprintf(stderr, "Error making subdirectory: %s",
                            directory);
                    fprintf(stderr, "Error making directory: %s", path);
                    free(dname);
                    return ierr;
                }
            } // End check on if directory exists
        }
        // Add directory delimiter
#if defined WINNT || defined WIN32 || defined WIN64
        strcat(directory, "//\0");
#else
        strcat(directory, "/\0");
#endif
        dname = dname + indx + 1;
        where = strpbrk(dname, find);
    } // End while
    // Add name of final subdirectory and make it
    strcat(directory, dname);
    if (!fteik_os_path_isdir(directory))
    {
         ierr = fteik_os_mkdir(directory);
         if (ierr != 0) //ISCL_SUCCESS)
         {
             fprintf(stderr, "Error making directory: %s\n", directory);
             free(work);
             return ierr;
         }
    }
    // Free space
    dname = NULL;
    free(work);
    return ierr;
}
