#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <sys/stat.h>
#include <errno.h>
#include "fteik_os.h"
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
