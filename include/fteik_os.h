#ifndef FTEIK_OS_H__
#define FTEIK_OS_H__ 1
#include <limits.h>
#ifdef __cplusplus
extern "C"
{
#endif

bool fteik_os_isfile(const char *filenm);
int fteik_os_basename_work(const char *path, char baseName[PATH_MAX]);
bool fteik_os_path_isdir(const char *dirnm);
int fteik_os_dirname_work(const char *path, char dirName[PATH_MAX]);
int fteik_os_mkdir(const char *dirnm);
int fteik_os_makedirs(const char *path);

#ifdef __cplusplus
}
#endif

#endif
