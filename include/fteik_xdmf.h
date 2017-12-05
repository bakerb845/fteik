#ifndef FTEIK_XDMF_H__
#define FTEIK_XDMF_H__ 1

enum fteikDataType_enum
{
    FTEIK_DTYPE_UNKNOWN = 0, /*!< undefined */
    FTEIK_INT16_T = 1,       /*!< int16_t */
    FTEIK_INT32_T = 2,       /*!< int32_t */
    FTEIK_FLOAT32_T = 3,     /*!< float */
    FTEIK_FLOAT64_T = 4      /*!< double */
};

struct fteikXDMFGrid_struct
{
    char **h5flNames; /*!< HDF5 file names containing the datasets.  This is
                           an array of dimension [mDataSets]. */
    char **dataSets;  /*!< HDF5 dataset names.  This is an array of dimension
                           [mDataSets]. */
    enum fteikDataType_enum
         *precision;  /*!< Precision of datasets. */
    bool *lcell;      /*!< If true then the dataset is cell based.  Otherwise
                           it is node based.  This is an array of dimension
                           [mDataSets]. */
    char xdmfFile[PATH_MAX]; /*!< Name of xdmf file. */
    double x0;        /*!< x origin (meters). */
    double y0;        /*!< y origin (meters). */
    double z0;        /*!< z origin (meters). */
    double dx;        /*!< Grid spacing in x (meters). */
    double dy;        /*!< Grid spacing in y (meters). */
    double dz;        /*!< Grid spacing in z (meters). */
    int nx;           /*!< Number of grid points in x. */
    int ny;           /*!< Number of grid points in y. */
    int nz;           /*!< Number of grid points in z. */
    size_t mDataSets; /*!< Max number of datasets. */
    size_t nDataSets; /*!< Number of datasets in collection. */
    bool linit;       /*!< Flag indicating the struture has been initialized. */
    char pad[3];
};

#ifdef __cplusplus
extern "C"
{
#endif

/* Frees the xdmf grid structure. */
int fteik_xdmfGrid_free(struct fteikXDMFGrid_struct *xdmf);
/* Initializes the xdmf grid structure. */ 
int fteik_xdmfGrid_initialize(
    const char *xdmfFile,
    const int nx, const int ny, const int nz, 
    const double dx, const double dy, const double dz, 
    const double x0, const double y0, const double z0, 
    struct fteikXDMFGrid_struct *xdmf);
/* Adds a dataset to the xdmf grid structure. */
int fteik_xdmfGrid_add(const char *h5flName, const char *dataSet,
                       const bool lcell, const enum fteikDataType_enum prec,
                       struct fteikXDMFGrid_struct *xdmf);
/* Writes the xdmf grid structure. */
int fteik_xdmfGrid_write(const struct fteikXDMFGrid_struct xdmf);

/* TODO - deprecate */
int fteik_xdmf_writeVelocityModel(const int job, 
                                  const char *h5flName,
                                  const char *velModel,
                                  const bool lflip,
                                  const int nz, const int nx, const int ny, 
                                  const double dz, const double dx, 
                                  const double dy, 
                                  const double z0, const double x0, 
                                  const double y0);

#ifdef __cplusplus
}
#endif

#endif
