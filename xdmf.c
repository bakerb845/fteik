#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <limits.h>
#include "fteik_h5io.h"

int fteik_xdmf_initialize(const char *h5dir,
                          const char *h5flName, const char *projnm,
                          const int nz, const int nx, const int ny,
                          struct xdmf_struct *xdmf)
{
    char fname[PATH_MAX];
    size_t lenos;
    memset(xdmf, 0, sizeof(struct xdmf_struct));
    strcpy(xdmf->h5fl, h5flName);
    strcpy(xdmf->projnm, projnm);
    xdmf->nz = nz;
    xdmf->nx = nx;
    xdmf->ny = ny;
    xdmf->nnodes = nx*ny*nz;
    xdmf->nelem = (nx - 1)*(ny - 1)*(nz - 1);
    xdmf->ipart = 1;
    memset(fname, 0, PATH_MAX*sizeof(char));
    if (h5dir != NULL)
    {
        lenos = strlen(h5dir); 
        if (lenos > 0)
        {
            strcpy(fname, h5dir);
            if (fname[lenos-1] != '/'){fname[lenos] = '/';}
        }
        else
        {
            strcpy(fname, "./\0");
        }
    }
    else
    {
        strcpy(fname, "./\0");
    }
    strcat(fname, projnm);
    strcat(fname, ".xdmf\0");
    xdmf->xdmfl = fopen(fname, "w");
    fprintf(xdmf->xdmfl, "<?xml version=\"1.0\" ?>\n");
    fprintf(xdmf->xdmfl, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
    fprintf(xdmf->xdmfl, "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.0\">\n");
    fprintf(xdmf->xdmfl, "  <Domain>\n");
    //fprintf(xdmf->xdmfl, "    <Grid Collection=\"%s\" Name=\"Partition %d\">\n",
    //        xdmf->projnm, xdmf->ipart);
    xdmf->lfileOpen = true;
    xdmf->linit = true;
    return 0;
}

int fteik_xdmf_addLevelSet(const struct xdmf_struct xdmf,
                           const int sweep, const char *levelName)
{
    const char *fcnm = "fteik_xdmf_addLevelSet\0";
    FILE *fl;
    if (!xdmf.linit || !xdmf.lfileOpen)
    {
        printf("%s: Error xdmf structure not initialized\n", fcnm);
        return -1;
    }
    if (sweep < 1 || sweep > 8)
    {
        printf("%s: Error sweep %d must be between 1 and 8\n", fcnm, sweep);
        return -1;
    }
    fl = xdmf.xdmfl;
    if (sweep == 1)
    {
        fprintf(fl, "    <Grid Name=\"Level Sets\" GridType=\"Uniform\">\n");
    }
    fprintf(fl, "      <Topology TopologyType=\"Hexahedron\" NumberOfElements=\"%d\">\n", xdmf.nelem);
    fprintf(fl, "        <DataItem Dimensions=\"%d\" Numbertype=\"Int\" Precision=\"4\" Format=\"HDF\">\n",
                         8*xdmf.nelem);
    fprintf(fl, "          %s:/Geometry/Connectivity\n", xdmf.h5fl);
    fprintf(fl, "        </DataItem>\n");
    fprintf(fl, "      </Topology>\n");
    fprintf(fl, "      <Geometry GeometryType=\"X_Y_Z\">\n");
    fprintf(fl, "        <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n",
                         xdmf.nnodes);
    fprintf(fl, "          %s:/Geometry/xLocations\n", xdmf.h5fl);
    fprintf(fl, "        </DataItem>\n");
    fprintf(fl, "        <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n",
                         xdmf.nnodes);
    fprintf(fl, "          %s:/Geometry/yLocations\n", xdmf.h5fl);
    fprintf(fl, "        </DataItem>\n");
    fprintf(fl, "        <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n",
                         xdmf.nnodes);
    fprintf(fl, "          %s:/Geometry/zLocations\n", xdmf.h5fl);
    fprintf(fl, "        </DataItem>\n");
    fprintf(fl, "      </Geometry>\n");
    fprintf(fl, "      <Attribute Name=\"Sweep %d\" AttributeType=\"Scalar\" Center=\"Node\">\n", sweep);
    fprintf(fl, "        <DataItem Dimensions=\"%d\" NumberType=\"Int\" Precision=\"4\" Format=\"HDF\">\n",
                         xdmf.nnodes);
    fprintf(fl, "          %s:/LevelSet/%s\n", xdmf.h5fl, levelName);
    fprintf(fl, "        </DataItem>\n");
    fprintf(fl, "      </Attribute>\n");
    if (sweep == 8)
    {
        fprintf(fl, "    </Grid>\n");
    }
    fl = NULL;
    return 0;
}

int fteik_xdmf_addVelocityModel(const struct xdmf_struct xdmf,
                                const char *vmodelName)
{
    const char *fcnm = "fteik_xdmf_addVelocityModel\0";
    FILE *fl;
    if (!xdmf.linit || !xdmf.lfileOpen)
    {   
        printf("%s: Error xdmf structure not initialized\n", fcnm);
        return -1; 
    }   
    fl = xdmf.xdmfl;
    fprintf(fl, "    <Grid Name=\"%s\" GridType=\"Uniform\">\n", vmodelName);
    fprintf(fl, "      <Topology TopologyType=\"Hexahedron\" NumberOfElements=\"%d\">\n", xdmf.nelem);
    fprintf(fl, "        <DataItem Dimensions=\"%d\" Numbertype=\"Int\" Precision=\"4\" Format=\"HDF\">\n",
                         8*xdmf.nelem);
    fprintf(fl, "          %s:/Geometry/Connectivity\n", xdmf.h5fl);
    fprintf(fl, "        </DataItem>\n");
    fprintf(fl, "      </Topology>\n");
    fprintf(fl, "      <Geometry GeometryType=\"X_Y_Z\">\n");
    fprintf(fl, "        <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n",
                         xdmf.nnodes);
    fprintf(fl, "          %s:/Geometry/xLocations\n", xdmf.h5fl);
    fprintf(fl, "        </DataItem>\n");
    fprintf(fl, "        <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n",
                         xdmf.nnodes);
    fprintf(fl, "          %s:/Geometry/yLocations\n", xdmf.h5fl);
    fprintf(fl, "        </DataItem>\n");
    fprintf(fl, "        <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n",
                         xdmf.nnodes);
    fprintf(fl, "          %s:/Geometry/zLocations\n", xdmf.h5fl);
    fprintf(fl, "        </DataItem>\n");
    fprintf(fl, "      </Geometry>\n");
    fprintf(fl, "      <Attribute Name=\"Velocity Model (m/s)\" AttributeType=\"Scalar\" Center=\"Cell\">\n");
    fprintf(fl, "        <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n",
                         xdmf.nelem);
    fprintf(fl, "          %s:/VelocityModel/%s\n", xdmf.h5fl, vmodelName);
    fprintf(fl, "        </DataItem>\n");
    fprintf(fl, "      </Attribute>\n");
    fprintf(fl, "    </Grid>\n");
    fl = NULL;
    return 0;
}

int fteik_xdmf_addTravelTimes(const struct xdmf_struct xdmf,
                              const char *ttName)
{
    const char *fcnm = "fteik_xdmf_addVelocityModel\0";
    FILE *fl;
    if (!xdmf.linit || !xdmf.lfileOpen)
    {   
        printf("%s: Error xdmf structure not initialized\n", fcnm);
        return -1; 
    }   
    fl = xdmf.xdmfl;
    fprintf(fl, "    <Grid Name=\"%s\" GridType=\"Uniform\">\n", ttName);
    fprintf(fl, "      <Topology TopologyType=\"Hexahedron\" NumberOfElements=\"%d\">\n", xdmf.nelem);
    fprintf(fl, "        <DataItem Dimensions=\"%d\" Numbertype=\"Int\" Precision=\"4\" Format=\"HDF\">\n",
                         8*xdmf.nelem);
    fprintf(fl, "          %s:/Geometry/Connectivity\n", xdmf.h5fl);
    fprintf(fl, "        </DataItem>\n");
    fprintf(fl, "      </Topology>\n");
    fprintf(fl, "      <Geometry GeometryType=\"X_Y_Z\">\n");
    fprintf(fl, "        <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n",
                         xdmf.nnodes);
    fprintf(fl, "          %s:/Geometry/xLocations\n", xdmf.h5fl);
    fprintf(fl, "        </DataItem>\n");
    fprintf(fl, "        <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n",
                         xdmf.nnodes);
    fprintf(fl, "          %s:/Geometry/yLocations\n", xdmf.h5fl);
    fprintf(fl, "        </DataItem>\n");
    fprintf(fl, "        <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n",
                         xdmf.nnodes);
    fprintf(fl, "          %s:/Geometry/zLocations\n", xdmf.h5fl);
    fprintf(fl, "        </DataItem>\n");
    fprintf(fl, "      </Geometry>\n");
    fprintf(fl, "      <Attribute Name=\"%s Travel Times (s)\" AttributeType=\"Scalar\" Center=\"Node\">\n", ttName);
    fprintf(fl, "        <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n",
                         xdmf.nnodes);
    fprintf(fl, "          %s:/TravelTimes/%s\n", xdmf.h5fl, ttName);
    fprintf(fl, "        </DataItem>\n");
    fprintf(fl, "      </Attribute>\n");
    fprintf(fl, "    </Grid>\n");
    fl = NULL;
    return 0;
}

int fteik_xdmf_finalize(struct xdmf_struct *xdmf)
{
    const char *fcnm = "fteik_xdmf_finalize\0";
    if (!xdmf->linit || !xdmf->lfileOpen)
    {
        printf("%s: Error xdmf structure not initialized\n", fcnm);
        return -1;
    } 
    //fprintf(xdmf->xdmfl, "    </Grid>\n");
    fprintf(xdmf->xdmfl, "  </Domain>\n");
    fprintf(xdmf->xdmfl, "</Xdmf>\n");
    fclose(xdmf->xdmfl);
    memset(xdmf, 0, sizeof(struct xdmf_struct));
    return 0;
}
