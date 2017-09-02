#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <libgen.h>
#include <limits.h>
#include <math.h>
#include <libxml/parser.h>
#include <libxml/xpath.h>
#include <libxml/tree.h>
#include <libxml/encoding.h>
#include <libxml/xmlwriter.h>
#include <libxml/parser.h>


#define XML_ENCODING "utf-8"
#define XDMF_VERSION "3"
#define TOPOLOGY_TYPE "3DCoRectMesh"
#define GEOMETRY_TYPE "ORIGIN_DXDYDZ"


static int writeOrigin(const double z0, const double x0, const double y0, 
                       xmlTextWriterPtr *writer);
static int writeGridSpacing(const double dz, const double dx, const double dy, 
                            xmlTextWriterPtr *writer);
static int writeTopology(const int nz, const int nx, const int ny,
                         xmlTextWriterPtr *writer);

/*!
 * @brief Writes the XDMF 
 *
 * @param[in] job          If job is 1 then the file file is initialized. \n
 *                         If job is 2 then the velocity model is appended. \n
 *                         If job is 3 then the XDMF file is closed.
 * @param[in] h5flName     Name of the HDF5 file.  This is only accessed when
 *                         job = 1.
 * @param[in] velModel     Name of the velocity model in the HDF5 file. 
 *                         This is only accessed when job = 2.
 * @param[in] lflip        If true then the velocity model in the solver
 *                         in the solver is z positive down.  In this
 *                         case the z direction must be reversed for plotting
 *                         in ParaView. \n
 *                         Otherwise, dz will not be negated.
 *                         This is only accessed when job = 1.
 * @param[in] nz           Number of z grid points.  There are nz - 1 cells
 *                         in the velocity model.  This is only accessed when
 *                         job = 1.
 * @param[in] nx           Number of x grid points.  There are nx - 1 cells
 *                         in the velocity model.  This is only accessed when
 *                         job = 1.
 * @param[in] ny           Number of y grid points.  There are ny - 1 cells
 *                         in the velocity model.  This is only accessed when
 *                         job = 1.
 * @param[in] z0           z model origin (meters).  This is only accessed when
 *                         job = 1.
 * @param[in] x0           x model origin (meters).  This is only accessed when
 *                         job = 1.
 * @param[in] y0           y model origin (meters).  This is only accessed when
 *                         job = 1.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker
 *
 */
int fteik_xdmf_writeVelocityModel(const int job,
                                  const char *h5flName,
                                  const char *velModel,
                                  const bool lflip,
                                  const int nz, const int nx, const int ny, 
                                  const double dz, const double dx, 
                                  const double dy, 
                                  const double z0, const double x0, 
                                  const double y0)
{
    static xmlTextWriterPtr writer;
    static xmlBufferPtr buf;
    static char cbase[PATH_MAX];
    char ch5[PATH_MAX], cdim[128], clevel[128];
    static char *h5work, *h5fl;
    char *xmlmsg; 
    size_t msglen;
    int i, rc;
    static int ncellx, ncelly, ncellz;
    static double dzOut, dxOut, dyOut;
    static bool linit = false;

    if (job == 1)
    {
        if (linit)
        {
            fprintf(stderr, "%s: Error must close xdmf file\n", __func__);
            return -1;
        }
        ncellz = nz - 1;
        ncellx = nx - 1;
        ncelly = ny - 1;

        dzOut = dz;
        dxOut = dx;
        dyOut = dy;
        if (lflip)
        {    dzOut =-dzOut;
             dyOut =-dyOut;
        }
        memset(cbase, 0, PATH_MAX*sizeof(char));
        h5work = strdup(h5flName);
        h5fl = basename(h5work); //h5flName);
        sprintf(cbase, "%s:/VelocityModels/", h5fl); //, velModel);
        // Create a new XML buffer to which the XML document will be written
        buf = xmlBufferCreate();
        if (buf == NULL)
        {   
            fprintf(stderr, "%s: Error creating XML buffer!", __func__);
            return -1;
        }
        // Create a new xmlWriter for uri with no compression
        writer = xmlNewTextWriterMemory(buf, 0);
        if (writer == NULL)
        {
            fprintf(stderr, "%s: Error creating xml writer", __func__);
            return -1;
        }
        xmlTextWriterSetIndentString(writer, "  "); //"\t");
        // Start the document with default xml version
        rc = xmlTextWriterStartDocument(writer, NULL, XML_ENCODING, NULL);
        if (rc < 0)
        {
            fprintf(stderr, "%s: Error starting writer", __func__);
            return -1;
        }
        // <Xdmf>
        xmlTextWriterSetIndent(writer, 1);
        rc = xmlTextWriterStartElement(writer, BAD_CAST "Xdmf\0");
        rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "xmlns:xi\0",
                                  BAD_CAST "http://www.w3.org/2003/XInclude\0");
        rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Version\0",
                                         BAD_CAST XDMF_VERSION); 
        // <Domain>
        rc = xmlTextWriterStartElement(writer, BAD_CAST "Domain\0"); 
        // <Grid>
        rc = xmlTextWriterStartElement(writer, BAD_CAST "Grid\0");
        rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Name\0",
                                         BAD_CAST "Level Sets\0");
        rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "GridType\0",
                                         BAD_CAST "Uniform\0");
        // <Topology>
        writeTopology(ncellz+1, ncellx+1, ncelly+1, &writer);
        // <Geometry>
        rc = xmlTextWriterStartElement(writer, BAD_CAST "Geometry\0");
        rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "GeometryType\0",
                                         BAD_CAST GEOMETRY_TYPE);
        // Origin
        writeOrigin(z0, x0, y0, &writer);
        // Grid spacing
        writeGridSpacing(dzOut, dxOut, dyOut, &writer);
        // </Geometry> 
        rc = xmlTextWriterEndElement(writer); //</Geometry>

        // Write level schedules
        for (i=1; i<=8; i++)
        {
            memset(clevel, 0, 128*sizeof(char));
            memset(ch5, 0, PATH_MAX*sizeof(char));
            sprintf(clevel, "Level Schedule %d", i); 
            sprintf(ch5, "%s:/LevelSchedules/Sweep_%d", h5fl, i); 
            // <Attribute>
            rc = xmlTextWriterStartElement(writer, BAD_CAST "Attribute\0");
            rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Name\0",
                                             BAD_CAST clevel);
            rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Type\0",
                                             BAD_CAST "Scalar\0");
            rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Center\0",
                                             BAD_CAST "Node\0");
            // <DataItem>
            memset(cdim, 0, 128*sizeof(char));
            sprintf(cdim, "%d %d %d", ncellz+1, ncelly+1, ncellx+1);
            rc = xmlTextWriterStartElement(writer, BAD_CAST "DataItem\0");
            rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "DataType\0",
                                             BAD_CAST "Short\0");
            rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Dimensions\0",
                                             BAD_CAST cdim);
            rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Format\0",
                                             BAD_CAST "HDF\0");
            rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Precision\0",
                                             BAD_CAST "2\0");
            rc = xmlTextWriterWriteString(writer, BAD_CAST ch5);
            rc = xmlTextWriterEndElement(writer); //</DataItem>
            // </DataItem>
            rc = xmlTextWriterEndElement(writer); //</Attribute>
            if (rc < 0)
            {
                fprintf(stderr, "%s: Error writing model\n", __func__);
            }
        }
        if (rc < 0)
        {
            fprintf(stderr, "%s: Error ending geometry\n", __func__);
            return -1;
        }
/*
        rc = xmlTextWriterEndElement(writer); //</Grid>
        //--------------------------------------------------------------------//
        // <Grid>
        rc = xmlTextWriterStartElement(writer, BAD_CAST "Grid\0");
        rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Name\0",
                                         BAD_CAST "Structured Grid\0");
        rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "GridType\0",
                                         BAD_CAST "Uniform\0");
        // <Topology>
        writeTopology(ncellz, ncellx, ncelly, &writer);
        // <Geometry>
        rc = xmlTextWriterStartElement(writer, BAD_CAST "Geometry\0");
        rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "GeometryType\0",
                                         BAD_CAST GEOMETRY_TYPE);
        // Origin
        writeOrigin(z0, x0, y0, &writer);
        // Grid spacing
        writeGridSpacing(dzOut, dxOut, dyOut, &writer);
        // </Geometry> 
        rc = xmlTextWriterEndElement(writer); //</Geometry>
        if (rc < 0)
        {
            fprintf(stderr, "%s: Error ending geometry\n", __func__);
            return -1;
        }
        free(h5work);
        h5fl = NULL;
*/
        linit = true;
    }
    if (job == 2)
    {   
/*
        if (!linit)
        {
            fprintf(stderr, "%s: xdmf file not initialized\n", __func__); 
            return -1; 
        }
        memset(ch5, 0, PATH_MAX*sizeof(char));
        sprintf(ch5, "%s%s", cbase, velModel);
        // <Attribute>
        rc = xmlTextWriterStartElement(writer, BAD_CAST "Attribute\0");
        rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Name\0",
                                         BAD_CAST velModel);
        rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Type\0",
                                         BAD_CAST "Scalar\0");
        rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Center\0",
                                         BAD_CAST "Cell\0");
        // <DataItem>
        memset(cdim, 0, 128*sizeof(char));
        sprintf(cdim, "%d %d %d", ncellz, ncelly, ncellx);
        rc = xmlTextWriterStartElement(writer, BAD_CAST "DataItem\0");
        rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "DataType\0",
                                         BAD_CAST "Short\0");
        rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Dimensions\0",
                                         BAD_CAST cdim);
        rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Format\0",
                                         BAD_CAST "HDF\0");
        rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Precision\0",
                                         BAD_CAST "2\0");
        rc = xmlTextWriterWriteString(writer, BAD_CAST ch5);
        rc = xmlTextWriterEndElement(writer); //</DataItem>
        // </DataItem>
        rc = xmlTextWriterEndElement(writer); //</Attribute>
        if (rc < 0)
        {
            fprintf(stderr, "%s: Error writing model\n", __func__);
        }
*/
    }
    if (job == 3)
    {   
        if (!linit)
        {
            fprintf(stderr, "%s: xdmf file not initialized\n", __func__); 
            return -1; 
        }
        memset(ch5, 0, PATH_MAX*sizeof(char));
        sprintf(ch5, "%s:/TravelTimeFields/%s", h5fl, velModel);
        // <Attribute>
        rc = xmlTextWriterStartElement(writer, BAD_CAST "Attribute\0");
        rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Name\0",
                                         BAD_CAST velModel);
        rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Type\0",
                                         BAD_CAST "Scalar\0");
        rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Center\0",
                                         BAD_CAST "Node\0");
        // <DataItem>
        memset(cdim, 0, 128*sizeof(char));
        sprintf(cdim, "%d %d %d", ncellz+1, ncelly+1, ncellx+1);
        rc = xmlTextWriterStartElement(writer, BAD_CAST "DataItem\0");
        rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "DataType\0",
                                         BAD_CAST "Float\0");
        rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Dimensions\0",
                                         BAD_CAST cdim);
        rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Format\0",
                                         BAD_CAST "HDF\0");
        rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Precision\0",
                                         BAD_CAST "4\0");
        rc = xmlTextWriterWriteString(writer, BAD_CAST ch5);
        rc = xmlTextWriterEndElement(writer); //</DataItem>
        // </DataItem>
        rc = xmlTextWriterEndElement(writer); //</Attribute>
        if (rc < 0)
        {
            fprintf(stderr, "%s: Error writing model\n", __func__);
        }
    }
    else if (job == 4)
    {
        if (!linit)
        {
            fprintf(stderr, "%s: xdmf file not initialized\n", __func__); 
            return -1;
        }
        memset(ch5, 0, PATH_MAX*sizeof(char));
        sprintf(ch5, "%s%s", cbase, velModel);
        // <Attribute>
        rc = xmlTextWriterStartElement(writer, BAD_CAST "Attribute\0");
        rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Name\0",
                                         BAD_CAST velModel);
        rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Type\0",
                                         BAD_CAST "Scalar\0");
        rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Center\0",
                                         BAD_CAST "Cell\0");
        // <DataItem>
        memset(cdim, 0, 128*sizeof(char));
        sprintf(cdim, "%d %d %d", ncellz, ncelly, ncellx);
        rc = xmlTextWriterStartElement(writer, BAD_CAST "DataItem\0");
        rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "DataType\0",
                                         BAD_CAST "Short\0");
        rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Dimensions\0",
                                         BAD_CAST cdim);
        rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Format\0",
                                         BAD_CAST "HDF\0");
        rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Precision\0",
                                         BAD_CAST "2\0");
        rc = xmlTextWriterWriteString(writer, BAD_CAST ch5);
        rc = xmlTextWriterEndElement(writer); //</DataItem>
        // </DataItem>
        rc = xmlTextWriterEndElement(writer); //</Attribute>
        if (rc < 0)
        {
            fprintf(stderr, "%s: Error writing model\n", __func__);
        }
    }
    else if (job == 5)
    {
        rc = xmlTextWriterEndElement(writer); //</Grid>
        rc = xmlTextWriterEndElement(writer); //</Domain>
        rc = xmlTextWriterEndElement(writer); //</Xdmf>
        if (rc < 0)
        {
            fprintf(stderr, "%s: Error closing Xdmf\n", __func__);
            return -1;
        }
        // Clean up
        xmlFreeTextWriter(writer);
        xmlCleanupCharEncodingHandlers();
        // Finally copy the char * XML message
        msglen = xmlStrlen(buf->content); //strlen((const char *)buf->content);
        xmlmsg = (char *) calloc((size_t) (msglen+1), sizeof(char));
        strncpy(xmlmsg, (const char *)buf->content, msglen);

        printf("%s\n", xmlmsg);
        FILE *ofl = fopen("debug.xdmf", "w");
        fprintf(ofl, xmlmsg);
        fclose(ofl);
        free(xmlmsg);

        xmlCleanupParser();
        xmlBufferFree(buf);
        xmlDictCleanup();
        xmlCleanupThreads();

        free(h5work);
        h5fl = NULL;
        linit = false;
    }
    return 0;
}
//============================================================================//

//============================================================================//
/*!
 * @brief Convenience function to write the origin.
 */
static int writeOrigin(const double z0, const double x0, const double y0,
                       xmlTextWriterPtr *writer)
{
    char cori[128];
    int rc;
    // Origin
    rc = xmlTextWriterStartElement(*writer, BAD_CAST "DataItem\0");
    rc = xmlTextWriterWriteAttribute(*writer, BAD_CAST "Name\0",
                                     BAD_CAST "Origin\0");
    rc = xmlTextWriterWriteAttribute(*writer, BAD_CAST "Dimensions\0",
                                     BAD_CAST "3\0");
    rc = xmlTextWriterWriteAttribute(*writer, BAD_CAST "DataType\0",
                                     BAD_CAST "Float\0");
    rc = xmlTextWriterWriteAttribute(*writer, BAD_CAST "Precision\0",
                                     BAD_CAST "8\0");
    rc = xmlTextWriterWriteAttribute(*writer, BAD_CAST "Format\0",
                                     BAD_CAST "XML\0");
    memset(cori, 0, 128*sizeof(char));
    sprintf(cori, "%e %e %e", z0, y0, x0); 
    rc = xmlTextWriterWriteString(*writer, BAD_CAST cori);
    rc = xmlTextWriterEndElement(*writer); //</DataItem>
    return rc;
}
//============================================================================//
/*!
 * @brief Convenience function to write the grid spacing.
 */
static int writeGridSpacing(const double dz, const double dx, const double dy,
                            xmlTextWriterPtr *writer)
{
    char cspace[128];
    int rc;
    memset(cspace, 0, 128*sizeof(char));
    sprintf(cspace, "%e %e %e", dz, dy, dx);
    // Grid spacing
    rc = xmlTextWriterStartElement(*writer, BAD_CAST "DataItem\0");
    rc = xmlTextWriterWriteAttribute(*writer, BAD_CAST "Name\0",
                                     BAD_CAST "Spacing\0");
    rc = xmlTextWriterWriteAttribute(*writer, BAD_CAST "Dimensions\0",
                                     BAD_CAST "3\0");
    rc = xmlTextWriterWriteAttribute(*writer, BAD_CAST "DataType\0",
                                     BAD_CAST "Float\0");
    rc = xmlTextWriterWriteAttribute(*writer, BAD_CAST "Precision\0",
                                     BAD_CAST "8\0");
    rc = xmlTextWriterWriteAttribute(*writer, BAD_CAST "Format\0",
                                     BAD_CAST "XML\0");
    rc = xmlTextWriterWriteString(*writer, BAD_CAST cspace);
    rc = xmlTextWriterEndElement(*writer); //</DataItem>
    return rc;
}
//============================================================================//
/*!
 * @brief Convenience function to write the topology.
 */
static int writeTopology(const int nz, const int nx, const int ny,
                         xmlTextWriterPtr *writer)
{
    char cdim[128];
    int rc;
    memset(cdim, 0, 128*sizeof(char));
    sprintf(cdim, "%d %d %d", nz, ny, nx);
    // Topology
    rc = xmlTextWriterStartElement(*writer, BAD_CAST "Topology\0"); 
    rc = xmlTextWriterWriteAttribute(*writer, BAD_CAST "TopologyType\0",
                                     BAD_CAST TOPOLOGY_TYPE);
    rc = xmlTextWriterWriteAttribute(*writer, BAD_CAST "Dimensions\0",
                                     BAD_CAST cdim);
    rc = xmlTextWriterEndElement(*writer); //</Topology>
    return rc;
}


//============================================================================//
/*
extern "C"
int fteik_xdmf_writeVelocityModel2(const int nz, const int nx, const int ny,
                                  const double dz, const double dx,
                                  const double dy,
                                  const double z0, const double x0,
                                  const double y0)
{
    int ncellx, ncelly, ncellz;
    ncellz = nz - 1;
    ncellx = nx - 1;
    ncelly = ny - 1;
 
    shared_ptr<XdmfDomain> root = XdmfDomain::New();
    // Information
    shared_ptr<XdmfInformation> info = XdmfInformation::New("Name", "Model");
    //info->setName("Model");
    root->insert(info); // XdmfDomain is the root of the tree - this adds info
    // Dimensions
    shared_ptr<XdmfArray> dimensions = XdmfArray::New();
    dimensions->pushBack((unsigned int) ncellz);
    dimensions->pushBack((unsigned int) ncelly);
    dimensions->pushBack((unsigned int) ncellx);
    // Orgin
    shared_ptr<XdmfArray> origin = XdmfArray::New();
    origin->setName("Origin");
    origin->pushBack(z0);
    origin->pushBack(y0);
    origin->pushBack(x0);
    // Brick size
    shared_ptr<XdmfArray> brick = XdmfArray::New();
    brick->setName("Spacing");
    brick->pushBack(dz);
    brick->pushBack(dy);
    brick->pushBack(dx);
    // Grid
    shared_ptr<XdmfRegularGrid> grid = XdmfRegularGrid::New(brick, dimensions, origin);
    grid->setName("Structured Grid");
    root->insert(grid);
    // Attributes
    shared_ptr<XdmfAttribute> attr = XdmfAttribute::New();
    attr->setName("Velocity (m/s)");
    attr->setCenter(XdmfAttributeCenter::Cell()); // Could be ::Node()
    attr->setType(XdmfAttributeType::Scalar());
    // Describe the HDF5 data item
    std::string hdf5FilePath = "dum file path"; //./debug.h5";
    std::string dataSetPath = "/VelocityModels/debugModel";

    //shared_ptr<const XdmfArrayType> p(Int16() atype;//const XdmfArrayType> Int16();
    //shared_ptr<const XdmfArrayType> Int16 p; //("Short", 2, type); //Signed);
    shared_ptr<const XdmfArrayType> dataType = XdmfArrayType::Int16();
    std::vector<unsigned int> start;
    start.push_back((unsigned int) 0);
    start.push_back((unsigned int) 0);
    start.push_back((unsigned int) 0);
    std::vector<unsigned int> stride;
    stride.push_back((unsigned int) 1);
    stride.push_back((unsigned int) 1);
    stride.push_back((unsigned int) 1);
    std::vector<unsigned int> dims;
    dims.push_back((unsigned int) ncellz);
    dims.push_back((unsigned int) ncelly);
    dims.push_back((unsigned int) ncellx);
    std::vector<unsigned int> dataSpaceDimensions;
    dataSpaceDimensions.push_back(ncellz);
    dataSpaceDimensions.push_back(ncelly);
    dataSpaceDimensions.push_back(ncellx);
    shared_ptr<XdmfPlaceholder> vmodelH5
         = XdmfPlaceholder::New("/VelocityModels/debugModel", dataType, //dataSetPath, dataType,
                                   start, stride, dims, dataSpaceDimensions);
    std::cout << vmodelH5->getName() << std::endl;
    attr->insert(vmodelH5);

    start.clear();
    stride.clear();
    dims.clear();
    dataSpaceDimensions.clear();


    grid->insert(attr);

    // Generate a writer
    shared_ptr<XdmfWriter> writer = XdmfWriter::New("debug.xmf"); 
    root->accept(writer);
 
    return 0;
}
*/
