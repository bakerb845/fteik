#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
/*
#include <iostream>
#include "Xdmf.hpp"
#include "XdmfArray.hpp"
#include "XdmfArrayType.hpp"
#include "XdmfAttribute.hpp"
#include "XdmfAttributeType.hpp"
#include "XdmfAttributeCenter.hpp"
#include "XdmfDomain.hpp"
#include "XdmfHeavyDataDescription.hpp"
#include "XdmfHDF5Controller.hpp"
#include "XdmfInformation.hpp"
#include "XdmfRegularGrid.hpp"
#include "XdmfPlaceholder.hpp"
#include "XdmfWriter.hpp"
*/

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

int fteik_xdmf_writeVelocityModel(const int nz, const int nx, const int ny, 
                                  const double dz, const double dx, 
                                  const double dy, 
                                  const double z0, const double x0, 
                                  const double y0)
{
    xmlTextWriterPtr writer;
    xmlBufferPtr buf;
    char cdim[128], cori[128], cspace[128];
    char *xmlmsg; 
    size_t msglen;
    int rc;
    int ncellx, ncelly, ncellz;
    ncellz = nz - 1;
    ncellx = nx - 1;
    ncelly = ny - 1;
    memset(cdim, 0, 128*sizeof(char));
    sprintf(cdim, "%d %d %d", ncellz, ncelly, ncellx);

    memset(cspace, 0, 128*sizeof(char));
    sprintf(cspace, "%e %e %e", -dz, dy, dx); // TODO - make this an option

const char *ch5 = "debug.h5:/VelocityModels/debugModel\0"; // TODO - change this

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
                                     BAD_CAST "Structured Grid\0");
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "GridType\0",
                                     BAD_CAST "Uniform\0");
    // <Topology>
    rc = xmlTextWriterStartElement(writer, BAD_CAST "Topology\0"); 
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "TopologyType\0",
                                     BAD_CAST TOPOLOGY_TYPE);
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Dimensions\0",
                                     BAD_CAST cdim);
    rc = xmlTextWriterEndElement(writer); //</Topology>
    // <Geometry>
    rc = xmlTextWriterStartElement(writer, BAD_CAST "Geometry\0");
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "GeometryType\0",
                                     BAD_CAST GEOMETRY_TYPE);
    // Origin
    rc = xmlTextWriterStartElement(writer, BAD_CAST "DataItem\0");
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Name\0",
                                     BAD_CAST "Origin\0");
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Dimensions\0",
                                     BAD_CAST "3\0");
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "DataType\0",
                                     BAD_CAST "Float\0");
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Precision\0",
                                     BAD_CAST "8\0");
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Format\0",
                                     BAD_CAST "XML\0");
    memset(cori, 0, 128*sizeof(char));
    sprintf(cori, "%e %e %e", z0, y0, x0); 
    rc = xmlTextWriterWriteString(writer, BAD_CAST cori);
    rc = xmlTextWriterEndElement(writer); //</DataItem>
    // Grid spacing
    rc = xmlTextWriterStartElement(writer, BAD_CAST "DataItem\0");
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Name\0",
                                     BAD_CAST "Spacing\0");
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Dimensions\0",
                                     BAD_CAST "3\0");
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "DataType\0",
                                     BAD_CAST "Float\0");
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Precision\0",
                                     BAD_CAST "8\0");
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Format\0",
                                     BAD_CAST "XML\0");
    rc = xmlTextWriterWriteString(writer, BAD_CAST cspace);
    rc = xmlTextWriterEndElement(writer); //</DataItem>
    // </Geometry> 
    rc = xmlTextWriterEndElement(writer); //</Geometry>
   
    // <Attribute>
    rc = xmlTextWriterStartElement(writer, BAD_CAST "Attribute\0");
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Name\0",
                                     BAD_CAST "Compressional Velocity\0");
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Type\0",
                                     BAD_CAST "Scalar\0");
    rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "Center\0",
                                     BAD_CAST "Cell\0");
    // <DataItem>
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
    return 0;
}

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
