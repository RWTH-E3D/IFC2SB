#include "eaf_interface.h"

int OFCInterfaceEAF::inputGeometry(OFCGeometry &geom, const char *input_filename) {

    // clear all old geometry
    if (geom.clear() != 0) return 1;

    // set reference to maps (i.e. shortcuts)
    std::map<int, OFCNode> &nds = geom.getNdsMap();
    std::map<int, OFCTri> &trs = geom.getTrsMap();

    stl_reader::StlMesh<float, unsigned int> mesh(input_filename);

    for (size_t idx_v = 0; idx_v < mesh.num_vrts(); ++idx_v) {
        const float *pnt = mesh.vrt_coords(idx_v);

        //std::cout << idx_v << " V " << pnt[0] << " " << pnt[1] << " " << pnt[2] << std::endl;

        // nodes
        OFCNode nd;
        nd.id = idx_v;
        nd.coord[0] = pnt[0];
        nd.coord[1] = pnt[1];
        nd.coord[2] = pnt[2];

        // add node to the map
        nds[nd.id] = nd;

        //std::cout << nd.id << " " << nd.coord[0] << " " << nd.coord[1] << " " << nd.coord[2] << std::endl;
    }

    for (size_t idx_f = 0; idx_f < mesh.num_tris(); ++idx_f) {
        const unsigned int *face = mesh.tri_corner_inds(idx_f);
        // std::cout << i << " F " << face[0] << " " << face[1] << " " << face[2] << std::endl;

        // tri elements
        OFCTri tr;
        tr.id = idx_f;
        tr.nd[0] = &(nds[face[0]]);
        tr.nd[1] = &(nds[face[1]]);
        tr.nd[2] = &(nds[face[2]]);

        // add attributes
        OFCAttrib *attr = nullptr;
        attr = new OFCAttrib(ATTR_STRING);
        std::string val = "JO";
        std::string key = "GUID";
        attr->setAttrib(&val);
        attr->key = key;
        tr.attrib.push_back(attr);

        // add tri to list
        trs[tr.id] = tr;

    }

    return 0;
}

int OFCInterfaceEAF::inputGeometry(OFCGeometry &geom, const std::vector<std::tuple<double, double, double>> &vertices, const std::vector<std::tuple<int, int, int>> &faces, const std::vector<std::string> &attrs) {

    // clear all old geometry
    if (geom.clear() != 0) return 1;

    // set reference to maps (i.e. shortcuts)
    std::map<int, OFCNode> &nds = geom.getNdsMap();
    std::map<int, OFCTri> &trs = geom.getTrsMap();

    for (size_t i = 0; i < vertices.size(); i++) {

        auto P = vertices[i];

        // nodes
        OFCNode nd;
        nd.id = i;
        nd.coord[0] = std::get<0>(P);
        nd.coord[1] = std::get<1>(P);
        nd.coord[2] = std::get<2>(P);

        // add node to the map
        nds[nd.id] = nd;
    }

    for (size_t i = 0; i < faces.size(); i++) {

        auto F = faces[i];

        // tri elements
        OFCTri tr;
        tr.id = i;
        tr.nd[0] = &(nds[std::get<0>(F)]);
        tr.nd[1] = &(nds[std::get<1>(F)]);
        tr.nd[2] = &(nds[std::get<2>(F)]);

        // add attributes
        OFCAttrib *attr = nullptr;
        attr = new OFCAttrib(ATTR_STRING);
        std::string val = attrs[i];
        std::string key = "GUID";
        attr->setAttrib(&val);
        attr->key = key;
        tr.attrib.push_back(attr);

        // add tri to list
        trs[tr.id] = tr;
    }

    return 0;
}