#include "octree_includes.h"

OFCAttrib::OFCAttrib(const int _type)
        : type(_type), valid(false) {
    // set key to standard value
    key = "";

    val_integer = 0;
    val_double = 0;
    val_float = 0;
    val_string = "";

    switch (type) {
        case (ATTR_INTEGER) :
        case (ATTR_DOUBLE)  :
        case (ATTR_FLOAT)   :
        case (ATTR_STRING)  :
            // everything ok
            valid = true;
            break;
        default :
            std::cerr << "ERROR in OFCAttrib::OFCAttrib(..) : Type " << type << " not defined as attribute at address " << this << "!\n\n";
    }

}

OFCAttrib::~OFCAttrib() = default;

int OFCAttrib::setAttrib(void *_val) {
    if (!valid) return 1;

    switch (type) {
        case (ATTR_INTEGER) :
            val_integer = *((int *) _val);
            break;
        case (ATTR_DOUBLE)  :
            val_double = *((double *) _val);
            break;
        case (ATTR_FLOAT)   :
            val_float = *((float *) _val);
            break;
        case (ATTR_STRING)  :
            val_string = *((std::string *) _val);
            break;
        default :
            return 1;
    }
    return 0;
}

void *OFCAttrib::getAttrib() {
    switch (type) {
        case (ATTR_INTEGER) :
            return &val_integer;
        case (ATTR_DOUBLE)  :
            return &val_double;
        case (ATTR_FLOAT)   :
            return &val_float;
        case (ATTR_STRING)  :
            return &val_string;
        default :
            return nullptr;
    }
}

OFCNode::OFCNode() {
    id = -1;
    coord[0] = 0;
    coord[1] = 0;
    coord[2] = 0;
}

OFCNode::~OFCNode() = default;

OFCTri::OFCTri() {
    id = -1;
    nd[0] = nullptr;
    nd[1] = nullptr;
    nd[2] = nullptr;
}

OFCTri::~OFCTri() = default;

OFCGeometry::OFCGeometry() {
}

OFCGeometry::~OFCGeometry() {
    clear();
}

// clear all geometry
int OFCGeometry::clear() {
    // delete node attributes
    std::map<int, OFCNode>::iterator nit;
    for (nit = nds.begin(); nit != nds.end(); ++nit) {
        OFCNode &nd = nit->second;

        if (!nd.attrib.empty()) {
            std::list<OFCAttrib *>::iterator lit;
            for (lit = nd.attrib.begin(); lit != nd.attrib.end(); ++lit) {
                if (*lit != 0) delete *lit;
                *lit = 0;
            }
            nd.attrib.clear();
        }
    }

    // delete tri attributes
    std::map<int, OFCTri>::iterator tit;
    for (tit = trs.begin(); tit != trs.end(); ++tit) {
        OFCTri &tr = tit->second;

        if (!tr.attrib.empty()) {
            std::list<OFCAttrib *>::iterator lit;
            for (lit = tr.attrib.begin(); lit != tr.attrib.end(); ++lit) {
                if (*lit != 0) delete *lit;
                *lit = 0;
            }
            tr.attrib.clear();
        }
    }

    trs.clear();
    nds.clear();
    return 0;
}

// get references to maps
std::map<int, OFCNode> &OFCGeometry::getNdsMap() {
    return nds;
}

std::map<int, OFCTri> &OFCGeometry::getTrsMap() {
    return trs;
}

const std::map<int, OFCNode> &OFCGeometry::getConstNdsMap() const {
    return nds;
}

const std::map<int, OFCTri> &OFCGeometry::getConstTrsMap() const {
    return trs;
}

int OFCGeometry::translate(const double &tl_x, const double &tl_y, const double &tl_z) {
    std::map<int, OFCNode>::iterator nit;
    for (nit = nds.begin(); nit != nds.end(); ++nit) {
        OFCNode &nd = nit->second;
        nd.coord[0] += tl_x;
        nd.coord[1] += tl_y;
        nd.coord[2] += tl_z;
    }
    return 0;
}

int OFCGeometry::scale(const double &sc_x, const double &sc_y, const double &sc_z) {
    std::map<int, OFCNode>::iterator nit;
    for (nit = nds.begin(); nit != nds.end(); ++nit) {
        OFCNode &nd = nit->second;
        nd.coord[0] *= sc_x;
        nd.coord[1] *= sc_y;
        nd.coord[2] *= sc_z;
    }
    return 0;
}

int OFCGeometry::scale(const double &sc) {
    return scale(sc, sc, sc);
}

