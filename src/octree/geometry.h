#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "octree_includes.h"

// --- definition of attribute types -------------------------------------------
#define ATTR_UNDEF    0
#define ATTR_INTEGER  1
#define ATTR_DOUBLE   2
#define ATTR_FLOAT    3
#define ATTR_STRING   4

// --- base class template containing attributes -------------------------------
class OFCAttrib {
private:
    // diffent types of attributes
    int val_integer;
    double val_double;
    float val_float;
    std::string val_string;

    // validity of attribute
    bool valid;
public:
    const int type;    // attribute type (according to defines)
    std::string key;        // string containing key to attribute (specifying which attribute is named)

    // get and set
    void *getAttrib();

    int setAttrib(void *_val);

    // constructor and destructor
    OFCAttrib(const int _type);

    ~OFCAttrib();
};

// --- struct containing nodal information --------------------------------------
class OFCNode {
public:
    int id;             // node ID
    double coord[3];    // global coordinates of the node (x,y,z)
    std::list<OFCAttrib *> attrib; // list containing attributes associated to nodes

    OFCNode();

    ~OFCNode();
};

// --- struct containing tri element information -------------------------------
class OFCTri {
public:
    int id;             // tri ID
    OFCNode *nd[3];     // pointer to nodes (nd1, nd2, nd3)
    std::list<OFCAttrib *> attrib; // list containing attributes associated to tri

    OFCTri();

    ~OFCTri();
};

// --- class containing all geometry information -------------------------------
class OFCGeometry {
public:
    // constructor and destructor
    OFCGeometry();

    ~OFCGeometry();

    // clear all geometry and attributes data
    int clear();

    // get references to maps
    std::map<int, OFCNode> &getNdsMap();

    std::map<int, OFCTri> &getTrsMap();

    const std::map<int, OFCNode> &getConstNdsMap() const;

    const std::map<int, OFCTri> &getConstTrsMap() const;

    // change geometry
    int translate(const double &tl_x, const double &tl_y, const double &tl_z);

    int scale(const double &sc_x, const double &sc_y, const double &sc_z);

    int scale(const double &sc);


private:
    // map containing all the nodes
    std::map<int, OFCNode> nds;

    // map containing all the tri elements
    std::map<int, OFCTri> trs;

};


// some nice and usefull properties --------------------------------------------

inline OFCNode diff(OFCNode nd1, const OFCNode &nd2) {
    for (int i = 0; i < 3; ++i) nd1.coord[i] -= nd2.coord[i];
    return nd1;
}

inline OFCNode cross(const OFCNode &v1, const OFCNode &v2) {
    OFCNode v;
    v.coord[0] = v1.coord[1] * v2.coord[2] - v1.coord[2] * v2.coord[1];
    v.coord[1] = v1.coord[2] * v2.coord[0] - v1.coord[0] * v2.coord[2];
    v.coord[2] = v1.coord[0] * v2.coord[1] - v1.coord[1] * v2.coord[0];
    return v;
}

inline double scal(const OFCNode &v1, const OFCNode &v2) {
    return v1.coord[0] * v2.coord[0] + v1.coord[1] * v2.coord[1] + v1.coord[2] * v2.coord[2];
}

inline double norm(const OFCNode &v1) {
    return sqrt(v1.coord[0] * v1.coord[0] + v1.coord[1] * v1.coord[1] + v1.coord[2] * v1.coord[2]);
}

#endif //GEOMETRY_H
