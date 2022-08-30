#ifndef EAF_INTERFACE_H
#define EAF_INTERFACE_H

#include "octree_includes.h"

// base class for interface functions
class OFCInterfaceEAF : public OFCInterface {
public:
    // input of geometry (pure virtual function)
    int inputGeometry(OFCGeometry &geom, const char *input_filename);

    int inputGeometry(OFCGeometry &geom, const std::vector<std::tuple<double, double, double>> &vertices, const std::vector<std::tuple<int, int, int>> &faces, const std::vector<std::string> &attrs);

};

#endif // EAF_INTERFACE_H