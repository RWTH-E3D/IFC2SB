#ifndef BASE_INTERFACE_H
#define BASE_INTERFACE_H

#include "octree_includes.h"

// base class for interface functions
class OFCInterface {

public:
    OFCInterface();

    virtual ~OFCInterface();

    // input of geometry (pure virtual function)
    virtual int inputGeometry(OFCGeometry &geom, const char *input_filename) = 0;

};

#endif //BASE_INTERFACE_H