// Copyright 2022 Eric Fichter
#ifndef PRODUCT_H
#define PRODUCT_H

#include "headers.h"

class Product {
public:

    Product(IfcUtil::IfcBaseEntity *_ifcproduct, std::string _guid, TopoDS_Shape _shape);

    bool operator==(const Product &P) const { return guid == P.guid; }

    IfcUtil::IfcBaseEntity *ifcproduct;
    TopoDS_Shape shape;
    const std::string guid;
    bool hasOnlyPolygons;
    bool valid;
    std::list<IfcUtil::IfcBaseClass *> IfcOpeningElements;
    std::list<std::pair<IfcUtil::IfcBaseClass *, IfcUtil::IfcBaseClass *>> VoidFillingElements; // IfcOpeningElement in pair with IfcWindow/IfcDoor
    std::list<oFace *> orig_faces; // only valid before adding window/door openings

    //! Print attributes.
    void Log(unsigned int i) const;

    //! Returns true if shell(s) are closed.
    bool CheckShells() const;

    //! Returns name of the IfcClass.
    std::string IfcClass() const;

    //! Returns true if product is of specified IfcClass.
    bool IsIfcClass(const std::string& s) const;

    //! String with basic info.
    std::string Info() const;
};

#endif //PRODUCT_H