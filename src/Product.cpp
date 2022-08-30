// Copyright 2022 Eric Fichter
#include "Product.h"

Product::Product(IfcUtil::IfcBaseEntity *_ifcproduct, std::string _guid, TopoDS_Shape _shape) : ifcproduct(_ifcproduct), guid(std::move(_guid)), shape(std::move(_shape)) {

    hasOnlyPolygons = true;
    valid = false;
    IfcOpeningElements.clear();
    VoidFillingElements.clear();
    orig_faces.clear();
}

void Product::Log(unsigned int i) const {
    std::string s = ifcproduct->data().toString();
    s.resize(55);
    std::cout << std::setfill('0') << std::setw(6) << i;
    std::cout << "\t" << s;
    std::cout << "\t" << std::boolalpha << hasOnlyPolygons;
    std::cout << "\t" << std::boolalpha << valid;
    std::cout << "\t" << std::setfill('0') << std::setw(11) << Kernel::hash(shape);
    std::cout << "\t" << std::boolalpha << shape.IsNull();
    std::cout << "\t" << std::setfill('0') << std::setw(3) << VoidFillingElements.size();
    std::cout << "\t" << std::setfill('0') << std::setw(3) << IfcOpeningElements.size();
    std::cout << "\t" << std::setfill('0') << std::setw(4) << orig_faces.size();
    if (!shape.IsNull()) {
        Bnd_Box box = Kernel::aabb(shape, 1.0e-5);
        gp_Pnt P1 = box.CornerMin();
        gp_Pnt P2 = box.CornerMax();
        std::cout << "\t" << std::setfill('0') << std::setw(6) << P1.X() << "|" << std::setfill('0') << std::setw(6) << P1.Y() << "|" << std::setfill('0') << std::setw(6) << P1.Z();
        std::cout << "\t" << std::setfill('0') << std::setw(6) << P2.X() << "|" << std::setfill('0') << std::setw(6) << P2.Y() << "|" << std::setfill('0') << std::setw(6) << P2.Z();
    }
    std::cout << "\n";
}

bool Product::CheckShells() const {

    if (shape.IsNull()) {
        std::cerr << "[Warning] Shape is Null " << ifcproduct->data().toString() << std::endl;
        return false;
    }

    // First check if all faces belong to a shell
    TopTools_IndexedDataMapOfShapeListOfShape faceShellMap;
    TopExp::MapShapesAndAncestors(shape, TopAbs_FACE, TopAbs_SHELL, faceShellMap);

    TopExp_Explorer Ex;
    for (Ex.Init(shape, TopAbs_FACE); Ex.More(); Ex.Next()) {
        const TopTools_ListOfShape &shells = faceShellMap.FindFromKey(Ex.Current());
        if (shells.Size() == 0) {
            std::cerr << "[Warning] Face does not belong to a shell in product " << ifcproduct->data().toString() << std::endl;
            return false;
        }
    }

    // Check if all shells are closed
    for (Ex.Init(shape, TopAbs_SHELL); Ex.More(); Ex.Next()) {

        BRepCheck_Shell c = BRepCheck_Shell(TopoDS::Shell(Ex.Current()));

        BRepCheck_Status Closed = c.Closed();
        BRepCheck_Status Orientation = c.Orientation();

        if (Closed != BRepCheck_NoError) {
            std::cerr << "[Warning] Shell is not closed: " << ifcproduct->data().toString() << std::endl;
            return false;
        }

        if (Orientation != BRepCheck_NoError) {
            std::cerr << "[Warning] Shell is not consistently oriented: " << ifcproduct->data().toString() << std::endl;
            return false;
        }

        if (Kernel::volume(Ex.Current()) <= 0) {
            std::cerr << "[Warning] Shell normals pointing inwards: " << ifcproduct->data().toString() << std::endl;
            return false;
        }

        ShapeAnalysis_Shell a;
        a.LoadShells(TopoDS::Shell(Ex.Current()));
        a.CheckOrientedShells(TopoDS::Shell(Ex.Current()), true, true);
        if (a.HasBadEdges()) {
            std::cerr << "[Warning] Shell has bad edges: " << ifcproduct->data().toString() << std::endl;
            return false;
        }
        if (a.HasFreeEdges()) {
            std::cerr << "[Warning] Shell has free edges: " << ifcproduct->data().toString() << std::endl;
            return false;
        }
    }

    return true;
}

std::string Product::IfcClass() const {

    if (ifcproduct == nullptr) {
        std::cerr << "[Warning] No IFC product." << std::endl;
        return "No IFC product";
    } else return ifcproduct->declaration().name();
}

bool Product::IsIfcClass(const std::string& s) const {

    if (ifcproduct == nullptr) {
        std::cerr << "[Warning] No IFC product." << std::endl;
        return false;
    } else return ifcproduct->declaration().is(s);
}

std::string Product::Info() const { return IfcClass() + ", " + guid; }

