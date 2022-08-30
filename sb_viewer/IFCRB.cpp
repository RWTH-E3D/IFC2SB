#define IFC2X3      0
#define IFC4        1
#define IFC4x1      2
#define IFC4x2      3
#define IFC4x3_RC1  4

/* define IfcVersion here */
#define IFCVERSION  IFC4

#if (IFCVERSION == IFC2X3)
#define IfcSchema Ifc2x3
#elif (IFCVERSION == IFC4)
#define IfcSchema Ifc4
#elif (IFCVERSION == IFC4x1)
#define IfcSchema Ifc4x1
#elif (IFCVERSION == IFC4x2)
#define IfcSchema Ifc4x2
#elif (IFCVERSION == IFC4x3_RC1)
#define IfcSchema Ifc4x3_rc1
#endif

#define TOSTRING(x) STRINGIFY(x)

#include "../src/viewer/ViewerHelper.h"
#include "../src/viewer/ViewerMain.h"
#include <BRepLProp_SLProps.hxx>
#include <IfcGeom.h>
#include <IfcGeomTree.h>
#include <memory>
#include <thread>
#include <utility>
#include <BRepTools_WireExplorer.hxx>

enum VIEW_MODE {
    FirstLevel, SecondLevel, RelSB
};

/* function to calculate face area */
double face_area(const TopoDS_Shape &shape) {
    GProp_GProps gprop;
    BRepGProp::SurfaceProperties(shape, gprop);
    return gprop.Mass();
}

/* function to calculate face normal */
gp_Dir face_normal(const TopoDS_Face &face) {

    Standard_Real UMin, UMax, VMin, VMax;
    BRepTools::UVBounds(face, UMin, UMax, VMin, VMax);
    BRepAdaptor_Surface surface = BRepAdaptor_Surface(face); // create surface to access the geometry of face
    gp_Dir normal = BRepLProp_SLProps(surface, UMin, VMin, 1, 1.e-3).Normal();

    if (face.Orientation() == TopAbs_FORWARD)
        return normal;
    else if (face.Orientation() == TopAbs_REVERSED)
        return normal.Reversed();
    else {
        std::cerr << "[Error] in normal calculation: " << normal.X() << " " << normal.Y() << " " << normal.Z() << std::endl;
        return normal;
    }
}

/* function to calculate face center */
gp_Pnt face_center(const TopoDS_Face &face) {
    GProp_GProps p;
    BRepGProp::SurfaceProperties(face, p);
    return p.CentreOfMass();
}

double volume(const TopoDS_Shape &shape) {
    GProp_GProps gprop;
    BRepGProp::VolumeProperties(shape, gprop);
    return gprop.Mass();
}

long projection_axis(const TopoDS_Face &face) {
    gp_Dir n = face_normal(face);
    std::vector<double> temp = {fabs(n.X()), fabs(n.Y()), fabs(n.Z())}; // determine projection axis (0,1 or 2), such that v[proj_axis] != 0
    return std::distance(temp.begin(), std::max_element(temp.begin(), temp.end()));
}

Bnd_Box aabb(const TopoDS_Shape &shape, double gap) {
    Bnd_Box bnd;
    bnd.SetGap(gap);
    BRepBndLib::Add(shape, bnd);
    return bnd;
}

TopoDS_ListOfShape loop_topo(const TopoDS_Shape &shape, const TopAbs_ShapeEnum &topologyType) {

    TopExp_Explorer Ex;
    TopoDS_ListOfShape elements;

    for (Ex.Init(shape, topologyType); Ex.More(); Ex.Next())
        elements.Append(Ex.Current());

    return elements;
}

TopoDS_ListOfShape get_ordered_vertices_from_wire(const TopoDS_Wire &W) {

    BRepTools_WireExplorer Exp(W);
    TopoDS_ListOfShape vertices;

    while (Exp.More()) {
        vertices.Append(Exp.CurrentVertex());
        Exp.Next();
    }
    return vertices;
}

std::vector<gp_Pnt> vertex_to_point_list(const TopoDS_ListOfShape &V) {
    std::vector<gp_Pnt> P;
    P.reserve(V.Size());
    for (const auto &t: V)
        P.push_back(BRep_Tool::Pnt(TopoDS::Vertex(t)));
    return P;
}

std::vector<gp_Pnt2d> project_wire_to_2D(long a, const TopoDS_Wire &W) {

    // a ... Projection Axis
    TopoDS_ListOfShape V = get_ordered_vertices_from_wire(W);
    std::vector<gp_Pnt> Points = vertex_to_point_list(V);

    std::vector<gp_Pnt2d> L;
    L.reserve(Points.size());

    if (a == 0) {
        for (const auto &P: Points)
            L.emplace_back(P.Y(), P.Z());
    } else if (a == 1) {
        for (const auto &P: Points)
            L.emplace_back(P.X(), P.Z());
    } else if (a == 2) {
        for (const auto &P: Points)
            L.emplace_back(P.X(), P.Y());
    } else
        std::cerr << "[Error] Projection axis is wrong! " << a << std::endl;

    return L;
}

double determinant_three_points(gp_Pnt2d V1, gp_Pnt2d V2, gp_Pnt2d V3) {
    gp_Vec2d v1(V2.X() - V1.X(), V2.Y() - V1.Y());
    gp_Vec2d v2(V3.X() - V2.X(), V3.Y() - V2.Y());
    return v1.X() * v2.Y() - v1.Y() * v2.X();
}

double unsigned_angle_between_points(gp_Pnt p1, gp_Pnt mid, gp_Pnt p3) {
    gp_Vec v1(mid, p1);
    gp_Vec v2(mid, p3);
    return v1.Angle(v2);
}

bool are_points_colinear(gp_Pnt p1, gp_Pnt p2, gp_Pnt p3, double tol) {

    // calculate the smallest (perpendicular) distance of p2 to the vector between p1 and p3
    gp_Vec v1(p1, p2);
    gp_Vec v2(p1, p3);
    gp_Vec v3(p2, p3);

    if (v1.Magnitude() < 1e-12 || v2.Magnitude() < 1e-12 || v3.Magnitude() < 1e-12) {
        std::cerr << "[Warning] Vector length is zero. (Self-intersection?) P1: " << p1.X() << " " << p1.Y() << " P2: " << p2.X() << " " << p2.Y() << " P3: " << p3.X() << " " << p3.Y() << std::endl;
        return true;
    }

    double dist = v1.Crossed(v2).Magnitude() / v2.Magnitude();
    if (dist > tol) return false;

    // calculate angle between edges around p2
    double a = unsigned_angle_between_points(p1, p2, p3);
    return a > 3.1398473;
}

bool are_points_coincident(gp_Pnt p1, gp_Pnt p2, double tol) {
    return p1.Distance(p2) < tol;
}

bool has_collinear_points(std::vector<gp_Pnt> Pnts, double tol) {

    gp_Pnt F = Pnts.front();
    gp_Pnt B = Pnts.back();

    Pnts.insert(Pnts.begin(), B);
    Pnts.push_back(F);

    for (unsigned int i = 1; i < Pnts.size() - 1; i++)
        if (are_points_colinear(Pnts[i - 1], Pnts[i], Pnts[i + 1], tol))
            return true;

    return false;
}

bool has_coincident_points(std::vector<gp_Pnt> Pnts, double tol) {

    Pnts.push_back(Pnts.front());

    for (unsigned int i = 0; i < Pnts.size() - 1; i++)
        if (are_points_coincident(Pnts[i], Pnts[i + 1], tol))
            return true;

    return false;
}

bool is_convex(const TopoDS_Face &face) {

    std::vector<gp_Pnt2d> l = project_wire_to_2D(projection_axis(face), BRepTools::OuterWire(face));
    if (l.size() < 3) return false;

    // first point
    const bool b = (determinant_three_points(l.back(), l[0], l[1]) >= 0); // true if det > 0

    // in between
    for (unsigned int i = 1; i < l.size() - 1; i++)
        if (b != (determinant_three_points(l[i - 1], l[i], l[i + 1]) >= 0))
            return false;

    // last
    if (b != (determinant_three_points(l[l.size() - 2], l.back(), l[0]) >= 0))
        return false;

    return true;
}


class SpaceBoundary {
public:
    TopoDS_Shape Shape;
    std::string IfcInternalOrExternalEnum;
    std::string IfcPhysicalOrVirtualEnum;
    std::string Type;
    std::string RelatedSpaceString;
    IfcUtil::IfcBaseClass *RelatedSpace;
    IfcUtil::IfcBaseClass *SB;
    unsigned int id;
    unsigned int space_id;

    SpaceBoundary(IfcUtil::IfcBaseClass *SBt, TopoDS_Shape S, std::string PhysicalOrVirtualEnum, std::string InternalOrExternalEnum, std::string Description, IfcUtil::IfcBaseClass *_space) {
        SB = SBt;
        id = SB->data().id();
        Shape = std::move(S);
        IfcPhysicalOrVirtualEnum = std::move(PhysicalOrVirtualEnum);
        IfcInternalOrExternalEnum = std::move(InternalOrExternalEnum);
        Type = std::move(Description);
        RelatedSpace = _space;
        RelatedSpaceString = RelatedSpace->declaration().name();
        space_id = RelatedSpace->data().id();
    }

    std::set<std::string> categories(std::map<IfcUtil::IfcBaseClass *, SpaceBoundary *> &M) const {

        std::set<std::string> s;
        s.insert(IfcInternalOrExternalEnum);
        s.insert(Type);

        std::set<std::string> t = s;
        for (const auto &i: t)
            for (const auto &j: t)
                if (i != j && s.find(j + "_" + i) == s.end()) {
                    s.insert(i + "_" + j);
                    s.insert(i + "_" + j + "_" + RelatedSpaceString.substr(0, 9));
                }

        s.insert(IfcInternalOrExternalEnum + "_" + RelatedSpaceString.substr(0, 9));
        s.insert("Rel_" + RelatedSpaceString.substr(0, 9));
        s.insert(IfcPhysicalOrVirtualEnum);
        if (IfcPhysicalOrVirtualEnum == "VIRTUAL" && Type == "2b") s.insert("VIRTUAL_2b");
        if (IfcPhysicalOrVirtualEnum == "PHYSICAL" && Type == "2b") s.insert("PHYSICAL_2b");
        s.insert("X_" + RelatedSpace->data().getArgument(0)->toString());
        s.insert("Z_" + SB->data().getArgument(0)->toString());
        s.insert("ZZ_#" + std::to_string(id) + "=" + "_#" + std::to_string(space_id) + "=");
        s.insert("All");

        bool parent_no = false;
        bool inner_bound = false;

#if (IFCVERSION != IFC2X3)
        // corresponding boundary
        if (SB->declaration().is("IfcRelSpaceBoundary2ndLevel")) {
            auto IfcRelSpaceBoundary2ndLevel = SB->as<IfcSchema::IfcRelSpaceBoundary2ndLevel>();
            if (IfcRelSpaceBoundary2ndLevel->hasCorrespondingBoundary()) {
                s.insert("Corresponding_Yes");
                if (M.find(IfcRelSpaceBoundary2ndLevel->CorrespondingBoundary()) != M.end()) {

                    const auto & Shape2 = M[IfcRelSpaceBoundary2ndLevel->CorrespondingBoundary()]->Shape;
                    double A1 = face_area(Shape);
                    double A2 = face_area(Shape2);
                    double a = face_normal(TopoDS::Face(Shape)).Angle(face_normal(TopoDS::Face(M[IfcRelSpaceBoundary2ndLevel->CorrespondingBoundary()]->Shape)));

                    if (fabs(A1 - A2) / std::min(A1, A2) > 0.001 || (a > 0.0174533 && a < 3.12414))
                        s.insert("Corresponding_Bad");
                    else
                        s.insert("Corresponding_Good");

                    if (loop_topo(Shape, TopAbs_VERTEX).Size() == loop_topo(Shape2, TopAbs_VERTEX).Size())
                        s.insert("Corresponding_NumberVertices_Same");
                    else
                        s.insert("Corresponding_NumberVertices_Diff");
                }
            } else
                s.insert("Corresponding_No");
        }

        // parent boundary
        if (SB->declaration().is("IfcRelSpaceBoundary2ndLevel")) {
            auto IfcRelSpaceBoundary2ndLevel = SB->as<IfcSchema::IfcRelSpaceBoundary2ndLevel>();
            if (IfcRelSpaceBoundary2ndLevel->hasParentBoundary())
                s.insert("Parent_Yes");
            else {
                s.insert("Parent_No");
                parent_no = true;
            }
        } else if (SB->declaration().is("IfcRelSpaceBoundary1stLevel")) {
            auto IfcRelSpaceBoundary1stLevel = SB->as<IfcSchema::IfcRelSpaceBoundary1stLevel>();
            if (IfcRelSpaceBoundary1stLevel->hasParentBoundary())
                s.insert("Parent_Yes");
            else {
                s.insert("Parent_No");
                parent_no = true;
            }
        } else {
            s.insert("Parent_No");
            parent_no = true;
        }

        // inner boundaries
        if (SB->declaration().is("IfcRelSpaceBoundary2ndLevel")) {
            auto IfcRelSpaceBoundary2ndLevel = SB->as<IfcSchema::IfcRelSpaceBoundary2ndLevel>();
            if ((*IfcRelSpaceBoundary2ndLevel->InnerBoundaries()).size() != 0) {
                s.insert("InnerBoundaries_Yes");
                inner_bound = true;
            }
            else
                s.insert("InnerBoundaries_No");
        } else if (SB->declaration().is("IfcRelSpaceBoundary1stLevel")) {
            auto IfcRelSpaceBoundary1stLevel = SB->as<IfcSchema::IfcRelSpaceBoundary1stLevel>();
            if ((*IfcRelSpaceBoundary1stLevel->InnerBoundaries()).size() != 0) {
                s.insert("InnerBoundaries_Yes");
                inner_bound = true;
            }
            else
                s.insert("InnerBoundaries_No");
        } else
            s.insert("InnerBoundaries_NA");

        // inner loops
        auto IfcRelSpaceBoundary = SB->as<IfcSchema::IfcRelSpaceBoundary>();
        if (IfcRelSpaceBoundary->hasConnectionGeometry()) {
            auto ConnectionGeometry = IfcRelSpaceBoundary->ConnectionGeometry()->as<IfcSchema::IfcConnectionGeometry>();

            if (ConnectionGeometry->declaration().is("IfcConnectionSurfaceGeometry")) {
                auto IfcConnectionSurfaceGeometry = IfcRelSpaceBoundary->ConnectionGeometry()->as<IfcSchema::IfcConnectionSurfaceGeometry>();
                auto IfcCurveBoundedPlane = IfcConnectionSurfaceGeometry->SurfaceOnRelatingElement()->as<IfcSchema::IfcCurveBoundedPlane>();

                if ((*IfcCurveBoundedPlane->InnerBoundaries()).size() != 0)
                    s.insert("InnerLoops_Yes");
                else
                    s.insert("InnerLoops_No");
            }
        }

        // has opening parent
        if (SB->declaration().is("IfcRelSpaceBoundary2ndLevel")) {
            auto IfcRelSpaceBoundary2ndLevel = SB->as<IfcSchema::IfcRelSpaceBoundary2ndLevel>();
            if (IfcRelSpaceBoundary2ndLevel->RelatedBuildingElement() != nullptr)
                if (IfcRelSpaceBoundary2ndLevel->RelatedBuildingElement()->declaration().is("IfcWindow") || IfcRelSpaceBoundary2ndLevel->RelatedBuildingElement()->declaration().is("IfcDoor"))
                    if (!IfcRelSpaceBoundary2ndLevel->hasParentBoundary())
                        s.insert("Opening_No_Parent");
        } else if (SB->declaration().is("IfcRelSpaceBoundary1stLevel")) {
            auto IfcRelSpaceBoundary1stLevel = SB->as<IfcSchema::IfcRelSpaceBoundary1stLevel>();
            if (IfcRelSpaceBoundary1stLevel->RelatedBuildingElement() != nullptr)
                if (IfcRelSpaceBoundary1stLevel->RelatedBuildingElement()->declaration().is("IfcWindow") || IfcRelSpaceBoundary1stLevel->RelatedBuildingElement()->declaration().is("IfcDoor"))
                    if (!IfcRelSpaceBoundary1stLevel->hasParentBoundary())
                        s.insert("Opening_No_Parent");
        }
#endif

        // RelatingSpace
        if (IfcRelSpaceBoundary->data().getArgument(4)->isNull()) s.insert("V_NoRelatingSpace");
        else {

            IfcSchema::IfcObjectDefinition *SE = nullptr;
            if (IfcRelSpaceBoundary->RelatingSpace()->declaration().is("IfcObjectDefinition")) {
                auto IfcObjectDefinition = IfcRelSpaceBoundary->RelatingSpace()->as<IfcSchema::IfcObjectDefinition>();
                if (IfcObjectDefinition->Decomposes()->size() > 0) {
                    try { SE = (*IfcObjectDefinition->Decomposes()->begin())->RelatingObject(); }
                    catch (...) { std::cout << "[Info] No RelatingStructure found for space boundary by RelatingSpace! " << IfcObjectDefinition->data().toString() << std::endl; }
                }
            }
            if (SE != nullptr) s.insert("RelStructureRS_" + SE->Name() + "_" + SE->declaration().name());
            else s.insert("RelStructureRS_None");
        }

        // storey by space
        if (IfcRelSpaceBoundary->data().getArgument(5)->isNull()) s.insert("V_NoRelatedBuildingElement");
        else {
            s.insert("V_" + IfcRelSpaceBoundary->RelatedBuildingElement()->declaration().name().substr(0, 20));
            s.insert("Y_" + IfcRelSpaceBoundary->RelatedBuildingElement()->GlobalId());
            if (parent_no && (IfcRelSpaceBoundary->RelatedBuildingElement()->declaration().is("IfcWindow") || IfcRelSpaceBoundary->RelatedBuildingElement()->declaration().is("IfcDoor")))
                s.insert("Parent_No_WinDoor");

            IfcSchema::IfcProduct *SE = nullptr;
            if (IfcRelSpaceBoundary->RelatedBuildingElement()->declaration().is("IfcElement")) {
                auto IfcElement = IfcRelSpaceBoundary->RelatedBuildingElement()->as<IfcSchema::IfcElement>();
                if (IfcElement->ContainedInStructure()->size() > 0) {
                    try { SE = (*IfcElement->ContainedInStructure()->begin())->RelatingStructure(); }
                    catch (...) { std::cout << "[Info] No RelatingStructure found for space boundary by RelatedBuildingElement! " << IfcElement->data().toString() << std::endl; }
                }
            }
            //if (SE != nullptr) s.insert("RelStructureRBE_" + SE->Name() + "_" + SE->declaration().name());
            //else s.insert("RelStructureRBE_None");
        }

        // area
        if (face_area(Shape) < 0.0025)
            s.insert("W_SmallArea");
        else
            s.insert("W_NormalArea");

        // concavity
        if (is_convex(TopoDS::Face(Shape)))
            s.insert("Convex");
        else {
            s.insert("Concave");
            if (inner_bound) s.insert("Concave_InnerBoundaries_Yes");
            else s.insert("Concave_InnerBoundaries_No");
        }

        // collinearity and coincidence
        bool has_collinearities = false;
        bool has_coincidence = false;
        double tol = 1.0e-3;

        for (const auto &wire: loop_topo(Shape, TopAbs_WIRE)) {
            auto vertices = get_ordered_vertices_from_wire(TopoDS::Wire(wire));
            auto points = vertex_to_point_list(vertices);
            if (has_collinear_points(points, tol))
                has_collinearities = true;
            if (has_coincident_points(points, tol))
                has_coincidence = true;
            if(has_coincidence && has_collinearities) break;
        }
        has_collinearities ? s.insert("Collinearity_Yes") : s.insert("Collinearity_No");
        has_coincidence ? s.insert("Coincidence_Yes") : s.insert("Coincidence_No");

        return s;
    };

    ~SpaceBoundary() = default;
};

class SpaceBoundarySelect {
public:
    TopoDS_Shape Shape;
    std::string IfcClass;

    SpaceBoundarySelect(TopoDS_Shape S, std::string _IfcClass) {
        Shape = std::move(S);
        IfcClass = std::move(_IfcClass);
    }

    std::set<std::string> categories() const {
        std::set<std::string> s;
        s.insert(IfcClass);
        if (volume(Shape) < 5) s.insert(IfcClass + "_<5m3");
        else if (volume(Shape) < 100) s.insert(IfcClass + "_5to100m3");
        else if (volume(Shape) < 500) s.insert(IfcClass + "_100to500m3");
        else s.insert(IfcClass + "_>500m3");
        return s;
    };

    ~SpaceBoundarySelect() = default;
};

void generate_shapes_from_ifc_classes(IfcParse::IfcFile *model, const std::set<std::string> &include_entities, std::list<SpaceBoundarySelect> &Spaces) {

    IfcGeom::IteratorSettings settings;
    settings.set(IfcGeom::IteratorSettings::FASTER_BOOLEANS, true);
    settings.set(IfcGeom::IteratorSettings::SEW_SHELLS, true);
    settings.set(IfcGeom::IteratorSettings::USE_WORLD_COORDS, false);
    settings.set(IfcGeom::IteratorSettings::DISABLE_TRIANGULATION, true);
    settings.set(IfcGeom::IteratorSettings::DISABLE_OPENING_SUBTRACTIONS, false);

    IfcGeom::entity_filter entity_filter;
    entity_filter.include = true;
    entity_filter.traverse = false;
    entity_filter.entity_names = include_entities;
    std::vector<IfcGeom::filter_t> filter_funcs;
    filter_funcs.emplace_back(boost::ref(entity_filter));

    IfcGeom::Iterator<double> geom_iterator(settings, model, filter_funcs, std::thread::hardware_concurrency());

    if (!geom_iterator.initialize()) {
        std::cerr << "[Info] No spaces found." << std::endl;
        return;
    }

    do {
        IfcGeom::Element<double> *geom_object = geom_iterator.get();
        const auto *o = dynamic_cast<const IfcGeom::BRepElement<double> *>(geom_object);
        TopoDS_Shape shape = o->geometry().as_compound();
        gp_Trsf trsf = o->transformation().data();
        shape.Move(trsf);
        Spaces.emplace_back(shape, o->product()->declaration().name());

    } while (geom_iterator.next());
}

int main(int argc, char **argv) {

    if (argc != 3) {
        std::cout << "usage: IfcParseExamples <filename.ifc> <\"1\" or \"2\" for first or second lvl SB or \"3\" for IfcRelSpaceBoundary>" << std::endl;
        return EXIT_FAILURE;
    }

    VIEW_MODE mode = SecondLevel;

    std::string input_mode = argv[2];
    if (input_mode == "1")
        mode = FirstLevel;
    else if (input_mode == "3")
        mode = RelSB;

    // Redirect the output (both progress and log) to stdout
    Logger::SetOutput(&std::cout, &std::cout);

    // Parse the IFC file provided in argv[1]
    auto file = std::make_unique<IfcParse::IfcFile>(argv[1]);
    // IfcParse::IfcFile file(argv[1]);
    if (!file->good()) {
        std::cout << "[Error] Unable to parse .ifc file.get()-> " << file->good().value() << std::endl;
        return EXIT_FAILURE;
    }

    std::string IfcVersionProg = TOSTRING(IfcSchema);
    std::string IfcVersionFile = file->schema()->name();
    std::cout << "[Info] Ifc file version: " << IfcVersionFile << std::endl;
    transform(IfcVersionProg.begin(), IfcVersionProg.end(), IfcVersionProg.begin(), ::toupper);
    if (IfcVersionProg != IfcVersionFile) {
        std::cout << "[Error] IfcVersion mismatch!" << std::endl;
        return EXIT_FAILURE;
    }

    /* Return if Ifcfile doesn't have any IfcRelSpaceBoundaries */
    if (file->instances_by_type("IfcRelSpaceBoundary") == nullptr) {
        std::cout << "[Error] Ifc file doesn't have IfcRelSpaceBoundary." << std::endl;
        return EXIT_SUCCESS;
    }

    IfcGeom::Kernel kernel(file.get());

    // Space Boundaries ********************************************************
    std::list<SpaceBoundary> SB;
    std::set<std::string> errors;

    std::map<IfcUtil::IfcBaseClass *, SpaceBoundary *> M;

    for (const auto &entity: *file->instances_by_type("IfcRelSpaceBoundary")) {
        auto *IfcRelSpaceBoundary = entity->as<typename IfcSchema::IfcRelSpaceBoundary>();

        if (mode == RelSB) {
            if (IfcRelSpaceBoundary->declaration().is("IfcRelSpaceBoundary1stLevel")) {
                std::cout << "[Info] SpaceBoundary (" << IfcRelSpaceBoundary->declaration().name() << ") " << IfcRelSpaceBoundary->GlobalId() << " will be skipped" << std::endl;
                continue;
            }
        } else if (mode == FirstLevel) {
            if (IfcRelSpaceBoundary->declaration().name() == "IfcRelSpaceBoundary" || IfcRelSpaceBoundary->declaration().name() == "IfcRelSpaceBoundary2ndLevel") {
                std::cout << "[Info] SpaceBoundary (" << IfcRelSpaceBoundary->declaration().name() << ") " << IfcRelSpaceBoundary->GlobalId() << " will be skipped" << std::endl;
                continue;
            }
        } else if (mode == SecondLevel) {
            if (IfcRelSpaceBoundary->declaration().name() == "IfcRelSpaceBoundary" || IfcRelSpaceBoundary->declaration().name() == "IfcRelSpaceBoundary1stLevel") {
                std::cout << "[Info] SpaceBoundary (" << IfcRelSpaceBoundary->declaration().name() << ") " << IfcRelSpaceBoundary->GlobalId() << " will be skipped" << std::endl;
                continue;
            }
        }

//        } else if (IfcRelSpaceBoundary->declaration().name() == "IfcRelSpaceBoundary2ndLevel" || IfcRelSpaceBoundary->Name() == "2ndLevel") {
//            std::cout << "[Info] Second level SpaceBoundary " << IfcRelSpaceBoundary->GlobalId() << " will be skipped" << std::endl;
//            continue;
//        }
//        else if (IfcRelSpaceBoundary->declaration().name() == "IfcRelSpaceBoundary2ndLevel" || IfcRelSpaceBoundary->Name() == "2ndLevel") {
//            std::cout << "[Info] Second level SpaceBoundary " << IfcRelSpaceBoundary->GlobalId() << " will be skipped" << std::endl;
//            continue;
//        }

        /* get loc */
        gp_Trsf trsf;

        try { auto op = IfcRelSpaceBoundary->RelatingSpace()->as<IfcSchema::IfcProduct>()->ObjectPlacement(); }
        catch (...) {
            std::cout << "[Warning] Entity linked via RelatingSpace attribute has no ObjectPlacement. SpaceBoundary " << IfcRelSpaceBoundary->GlobalId() << " will be skipped" << std::endl;
            continue;
        }

        if (IfcRelSpaceBoundary->RelatingSpace()->declaration().is("IfcSpace"))
            kernel.convert_placement(IfcRelSpaceBoundary->RelatingSpace()->as<IfcSchema::IfcSpace>()->ObjectPlacement(), trsf);
#if (IFCVERSION != IFC2X3)
        else if (IfcRelSpaceBoundary->RelatingSpace()->declaration().is("IfcExternalSpatialElement"))
            kernel.convert_placement(IfcRelSpaceBoundary->RelatingSpace()->as<IfcSchema::IfcExternalSpatialElement>()->ObjectPlacement(), trsf);
#endif
        else {
            std::cout << "[Warning] Relating Space is no IfcSpaceBoundarySelect. SpaceBoundary " << IfcRelSpaceBoundary->GlobalId() << " will be skipped" << std::endl;
            continue;
        }

        /* get shape */
        auto ConnectionGeometry = IfcRelSpaceBoundary->ConnectionGeometry()->as<IfcSchema::IfcConnectionGeometry>();
        auto ConnectionSurfaceGeometry = ConnectionGeometry->as<IfcSchema::IfcConnectionSurfaceGeometry>();
        if (ConnectionSurfaceGeometry->hasSurfaceOnRelatedElement()) {
            std::cout << "[Warning] ConnectionSurfaceGeometry has SurfaceOnRelatedElement." << std::endl;
        }

        std::vector<IfcGeom::IfcRepresentationShapeItem> Shape_List;
        try { Shape_List = kernel.convert(ConnectionSurfaceGeometry->SurfaceOnRelatingElement()); }
        catch (...) {
            std::cerr << "[Warning] Error in geometry. SpaceBoundary " << IfcRelSpaceBoundary->GlobalId() << " will be skipped" << std::endl;
            errors.insert("Error in geometry of SpaceBoundary " + IfcRelSpaceBoundary->GlobalId());
            continue;
        }

        TopoDS_Shape S = BRepBuilderAPI_Transform(Shape_List.begin()->Shape(), trsf, false);

        std::string Description;
        if (IfcRelSpaceBoundary->hasDescription())
            Description = IfcRelSpaceBoundary->Description();
        else
            Description = "NoDescription";

        /* save */
        SB.emplace_back(entity, S, IfcSchema::IfcPhysicalOrVirtualEnum::ToString(IfcRelSpaceBoundary->PhysicalOrVirtualBoundary()), IfcSchema::IfcInternalOrExternalEnum::ToString(IfcRelSpaceBoundary->InternalOrExternalBoundary()), Description,
                        IfcRelSpaceBoundary->RelatingSpace());

        M[entity] = &SB.back();
    }

    // Spaces ***************************************************************************
    std::list<SpaceBoundarySelect> Spaces;
    std::set<std::string> include_entities = {"IfcSpace", "IfcExternalSpatialElement"};
    generate_shapes_from_ifc_classes(file.get(), include_entities, Spaces);


    //***********************************************************************************
    std::cout << "\n\nVisualize shapes\n\n\n";
    std::list<viewerHelper::DisplayShapes_SB> ds;
    TopExp_Explorer Ex;

    std::map<std::string, unsigned int> count;

    for (auto &b: SB)
        for (Ex.Init(b.Shape, TopAbs_FACE); Ex.More(); Ex.Next()) {
            ds.emplace_back();
            ds.back().shape = TopoDS::Face(Ex.Current());
            ds.back().show_adaptor_face_normal = true;
            ds.back().adaptor_face_normal = face_normal(TopoDS::Face(ds.back().shape));
            ds.back().center_of_face = face_center(TopoDS::Face(ds.back().shape));
            ds.back().categories = b.categories(M);
            ds.back().transparency = 0.7;

            for (const auto &c: ds.back().categories)
                if (count.find(c) == count.end())
                    count[c] = 1;
                else
                    count[c] += 1;
        }

    for (auto &b: Spaces) {
        ds.emplace_back();
        ds.back().clr_random = true;
        ds.back().shape = b.Shape;
        ds.back().show_adaptor_face_normal = true;

        // show at least one vector of random face
        for (Ex.Init(b.Shape, TopAbs_FACE); Ex.More(); Ex.Next()) {
            ds.back().adaptor_face_normal = face_normal(TopoDS::Face(Ex.Current()));
            ds.back().center_of_face = face_center(TopoDS::Face(Ex.Current()));
            break;
        }

        ds.back().categories = b.categories();
        ds.back().transparency = 1.0;
    }

    // brute-force virtual sbs
    bool brute_force_virtuals = false;

    if (brute_force_virtuals) {
        std::cout << "Calculate virtuals by space x space operation" << "\n";
        std::map<SpaceBoundarySelect *, Bnd_Box> boxes;
        for (auto &b1: Spaces)
            boxes[&b1] = aabb(b1.Shape, 0.01);
        unsigned int i1 = 0;
        for (auto &s1: Spaces) {
            auto b1 = boxes[&s1];
            unsigned int i2 = 0;

            for (auto &s2: Spaces) {
                if (i1 < i2) {
                    auto b2 = boxes[&s2];
                    if (!b1.IsOut(b2)) {
                        TopoDS_Shape shell1 = loop_topo(s1.Shape, TopAbs_SHELL).First();
                        TopoDS_Shape shell2 = loop_topo(s2.Shape, TopAbs_SHELL).First();
                        TopoDS_ListOfShape L = loop_topo(BRepAlgoAPI_Common(shell1, shell2).Shape(), TopAbs_FACE);

                        for (const auto &l: L) {
                            ds.emplace_back();
                            ds.back().shape = TopoDS::Face(l);
                            ds.back().transparency = 0.7;
                            std::set<std::string> s;
                            s.insert("SpacexSpaceVirtuals");
                            ds.back().categories = s;

                            for (const auto &c: ds.back().categories)
                                if (count.find(c) == count.end())
                                    count[c] = 1;
                                else
                                    count[c] += 1;
                        }
                    }
                }
                i2++;
            }
            i1++;
        }
    }

    std::cout << "SB, total" << "\t: " << SB.size() << "\n";
    std::cout << "SpaceBoundarySelect, total" << "\t: " << Spaces.size() << "\n";
    for (const auto &c: count)
        if (c.first.substr(0, 2) != "W_" && c.first.substr(0, 2) != "X_" && c.first.substr(0, 2) != "Y_" && c.first.substr(0, 2) != "Z_" && c.first.substr(0, 3) != "ZZ_")
            std::cout << c.first << "\t: " << c.second << "\n";

    std::cout << "\n";
    for (const auto &e: errors)
        std::cout << e << "\n";

    ViewerMain::start_viewer(ds);

    return EXIT_SUCCESS;
}