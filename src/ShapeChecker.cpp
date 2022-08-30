// Copyright 2022 Eric Fichter
#include "ShapeChecker.h"

ShapeChecker::ShapeChecker(const TopoDS_Shape &_shape) : shape(_shape), tol(1.0e-7) {
    errors.clear();
}

ShapeChecker::~ShapeChecker() {}

void ShapeChecker::full_check() {

    vertices = get_vertices();
    edges = get_edges();
    wires = get_wires();
    faces = get_faces();
    shells = get_shells();
    solids = get_solids();
    compsolids = get_compsolids();
    compounds = get_compounds();

    BRepCheck_Analyzer A(shape);
    analyzer = &A;
    isValid = analyzer->IsValid();

    std::cerr << "Shape is Valid?: " << std::boolalpha << isValid << "\n";

    std::cerr << "Number of vertices (redundant): " << vertices.Size() << "\n";
    std::cerr << "Number of edges (redundant): " << edges.Size() << "\n";
    std::cerr << "Number of wires: " << wires.Size() << "\n";
    std::cerr << "Number of faces: " << faces.Size() << "\n";
    std::cerr << "Number of shells: " << shells.Size() << "\n";
    std::cerr << "Number of solids: " << solids.Size() << "\n";
    std::cerr << "Number of compsolids: " << compsolids.Size() << "\n";
    std::cerr << "Number of compounds: " << compounds.Size() << "\n";

    check_vertices();
    check_edges();
    check_wires();
    check_faces();
    check_shells();

    print_errors();

    //delete analyzer;
}

void ShapeChecker::dump_topology(int level) {
    dump_topology(shape, level);
}

void ShapeChecker::dump_topology(const TopoDS_Shape &shape, int level) {

    if (level == 0) {
        ShapeAnalysis_ShapeContents a;
        a.Perform(shape);
        std::cerr << "Number of vertices (redundant): " << a.NbVertices() << "\n";
        std::cerr << "Number of edges (redundant): " << a.NbEdges() << "\n";
        std::cerr << "Number of shared edges (not redundant): " << a.NbSharedEdges() << "\n"; // Unique Edges?
        std::cerr << "Number of wires: " << Topo(shape).wires().Size() << "\n";
        std::cerr << "Number of faces: " << Topo(shape).faces().Size() << "\n";
        std::cerr << "Number of shells: " << Topo(shape).shells().Size() << "\n";
        std::cerr << "Number of solids: " << Topo(shape).solids().Size() << "\n";
    }

    if (shape.IsNull()) return;

    TopAbs_ShapeEnum s = shape.ShapeType();
    if (s == TopAbs_VERTEX) {
        gp_Pnt pnt = BRep_Tool::Pnt(TopoDS::Vertex(shape));
        std::cerr << std::string(level * 2, ' ') << shapeTypeString(shape.ShapeType()) << " " << shape.HashCode(INT_MAX) << "\t\tOr: " << shape.Orientation() << "\t\tLoc: " << shape.Location().HashCode(INT_MAX) << "\t\t" << pnt.X() << " " << pnt.Y() << " " << pnt.Z() << "\n";
    } else if (s == TopAbs_EDGE) {
        std::cerr << std::string(level * 2, ' ') << shapeTypeString(shape.ShapeType()) << " " << shape.HashCode(INT_MAX) << "\t\tOr: " << shape.Orientation() << "\t\tLoc: " << shape.Location().HashCode(INT_MAX) << "\t\tDgn: " << BRep_Tool::Degenerated(TopoDS::Edge(shape)) << "\t\tLin: "
                  << Kernel::is_edge_line(TopoDS::Edge(shape)) << "\t\tLength: " << Kernel::length(TopoDS::Edge(shape)) << "\n";
    } else if (s == TopAbs_FACE) {
        gp_Dir n = Kernel::face_normal(TopoDS::Face(shape));
        std::cerr << std::string(level * 2, ' ') << shapeTypeString(shape.ShapeType()) << " " << shape.HashCode(INT_MAX) << "\t\tOr: " << shape.Orientation() << "\t\tLoc: " << shape.Location().HashCode(INT_MAX) << "\t\tPln: " << Kernel::face_is_planar(TopoDS::Face(shape)) << "\t\tN: "
                  << Kernel::round_double_two_digits(n.X()) << " " << Kernel::round_double_two_digits(n.Y()) << " " << Kernel::round_double_two_digits(n.Z()) << "\t\tOuterWire: " << Kernel::hash(BRepTools::OuterWire(TopoDS::Face(shape))) << "\n";
    } else if (s == TopAbs_WIRE) {
        std::set<unsigned int> hashes;
        for (auto &v: Topo(shape).vertices()) hashes.insert(Kernel::hash(v));
        auto nv1 = hashes.size();
        auto nv2 = Topo(shape).ordered_vertices_of_wire().Size();
        std::cerr << std::string(level * 2, ' ') << shapeTypeString(shape.ShapeType()) << " " << shape.HashCode(INT_MAX) << "\t\tOr: " << shape.Orientation() << "\t\tLoc: " << shape.Location().HashCode(INT_MAX) << "\t\tnVertices: " << nv1 << "\t\tnVertices: " << nv2 << "\n";
    } else if (s == TopAbs_SHELL) {
        std::cerr << std::string(level * 2, ' ') << shapeTypeString(shape.ShapeType()) << " " << shape.HashCode(INT_MAX) << "\t\tOr: " << shape.Orientation() << "\t\tLoc: " << shape.Location().HashCode(INT_MAX) << "\t\tVol: " << Kernel::volume(shape) << "\n";
    } else
        std::cerr << std::string(level * 2, ' ') << shapeTypeString(shape.ShapeType()) << " " << shape.HashCode(INT_MAX) << "\t\tOr: " << shape.Orientation() << "\t\tLoc: " << shape.Location().HashCode(INT_MAX) << "\n";

    TopoDS_Iterator it = TopoDS_Iterator(shape);

    while (it.More()) {
        TopoDS_Shape shp = it.Value();
        dump_topology(shp, level + 1);
        it.Next();
    }
}

std::string ShapeChecker::shapeTypeString(const TopAbs_ShapeEnum &st) {

    std::string s;

    if (st == TopAbs_VERTEX)
        s = "Vertex";
    else if (st == TopAbs_SOLID)
        s = "Solid";
    else if (st == TopAbs_EDGE)
        s = "Edge";
    else if (st == TopAbs_FACE)
        s = "Face";
    else if (st == TopAbs_SHELL)
        s = "Shell";
    else if (st == TopAbs_WIRE)
        s = "Wire";
    else if (st == TopAbs_COMPOUND)
        s = "Compound.";
    else if (st == TopAbs_COMPSOLID)
        s = "Compsolid.";
    else
        s = "?";

    return s;
}

TopoDS_ListOfShape ShapeChecker::loop_topo(const TopAbs_ShapeEnum &topologyType, bool ignore_same_hashes) {

    TopExp_Explorer Ex;
    TopoDS_ListOfShape elements;

    if (ignore_same_hashes) {
        std::set<Standard_Integer> hashes;
        for (Ex.Init(shape, topologyType); Ex.More(); Ex.Next()) {
            auto search = hashes.find(Ex.Current().HashCode(INT_MAX));
            if (search == hashes.end()) {
                hashes.insert(Ex.Current().HashCode(INT_MAX));
                elements.Append(Ex.Current());
            }
        }
    } else
        for (Ex.Init(shape, topologyType); Ex.More(); Ex.Next())
            elements.Append(Ex.Current());

    return elements;
}

void ShapeChecker::get_results_from_analyzer(const TopoDS_Shape &shape) {

    auto res = analyzer->Result(shape);
    auto status = res->Status();
    auto iterator = BRepCheck_ListIteratorOfListOfStatus(status);

    while (iterator.More()) {
        BRepCheck_Status v = iterator.Value();
        std::cerr << "\t" << int_to_status(v) << std::endl;
        add_to_errors(shape, v);
        iterator.Next();
    }

}

std::string ShapeChecker::int_to_status(int i) {

    switch (i) {
        case BRepCheck_NoError:
            return "BRepCheck_NoError";
            break;
        case BRepCheck_InvalidPointOnCurve:
            return "BRepCheck_InvalidPointOnCurve";
            break;
        case BRepCheck_InvalidPointOnCurveOnSurface:
            return "BRepCheck_InvalidPointOnCurveOnSurface";
            break;
        case BRepCheck_InvalidPointOnSurface:
            return "BRepCheck_InvalidPointOnSurface";
            break;
        case BRepCheck_No3DCurve:
            return "BRepCheck_No3DCurve";
            break;
        case BRepCheck_Multiple3DCurve:
            return "BRepCheck_Multiple3DCurve";
            break;
        case BRepCheck_Invalid3DCurve:
            return "BRepCheck_Invalid3DCurve";
            break;
        case BRepCheck_NoCurveOnSurface:
            return "BRepCheck_NoCurveOnSurface";
            break;
        case BRepCheck_InvalidCurveOnSurface:
            return "BRepCheck_InvalidCurveOnSurface";
            break;
        case BRepCheck_InvalidCurveOnClosedSurface:
            return "BRepCheck_InvalidCurveOnClosedSurface";
            break;
        case BRepCheck_InvalidSameRangeFlag:
            return "BRepCheck_InvalidSameRangeFlag";
            break;
        case BRepCheck_InvalidSameParameterFlag:
            return "BRepCheck_InvalidSameParameterFlag";
            break;
        case BRepCheck_InvalidDegeneratedFlag:
            return "BRepCheck_InvalidDegeneratedFlag";
            break;
        case BRepCheck_FreeEdge:
            return "BRepCheck_FreeEdge";
            break;
        case BRepCheck_InvalidMultiConnexity:
            return "BRepCheck_InvalidMultiConnexity";
            break;
        case BRepCheck_InvalidRange:
            return "BRepCheck_InvalidRange";
            break;
        case BRepCheck_EmptyWire:
            return "BRepCheck_EmptyWire";
            break;
        case BRepCheck_RedundantEdge:
            return "BRepCheck_RedundantEdge";
            break;
        case BRepCheck_SelfIntersectingWire:
            return "BRepCheck_SelfIntersectingWire";
            break;
        case BRepCheck_NoSurface:
            return "BRepCheck_NoSurface";
            break;
        case BRepCheck_InvalidWire:
            return "BRepCheck_InvalidWire";
            break;
        case BRepCheck_RedundantWire:
            return "BRepCheck_RedundantWire";
            break;
        case BRepCheck_IntersectingWires:
            return "BRepCheck_IntersectingWires";
            break;
        case BRepCheck_InvalidImbricationOfWires:
            return "BRepCheck_InvalidImbricationOfWires";
            break;
        case BRepCheck_EmptyShell:
            return "BRepCheck_EmptyShell";
            break;
        case BRepCheck_RedundantFace:
            return "BRepCheck_RedundantFace";
            break;
        case BRepCheck_InvalidImbricationOfShells:
            return "BRepCheck_InvalidImbricationOfShells";
            break;
        case BRepCheck_UnorientableShape:
            return "BRepCheck_UnorientableShape";
            break;
        case BRepCheck_NotClosed:
            return "BRepCheck_NotClosed";
            break;
        case BRepCheck_NotConnected:
            return "BRepCheck_NotConnected";
            break;
        case BRepCheck_SubshapeNotInShape:
            return "BRepCheck_SubshapeNotInShape";
            break;
        case BRepCheck_BadOrientation:
            return "BRepCheck_BadOrientation";
            break;
        case BRepCheck_BadOrientationOfSubshape:
            return "BRepCheck_BadOrientationOfSubshape";
            break;
        case BRepCheck_InvalidPolygonOnTriangulation:
            return "BRepCheck_InvalidPolygonOnTriangulation";
            break;
        case BRepCheck_InvalidToleranceValue:
            return "BRepCheck_InvalidToleranceValue";
            break;
        case BRepCheck_EnclosedRegion:
            return "BRepCheck_EnclosedRegion";
            break;
        case BRepCheck_CheckFail:
            return "BRepCheck_CheckFail";
            break;
        default:
            return "Unknown";
    }
}

TopoDS_ListOfShape ShapeChecker::get_vertices() { return loop_topo(TopAbs_VERTEX); }

TopoDS_ListOfShape ShapeChecker::get_edges() { return loop_topo(TopAbs_EDGE); }

TopoDS_ListOfShape ShapeChecker::get_wires() { return loop_topo(TopAbs_WIRE); }

TopoDS_ListOfShape ShapeChecker::get_faces() { return loop_topo(TopAbs_FACE); }

TopoDS_ListOfShape ShapeChecker::get_shells() { return loop_topo(TopAbs_SHELL); }

TopoDS_ListOfShape ShapeChecker::get_solids() { return loop_topo(TopAbs_SOLID); }

TopoDS_ListOfShape ShapeChecker::get_compsolids() { return loop_topo(TopAbs_COMPSOLID); }

TopoDS_ListOfShape ShapeChecker::get_compounds() { return loop_topo(TopAbs_COMPOUND); }

void ShapeChecker::check_vertices() {
    for (auto it = vertices.begin(); it != vertices.end(); ++it) {
        std::cerr << "\n=== Vertex " << (*it).HashCode(INT_MAX) << " ===\n";
        check_vertex(TopoDS::Vertex(*it));
    }
}

void ShapeChecker::check_edges() {
    for (auto it = edges.begin(); it != edges.end(); ++it) {
        std::cerr << "\n=== Edge " << (*it).HashCode(INT_MAX) << " ===\n";
        check_edge(TopoDS::Edge(*it));
    }
}

void ShapeChecker::check_wires() {
    for (auto it = wires.begin(); it != wires.end(); ++it) {
        std::cerr << "\n=== Wire " << (*it).HashCode(INT_MAX) << " ===\n";
        check_wire(TopoDS::Wire(*it));
    }
}

void ShapeChecker::check_faces() {
    for (auto it = faces.begin(); it != faces.end(); ++it) {
        std::cerr << "\n=== Face " << (*it).HashCode(INT_MAX) << " ===\n";
        check_face(TopoDS::Face(*it));
    }
}

void ShapeChecker::check_shells() {
    for (auto it = shells.begin(); it != shells.end(); ++it) {
        std::cerr << "\n=== Shell " << (*it).HashCode(INT_MAX) << " ===\n";
        check_shell(TopoDS::Shell(*it));
    }
}

void ShapeChecker::check_vertex(const TopoDS_Vertex &V) {

    get_results_from_analyzer(V);

}

void ShapeChecker::check_edge(const TopoDS_Edge &E) {

    get_results_from_analyzer(E);

    BRepCheck_Edge c = BRepCheck_Edge(E);
    BRepCheck_Status CheckPolygonOnTriangulation = c.CheckPolygonOnTriangulation(E);
    add_to_errors(E, CheckPolygonOnTriangulation);

    ShapeAnalysis_Edge a = ShapeAnalysis_Edge();
    Standard_Boolean IsClosed3d = a.IsClosed3d(E);

    std::cerr << "\tCheckPolygonOnTriangulation: " << int_to_status(CheckPolygonOnTriangulation) << "\n";
    std::cerr << "\tIsClosed3d: " << std::boolalpha << IsClosed3d << "\n";
}

void ShapeChecker::check_wire(const TopoDS_Wire &W) {

    get_results_from_analyzer(W);

    BRepCheck_Wire c = BRepCheck_Wire(W);

    Standard_Real ContourArea = ShapeAnalysis::ContourArea(W);

    std::cerr << "\tContourArea: " << ContourArea << "\n";
    std::cerr << "\tClosed: " << std::boolalpha << W.Closed() << "\n";
    std::cerr << "\tInfinite: " << std::boolalpha << W.Infinite() << "\n";

    TopoDS_ListOfShape faces = faces_from_wire(W);

    for (auto it = faces.begin(); it != faces.end(); ++it) {

        BRepCheck_Status Closed = c.Closed();
        BRepCheck_Status Closed2d = c.Closed2d(TopoDS::Face(*it), false);
        BRepCheck_Status Orientation = c.Orientation(TopoDS::Face(*it), false);
        add_to_errors(W, Closed);
        add_to_errors(W, Closed2d);
        add_to_errors(W, Orientation);

        std::cerr << "\tClosed: " << int_to_status(Closed) << "\n";
        std::cerr << "\tClosed2d: " << int_to_status(Closed2d) << "\n";
        std::cerr << "\tOrientation: " << int_to_status(Orientation) << "\n";

        ShapeAnalysis_Wire a = ShapeAnalysis_Wire(W, TopoDS::Face(*it), tol);
        Standard_Boolean Perform = a.Perform(); // False, does most of the checks below
        Standard_Boolean CheckOrder = a.CheckOrder(); // False <- expectation
        Standard_Boolean CheckConnected = a.CheckConnected(tol); // False
        Standard_Boolean CheckSmall = a.CheckSmall(tol); // False
        Standard_Boolean CheckEdgeCurves = a.CheckEdgeCurves(); // False
        Standard_Boolean CheckDegenerated = a.CheckDegenerated(); // False
        Standard_Boolean CheckClosed = a.CheckClosed(tol); // False
        Standard_Boolean CheckSelfIntersection = a.CheckSelfIntersection(); // False, only adjacent edges are check
        Standard_Boolean CheckLacking = a.CheckLacking(); // False
        Standard_Boolean CheckIntersectingEdges = a.CheckIntersectingEdges(0, 1);
        Standard_Real MinDistance3d = a.MinDistance3d();
        Standard_Real MaxDistance3d = a.MaxDistance3d();

        ShapeExtend_WireData se(W);
        se.Add(a.WireData());
        Standard_Integer NbNonManifoldEdges = se.NbNonManifoldEdges();

        std::cerr << "\tPerform: " << std::boolalpha << Perform << "\n";
        std::cerr << "\tCheckOrder: " << std::boolalpha << CheckOrder << "\n";
        std::cerr << "\tCheckConnected: " << std::boolalpha << CheckConnected << "\n";
        std::cerr << "\tCheckSmall: " << std::boolalpha << CheckSmall << "\n";
        std::cerr << "\tCheckEdgeCurves: " << std::boolalpha << CheckEdgeCurves << "\n";
        std::cerr << "\tCheckDegenerated: " << std::boolalpha << CheckDegenerated << "\n";
        std::cerr << "\tCheckClosed: " << std::boolalpha << CheckClosed << "\n";
        std::cerr << "\tCheckSelfIntersection: " << std::boolalpha << CheckSelfIntersection << "\n";
        std::cerr << "\tCheckLacking: " << std::boolalpha << CheckLacking << "\n";
        std::cerr << "\tCheckIntersectingEdges: " << std::boolalpha << CheckIntersectingEdges << "\n";
        std::cerr << "\tMinDistance3d: " << MinDistance3d << "\n";
        std::cerr << "\tMaxDistance3d: " << MaxDistance3d << "\n";
        std::cerr << "\tNbNonManifoldEdges: " << NbNonManifoldEdges << "\n";
    }

}

void ShapeChecker::check_face(const TopoDS_Face &F) {

    get_results_from_analyzer(F);

    BRepCheck_Face c = BRepCheck_Face(F);

    BRepCheck_Status IntersectWires = c.IntersectWires();
    BRepCheck_Status ClassifyWires = c.ClassifyWires();
    BRepCheck_Status OrientationOfWires = c.OrientationOfWires();
    Standard_Boolean IsUnorientable = c.IsUnorientable();

    add_to_errors(F, IntersectWires);
    add_to_errors(F, ClassifyWires);
    add_to_errors(F, OrientationOfWires);

    std::cerr << "\tIntersectWires: " << int_to_status(IntersectWires) << "\n";
    std::cerr << "\tClassifyWires: " << int_to_status(ClassifyWires) << "\n";
    std::cerr << "\tOrientationOfWires: " << int_to_status(OrientationOfWires) << "\n";
    std::cerr << "\tIsUnorientable: " << std::boolalpha << IsUnorientable << "\n";
    std::cerr << "\tInfinite: " << std::boolalpha << F.Infinite() << "\n";
    std::cerr << "\tOrientation: " << F.Orientation() << "\n";

}

void ShapeChecker::check_shell(const TopoDS_Shell &S) {

    get_results_from_analyzer(S);

    BRepCheck_Shell c = BRepCheck_Shell(S);

    BRepCheck_Status Closed = c.Closed();
    BRepCheck_Status Orientation = c.Orientation();
    Standard_Boolean IsUnorientable = c.IsUnorientable();

    add_to_errors(S, Closed);
    add_to_errors(S, Orientation);

    std::cerr << "\tClosed: " << int_to_status(Closed) << "\n";
    std::cerr << "\tOrientation: " << int_to_status(Orientation) << "\n";
    std::cerr << "\tIsUnorientable: " << std::boolalpha << IsUnorientable << "\n";

    ShapeAnalysis_Shell a;
    a.LoadShells(S);
    Standard_Boolean CheckOrientedShells = a.CheckOrientedShells(S, true, true);
    Standard_Boolean HasBadEdges = a.HasBadEdges();
    Standard_Boolean HasFreeEdges = a.HasFreeEdges();
    std::cerr << "\tCheckOrientedShells: " << std::boolalpha << CheckOrientedShells << "\n";
    std::cerr << "\tHasBadEdges: " << std::boolalpha << HasBadEdges << "\n";
    std::cerr << "\tHasFreeEdges: " << std::boolalpha << HasFreeEdges << "\n";
    std::cerr << "\tInfinite: " << std::boolalpha << S.Infinite() << "\n";
    int n = 0;
    TopExp_Explorer Ex;
    for (Ex.Init(S, TopAbs_FACE); Ex.More(); Ex.Next())
        n++;
    std::cerr << "\tNumber of faces: " << n << "\n";
}

TopoDS_ListOfShape ShapeChecker::faces_from_wire(const TopoDS_Wire &W) { return map_shapes_and_ancestors(TopAbs_WIRE, TopAbs_FACE, W); }

TopoDS_ListOfShape ShapeChecker::map_shapes_and_ancestors(TopAbs_ShapeEnum topoTypeA, TopAbs_ShapeEnum topoTypeB, TopoDS_Shape topologicalEntity) {

    TopTools_IndexedDataMapOfShapeListOfShape _map;
    TopExp::MapShapesAndAncestors(shape, topoTypeA, topoTypeB, _map);
    return _map.FindFromKey(topologicalEntity);
}

void ShapeChecker::add_to_errors(const TopoDS_Shape &shape, BRepCheck_Status status) {

    if (status == BRepCheck_NoError)
        return;

    auto it = errors.find(&shape);
    if (it != errors.end())
        it->second.push_back(status);
    else {
        std::list<BRepCheck_Status> l{status};
        errors[&shape] = l;
    }
}

void ShapeChecker::print_errors() {

    std::cerr << "\n";
    for (auto &it: errors) {
        std::cerr << shapeTypeString(it.first->ShapeType()) << " " << it.first->HashCode(INT_MAX) << ":";
        for (auto it2 = it.second.cbegin(); it2 != it.second.cend(); ++it2)
            std::cerr << " " << int_to_status(*it2);
        std::cerr << "\n";
    }

}

bool ShapeChecker::face_is_planar(const TopoDS_Face &face) const {
    auto surface = BRep_Tool::Surface(face);
    GeomLib_IsPlanarSurface is_planar_surface(surface, tol);
    return is_planar_surface.IsPlanar();

}

TopoDS_ListOfShape ShapeChecker::free_edges() {

    std::map<unsigned int, TopoDS_Edge> E;
    std::map<unsigned int, std::set<TopAbs_Orientation>> M;

    for (auto &edge: Topo(shape).edges()) {
        unsigned int h = Kernel::hash(edge);
        M[h].insert(edge.Orientation());
        E[h] = TopoDS::Edge(edge);
    }

    TopoDS_ListOfShape free;

    for (auto &m: M) {
        if (m.second.find(TopAbs_FORWARD) != m.second.end() && m.second.find(TopAbs_REVERSED) != m.second.end()) continue;
        free.Append(E[m.first]);
    }

    return free;
}
