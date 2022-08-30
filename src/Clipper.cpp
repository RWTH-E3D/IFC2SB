// Copyright 2022 Eric Fichter
#include "Clipper.h"

Clipper::Clipper(TopoDS_Face tF1, gp_Dir tn1, bool _dimension) :
        F1(std::move(tF1)),
        n1(tn1),
        d1(gp_Vec(tn1)),
        dimension(_dimension),
        f(1e6),
        P_plane(BRep_Tool::Pnt(TopoDS::Vertex(Topo(F1).vertices().First()))) {

    success = false;
    same = false; // true if projected and clipped face is same as second face
    Result.Clear();
    Result2DWires.clear();
    F1_wire2D.clear();
    F1_holes2D.clear();
    P1_2D.clear();
    P1_holes_2D.clear();

    // Geometric values
    calculate_algebraic_values();

    // Tolerances
    f_inv = 1.0 / f;
    tol_area = 1.0e-6 * f * f; // 1e-6m² is 1mm*1mm
    tol_length = 1.0e-3 * f;

    TopoDS_Wire W1 = BRepTools::OuterWire(F1);
    if (W1.IsNull())
        std::cerr << "[Error] Outer wire is faulty!" << std::endl;

    // 2D wire
    P1_2D = wire_to_projected_points_1(W1);

    // 2D holes
    if (Topo(F1).wires().Size() > 1)
        P1_holes_2D = face_to_projected_holes_1(F1, W1);

}

void Clipper::clip(const TopoDS_Face &F2, gp_Dir n2) {

    success = false;
    same = false; // true if projected and clipped face is same as second face
    Result.Clear();
    Result2DWires.clear();
    F1_wire2D.clear();
    F1_holes2D.clear();

    if (fabs(n2.Dot(n1)) < 1.0e-12) return;

    TopoDS_Wire W2 = BRepTools::OuterWire(F2);

    if (W2.IsNull()) {
        std::cerr << "[Error] Outer wire is faulty!" << std::endl;
        return;
    }

    auto P2_2D = wire_to_projected_points_2(W2);

    // aabb intersection test
    if (!do_boxes_intersect(aabb_polygon(P1_2D), aabb_polygon(P2_2D))) return;

    // Holes
    std::vector<PointChain2D> P3_holes_2D;

    if (Topo(F2).wires().Size() > 1)
        P3_holes_2D = face_to_projected_holes_2(F2, W2);

    // Clipping of 2D polygons
    std::list<PointChain2D> polygons_2D, holes_2D;
    clipping(P1_2D, P2_2D, P1_holes_2D, P3_holes_2D, polygons_2D, holes_2D);

    if (polygons_2D.empty()) return;

    if (dimension) {

        // Reproject 2D polygons to have 3D points
        auto polygons_3D = reproject_2D_to_3D(polygons_2D);
        TopoDS_ListOfShape Holes;
        if (!holes_2D.empty()) {
            auto holes_3D = reproject_2D_to_3D(holes_2D);
            for (const auto &hole_3D: holes_3D)
                Holes.Append(polygon_to_wire(hole_3D));
        }

        // Create TopoDS_Faces
        for (const auto &polygon_3D: polygons_3D) {
            TopoDS_Face F = polygon_to_face(polygon_3D, Holes);
            if (!F.IsNull())
                Result.Append(F);
        }

        if (!Result.IsEmpty()) success = true;

    } else {

        chain2D L;

        for (const auto &poly: polygons_2D) {
            L.clear();
            for (const auto &P: poly)
                L.emplace_back(P.X, P.Y);
            Result2DWires.push_back(L);
        }

        for (const auto &poly: holes_2D) {
            L.clear();
            for (const auto &P: poly)
                L.emplace_back(P.X, P.Y);
            Result2DWires.push_back(L);
        }

        for (const auto &P: P1_2D)
            F1_wire2D.emplace_back(P.X, P.Y);

        for (const auto &poly: P1_holes_2D) {
            L.clear();
            for (const auto &P: poly)
                L.emplace_back(P.X, P.Y);
            F1_holes2D.push_back(L);
        }

        success = true;
    }
}

Clipper::PointChain2D Clipper::wire_to_projected_points_1(const TopoDS_Wire &W) const {
    TopoDS_ListOfShape V = Topo(W).ordered_vertices_of_wire();
    auto P = vertex_to_point_list(V);
    return project_points_to_2D(P);
}

Clipper::PointChain2D Clipper::wire_to_projected_points_2(const TopoDS_Wire &W) const {
    TopoDS_ListOfShape V = Topo(W).ordered_vertices_of_wire(); // Vertex list
    auto P2 = vertex_to_point_list(V); // Point list
    auto P3 = projects_points_onto_face_plane(P2);     // Project points of F2 onto F1 in the direction of F2 normal n2
    return project_points_to_2D(P3); // Project 3D points to 2D
}

std::vector<Clipper::PointChain2D> Clipper::face_to_projected_holes_1(const TopoDS_Face &F, const TopoDS_Wire &W) const {
    auto V_holes = ordered_vertices_of_hole_wires(F, W);
    auto P_holes = vertex_to_point_list_holes(V_holes);
    return project_points_to_2D_holes(P_holes);
}

std::vector<Clipper::PointChain2D> Clipper::face_to_projected_holes_2(const TopoDS_Face &F, const TopoDS_Wire &W) const {
    auto V_holes = ordered_vertices_of_hole_wires(F, W);
    auto P2_holes = vertex_to_point_list_holes(V_holes);
    auto P3_holes = projects_points_onto_face_plane_holes(P2_holes);
    return project_points_to_2D_holes(P3_holes);
}

std::vector<gp_Pnt> Clipper::vertex_to_point_list(const TopoDS_ListOfShape &V) {
    std::vector<gp_Pnt> P;
    for (const auto &t: V)
        P.push_back(BRep_Tool::Pnt(TopoDS::Vertex(t)));
    return P;
}

void Clipper::calculate_algebraic_values() { calculate_algebraic_values(P_plane, n1, a, proj_axis); }

void Clipper::calculate_algebraic_values(gp_Pnt P, gp_Dir n, double &temp_a, unsigned int &temp_proj_axis) {
    temp_a = P.X() * n.X() + P.Y() * n.Y() + P.Z() * n.Z();
    std::vector<double> temp = {fabs(n.X()), fabs(n.Y()), fabs(n.Z())}; // determine projection axis (0,1 or 2), such that v[proj_axis] != 0
    temp_proj_axis = std::distance(temp.begin(), std::max_element(temp.begin(), temp.end()));
}

gp_Pnt Clipper::project_point_onto_plane(gp_Vec P) const {
    // P ... Point to be projected

    // q_proj = P - dot(P - P_plane, n1) * n1
    const double t1 = P.Subtracted(P_plane.XYZ()).Dot(n1);
    gp_Vec t2 = d1.Multiplied(t1);
    return {P.Subtracted(t2).XYZ()};

    /*
     // d2 is the vector along the point is projected
     // if d2 = n1, implementation on top is ok
    const double t1 = P.Subtracted(P_plane.XYZ()).Dot(d1);
    const double t2 = d2.Dot(d1);
    gp_Vec t3 = d2.Multiplied(t1 / t2); // this will be division by zero, if faces have 90° angle and the projected face would be a line
    return {P.Subtracted(t3).XYZ()};
    */
}

std::vector<gp_Pnt> Clipper::projects_points_onto_face_plane(const std::vector<gp_Pnt> &R) const {
    std::vector<gp_Pnt> Q;
    for (const auto &P: R) Q.push_back(project_point_onto_plane(gp_Vec(P.XYZ())));
    return Q;
}

Clipper::Point2D Clipper::project(gp_Pnt P) const {

    // Project onto either the xy, yz, or xz plane, depending on proj_axis
    Point2D L{};

    if (proj_axis == 0) {
        L.X = P.Y();
        L.Y = P.Z();
    } else if (proj_axis == 1) {
        L.X = P.X();
        L.Y = P.Z();
    } else {
        L.X = P.X();
        L.Y = P.Y();
    }

    return L;
}

Clipper::PointChain2D Clipper::project_points_to_2D(const std::vector<gp_Pnt> &R) const {

    PointChain2D R_2D;
    for (const auto &P: R) R_2D.push_back(project(P));
    return R_2D;
}

std::list<TopoDS_ListOfShape> Clipper::ordered_vertices_of_hole_wires(const TopoDS_Face &F, const TopoDS_Wire &outerWire) {
    std::list<TopoDS_ListOfShape> L;
    for (const auto &W: Topo(F).wires())
        if (!W.IsSame(outerWire))
            L.push_back(Topo(W).ordered_vertices_of_wire());
    return L;
}

std::list<std::vector<gp_Pnt>> Clipper::vertex_to_point_list_holes(const std::list<TopoDS_ListOfShape> &V) {
    std::list<std::vector<gp_Pnt>> P;
    for (const auto &h: V) {
        P.emplace_back();
        for (const auto &t: h) P.back().push_back(BRep_Tool::Pnt(TopoDS::Vertex(t)));
    }
    return P;
}

std::list<std::vector<gp_Pnt>> Clipper::projects_points_onto_face_plane_holes(const std::list<std::vector<gp_Pnt>> &R) const {

    std::list<std::vector<gp_Pnt>> Q;
    for (const auto &h: R) {
        Q.emplace_back();
        for (const auto &P: h) Q.back().push_back(project_point_onto_plane(gp_Vec(P.XYZ())));
    }
    return Q;
}

std::vector<Clipper::PointChain2D> Clipper::project_points_to_2D_holes(const std::list<std::vector<gp_Pnt>> &R) const {

    std::vector<PointChain2D> Q;
    Q.reserve(R.size());

    for (const auto &h: R) {
        Q.emplace_back();
        for (const auto &P: h) Q.back().push_back(project(P));
    }
    return Q;
}

void Clipper::clipping(const PointChain2D &PolyPoints1, const PointChain2D &PolyPoints2, const std::vector<PointChain2D> &HolePoints1, const std::vector<PointChain2D> &HolePoints2,
                       std::list<PointChain2D> &common_2D_polygons, std::list<PointChain2D> &common_2D_holes) {

    ClipperLib::Clipper c;
    ClipperLib::Paths poly1(1 + HolePoints1.size()), poly2(1 + HolePoints2.size()), solution;

    for (const auto &P: PolyPoints1)
        poly1[0] << ClipperLib::IntPoint((unsigned long long) (f * P.X), (unsigned long long) (f * P.Y));

    for (const auto &P: PolyPoints2)
        poly2[0] << ClipperLib::IntPoint((unsigned long long) (f * P.X), (unsigned long long) (f * P.Y));

    if (!HolePoints1.empty())
        for (unsigned int i = 0; i < HolePoints1.size(); i++)
            for (const auto &P: HolePoints1[i])
                poly1[i + 1] << ClipperLib::IntPoint((unsigned long long) (f * P.X), (unsigned long long) (f * P.Y));

    if (!HolePoints2.empty())
        for (unsigned int i = 0; i < HolePoints2.size(); i++) {
            for (const auto &P: HolePoints2[i])
                poly2[i + 1] << ClipperLib::IntPoint((unsigned long long) (f * P.X), (unsigned long long) (f * P.Y));
        }

    c.AddPaths(poly1, ClipperLib::ptSubject, true);
    c.AddPaths(poly2, ClipperLib::ptClip, true);
    c.Execute(ClipperLib::ctIntersection, solution, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
    ClipperLib::CleanPolygons(solution, tol_length);
    ClipperLib::SimplifyPolygons(solution);

    // count non-hole polygons, to return if no result or if clipped polygon has same area than original polygon
    long count = std::count_if(solution.begin(), solution.end(), [](ClipperLib::Path &path) { return ClipperLib::Orientation(path); });

    if (count == 0) return;

    else if (count == 1) {
        double A1 = fabs(ClipperLib::Area(poly1[0]));
        double A2 = 0;

        // For closed paths (polygons) in the solution returned by Clipper's Execute method, their orientations will always be true for outer polygons and false for hole polygons (unless the ReverseSolution property has been enabled).
        for (const auto &s: solution)
            if (ClipperLib::Orientation(s))
                A2 += fabs(ClipperLib::Area(s));
            else
                A2 -= fabs(ClipperLib::Area(s));

        if (fabs(A1 - A2) < tol_area) {
            same = true;
            return;
        }
    }

    // outer wire will be returned by clipper with "true" orientation, holes with "false"
    for (const auto &s: solution) {

        if (fabs(ClipperLib::Area(s)) < tol_area) continue; // skip tiny polygons

        if (ClipperLib::Orientation(s)) {
            PointChain2D common_2D_polygon;
            for (const auto &clipperPoint: s)
                common_2D_polygon.push_back({clipperPoint.X * f_inv, clipperPoint.Y * f_inv});
            common_2D_polygons.push_back(common_2D_polygon);
        } else {
            PointChain2D common_2D_hole;
            for (const auto &clipperPoint: s)
                common_2D_hole.push_back({clipperPoint.X * f_inv, clipperPoint.Y * f_inv});
            common_2D_holes.push_back(common_2D_hole);
        }
    }
}

gp_Pnt Clipper::project_inv(Point2D P_2D) const {

    std::vector<double> P_3D(3);
    double t;

    if (proj_axis == 0) {
        t = n1.X();
        P_3D[0] = 0.0;
        P_3D[1] = P_2D.X;
        P_3D[2] = P_2D.Y;
    } else if (proj_axis == 1) {
        t = n1.Y();
        P_3D[0] = P_2D.X;
        P_3D[1] = 0.0;
        P_3D[2] = P_2D.Y;
    } else {
        t = n1.Z();
        P_3D[0] = P_2D.X;
        P_3D[1] = P_2D.Y;
        P_3D[2] = 0.0;
    }

    double z = a;
    z -= P_3D[0] * n1.X();
    z -= P_3D[1] * n1.Y();
    z -= P_3D[2] * n1.Z();
    z /= t;
    P_3D[proj_axis] = z;
    return {P_3D[0], P_3D[1], P_3D[2]};
}

std::vector<std::vector<gp_Pnt>> Clipper::reproject_2D_to_3D(const std::list<PointChain2D> &R) const {

    std::vector<std::vector<gp_Pnt>> Q;

    for (const auto &poly: R) {
        std::vector<gp_Pnt> L;
        for (const auto &P_2D: poly)
            L.emplace_back(project_inv(P_2D));
        Q.push_back(L);
    }

    return Q;
}

TopoDS_Wire Clipper::polygon_to_wire(const std::vector<gp_Pnt> &polygon_3D) {

    TopoDS_ListOfShape E;
    for (int i = 0; i < polygon_3D.size() - 1; i++)
        E.Append(BRepBuilderAPI_MakeEdge(polygon_3D[i], polygon_3D[i + 1]).Edge());
    E.Append(BRepBuilderAPI_MakeEdge(polygon_3D[polygon_3D.size() - 1], polygon_3D[0]).Edge());

    BRepBuilderAPI_MakeWire W;

    for (const auto &e: E)
        W.Add(TopoDS::Edge(e));

    if (W.IsDone())
        return W.Wire();
    else {
        std::cerr << "[Error] No wire created!" << std::endl;
        return {};
    }
}

TopoDS_Face Clipper::polygon_to_face(const std::vector<gp_Pnt> &polygon_3D, const TopoDS_ListOfShape &Holes) {

    TopoDS_ListOfShape E;
    for (int i = 0; i < polygon_3D.size() - 1; i++)
        E.Append(BRepBuilderAPI_MakeEdge(polygon_3D[i], polygon_3D[i + 1]).Edge());
    E.Append(BRepBuilderAPI_MakeEdge(polygon_3D[polygon_3D.size() - 1], polygon_3D[0]).Edge());

    BRepBuilderAPI_MakeWire W;

    for (const auto &e: E)
        W.Add(TopoDS::Edge(e));

    if (Holes.IsEmpty()) {
        if (W.IsDone())
            return BRepBuilderAPI_MakeFace(W.Wire()).Face();
        else {
            std::cerr << "[Error] No face created!" << std::endl;
            return {};
        }
    } else {
        if (W.IsDone()) {
            BRepBuilderAPI_MakeFace mFace(W.Wire());
            for (const auto &hole: Holes) mFace.Add(TopoDS::Wire(hole));
            return mFace.Face();
        } else {
            std::cerr << "[Error] No face created!" << std::endl;
            return {};
        }
    }
}

std::vector<double> Clipper::aabb_polygon(const PointChain2D &poly) {

    std::vector<double> aabb(4);

    std::vector<double> X, Y;
    for (const auto &P: poly) {
        X.push_back(P.X);
        Y.push_back(P.Y);
    }

    aabb[0] = *min_element(X.begin(), X.end()); // xmin
    aabb[1] = *max_element(X.begin(), X.end()); // xmax
    aabb[2] = *min_element(Y.begin(), Y.end()); // ymin
    aabb[3] = *max_element(Y.begin(), Y.end()); // ymax

    return aabb;
}

bool Clipper::do_boxes_intersect(const std::vector<double> &A, const std::vector<double> &B) { return !(A[0] > B[1] || B[0] > A[1] || A[3] < B[2] || B[3] < A[2]); }