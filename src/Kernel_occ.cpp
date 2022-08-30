// Copyright 2022 Eric Fichter
#include "Kernel.h"

// Functions dealing with OCC stuff or small geometrical calculations.

unsigned int Kernel::hash(const TopoDS_Shape &s) { return s.HashCode(INT_MAX); } // It is computed from the TShape and the Location of the shape. The Orientation is not used. So coplanar faces can have same ID, even if they are oppositely oriented.

double Kernel::volume(const TopoDS_Shape &shape) {
    GProp_GProps gprop;
    BRepGProp::VolumeProperties(shape, gprop);
    return gprop.Mass();
}

double Kernel::length(const TopoDS_Edge &e) {
    GProp_GProps gprop;
    BRepGProp::LinearProperties(e, gprop);
    return gprop.Mass();
}

TopoDS_Shape Kernel::shape_copy(const TopoDS_Shape &S) {
    BRepBuilderAPI_Copy A;
    A.Perform(S);
    return A.Shape();
}

void Kernel::move(TopoDS_Shape &S, gp_Vec v) {
    gp_Trsf trsf;
    trsf.SetTranslation(v);
    S.Move(TopLoc_Location(trsf));
}

TopoDS_Shape Kernel::moved(TopoDS_Shape &S, gp_Vec v) {
    gp_Trsf trsf;
    trsf.SetTranslation(v);
    return S.Moved(TopLoc_Location(trsf));
}

void Kernel::rotate(TopoDS_Shape &S, gp_Pnt P, gp_Dir n, double angle) {
    gp_Trsf trsf;
    trsf.SetRotation(gp_Ax1(P, n), angle);
    S.Move(TopLoc_Location(trsf));
}

bool Kernel::face_is_polygon(const TopoDS_Face &face) {
    return (face_is_planar(face) && Topo(face).all_edges_non_curved());
}

bool Kernel::face_is_planar(const TopoDS_Face &face) {
    auto surface = BRep_Tool::Surface(face);
    GeomLib_IsPlanarSurface is_planar_surface(surface, 1.0e-6);
    return is_planar_surface.IsPlanar();
}

gp_XY Kernel::surface_uv_coordinates_at_point(const Handle(Geom_Surface) &surface, gp_Pnt P) {
    ShapeAnalysis_Surface a(surface);
    gp_Pnt2d uv = a.ValueOfUV(P, 1.0e-5);
    return uv.Coord();
}

gp_Dir Kernel::face_normal_at_point(const TopoDS_Face &face, gp_Pnt P, bool &success) {

    // https://github.com/tpaviot/pythonocc-utils/blob/master/OCCUtils/face.py
    success = false;
    Handle(Geom_Surface) surface = BRep_Tool::Surface(face);
    auto uv = surface_uv_coordinates_at_point(surface, P);
    GeomLProp_SLProps curv(surface, uv.X(), uv.Y(), 2, 1e-7);

    curv.SetParameters(uv.X(), uv.Y());

    if (curv.IsNormalDefined()) {
        gp_Dir normal = curv.Normal();
        if (face.Orientation() == TopAbs_FORWARD) {
            success = true;
            return normal;
        } else if (face.Orientation() == TopAbs_REVERSED) {
            success = true;
            return normal.Reversed();
        } else {
            std::cerr << "[Error] Wrong Orientation in normal calculation: " << normal.X() << "\t" << normal.Y() << "\t" << normal.Z() << "\t" << face.Orientation() << "\t" << hash(face) << std::endl;
            return normal;
        }
    } else {
        std::cerr << "[Error] Undefined normal! " << face.Orientation() << std::endl;
        return {1, 0, 0};
    }
}

gp_Dir Kernel::face_normal(const TopoDS_Face &face) {

    Standard_Real UMin, UMax, VMin, VMax;
    BRepTools::UVBounds(face, UMin, UMax, VMin, VMax);
    BRepAdaptor_Surface surface = BRepAdaptor_Surface(face); // create surface to access the geometry of face
    gp_Dir normal = BRepLProp_SLProps(surface, UMin, VMin, 1, 1.e-3).Normal();

    if (face.Orientation() == TopAbs_FORWARD)
        return normal;
    else if (face.Orientation() == TopAbs_REVERSED)
        return normal.Reversed();
    else {
        std::cerr << "[Error] Wrong Orientation in normal calculation: " << normal.X() << "\t" << normal.Y() << "\t" << normal.Z() << "\t" << face.Orientation() << "\t" << hash(face) << std::endl;
        return normal;
    }
}

bool Kernel::offset_faces_outer_wire_cut(const TopoDS_Face &face, TopoDS_Face &offset_face, double offset) {

    TopoDS_Wire outerWire = BRepTools::OuterWire(face);

    if (outerWire.IsNull()) {
        std::cerr << "Error: No outer wire!\n";
        return false;
    }

    TopoDS_Wire offsetWire;

    if (!offset_wire_cut(outerWire, offsetWire, offset, GeomAbs_Intersection)) { // Intersection offset didn't succeed
        if (!offset_wire_cut(outerWire, offsetWire, offset, GeomAbs_Arc)) {// Intersection and Arc offset didn't succeed
            std::cerr << "Error: Both offset functions failed!. Return original face.\n";
            return false;
        }
    }

    // Intersection or Arc offset succeeded
    BRepBuilderAPI_MakeFace maker(offsetWire);
    //maker.Add(TopoDS::Wire(outerWire.Complemented())); // Complementing or reversing wire seems to fail sometimes, yielding to a constructed face with an island outerWire and empty offsetWire
    maker.Add(BRepTools::OuterWire(TopoDS::Face(BRepBuilderAPI_MakeFace(outerWire).Face().Complemented()))); // that seems not to happen when building a face first, that is reversed

    if (!maker.IsDone()) {
        std::cerr << "construct_face: Is not Done!\n";
        return false;
    }

    offset_face = maker.Face();

    if (offset_face.IsNull()) {
        std::cerr << "Error: Offsetted face is Null!\n";
        return false;
    }

    /*
     * // ShapeFixing will fix Orientation problems with islands, but are very time consuming
     * ShapeFix_Shape sfs;
     * sfs.Init(offset_face);
     * sfs.Perform();
     * offset_face = TopoDS::Face(sfs.Shape());
     * */

    return true;
}

bool Kernel::offset_wire_cut(const TopoDS_Wire &wire, TopoDS_Wire &offsetWire, double offset, GeomAbs_JoinType joinType) {

    if (wire.IsNull()) {
        std::cerr << "\t[Warning] Null Wire!" << std::endl;
        return false;
    }

    BRepOffsetAPI_MakeOffset o(wire, joinType, Standard_False);// if true, it sometimes returns unconnected edges
    o.Perform(offset, 0);

    // Check 1: Is offset done at all
    // o.Check();
    if (!o.IsDone()) { // This seems to happen for very long thin wires with one really small angle
        ShapeChecker(wire).dump_topology();
        std::cerr << "\t[Warning] Make Wire Offset " << joinType << " is not Done!" << std::endl;
        return false;
    }

    TopoDS_Shape offset_shape = o.Shape();

    // Check 2: Is there exactly on wire in the offset shape
    Topo top(offset_shape);
    unsigned int number_of_offset_wires = top.wires().Size();
    if (number_of_offset_wires != 1) {
        ShapeChecker(wire).dump_topology();
        ShapeChecker(offset_shape).dump_topology();
        std::cerr << "\t[Warning] No or more than 1 wire returned by MakeOffset " << joinType << std::endl;
        return false;
    }

    // Check 3: If Arc, create a Face from the single wire of the offset shape and triangulate face
    if (joinType == GeomAbs_Arc) {

        TopoDS_Face offset_face = BRepBuilderAPI_MakeFace(TopoDS::Wire(offset_shape)).Face();
        TopoDS_Shape triangulated_shape = polygonize_shape(offset_face, 12.0, false, 12.0, true);

        // get face
        Topo C(triangulated_shape);
        TopoDS_ListOfShape faces = C.faces();

        if (faces.Size() != 1) {
            std::cerr << "Offset: Triangulated shape of arc has less or more than one face!\n";
            return false;
        }

        offset_shape = BRepTools::OuterWire(TopoDS::Face(faces.First()));
    }

    // Check 4: Is offset shape Null
    if (offset_shape.IsNull()) {
        std::cerr << "\t[Warning] Offset wire " << joinType << " is Null." << std::endl;
        return false;
    }

    TopoDS_Wire offsetted_outerWire = TopoDS::Wire(offset_shape);

    // Check 5: Is offset wire Null
    if (offsetted_outerWire.IsNull()) {
        std::cerr << "\t[Warning] No offsetted " << joinType << " outer wire.\n";
        return false;
    }

    // Check 6: Is offset wire area big enough
    double area = ShapeAnalysis::ContourArea(offsetted_outerWire);
    if (area < 1.0e-9) {
        std::cerr << "\t[Warning] Offsetted " << joinType << " outer wire is faulty. Area: " << area << " (before: " << ShapeAnalysis::ContourArea(wire) << ")" << std::endl;
        return false;
    }

    // Finally return
    offsetWire = offsetted_outerWire;
    return true;

}

gp_Pnt Kernel::face_center(const TopoDS_Face &face) {
    GProp_GProps p;
    BRepGProp::SurfaceProperties(face, p);
    return p.CentreOfMass();
}

Bnd_Box Kernel::aabb(const TopoDS_Shape &shape, double gap) {
    Bnd_Box bnd;
    bnd.SetGap(gap);
    BRepBndLib::Add(shape, bnd);
    return bnd;
}

Bnd_OBB Kernel::obb(const TopoDS_Shape &shape) {
    Bnd_OBB obb;
    BRepBndLib::AddOBB(shape, obb);
    return obb;
}

TopoDS_Face Kernel::remove_seam_edges(const TopoDS_Face &f) {

    ShapeBuild_ReShape R;
    bool b = false;

    for (auto &e: Topo(f).edges())
        if (e.Orientation() == TopAbs_INTERNAL || e.Orientation() == TopAbs_EXTERNAL) {
            R.Remove(e);
            b = true;
        }

    if (!b) return f;

    TopoDS_Shape S = TopoDS::Face(R.Apply(f));
    R.Clear();

    for (auto &w: Topo(S).wires())
        if (Topo(w).edges().IsEmpty())
            R.Remove(w);

    S = TopoDS::Face(R.Apply(S));
    R.Clear();

    for (auto &w: Topo(S).faces())
        if (Topo(w).wires().IsEmpty())
            R.Remove(w);

    S = R.Apply(S);

    if (S.ShapeType() == TopAbs_FACE) return TopoDS::Face(S);
    else {
        std::cerr << "[Warning] Seam edge removal failed " << hash(f) << std::endl;
        return f;
    }
}

TopoDS_Shape Kernel::process_shape_for_export(TopoDS_Shape S) {

    // Remove seam edges
    TopExp_Explorer F;
    TopExp_Explorer E;
    ShapeBuild_ReShape R;

    for (F.Init(S, TopAbs_FACE); F.More(); F.Next())
        for (E.Init(F.Current(), TopAbs_EDGE); E.More(); E.Next())
            if (E.Current().Orientation() == TopAbs_INTERNAL || E.Current().Orientation() == TopAbs_EXTERNAL)
                R.Remove(E.Current());

    S = R.Apply(S);

    // Unify faces
    // ShapeUpgrade_UnifySameDomain U(S, true, true, false);
    ShapeUpgrade_UnifySameDomain U(S, false, true, false); // why does unifyEdges need be false, so that ifcopenshell doesn't throw segfault, when serializing shell?
    U.SetAngularTolerance(Precision::Angular() * 2);
    U.SetLinearTolerance(Precision::Confusion() * 2);
    U.Build();
    return TopoDS::Shell(U.Shape());
}

TopoDS_Shape Kernel::unify_shape(TopoDS_Shape &shp) {

    TopoDS_ListOfShape L = Topo(shp).solids();

    if (L.Size() > 1) {
        BRepAlgoAPI_Fuse f;
        TopTools_ListOfShape arg;
        TopTools_ListOfShape tools;

        arg.Append(L.First());
        auto it = L.begin();
        for (std::advance(it, 1); it != L.end(); ++it) // skip first shape
            tools.Append(*it);

        f.SetArguments(arg);
        f.SetTools(tools);
        f.SetFuzzyValue(1.0e-5);
        f.SetToFillHistory(false);
        f.Build();

        if (f.IsDone()) {
            ShapeUpgrade_UnifySameDomain shapeUpgrader(f.Shape(), true, true, false);
            shapeUpgrader.Build();
            return shapeUpgrader.Shape();
        } else
            return shp;

    } else if (L.Size() == 1)
        return shp;
    else {
        auto bnd = aabb(shp, 1.0e-5);
        return aabb_to_shape(bnd);
    }

}

Bnd_OBB Kernel::obb(const TopoDS_ListOfShape &L) {

    Bnd_OBB obb;
    TopoDS_Compound C = compound_from_shape_list(L);
    BRepBndLib::AddOBB(C, obb); // BRepBndLib::AddOBB(S, obb, true, true, true);

    return obb;
}

TopoDS_Shape Kernel::obb_to_shape(const Bnd_OBB &obb) {

    double tol = 1.0e-5;

    if (obb.XHSize() < tol || obb.YHSize() < tol || obb.ZHSize() < tol) {
        std::cerr << "[Error] OBB too flat for shape creation: (" << obb.XHSize() << ", " << obb.YHSize() << ", " << obb.ZHSize() << ")" << std::endl;
        return {};
    }

    gp_Pnt p(obb.Center());
    gp_Pnt m(p.XYZ() - obb.XDirection() * obb.XHSize() - obb.YDirection() * obb.YHSize() - obb.ZDirection() * obb.ZHSize());

    gp_Ax2 axis(p, gp_Dir(obb.ZDirection()), gp_Dir(obb.XDirection()));
    axis.SetLocation(m);

    return BRepPrimAPI_MakeBox(axis, 2.0 * obb.XHSize(), 2.0 * obb.YHSize(), 2.0 * obb.ZHSize()).Shape();
}

Bnd_Box Kernel::create_aabb_from_shape_list(const TopoDS_ListOfShape &L, double gap) {

    Bnd_Box bnd;
    bnd.SetGap(gap);
    for (const auto &S: L)
        BRepBndLib::Add(S, bnd);
    return bnd;
}

TopoDS_Shape Kernel::aabb_to_shape(const Bnd_Box &bnd) {

    if (bnd.IsVoid()) {
        std::cerr << "[Error] AABB to flat to create shape!" << std::endl;
        bnd.Dump();
        return {};
    }

    return BRepPrimAPI_MakeBox(bnd.CornerMin(), bnd.CornerMax()).Shape();
}

TopoDS_Compound Kernel::compound_from_shape_list(const TopoDS_ListOfShape &L) {

    TopoDS_Compound C;
    BRep_Builder B;
    B.MakeCompound(C);
    for (const auto &S: L)
        B.Add(C, S);

    return C;
}

void Kernel::calc_plane_uv_vectors(gp_Pnt O, gp_Pnt P, gp_Dir n, gp_Dir &u_plane, gp_Dir &v_plane) {

    // O ... Origin of face definition
    // P ... Point following O in face wire
    // n ... Face normal
    u_plane = gp_Dir(gp_Vec(O, P));
    v_plane = n.Crossed(u_plane);
}

void Kernel::calc_point_uv_parameters_on_plane(gp_Pnt P, gp_Pnt O, gp_Vec u_plane, gp_Vec v_plane, double &u, double &v) {

    // O ... Origin of face definition
    // P ... Any point on plane
    gp_Vec R(P.X() - O.X(), P.Y() - O.Y(), P.Z() - O.Z()); // point relative to origin point
    u = gp_Vec(u_plane).Dot(R);
    v = gp_Vec(v_plane).Dot(R);
}

gp_Pnt Kernel::move_point_along_scaled_unit_vector(gp_Pnt c, gp_Dir n, double s) {
    gp_Vec v(n);
    v.Scale(s);
    c.Translate(v);
    return c;
}

bool Kernel::is_point_in_face(const TopoDS_Face &face, gp_Pnt P, double tol) {

    BRepClass_FaceClassifier FB;
    FB.Perform(face, P, tol);
    return FB.State() == TopAbs_IN;
}

gp_Pnt Kernel::mid_uv_point_on_surface(const TopoDS_Face &face) {

    Standard_Real umin, umax, vmin, vmax;
    BRepTools::UVBounds(face, umin, umax, vmin, vmax);
    ShapeAnalysis_Surface A(BRep_Tool::Surface(face));
    Standard_Real u = 0.5 * (umax + umin);
    Standard_Real v = 0.5 * (vmax + vmin);
    return A.Value(u, v);

}

gp_Pnt Kernel::random_point_on_face_using_triangulation(const TopoDS_Shape &shp) {

    BRepMesh_IncrementalMesh(shp, 12.0, false, 12.0, false);

    TopExp_Explorer Ex;
    TopoDS_Face Face;
    int i1, i2, i3;
    TopLoc_Location L;

    for (Ex.Init(shp, TopAbs_FACE); Ex.More(); Ex.Next()) {
        Face = TopoDS::Face(Ex.Current());
        break;
    }

    Handle_Poly_Triangulation Poly_Triangulation = BRep_Tool::Triangulation(Face, L);
    if (Poly_Triangulation.IsNull()) {
        std::cerr << "[Warning] Null triangulation." << std::endl;
        return {0, 0, 0};
    }

    auto Poly_Triangle = Poly_Triangulation->Triangles().Value(1);
    Poly_Triangle.Get(i1, i2, i3);

    gp_Pnt p1 = Poly_Triangulation->Node(i1).Transformed(L.Transformation());
    gp_Pnt p2 = Poly_Triangulation->Node(i2).Transformed(L.Transformation());
    gp_Pnt p3 = Poly_Triangulation->Node(i3).Transformed(L.Transformation());

    double tol = 1.0e-6;
    if (p1.IsEqual(p2, tol) || p2.IsEqual(p3, tol) || p3.IsEqual(p1, tol)) {
        std::cerr << "[Warning] Overlapping points." << std::endl;
        return {0, 0, 0};
    }

    if (are_points_colinear(p1, p2, p3, 1.0e-5)) {
        std::cerr << "[Warning] Colinear points." << std::endl;
        return {0, 0, 0};
    }

    Standard_Real x = (p1.X() + p2.X() + p3.X()) / 3;
    Standard_Real y = (p1.Y() + p2.Y() + p3.Y()) / 3;
    Standard_Real z = (p1.Z() + p2.Z() + p3.Z()) / 3;

    return {x, y, z};

}

gp_Pnt Kernel::intersection_face_line(const TopoDS_Face &face, gp_Lin line, double tol) {

    IntCurvesFace_Intersector I(face, tol);
    I.Perform(line, RealFirst(), RealLast());
    return I.NbPnt() == 0 ? gp_Pnt(0, 0, 0) : I.Pnt(1);
}

gp_Pnt Kernel::point_on_face(TopoDS_Face &inputFace, const gp_Dir &n, std::string s) {

    TopoDS_Face face = TopoDS::Face(shape_copy(inputFace)); // create copy to not affect triangulation of inputFace

    gp_Pnt C;
    double tol = 1.0e-6;

    // 1) Centroid:                     Can bring wrong results for non-planar faces, faces with holes or concave faces.
    C = face_center(face);
    if (is_point_in_face(face, C, tol)) return C;

    // 2) Centroid-ray:                 Maybe helps in cases of curved faces
    C = intersection_face_line(face, gp_Lin(C, n), tol);
    if (is_point_in_face(face, C, tol)) return C;

    // 3) First triangle centroid:     Should lie in all cases in face, maybe for small faces it's detected as "on" edge
    C = random_point_on_face_using_triangulation(face);
    if (is_point_in_face(face, C, tol)) return C;

    // 4) First triangle centroid-ray:  Maybe helps in cases of curved faces
    C = intersection_face_line(face, gp_Lin(C, n), tol);
    if (is_point_in_face(face, C, tol)) return C;

    // 5) Mid of UV coordinates:        A lot of faces are trimmed surfaces. In that case bounds of u and v are unfortunately infinite
    C = mid_uv_point_on_surface(face);
    if (is_point_in_face(face, C, tol)) return C;

    std::cerr << "[Warning] No point on face was found. " << s << " " << area(face) << std::endl;

    // 6) Just take the vertex, better than nothing
    TopExp_Explorer Ex;
    for (Ex.Init(face, TopAbs_VERTEX); Ex.More(); Ex.Next()) {
        C = BRep_Tool::Pnt(TopoDS::Vertex(Ex.Current()));
        break;
    }

    return C;
}

double Kernel::area(const TopoDS_Shape &shape) {
    GProp_GProps gprop;
    BRepGProp::SurfaceProperties(shape, gprop);
    return gprop.Mass();
}

TopoDS_Shape Kernel::best_fitting_bbox(const TopoDS_ListOfShape &L) {

    TopoDS_Shape obb_shape = obb_to_shape(obb(L));
    TopoDS_Shape aabb_shape = aabb_to_shape(create_aabb_from_shape_list(L, 0));

    // only use obb if aabb is too big. Reason obb is not used by default is to keep right angles. obb in those cases can be a bit inclined
    return volume(aabb_shape) / volume(obb_shape) < 1.1 ? aabb_shape : obb_shape;
}

TopoDS_Shape Kernel::best_fitting_bbox(const TopoDS_Shape &S) {

    TopoDS_Shape obb_shape = obb_to_shape(obb(S));
    TopoDS_Shape aabb_shape = aabb_to_shape(aabb(S, 0));

    // only use obb if aabb is too big. Reason obb is not used by default is to keep right angles. obb in those cases can be a bit inclined
    return volume(aabb_shape) / volume(obb_shape) < 1.1 ? aabb_shape : obb_shape;
}

double Kernel::polygonize_shape_2a_curvature_distance(const TopoDS_Shape &F_tool, const TopoDS_Shape &common, const segmentation_prism &p) {

    IntCurvesFace_ShapeIntersector Intersector;
    Intersector.Load(common, 0.0001);
    Intersector.PerformNearest(gp_Lin(p.c, p.n), 0, RealLast());
    if (Intersector.NbPnt() == 0) return -1;

    IntCurvesFace_ShapeIntersector Intersector2;
    Intersector2.Load(F_tool, 0.0001);
    Intersector2.PerformNearest(gp_Lin(p.c, p.n.Reversed()), 0, RealLast());

    return Intersector2.NbPnt() == 0 ? Intersector.Pnt(1).Distance(p.c) : Intersector.Pnt(1).Distance(Intersector2.Pnt(1));
}

TopoDS_Shape Kernel::polygonize_shape_2a_curvature(const TopoDS_Shape &S, double linear_deflection, bool isRelative, double angular_deflection, bool isInParallel) {

    // if two curved faces opposite to each other exist

    TopoDS_ListOfShape f_planar, f_curved;

    for (auto &f: Topo(S).faces())
        face_is_planar(TopoDS::Face(f)) ? f_planar.Append(f) : f_curved.Append(f);

    if (f_curved.Size() != 2) return {};

    TopoDS_Shape F_tool, F_arg;

    if (area(f_curved.First()) > area(f_curved.Last())) { // bigger should be tool (>)
        F_tool = TopoDS::Face(f_curved.First());
        F_arg = TopoDS::Face(f_curved.Last());
    } else {
        F_tool = TopoDS::Face(f_curved.Last());
        F_arg = TopoDS::Face(f_curved.First());
    }

    TopoDS_Shape T_tool = polygonize_shape(F_tool, linear_deflection, isRelative, angular_deflection, isInParallel);
    if (T_tool.IsNull()) return {};

    //********************************************************
    std::list<segmentation_prism> prisms;
    double length = 10;

    for (auto &F: Topo(T_tool).faces()) {
        const auto &f = TopoDS::Face(F);
        gp_Dir n = face_normal(f).Reversed();
        TopoDS_Shape prism = BRepPrimAPI_MakePrism(f, gp_Vec(n).Scaled(length)).Shape();
        prisms.emplace_back(prism, f, n, face_center(f));
    }

    if (prisms.empty()) return {};

    //********************************************************
    TopoDS_ListOfShape planes;
    bool rotate_faces = true;

    for (auto &p: prisms) {
        TopoDS_Shape common = BRepAlgoAPI_Common(F_arg, p.S).Shape();

        if (common.IsNull() || Topo(common).faces().IsEmpty()) continue;

        if (rotate_faces) {
            double distance = polygonize_shape_2a_curvature_distance(F_tool, common, p);
            if (distance < 0) continue;
            TopoDS_Shape plane = moved(p.base, gp_Vec(p.n).Scaled(distance));
            planes.Append(plane);
        } else {
            BRepBuilderAPI_MakePolygon m;
            auto L = Topo(Topo(common).wires().First()).ordered_vertices_of_wire();
            for (auto &v: L) m.Add(TopoDS::Vertex(v));
            m.Close();
            TopoDS_Face plane = BRepBuilderAPI_MakeFace(m.Wire()).Face();
            planes.Append(plane);
        }
    }

    if (planes.IsEmpty()) return {};

    //********************************************************
    for (auto &p: prisms)
        F_arg = BRepAlgoAPI_Cut(F_arg, p.S).Shape();

    //********************************************************
    BOPAlgo_Builder fuse;
    for (auto &P: planes)
        fuse.AddArgument(P);
    if (!Topo(F_arg).faces().IsEmpty())
        fuse.AddArgument(polygonize_shape(F_arg, linear_deflection, isRelative, angular_deflection, isInParallel));
    fuse.Perform();

    f_planar.Append(fuse.Shape());
    f_planar.Append(T_tool);

    //********************************************************
    TopoDS_ListOfShape L;
    double offset = 0.1;
    for (auto &f_plan: f_planar)
        for (auto &f: Topo(f_plan).faces()) {
            L.Append(f);
/*            TopoDS_Face offset_face;
            if (offset_faces_outer_wire_cut(TopoDS::Face(f), offset_face, offset)) { // TODO OFFSET WIRE (BRepOffsetAPI_MakeOffset) SEGFAULTS IN PARALLEL
                if (!face_normal(TopoDS::Face(f)).IsEqual(face_normal(offset_face), 0.1))
                    offset_face.Reverse();
            }
            if (!offset_face.IsNull()) L.Append(offset_face);*/
        }

    //********************************************************
    TopoDS_Compound C = compound_from_shape_list(L);

    TopoDS_Shape R = ShapeHealing(C).heal_shape(true);
    return R.IsNull() ? TopoDS_Shape() : R;
}

TopoDS_Shape Kernel::polygonize_shape(const TopoDS_Shape &shp, double linear_deflection, bool isRelative, double angular_deflection, bool isInParallel) {

    BRepMesh_IncrementalMesh(shp, linear_deflection, isRelative, angular_deflection, isInParallel);
    TopoDS_ListOfShape Faces;
    double tol = 1.0e-6;

    TopExp_Explorer Ex;
    for (Ex.Init(shp, TopAbs_FACE); Ex.More(); Ex.Next()) {
        TopoDS_Face Face = TopoDS::Face(Ex.Current());
        TopLoc_Location L;
        Handle_Poly_Triangulation Poly_Triangulation = BRep_Tool::Triangulation(Face, L);

        if (Poly_Triangulation.IsNull()) {
            std::cerr << "[Warning] Null triangulation." << std::endl;
            continue;
        }

        const Poly_Array1OfTriangle &Poly_Array1OfTriangle = Poly_Triangulation->Triangles();

        for (int i = 1; i <= Poly_Triangulation->NbTriangles(); ++i) {

            auto Poly_Triangle = Poly_Array1OfTriangle.Value(i);

            int i1, i2, i3;
            if (Face.Orientation() == TopAbs_REVERSED)
                Poly_Triangle.Get(i1, i3, i2);
            else
                Poly_Triangle.Get(i1, i2, i3);

            gp_Pnt p1 = Poly_Triangulation->Node(i1).Transformed(L.Transformation()); // get gp_pnt of vertices
            gp_Pnt p2 = Poly_Triangulation->Node(i2).Transformed(L.Transformation());
            gp_Pnt p3 = Poly_Triangulation->Node(i3).Transformed(L.Transformation());

            if (p1.IsEqual(p2, tol) || p2.IsEqual(p3, tol) || p3.IsEqual(p1, tol)) {
                std::cerr << "[Warning] Overlapping points." << std::endl;
                continue;
            }

            if (are_points_colinear(p1, p2, p3, 1.0e-5)) {
                std::cerr << "[Warning] Colinear points." << std::endl;
                continue;
            }

            TopoDS_Edge E1 = BRepBuilderAPI_MakeEdge(p1, p2);
            TopoDS_Edge E2 = BRepBuilderAPI_MakeEdge(p2, p3);
            TopoDS_Edge E3 = BRepBuilderAPI_MakeEdge(p3, p1);

            TopoDS_Wire W = BRepBuilderAPI_MakeWire(E1, E2, E3);
            TopoDS_Face F = BRepBuilderAPI_MakeFace(W);
            Faces.Append(F);

        }
    }

    BRepBuilderAPI_Sewing sew;
    for (const auto &F: Faces)
        sew.Add(F);
    sew.Perform();
    TopoDS_Shape shells = sew.SewedShape();

    TopoDS_Shape solids = shells_to_solids(shells);

    ShapeUpgrade_UnifySameDomain shapeUpgrader(solids, true, true, false);
    shapeUpgrader.SetAngularTolerance(1.0e-10);
    //shapeUpgrader.SetLinearTolerance(1.0e-6);
    shapeUpgrader.Build();
    return shapeUpgrader.Shape();
}

TopoDS_Shape Kernel::shells_to_solids(const TopoDS_Shape &shp, bool complement) {

    TopoDS_ListOfShape L;
    for (auto &shell: Topo(shp).shells()) {
        TopoDS_Solid solid = BRepBuilderAPI_MakeSolid(TopoDS::Shell(shell)).Solid();

        if (complement && volume(solid) < 0)
            L.Append(solid.Complemented());
        else
            L.Append(solid);
    }
    return compound_from_shape_list(L);
}

bool Kernel::are_points_colinear(gp_Pnt p1, gp_Pnt p2, gp_Pnt p3, double tol) {

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

    //    gp_Vec v1(p1, p2);
    //    gp_Vec v2(p1, p3);
    //    gp_Vec n1 = gp_Vec(gp_Dir(v1));
    //    gp_Vec n2 = gp_Vec(gp_Dir(v2));
    //    return n1.Crossed(n2).Magnitude() < radiant_tol;
}

bool Kernel::are_points_colinear(gp_Pnt2d p1, gp_Pnt2d p2, gp_Pnt2d p3, double tol) {

    // calculate the smallest (perpendicular) distance of p2 to the vector between p1 and p3
    gp_Vec2d v1(p1, p2);
    gp_Vec2d v2(p1, p3);
    gp_Vec2d v3(p2, p3);

    if (v1.Magnitude() < 1e-12 || v2.Magnitude() < 1e-12 || v3.Magnitude() < 1e-12) {
        std::cerr << "[Warning] Vector length is zero. (Self-intersection?) P1: " << p1.X() << " " << p1.Y() << " P2: " << p2.X() << " " << p2.Y() << " P3: " << p3.X() << " " << p3.Y() << std::endl;
        return true;
    }

    double dist = v1.Crossed(v2) / v2.Magnitude();
    //double dist = fabs((v2.Y() * p2.X()) - (v2.X() * p2.Y()) - (v2.Y() * p1.X()) + (v2.X() * p1.Y())) / v2.Magnitude();

    if (dist > tol) return false;

    // calculate angle between edges around p2
    double a = unsigned_angle_between_points(p1, p2, p3);
    return a > 3.1398473;
}

bool Kernel::polygon_is_convex(const std::vector<gp_Pnt2d> &l) {

    if (l.size() < 3) return false;

    if (l.size() == 3) return true;

    // first point
    const bool b = (Kernel::determinant_three_points(l.back(), l[0], l[1]) >= 0); // true if det > 0

    // in between
    for (unsigned int i = 1; i < l.size() - 1; i++)
        if (b != (Kernel::determinant_three_points(l[i - 1], l[i], l[i + 1]) >= 0))
            return false;

    // last
    if (b != (Kernel::determinant_three_points(l[l.size() - 2], l.back(), l[0]) >= 0))
        return false;

    return true;
}

double Kernel::unsigned_angle_between_points(gp_Pnt p1, gp_Pnt mid, gp_Pnt p3) {
    gp_Vec v1(mid, p1);
    gp_Vec v2(mid, p3);
    return fabs(v1.Angle(v2));
}

double Kernel::unsigned_angle_between_points(gp_Pnt2d p1, gp_Pnt2d mid, gp_Pnt2d p3) {
    gp_Vec2d v1(mid, p1);
    gp_Vec2d v2(mid, p3);
    return fabs(v1.Angle(v2));
}

bool Kernel::are_points_coincident(gp_Pnt p1, gp_Pnt p2, double tol) { return p1.Distance(p2) < tol; }

bool Kernel::are_points_coincident(gp_Pnt2d p1, gp_Pnt2d p2, double tol) { return p1.Distance(p2) < tol; }

gp_Pnt Kernel::shape_center(const TopoDS_Shape &S) {
    GProp_GProps p;
    BRepGProp::VolumeProperties(S, p);
    return p.CentreOfMass();
}

TopoDS_Shape Kernel::fuse_shape(TopoDS_Shape &S, double tol, const std::string &s) {

    BOPAlgo_Builder builder;
    for (const auto &f: Topo(S).faces()) builder.AddArgument(f);
    builder.SetFuzzyValue(tol);
    builder.SetNonDestructive(true);
    builder.SetToFillHistory(false);
    builder.Perform();

    if (builder.HasWarnings()) {
        std::cerr << "Warnings: " << s << std::endl;
        builder.DumpWarnings(std::cerr);
    }
    if (builder.HasErrors()) {
        std::cerr << "Errors: " << s << std::endl;
        builder.DumpErrors(std::cerr);
    }
    if (!builder.HasModified())
        std::cerr << "Nothing was modified." << std::endl;

    BRepBuilderAPI_Sewing sew;
    for (auto &f: Topo(builder.Shape()).faces()) sew.Add(f);
    sew.Perform();
    return sew.SewedShape();
}

TopoDS_Shape Kernel::prism_from_face(const TopoDS_Face &F, gp_Dir n, double l, bool rm_tf, bool rm_bf) {

    gp_Vec v(n);
    v.Scale(l);
    TopoDS_Shape P = BRepPrimAPI_MakePrism(F, v).Shape();

    if (rm_tf || rm_bf) {

        ShapeBuild_ReShape R;
        if (rm_bf) R.Remove(F); // remove extruded face

        if (rm_tf) { // remove on other side of extruded face
            for (const auto &face: Topo(P).faces()) {
                if (F.IsSame(face)) continue;
                gp_Dir n2 = face_normal(TopoDS::Face(face));
                gp_Dir n3 = n2.Reversed();

                if (n.IsEqual(n2, 0.017) || n.IsEqual(n3, 0.017)) {
                    R.Remove(face);
                    break;
                }
            }
        }
        P = R.Apply(P);
    }
    return P;
}

std::vector<gp_Pnt> Kernel::remove_colinear_points_from_polygon(std::vector<gp_Pnt> Pnts, double tol) {

    std::vector<gp_Pnt> L;

    gp_Pnt F = Pnts.front();
    gp_Pnt B = Pnts.back();

    Pnts.insert(Pnts.begin(), B);
    Pnts.push_back(F);

    unsigned int start(0);

    for (unsigned int i = 1; i < Pnts.size() - 1; i++) {
        if (!are_points_colinear(Pnts[start], Pnts[i], Pnts[i + 1], tol)) {
            start = i;
            L.push_back(Pnts[i]);
        }
    }

    return L;
}

std::vector<gp_Pnt2d> Kernel::remove_colinear_points_from_polygon(std::vector<gp_Pnt2d> Pnts, double tol) {

    std::vector<gp_Pnt2d> L;

    gp_Pnt2d F = Pnts.front();
    gp_Pnt2d B = Pnts.back();

    Pnts.insert(Pnts.begin(), B);
    Pnts.push_back(F);

    unsigned int start(0);

    for (unsigned int i = 1; i < Pnts.size() - 1; i++) {
        if (!are_points_colinear(Pnts[start], Pnts[i], Pnts[i + 1], tol)) {
            start = i;
            L.push_back(Pnts[i]);
        }
    }

    return L;
}

std::vector<gp_Pnt> Kernel::sort_points_from_polygon_by_min(std::vector<gp_Pnt> Pnts, cFace *cface) {

    std::map<std::pair<double, double>, unsigned int> M;

    // find "minimal point" by map
    auto a = cface->ProjectionAxis();

    if (a == 0)
        for (unsigned int i = 0; i < Pnts.size(); i++)
            M[std::make_pair(Pnts[i].Y(), Pnts[i].Z())] = i;
    else if (a == 1)
        for (unsigned int i = 0; i < Pnts.size(); i++)
            M[std::make_pair(Pnts[i].X(), Pnts[i].Z())] = i;
    else
        for (unsigned int i = 0; i < Pnts.size(); i++)
            M[std::make_pair(Pnts[i].X(), Pnts[i].Y())] = i;

    unsigned int p = M.begin()->second;

    // start loop with minimal point
    std::vector<gp_Pnt> L;

    for (unsigned int i = p; i < Pnts.size(); i++)
        L.push_back(Pnts[i]);

    for (unsigned int i = 0; i < p; i++)
        L.push_back(Pnts[i]);

    return L;
}

std::vector<gp_Pnt> Kernel::sort_points_from_polygon_by_dist(std::vector<gp_Pnt> Pnts, cFace *cface) {

    std::map<double, unsigned int> M;

    // find "minimal point" by map
    auto a = cface->ProjectionAxis();
    gp_Pnt2d R(-1e6, -1e6);

    if (a == 0)
        for (unsigned int i = 0; i < Pnts.size(); i++)
            M[R.SquareDistance(gp_Pnt2d(Pnts[i].Y(), Pnts[i].Z()))] = i;
    else if (a == 1)
        for (unsigned int i = 0; i < Pnts.size(); i++)
            M[R.SquareDistance(gp_Pnt2d(Pnts[i].X(), Pnts[i].Z()))] = i;
    else
        for (unsigned int i = 0; i < Pnts.size(); i++)
            M[R.SquareDistance(gp_Pnt2d(Pnts[i].X(), Pnts[i].Y()))] = i;

    unsigned int p = M.begin()->second;

    // start loop with minimal point
    std::vector<gp_Pnt> L;

    for (unsigned int i = p; i < Pnts.size(); i++)
        L.push_back(Pnts[i]);

    for (unsigned int i = 0; i < p; i++)
        L.push_back(Pnts[i]);

    return L;
}

std::vector<gp_Pnt> Kernel::remove_coincident_points_from_polygon(std::vector<gp_Pnt> Pnts, double tol) {

    std::vector<gp_Pnt> L;

    gp_Pnt F = Pnts.front();
    gp_Pnt B = Pnts.back();

    Pnts.insert(Pnts.begin(), B);
    Pnts.push_back(F);

    unsigned int start(0);

    for (unsigned int i = 1; i < Pnts.size() - 1; i++) {
        if (!are_points_coincident(Pnts[start], Pnts[i], tol)) {
            start = i;
            L.push_back(Pnts[i]);
        }
    }

    return L;
}

std::vector<gp_Pnt2d> Kernel::remove_coincident_points_from_polygon(std::vector<gp_Pnt2d> Pnts, double tol) {

    std::vector<gp_Pnt2d> L;

    gp_Pnt2d F = Pnts.front();
    gp_Pnt2d B = Pnts.back();

    Pnts.insert(Pnts.begin(), B);
    Pnts.push_back(F);

    unsigned int start(0);

    for (unsigned int i = 1; i < Pnts.size() - 1; i++) {
        if (!are_points_coincident(Pnts[start], Pnts[i], tol)) {
            start = i;
            L.push_back(Pnts[i]);
        }
    }

    return L;
}

bool Kernel::is_edge_line(const TopoDS_Edge &E) {

    return BRepAdaptor_Curve(E).GetType() == GeomAbs_Line;

    auto curve = BRepAdaptor_Curve(E).Curve().Curve();
    if (curve.IsNull()) return false;

    if (curve->DynamicType() == STANDARD_TYPE(Geom_Line)) return true;
    else {
        if (curve->DynamicType() == STANDARD_TYPE(Geom_TrimmedCurve))
            return (Handle_Geom_TrimmedCurve::DownCast(curve)->BasisCurve()->DynamicType() == STANDARD_TYPE(Geom_Line));
        else
            return false;
    }
}

TopoDS_Face Kernel::rebuild_face(const TopoDS_Face &f, bool &success) {

    success = false;

    if (!face_is_polygon(TopoDS::Face(f))) {
        std::cerr << "[Info] Skip non-polygonal face." << std::endl;
        return f;
    }

    gp_Dir n = face_normal(f);
    gp_Pnt A = BRep_Tool::Pnt(TopoDS::Vertex(Topo(f).vertices().First()));

    TopoDS_Wire oW = BRepTools::OuterWire(f);
    if (oW.IsNull()) return f;

    TopoDS_Wire oW2 = rebuild_wire_projection(oW, n, A);
    if (oW2.IsNull()) return f;

    TopoDS_Face F2 = BRepBuilderAPI_MakeFace(oW2).Face();
    if (F2.IsNull()) return f;

    // face normal check
    if (!face_normal(F2).IsEqual(n, 0.1)) F2.Complement();

    // inner wires
    TopoDS_ListOfShape iWs = Topo(f).wires();
    iWs.Remove(oW);

    if (iWs.IsEmpty()) {
        success = true;
        return F2;
    }

    // add holes
    BRepBuilderAPI_MakeFace M(F2);
    for (const auto &iW: iWs) {
        TopoDS_Wire iW2 = rebuild_wire_projection(TopoDS::Wire(iW), n, A);
        if (iW2.IsNull()) continue;
        M.Add(iW2);
    }
    F2 = M.Face();
    if (F2.IsNull()) return f;

    success = true;
    return F2;
}

TopoDS_Wire Kernel::rebuild_wire(const TopoDS_Wire &w) {

    // points
    std::vector<gp_Pnt> P;
    for (auto &v: Topo(w).ordered_vertices_of_wire())
        P.push_back(BRep_Tool::Pnt(TopoDS::Vertex(v)));

    if (P.size() < 3) {
        std::cerr << "[Error] Wire has insufficient number of points!" << std::endl;
        return {};
    }

    // vertices
    std::vector<TopoDS_Vertex> V;
    for (auto &p: P)
        V.push_back(BRepBuilderAPI_MakeVertex(p).Vertex());

    // edges
    TopoDS_ListOfShape E;
    for (unsigned int i = 0; i < P.size() - 1; i++)
        E.Append(BRepBuilderAPI_MakeEdge(V[i], V[i + 1]).Edge());
    E.Append(BRepBuilderAPI_MakeEdge(V.back(), V[0]).Edge());

    // wire
    BRepBuilderAPI_MakeWire M;
    for (auto &e: E)
        M.Add(TopoDS::Edge(e));

    if (M.IsDone())
        return M.Wire();
    else {
        std::cerr << "[Error] Wire could not be build!" << std::endl;
        return {};
    }

}

TopoDS_Wire Kernel::rebuild_wire_projection(const TopoDS_Wire &w, gp_Dir n, gp_Pnt A) {

    // points
    std::vector<gp_Pnt> P;
    for (auto &v: Topo(w).ordered_vertices_of_wire()) {
        gp_Pnt p = BRep_Tool::Pnt(TopoDS::Vertex(v));
        gp_Pnt p2 = project_point_onto_plane(p, n, A);
        P.push_back(p2);
    }

    if (P.size() < 3) {
        std::cerr << "[Error] Wire has insufficient number of points!" << std::endl;
        return {};
    }

    // vertices
    std::vector<TopoDS_Vertex> V;
    for (auto &p: P)
        V.push_back(BRepBuilderAPI_MakeVertex(p).Vertex());

    // edges
    TopoDS_ListOfShape E;
    for (unsigned int i = 0; i < P.size() - 1; i++)
        E.Append(BRepBuilderAPI_MakeEdge(V[i], V[i + 1]).Edge());
    E.Append(BRepBuilderAPI_MakeEdge(V.back(), V[0]).Edge());

    // wire
    BRepBuilderAPI_MakeWire M;
    for (auto &e: E)
        M.Add(TopoDS::Edge(e));

    if (M.IsDone())
        return M.Wire();
    else {
        std::cerr << "[Error] Wire could not be build!" << std::endl;
        return {};
    }

}

double Kernel::area_by_triangulation(const TopoDS_Shape &shp) {

    BRepMesh_IncrementalMesh(shp, 0.5);
    double A = 0;
    double tol = 1.0e-6;

    TopExp_Explorer Ex;
    for (Ex.Init(shp, TopAbs_FACE); Ex.More(); Ex.Next()) {
        TopoDS_Face Face = TopoDS::Face(Ex.Current());
        TopLoc_Location L;
        Handle_Poly_Triangulation Poly_Triangulation = BRep_Tool::Triangulation(Face, L);

        if (Poly_Triangulation.IsNull()) {
            std::cerr << "[Warning] Null triangulation." << std::endl;
            continue;
        }

        const Poly_Array1OfTriangle &Poly_Array1OfTriangle = Poly_Triangulation->Triangles();

        for (int i = 1; i <= Poly_Triangulation->NbTriangles(); ++i) {

            auto Poly_Triangle = Poly_Array1OfTriangle.Value(i);

            int i1, i2, i3;
            Poly_Triangle.Get(i1, i2, i3);
            gp_Pnt p1 = Poly_Triangulation->Node(i1).Transformed(L.Transformation());
            gp_Pnt p2 = Poly_Triangulation->Node(i2).Transformed(L.Transformation());
            gp_Pnt p3 = Poly_Triangulation->Node(i3).Transformed(L.Transformation());

            if (p1.IsEqual(p2, tol) || p2.IsEqual(p3, tol) || p3.IsEqual(p1, tol)) {
                std::cerr << "[Warning] Overlapping points." << std::endl;
                continue;
            }

            if (are_points_colinear(p1, p2, p3, 1.0e-5)) {
                std::cerr << "[Warning] Colinear points." << std::endl;
                continue;
            }

            gp_Vec v12(p2.X() - p1.X(), p2.Y() - p1.Y(), p2.Z() - p1.Z());
            gp_Vec v13(p3.X() - p1.X(), p3.Y() - p1.Y(), p3.Z() - p1.Z());
            A += 0.5 * (v12.CrossMagnitude(v13));
        }
    }
    return A;
}

gp_Pnt Kernel::project_point_onto_plane(gp_Pnt P, gp_Dir n, gp_Pnt A) {
    // plane defined by normal and point A
    gp_Vec p(P.XYZ());
    gp_Vec a(A.XYZ());
    gp_Vec v(n.XYZ());

    v.Multiply(p.Subtracted(a).Dot(n));
    return {p.Subtracted(v).XYZ()};
}

std::vector<gp_Pnt> Kernel::vertex_to_point_list(const TopoDS_ListOfShape &V) {
    std::vector<gp_Pnt> P;
    P.reserve(V.Size());
    for (const auto &t: V)
        P.push_back(BRep_Tool::Pnt(TopoDS::Vertex(t)));
    return P;
}

double Kernel::determinant_three_points(gp_Pnt2d V1, gp_Pnt2d V2, gp_Pnt2d V3) {
    gp_Vec2d v1(V2.X() - V1.X(), V2.Y() - V1.Y());
    gp_Vec2d v2(V3.X() - V2.X(), V3.Y() - V2.Y());
    return v1.X() * v2.Y() - v1.Y() * v2.X();
}

std::vector<gp_Pnt2d> Kernel::project_wire_to_2D(long a, const TopoDS_Wire &W) {

    // a ... Projection Axis
    TopoDS_ListOfShape V = Topo(W).ordered_vertices_of_wire();
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

int Kernel::double_to_scaled_int(double d, int s) { return int(d * s); }

std::vector<std::string> Kernel::split_string_by_character(const std::string &s, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter))
        tokens.push_back(token);
    return tokens;
}

double Kernel::round_double_for_ifc_write(double d) const { return std::round(d * round_double_ifc_write) * round_ifc_write_inv; }

double Kernel::round_double_one_digit(double d) { return std::round(d * 10) * 0.1; }

double Kernel::round_double_two_digits(double d) { return std::round(d * 100) * 0.01; }

double Kernel::round_double_three_digits(double d) { return std::round(d * 1000) * 0.001; }

bool Kernel::do_aabbs_intersect(double A_Xmin, double A_Xmax, double A_Ymin, double A_Ymax, double B_Xmin, double B_Xmax, double B_Ymin, double B_Ymax) {
    return !(B_Xmin > A_Xmax || B_Xmax < A_Xmin || B_Ymax < A_Ymin || B_Ymin > A_Ymax);
}

double Kernel::overlapping_volume_aabbs(double A_Xmin, double A_Xmax, double A_Ymin, double A_Ymax, double A_Zmin, double A_Zmax, double B_Xmin, double B_Xmax, double B_Ymin, double B_Ymax, double B_Zmin, double B_Zmax) {

    double dx = std::min(A_Xmax, B_Xmax) - std::max(A_Xmin, B_Xmin);
    double dy = std::min(A_Ymax, B_Ymax) - std::max(A_Ymin, B_Ymin);
    double dz = std::min(A_Zmax, B_Zmax) - std::max(A_Zmin, B_Zmin);
    return dx * dy * dz;
}

double Kernel::round_double_to_n_decimal_places(double d, unsigned int n) {
    unsigned int f = pow(10, n);
    return std::round(d * f) / f;
}

gp_Pnt Kernel::project_inv(double X, double Y, long proj_axis, gp_Dir n, gp_Pnt R) {

    double a = R.X() * n.X() + R.Y() * n.Y() + R.Z() * n.Z();

    std::vector<double> P_3D(3);
    double t;

    if (proj_axis == 0) {
        t = n.X();
        P_3D[0] = 0.0;
        P_3D[1] = X;
        P_3D[2] = Y;
    } else if (proj_axis == 1) {
        t = n.Y();
        P_3D[0] = X;
        P_3D[1] = 0.0;
        P_3D[2] = Y;
    } else {
        t = n.Z();
        P_3D[0] = X;
        P_3D[1] = Y;
        P_3D[2] = 0.0;
    }

    double z = a;
    z -= P_3D[0] * n.X();
    z -= P_3D[1] * n.Y();
    z -= P_3D[2] * n.Z();
    z /= t;
    P_3D[proj_axis] = z;
    return {P_3D[0], P_3D[1], P_3D[2]};
}

std::list<std::array<gp_Pnt, 3>> Kernel::triangles_from_face(const TopoDS_Face &face) {

    double tol = 1.0e-5;
    unsigned int ct = 0;
    std::list<std::array<gp_Pnt, 3>> triangles;

    TopLoc_Location L;
    Handle_Poly_Triangulation Poly_Triangulation = BRep_Tool::Triangulation(face, L);

    if (Poly_Triangulation.IsNull()) {
        std::cerr << "[Warning] Null triangulation." << std::endl;
        return triangles;
    }

    const Poly_Array1OfTriangle &Poly_Array1OfTriangle = Poly_Triangulation->Triangles();

    for (int i = 1; i <= Poly_Triangulation->NbTriangles(); ++i) {

        auto Poly_Triangle = Poly_Array1OfTriangle.Value(i);

        int i1, i2, i3;
        if (face.Orientation() == TopAbs_REVERSED)
            Poly_Triangle.Get(i1, i3, i2);
        else
            Poly_Triangle.Get(i1, i2, i3);

        gp_Pnt p1 = Poly_Triangulation->Node(i1).Transformed(L.Transformation());
        gp_Pnt p2 = Poly_Triangulation->Node(i2).Transformed(L.Transformation());
        gp_Pnt p3 = Poly_Triangulation->Node(i3).Transformed(L.Transformation());

        if (p1.IsEqual(p2, tol) || p2.IsEqual(p3, tol) || p3.IsEqual(p1, tol)) {
            std::cerr << "[Warning] Overlapping points." << std::endl;
            continue;
        }
        if (are_points_colinear(p1, p2, p3, tol)) {
            std::cerr << "[Warning] Collinear points." << std::endl;
            continue;
        }

        std::array<gp_Pnt, 3> a = {p1, p2, p3};
        triangles.push_back(a);
    }
    return triangles;
}

TopoDS_Shape Kernel::sew_cfaces(const std::set<cFace *> &faces, double tol, bool nonmanifold) {
    // single faces are not stored in a shell but as single topological entity in the returned shape.
    BRepBuilderAPI_Sewing sew;

    for (const auto &cface: faces)
        sew.Add(cface->face);
    sew.SetTolerance(tol);
    sew.SetNonManifoldMode(nonmanifold);
    sew.Perform();
    return sew.SewedShape();
}

double Kernel::minimal_distance_between_shapes(const TopoDS_Shape &b1, const TopoDS_Shape &b2) {
    BRepExtrema_DistShapeShape d;
    d.LoadS1(b1);
    d.LoadS2(b2);
    d.SetFlag(Extrema_ExtFlag_MIN);
    d.Perform();
    return d.IsDone() ? d.Value() : 1e12;
}

gp_Pnt2d Kernel::centroid_2D_polygon(const std::vector<gp_Pnt2d> &poly) {

    std::pair<double, double> ans = {0, 0};

    auto n = poly.size();
    double signedArea = 0;

    for (int i = 0; i < poly.size(); i++) {

        double x0 = poly[i].X(), y0 = poly[i].Y();
        double x1 = poly[(i + 1) % n].X(), y1 = poly[(i + 1) % n].Y();

        // Calculate value of A using shoelace formula
        double A = (x0 * y1) - (x1 * y0);
        signedArea += A;

        // Calculating coordinates of  centroid of polygon
        ans.first += (x0 + x1) * A;
        ans.second += (y0 + y1) * A;
    }

    signedArea *= 0.5;
    ans.first = (ans.first) / (6 * signedArea);
    ans.second = (ans.second) / (6 * signedArea);

    return {ans.first, ans.second};
}

bool Kernel::is_point_inside_2D_polygon(const std::vector<gp_Pnt2d> &poly, const gp_Pnt2d &p) {

    double angle = 0;
    auto n = poly.size();

    for (unsigned int i = 0; i < n; i++)
        angle += Angle2D(
                poly[i].X() - p.X(),
                poly[i].Y() - p.Y(),
                poly[(i + 1) % n].X() - p.X(),
                poly[(i + 1) % n].Y() - p.Y()
        );

    return fabs(angle) > M_PI;
}

double Kernel::Angle2D(double x1, double y1, double x2, double y2) {

    double dtheta, theta1, theta2;

    theta1 = atan2(y1, x1);
    theta2 = atan2(y2, x2);
    dtheta = theta2 - theta1;
    while (dtheta > M_PI)
        dtheta -= 2 * M_PI;
    while (dtheta < -M_PI)
        dtheta += 2 * M_PI;

    return (dtheta);
}

gp_Pnt Kernel::mid_point(gp_Pnt P1, gp_Pnt P2) {
    gp_Vec vec = 0.5 * (gp_Vec(P1.XYZ()) + gp_Vec(P2.XYZ()));
    return {vec.XYZ()};
}

bool Kernel::check_face(const TopoDS_Face &F, unsigned int n, const std::string &g, double l_warn, double l_crit) {

    auto start = std::chrono::high_resolution_clock::now();

    ShapeAnalysis_ShapeContents a;
    bool good_face = true;
    unsigned int hF = hash(F);
    bool planar = face_is_planar(F);
    bool linear_edges = Topo(F).all_edges_non_curved();

    if (F.IsNull()) {
        std::cerr << "[Warning] Null face. " << n << "\t" << hF << "\t" << g << "\t" << planar << "\t" << linear_edges << std::endl;
        good_face = false;
    }

    double A = area(F);
    if (A < 1.0e-6)
        std::cerr << "[Warning] Small area. " << n << "\t" << hF << "\t" << g << "\t" << A << "\t" << planar << "\t" << linear_edges << std::endl;

    double A_crit = 1.0e-8;
    if (A < A_crit) good_face = false;

    ShapeAnalysis_CheckSmallFace saf;
    if (saf.CheckSpotFace(F, A_crit)) {
        std::cerr << "[Warning] Spot face. " << n << "\t" << hF << "\t" << g << "\t" << A << "\t" << planar << "\t" << linear_edges << std::endl;
        good_face = false;
    }

    a.Clear();
    a.Perform(F);
    if (linear_edges && a.NbVertices() < 3) std::cerr << "[Error] Critical vertex number. " << n << "\t" << hF << "\t" << g << "\t" << a.NbVertices() << "\t" << planar << "\t" << linear_edges << std::endl;

    if (linear_edges) {
        TopoDS_Wire outerWire = BRepTools::OuterWire(F);

        std::set<unsigned int> hashes;
        for (auto &v: Topo(outerWire).vertices()) hashes.insert(Kernel::hash(v));
        auto nv1 = hashes.size();
        auto nv2 = Topo(outerWire).ordered_vertices_of_wire().Size(); // will return less than correct number. maybe stops at self-intersection

        if (nv1 != nv2) {
            std::cerr << "[Warning] Seems like there is a self-intersection in outer wire (maybe in one point). " << n << "\t" << hF << "\t" << g << "\t" << nv1 << "\t" << nv2 << "\t" << planar << "\t" << linear_edges << std::endl;
            good_face = false;
        }

        if (nv2 < 3) {
            std::cerr << "[Error] Outer wire consists of less than three vertices. " << n << "\t" << hF << "\t" << g << "\t" << nv2 << "\t" << planar << "\t" << linear_edges << std::endl;
            good_face = false;
        }
    }

    if (linear_edges && a.NbEdges() < 3) {
        std::cerr << "[Error] Critical edge number. " << n << "\t" << hF << "\t" << g << "\t" << a.NbEdges() << "\t" << planar << "\t" << linear_edges << std::endl;
        good_face = false;
    }
    if (a.NbWires() == 0) {
        std::cerr << "[Error] Critical wire number. " << n << "\t" << hF << "\t" << g << "\t" << a.NbWires() << "\t" << planar << "\t" << linear_edges << std::endl;
        good_face = false;
    }

    //*****************************************************
    for (const auto &v: Topo(F).vertices()) {
        if (v.IsNull()) std::cerr << "[Error] Null vertex. " << n << "\t" << hF << "\t" << g << "\t" << hash(v) << "\t" << planar << "\t" << linear_edges << std::endl;

        if (v.Orientation() == TopAbs_INTERNAL || v.Orientation() == TopAbs_EXTERNAL)
            std::cerr << "[Warning] Bad orientation vertex. " << n << "\t" << hF << "\t" << g << "\t" << hash(v) << "\t" << v.Orientation() << "\t" << planar << "\t" << linear_edges << std::endl;
    }
    //*****************************************************

    //*****************************************************
    for (const auto &ed: Topo(F).edges()) {

        TopoDS_Edge e = TopoDS::Edge(ed);
        bool is_line = is_edge_line(e);

        if (e.IsNull()) std::cerr << "[Error] Null edge. " << n << "\t" << hF << "\t" << g << "\t" << hash(e) << "\t" << planar << "\t" << linear_edges << std::endl;

        if (BRepAdaptor_Curve(e).Curve().Curve().IsNull()) { // BRepAdaptor_Curve(e).Is3DCurve()
            if (is_line) {
                std::cerr << "[Error] Edge has no 3D curve. " << n << "\t" << hF << "\t" << g << "\t" << hash(e) << "\t" << planar << "\t" << linear_edges << std::endl;
                good_face = false;
            } else
                std::cerr << "[Warning] Edge has no 3D curve. " << n << "\t" << hF << "\t" << g << "\t" << hash(e) << "\t" << planar << "\t" << linear_edges << std::endl;
        }

        /*
         * Degenerated edges have no 3D-curve, so they have zero length. Degenerated edges have 2D-curve. The parametric range of the curve depends on surface to which degenerated edge belongs to.
         * Generally the existence or absence of degenerated edges depends on type of parametric representation of a surface and nothing more.
         * For e.g degenerated edges are nonsense for planes and degenerated edges must be definitely removed from planes.
         * On the other hand degenerated edges are important for cones, spheres, revolutions,
         */
        if (BRep_Tool::Degenerated(e)) { // BRepAdaptor_Curve(e).Is3DCurve()
            if (is_line) {
                std::cerr << "[Error] Edge is degenerated. " << n << "\t" << hF << "\t" << g << "\t" << hash(e) << "\t" << planar << "\t" << linear_edges << std::endl;
                good_face = false;
            } else std::cerr << "[Warning] Edge is degenerated. " << n << "\t" << hF << "\t" << g << "\t" << hash(e) << "\t" << planar << "\t" << linear_edges << std::endl;
        }

        if (e.Orientation() == TopAbs_INTERNAL)
            std::cerr << "[Warning] Seam edge. " << n << "\t" << hF << "\t" << g << "\t" << hash(e) << "\t" << e.Orientation() << "\t" << planar << "\t" << linear_edges << std::endl;
        else if (e.Orientation() == TopAbs_EXTERNAL) {
            std::cerr << "[Error] Bad orientation edge. " << n << "\t" << hF << "\t" << g << "\t" << hash(e) << "\t" << e.Orientation() << "\t" << planar << "\t" << linear_edges << std::endl;
            good_face = false;
        }

        double l = length(e);
        if (l < l_warn) {
            if (l < l_crit && !BRep_Tool::Degenerated(e)) {
                std::cerr << "[Error] Super short edge. " << n << "\t" << hF << "\t" << g << "\t" << hash(e) << "\t" << l << "\t" << planar << "\t" << linear_edges << std::endl;
                good_face = false;
            } else
                std::cerr << "[Warning] Short edge. " << n << "\t" << hF << "\t" << g << "\t" << hash(e) << "\t" << l << "\t" << planar << "\t" << linear_edges << std::endl;
        }

        auto vs = Topo(e).vertices();

        if (vs.Size() < 2) std::cerr << "[Error] No or one vertex on edge. " << n << "\t" << hF << "\t" << g << "\t" << hash(e) << "\t" << planar << "\t" << linear_edges << std::endl;
        else {

            TopoDS_Vertex v1;
            TopoDS_Vertex v2;

            if (e.Orientation() == TopAbs_INTERNAL) {
                v1 = TopoDS::Vertex(vs.First());
                v2 = TopoDS::Vertex(vs.Last());
            } else {
                v1 = TopExp::FirstVertex(e, Standard_True);
                v2 = TopExp::LastVertex(e, Standard_True);
            }

            if (is_line && v1.IsSame(v2)) {
                std::cerr << "[Error] Duplicate vertices on edge. " << n << "\t" << hF << "\t" << g << "\t" << hash(e) << "\t" << planar << "\t" << linear_edges << std::endl;
                good_face = false;
            }

            if (v1.IsNull()) {
                std::cerr << "[Error] Null vertex in edge. " << n << "\t" << hF << "\t" << g << "\t" << hash(e) << "\t" << planar << "\t" << linear_edges << std::endl;
                good_face = false;
            }
            if (v2.IsNull()) {
                std::cerr << "[Error] Null vertex in edge. " << n << "\t" << hF << "\t" << g << "\t" << hash(e) << "\t" << planar << "\t" << linear_edges << std::endl;
                good_face = false;
            }
        }
        for (const auto &v: vs)
            if (v.IsNull()) std::cerr << "[Error] Null vertex (v2). " << n << "\t" << hF << "\t" << g << "\t" << hash(v) << "\t" << planar << "\t" << linear_edges << std::endl;
    }
    //*****************************************************

    //*****************************************************
    for (const auto &w: Topo(F).wires()) {
        if (w.IsNull()) std::cerr << "[Error] Null wire. " << n << "\t" << hF << "\t" << g << "\t" << hash(w) << "\t" << planar << "\t" << linear_edges << std::endl;
        auto edges = Topo(w).edges();
        if (edges.Size() < 3) {
            if (edges.Size() == 1 && !is_edge_line(TopoDS::Edge(edges.First())))
                std::cerr << "[Warning] Less than three edges in wire. " << n << "\t" << hF << "\t" << g << "\t" << hash(w) << "\t" << planar << "\t" << linear_edges << std::endl;
            else if (edges.Size() == 1 && edges.First().Orientation() == TopAbs_INTERNAL)
                std::cerr << "[Warning] Less than three edges in wire. " << n << "\t" << hF << "\t" << g << "\t" << hash(w) << "\t" << planar << "\t" << linear_edges << std::endl;
            else if (edges.Size() == 2 && edges.First().Orientation() == TopAbs_INTERNAL && edges.Last().Orientation() == TopAbs_INTERNAL)
                std::cerr << "[Warning] Less than three edges in wire. " << n << "\t" << hF << "\t" << g << "\t" << hash(w) << "\t" << planar << "\t" << linear_edges << std::endl;
            else {
                std::cerr << "[Error] Less than three edges in wire. " << n << "\t" << hF << "\t" << g << "\t" << hash(w) << "\t" << planar << "\t" << linear_edges << std::endl;
                good_face = false;
            }
        }
    }
    //*****************************************************

    //*****************************************************
    if (F.Orientation() == TopAbs_INTERNAL || F.Orientation() == TopAbs_EXTERNAL) {
        std::cerr << "[Error] Bad orientation face. " << n << "\t" << hF << "\t" << g << "\t" << F.Orientation() << "\t" << planar << "\t" << linear_edges << std::endl;
        good_face = false;
    }
    //*****************************************************

    return good_face;
}