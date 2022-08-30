// Copyright 2022 Eric Fichter
#include "ShapeHealing.h"

ShapeHealing::ShapeHealing(TopoDS_Shape _Shp) : Shp(std::move(_Shp)) {

    component_id_counter = -1;
    fuse_tol = 0.0001;
}

TopoDS_Shape ShapeHealing::heal_shape(bool ignore_curved_edges) {

    if (Shp.IsNull())
        return {};

    remove_duplicated_hashes_from_shape();

    fuse();

    bool has_curved_edges, has_curved_faces;
    curvatures(has_curved_edges, has_curved_faces);

    if (has_curved_faces || (has_curved_edges && !ignore_curved_edges)) {
        std::cerr << "[Info] Shape to be healed will be triangulated." << std::endl;
        Shp = Kernel::polygonize_shape(Kernel::shape_copy(Shp), 0.8, false, 0.8, true);
        if (Shp.IsNull())
            return {};
        else
            fuse();
    }

    if (Shp.IsNull())
        return {};

    remove_bad_faces();

    if (Shp.IsNull())
        return {};

    TopExp::MapShapesAndAncestors(Shp, TopAbs_EDGE, TopAbs_FACE, edgeFaceMap);

    fill_face_list();

    if (sFaces.empty())
        return {};

    std::list<unsigned int> bad_hashes;
    if (!find_duplicated_hashes(bad_hashes))
        return {};

    //TODO eliminate duplicate hashes

    for (auto &sface: sFaces)
        sface.update_half_edges();

    setup_sFace_map();

    for (auto &sface: sFaces)
        sface.update_face_adjacencies(edgeFaceMap, id2sFace);

    identify_hanging_faces();

    remove_trash_and_face_adjacency();

    if (sFaces.empty())
        return {};

    find_components();

    TopoDS_ListOfShape Shells = components_to_solids();

    Shp = Kernel::compound_from_shape_list(Shells);

    return Shp;
}

void ShapeHealing::curvatures(bool &has_curved_edges, bool &has_curved_faces) {

    has_curved_edges = false;
    has_curved_faces = false;
    auto L = Topo(Shp).faces();

    for (auto &face: L)
        if (!Topo(face).all_edges_non_curved()) {
            has_curved_edges = true;
            break;
        }

    for (auto &face: L)
        if (!Kernel::face_is_planar(TopoDS::Face(face))) {
            has_curved_faces = true;
            break;
        }
}

void ShapeHealing::remove_duplicated_hashes_from_shape() {

    unsigned int i = 0;
    while (i < 3) { // try to remove duplicated hashes. In case a surface is used by e.g. two TopoDS_Faces with different orientations, let's ignore it, because trsf of face will alwyas be applied to both topo faces
        i++;
        std::list<unsigned int> bad_hashes;

        if (find_duplicated_hashes_in_shape(bad_hashes))
            break;

        std::cerr << "[Warning] Duplicated hashes in shape (" << bad_hashes.size() << ")" << std::endl;

        std::list<std::pair<TopoDS_Shape, TopoDS_Shape>> old_and_new_subshapes;

        TopExp_Explorer Ex;
        for (Ex.Init(Shp, TopAbs_FACE); Ex.More(); Ex.Next()) {
            unsigned int H = Kernel::hash(Ex.Current());
            if (std::find(bad_hashes.begin(), bad_hashes.end(), H) != bad_hashes.end()) {
                gp_Trsf trsf;
                trsf.SetTranslation(gp_Vec(0, 0, 0));
                old_and_new_subshapes.emplace_back(Ex.Current(), Ex.Current().Moved(TopLoc_Location(trsf)));
            }
        }
        Shp = Kernel::replace_subshapes_in_shape(old_and_new_subshapes, Shp);
    }
}

bool ShapeHealing::find_duplicated_hashes_in_shape(std::list<unsigned int> &bad_hashes) {

    std::set<unsigned int> A;

    TopExp_Explorer Ex;
    for (Ex.Init(Shp, TopAbs_FACE); Ex.More(); Ex.Next()) {
        unsigned int H = Kernel::hash(Ex.Current());
        if (A.find(H) != A.end())
            bad_hashes.push_back(H);
        A.insert(H);
    }

    return bad_hashes.empty();
}

void ShapeHealing::fuse() {

    BOPAlgo_Builder B;
    TopExp_Explorer Ex;

    for (Ex.Init(Shp, TopAbs_FACE); Ex.More(); Ex.Next())
        B.AddArgument(Ex.Current());

    B.SetFuzzyValue(fuse_tol);
    B.Perform();

    if (B.HasWarnings()) {
        std::cerr << "Warnings:" << std::endl;
        B.DumpWarnings(std::cerr);
    }
    if (B.HasErrors()) {
        std::cerr << "Errors:" << std::endl;
        B.DumpErrors(std::cerr);
    }
    //if (!B.HasModified())
    //    std::cerr << "Nothing was modified." << std::endl;

    Shp = B.Shape();

}

void ShapeHealing::fill_face_list() {

    for (auto &f: Topo(Shp).faces())
        sFaces.emplace_back(TopoDS::Face(f));
}

bool ShapeHealing::find_duplicated_hashes(std::list<unsigned int> &bad_hashes) {

    std::set<unsigned int> unique_ids;
    for (auto &sface: sFaces) {
        if (unique_ids.find(sface.face_id()) != unique_ids.end()) {
            std::cerr << "[Error] Same ID twice in cFaces! " << sface.face_id() << std::endl;
            return false;
        }
        unique_ids.insert(sface.face_id());
    }

    return true;
}

void ShapeHealing::setup_sFace_map() {

    for (auto &sface: sFaces) {
        if (id2sFace.find(sface.face_id()) != id2sFace.end())
            std::cerr << "[Error] Two faces with same ID in cFaces list! " << sface.face_id() << std::endl;
        id2sFace[sface.face_id()] = &(sface);
    }
}

void ShapeHealing::identify_hanging_faces() { while (identify_hanging_faces_while()) std::cout << "#"; }

bool ShapeHealing::identify_hanging_faces_while() {

    bool found_hanging_faces = false;
    for (auto &sface: sFaces)
        if (sface.set_to_hanging())
            found_hanging_faces = true;

    return found_hanging_faces;
}

void ShapeHealing::remove_trash_and_face_adjacency() {

    std::list<sFace *> del;

    for (auto &sface: sFaces) if (sface.isTrash) del.push_back(&sface);

    for (auto &sface: sFaces) {
        if (sface.isTrash) continue;
        sface.remove_deleted_adjacent_faces(del);
    }

    for (auto &it: del) sFaces.remove(*it);
}

void ShapeHealing::find_components() {

    while (true) {

        sFace *start_sface = nullptr;
        find_start_face_for_component_search(start_sface);

        if (start_sface == nullptr) break;

        component_id_counter++;

        std::stack<sFace *> S;
        S.push(start_sface);

        while (!S.empty()) {

            sFace *sface = S.top();
            S.pop();

            if (sface->component_id != -1) continue;

            sface->component_id = component_id_counter;

            for (auto &e: sface->edgeId2adjacentsFaces)
                for (const auto &adjacent_sface: e.second)
                    S.push(adjacent_sface);
        }
    }
}

void ShapeHealing::find_start_face_for_component_search(sFace *&start_sface) {

    for (auto &sface: sFaces)
        if (sface.component_id == -1) {
            start_sface = &sface;
            return;
        }
}

TopoDS_ListOfShape ShapeHealing::components_to_solids() {

    TopoDS_ListOfShape Solids;

    for (int i = 0; i <= component_id_counter; i++) {

        bool good = false;
        TopoDS_Solid solid = get_outer_hull_of_component(i, good);
        if (good) Solids.Append(solid);
    }

    return Solids;
}

TopoDS_Solid ShapeHealing::get_outer_hull_of_component(int i, bool &isGood) {

    // collect component faces
    std::list<sFace *> component;
    for (auto &sface: sFaces)
        if (sface.component_id == i)
            component.push_back(&sface);

    // find outer face
    IntCurvesFace_ShapeIntersector Intersector;
    Intersector.Load(Shp, 0.0001);

    sFace *start_sface = nullptr;

    for (auto &sface: component) {
        if (isOuterFace(sface->face, sface->c, sface->n, Intersector)) {
            sface->isOuter = true;
            start_sface = sface;
            break;
        } else if (isOuterFace(sface->face, sface->c, sface->n.Reversed(), Intersector)) {
            sface->complement();
            sface->isOuter = true;
            start_sface = sface;
            break;
        }
    }

    if (start_sface == nullptr) {
        std::cerr << "[Warning] No outer face found! " << std::endl;
        return {};
    }

    // collect outer faces using adjacency
    std::stack<sFace *> S;
    S.push(start_sface);

    while (!S.empty()) {

        sFace *sface = S.top();
        S.pop();

        if (!sface->isFlagged)
            topoHull(sface, S);
    }

    // sew faces
    BRepBuilderAPI_Sewing sew;
    for (auto &sface: component)
        if (sface->isOuter)
            sew.Add(Kernel::remove_seam_edges(sface->face));

    sew.Perform();
    TopoDS_Shell Shell = TopoDS::Shell(Topo(sew.SewedShape()).shells().First());

    TopoDS_Solid Solid = BRepBuilderAPI_MakeSolid(Shell).Solid();
    if (Kernel::volume(Solid) < 0) Solid.Complement();

    // unify same domain
    ShapeUpgrade_UnifySameDomain U(Solid, true, true, true);
    U.SetAngularTolerance(Precision::Angular() * 2);
    U.SetLinearTolerance(Precision::Confusion() * 2);
    U.Build();
    Solid = TopoDS::Solid(U.Shape());

    isGood = true;
    return Solid;
}

bool ShapeHealing::isOuterFace(const TopoDS_Face &face, gp_Pnt P, const gp_Dir &n, IntCurvesFace_ShapeIntersector &Intersector) {

    gp_Vec t = gp_Vec(n).Scaled(0.00001);
    P.Translate(t); // move point a little bit in direction of normal, to prevent self intersection
    gp_Lin Line(P, n);

    Intersector.PerformNearest(Line, 0, RealLast());
    if (Intersector.NbPnt() > 0) {
        if (Kernel::hash(Intersector.Face(1)) == Kernel::hash(face)) std::cerr << "[Error] Self hit! " << std::endl;
        return false;
    } else
        return true;
}

void ShapeHealing::topoHull(sFace *sface, std::stack<sFace *> &stack) {

    sface->isFlagged = true;

    for (auto &e: sface->edgeId2adjacentsFaces) {

        unsigned int e_id = e.first;
        std::list<sFace *> adjacent_sfaces = e.second;

        if (adjacent_sfaces.empty())
            std::cerr << "[Error] No adjacent faces!" << std::endl;

        else if (adjacent_sfaces.size() == 1) { // if only one adjacent face -> done

            // if the only adj face is already outer face, go on
            if (adjacent_sfaces.front()->isOuter) continue;

            if (sface->edgeId2halfEdge[e_id].Orientation() == adjacent_sfaces.front()->edgeId2halfEdge[e_id].Orientation())
                adjacent_sfaces.front()->complement(); // Complement if wrong face normal

            stack.push(adjacent_sfaces.front());
            adjacent_sfaces.front()->isOuter = true;

        } else { //  find adjacent face with smallest angle to face

            sFace *correct_sface = nullptr;
            gp_Dir n1 = sface->n; // normal of cface
            gp_Dir v1 = sface->calculate_edge_vector(e_id);  // vector of half-edge
            gp_Pnt c = sface->calculate_mid_point(e_id);  // midpoint half-edge

            double minimum_angle = 400;

            for (const auto &adjacent_sface: adjacent_sfaces) {

                if (adjacent_sface->isTrash) continue;

                if (!adjacent_sface->is_id_in_edgeId2halfEdge(e_id))
                    std::cerr << "[Error] Edge id " << e_id << " not found in map" << std::endl;

                if (sface->edgeId2halfEdge[e_id].Orientation() == adjacent_sface->edgeId2halfEdge[e_id].Orientation())
                    if (!adjacent_sface->isOuter) // if sface and adjacent_sface are both outer faces, normals should already be correct
                        adjacent_sface->complement(); // complement if wrong face normal
                    else continue; // if both faces are outer faces and therefore have correct normals, then, if orientation is not matching, skip

                gp_Dir n2 = adjacent_sface->n;  // normal of cface
                gp_Dir v2 = adjacent_sface->calculate_edge_vector(e_id);  // vector of half-edge

                double angle = cFace::angle_between_planar_adjacent_faces(n1, v1, n2, v2, c);

                if (angle < minimum_angle) {
                    minimum_angle = angle;
                    correct_sface = adjacent_sface;
                }
            }

            if (correct_sface == nullptr)
                std::cerr << "[Warning] No unvisited adjacent face found" << std::endl;
            else {
                if (!adjacent_sfaces.front()->isOuter) {
                    stack.push(correct_sface);
                    correct_sface->isOuter = true;
                }
            }
        }
    }
}

void ShapeHealing::remove_bad_faces() {

    ShapeBuild_ReShape R;
    for (auto & F : Topo(Shp).faces()) {
        BRepCheck_Analyzer A(F);
        if (!A.IsValid()) {
            R.Remove(F);
            std::cerr << "\tRemove\n";
            //ShapeChecker(F).dump_topology();
        }
    }
    Shp = R.Apply(Shp);
}

void ShapeHealing::free_and_bad_edges() {

    ShapeAnalysis_Shell a;
    a.LoadShells(Shp);
    Standard_Boolean CheckOrientedShells = a.CheckOrientedShells(Shp, true, true);
    Standard_Boolean HasBadEdges = a.HasBadEdges();
    Standard_Boolean HasFreeEdges = a.HasFreeEdges();
    std::cout << "\tCheckOrientedShells: " << std::boolalpha << CheckOrientedShells << "\n";
    std::cout << "\tHasBadEdges: " << std::boolalpha << HasBadEdges << "\n";
    std::cout << "\tHasFreeEdges: " << std::boolalpha << HasFreeEdges << "\n";
    for (const auto &E: Topo(a.BadEdges()).edges())
        std::cerr << "Bad Edge:  " << Kernel::hash(E) << "\t" << E.Orientation() << std::endl;
    for (const auto &E: Topo(a.FreeEdges()).edges())
        std::cerr << "Free Edge: " << Kernel::hash(E) << "\t" << E.Orientation() << std::endl;

//    std::list<viewerHelper::DisplayShapes> ds2;
//    ds2.emplace_back();
//    ds2.back().shape = a.BadEdges();
//    ds2.back().transparency = 0.7;
//    ds2.back().clr_by_string = true;
//    ds2.back().clr_string = "RED";
//
//    ds2.emplace_back();
//    ds2.back().shape = a.FreeEdges();
//    ds2.back().transparency = 0.7;
//    ds2.back().clr_by_string = true;
//    ds2.back().clr_string = "GREEN";
//
//    ViewerMain viewer2;
//    viewer2.start_viewer(ds2);


    TopExp::MapShapesAndAncestors(Shp, TopAbs_EDGE, TopAbs_FACE, edgeFaceMap);

    fill_face_list();

    if (sFaces.empty())
        return;

    std::list<unsigned int> bad_hashes;
    if (!find_duplicated_hashes(bad_hashes))
        return;

    for (auto &sface: sFaces)
        sface.update_half_edges();

    setup_sFace_map();

    for (auto &sface: sFaces)
        sface.update_face_adjacencies(edgeFaceMap, id2sFace);

    TopoDS_ListOfShape FreeEdges;
    TopoDS_ListOfShape BadEdges;
    std::set<sFace *> BadFaces;

    for (auto &sface: sFaces) {
        bool is_joined = true;

        for (auto &it: sface.edgeId2adjacentsFaces) {

            // check existence of adjacent faces on edge -> FREE Edges
            if (it.second.empty()) {
                is_joined = false;
                FreeEdges.Append(sface.edgeId2halfEdge[it.first]);
                break;
            }

            // check if orientation is same -> BAD Edges
            for (const auto &adjacent_sface: it.second)
                for (auto &adj: adjacent_sface->edgeId2halfEdge)
                    if (adj.first == it.first) {
                        if (adj.second.Orientation() == sface.edgeId2halfEdge[it.first].Orientation()) {
                            BadEdges.Append(sface.edgeId2halfEdge[it.first]);
                            BadFaces.insert(&sface);
                        }
                        break;
                    }
        }

        if (!is_joined)
            sface.isTrash = true;
    }

    for (const auto &E: BadEdges)
        std::cerr << "Bad Edge2:  " << Kernel::hash(E) << "\t" << E.Orientation() << std::endl;
    for (const auto &E: FreeEdges)
        std::cerr << "Free Edge2: " << Kernel::hash(E) << "\t" << E.Orientation() << std::endl;
    for (const auto &F: BadFaces)
        std::cerr << "Bad Face2: " << Kernel::hash(F->face) << "\t" << Kernel::area(F->face) << std::endl;

#ifdef VISUALIZATION
    std::list<viewerHelper::DisplayShapes> ds;
    for (const auto &sface: sFaces) {
        ds.emplace_back();
        ds.back().shape = sface.face;
        ds.back().transparency = 0.8;
        ds.back().clr_by_string = true;
//        ds.back().show_adaptor_face_normal = true;
//        ds.back().adaptor_face_normal = Kernel::face_normal(TopoDS::Face(ds.back().shape));
//        ds.back().center_of_face = Kernel::face_center(TopoDS::Face(ds.back().shape));

        if (sface.isTrash) {
            ds.back().clr_string = "RED";
            ds.back().transparency = 0.2;
            //std::cerr << "RED" << std::endl;
        } else ds.back().clr_string = "GREEN";

    }
    for (const auto &E: FreeEdges) {
        ds.emplace_back();
        ds.back().shape = E;
        ds.back().transparency = 0.0;
        ds.back().clr_by_string = true;
        ds.back().clr_string = "BLUE";
    }

    for (const auto &E: BadEdges) {
        ds.emplace_back();
        ds.back().shape = E;
        ds.back().transparency = 0.0;
        ds.back().clr_by_string = true;
        ds.back().clr_string = "YELLOW";
    }

    for (const auto &E: BadFaces) {
        ds.emplace_back();
        ds.back().shape = E->face;
        ds.back().transparency = 0.0;
        ds.back().clr_by_string = true;
        ds.back().clr_string = "ORANGE";
        ds.back().show_adaptor_face_normal = true;
        ds.back().adaptor_face_normal = Kernel::face_normal(TopoDS::Face(ds.back().shape));
        ds.back().center_of_face = Kernel::face_center(TopoDS::Face(ds.back().shape));
    }

    ViewerMain::start_viewer(ds);
#endif
}