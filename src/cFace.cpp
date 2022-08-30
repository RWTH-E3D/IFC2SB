// Copyright 2022 Eric Fichter
#include "cFace.h"

cFace::cFace(TopoDS_Face _face, oFace *_face_original, unsigned int _id, cFace *_superface) : face(_face), ancestor(_face_original), id(_id), superface(_superface) {

    halfedges.clear();
    adjacentfaces.clear();
    isTrash = false;
    isHanging = false;
    wasVisited = false;
    isInner = false;
    isCoplanar = false;
    space = nullptr;
    sb_type = SB_TYPE_2B;
    internalOrExternal = SB_IE_NOTDEFINED;
    physicalOrVirtual = SB_PV_NOTDEFINED;
    corresponding = nullptr;
    space_behind = nullptr;
    materials.clear();
    unifying_face = nullptr;
    parent = nullptr;
    subfaces.clear();
    ifcrelspaceboundary = nullptr;
    P = gp_Pnt();
    fixed_normal = gp_Dir();
    status_normal = FACE_NORMAL_UNKNOWN;
    distance = -1;
    CheckSetFaceNormal();
}

void cFace::CheckSetFaceNormal() {

    bool normal_was_checked;
    bool normal_is_ok = CheckFaceNormal(normal_was_checked);

    if (normal_was_checked) {
        if (!normal_is_ok) ComplementFace();
        SetNormalStatus(FACE_NORMAL_KNOWN);
    } else
        std::cout << "[Info] Correct face normal is unknown. " << Info() << "\n";
}

void cFace::SetNormalStatus(face_normal_status s) {
    status_normal = s;
    if (s == FACE_NORMAL_KNOWN)
        fixed_normal = FaceNormal();
}

bool cFace::check_face_normal(gp_Dir n, gp_Dir m, double AngularTolerance, const std::string &s) const {

    if (!m.IsEqual(n, AngularTolerance)) {
        double a = m.Angle(n);
        if (fabs(a - M_PI) < 0.01)
            std::cerr << "[Warning] Face was complemented compared to " << s << " " << Info() << "\t" << a << " | "
                      << Kernel::round_double_two_digits(n.X()) << " "
                      << Kernel::round_double_two_digits(n.Y()) << " "
                      << Kernel::round_double_two_digits(n.Z()) << " | "
                      << Kernel::round_double_two_digits(m.X()) << " "
                      << Kernel::round_double_two_digits(m.Y()) << " "
                      << Kernel::round_double_two_digits(m.Z()) << std::endl;
        return false;
    } else
        return true;
}

bool cFace::CheckFaceNormal(bool &normal_was_checked) const {

    normal_was_checked = false;
    double AngularTolerance = 1.57; // in radians (ca. 90°), // especially slim faces can be rotated while fusing, so we need a big tolerance
    gp_Dir nf = FaceNormal();

    if (NormalStatus() == FACE_NORMAL_KNOWN) {
        normal_was_checked = true;
        return check_face_normal(FaceNormal(), fixed_normal, AngularTolerance, "fixed normal");
    }

    if (superface != nullptr) {
        normal_was_checked = true;
        return check_face_normal(FaceNormal(), superface->FaceNormal(), AngularTolerance, "superface normal");
    }

    if (ancestor == nullptr) {
        std::cout << "[Info] Face has no ancestor. Skip test " << Info() << "\n";
        return false;
    }

    if (ancestor->NormalStatus() == FACE_NORMAL_KNOWN) {
        normal_was_checked = true;
        return check_face_normal(FaceNormal(), AncestorNormal(), AngularTolerance, "ifc face normal");
    }

    return false;
}

bool cFace::AngleToFaceSmaller90(const cFace &f) const {

    double AngularTolerance = 1.57;
    gp_Dir nf = FaceNormal();
    gp_Dir no = Kernel::face_normal(f.face);
    return nf.IsEqual(no, AngularTolerance);
}

void cFace::ComplementFace() { face.Complement(); }

unsigned int cFace::FaceID() const { return Kernel::hash(face); }

void cFace::UpdateHalfEdges() {

    halfedges.clear();

    TopExp_Explorer Ex;

    for (Ex.Init(face, TopAbs_EDGE); Ex.More(); Ex.Next()) { // skip seam edges
        if (Ex.Current().Orientation() == TopAbs_INTERNAL || Ex.Current().Orientation() == TopAbs_EXTERNAL)
            continue;
        halfedges[Kernel::hash(Ex.Current())] = Ex.Current();
    }
}

bool cFace::IsIdInHalfEdges(unsigned int theId) { return !(halfedges.find(theId) == halfedges.end()); }

double cFace::DistanceToSpaceBoundaryBehind() const { return distance; }

gp_Pnt cFace::MidPoint(unsigned int edgeID) {
    // computes the point that lies in the middle of P1 and P2
    TopoDS_Edge Edge = TopoDS::Edge(halfedges[edgeID]);
    gp_Pnt P1 = BRep_Tool::Pnt(TopExp::FirstVertex(Edge, Standard_True));
    gp_Pnt P2 = BRep_Tool::Pnt(TopExp::LastVertex(Edge, Standard_True));
    gp_Vec vec = 0.5 * (gp_Vec(P1.XYZ()) + gp_Vec(P2.XYZ()));
    return {vec.XYZ()};
}

gp_Dir cFace::vector_on_plane(gp_Dir n, gp_Dir v) { return n.Crossed(v); }

double cFace::degree_angle_from_unit_arc_edge(const TopoDS_Edge &Edge) {
    // length of edge, that is based on an arc on a unit circle is equal to the angle in radians
    GProp_GProps props;
    BRepGProp::LinearProperties(Edge, props); // this is angle in radians
    return props.Mass() * 180.0 / M_PI; // radian to degree
}

gp_Pnt cFace::point_on_plane(gp_Dir r, gp_Pnt c) {
    // point on face in distance of 1 to mid point of edge c
    gp_Vec v1(r);
    gp_Vec v2(c.XYZ());
    gp_Vec p = v2.Added(v1);
    return {p.XYZ()};
}

gp_Pnt cFace::Center() const { return Kernel::face_center(face); }

gp_Dir cFace::EdgeVector(unsigned int edgeID) {

    TopoDS_Edge Edge = TopoDS::Edge(halfedges[edgeID]);
    gp_Pnt P1 = BRep_Tool::Pnt(TopExp::FirstVertex(Edge, Standard_True));
    gp_Pnt P2 = BRep_Tool::Pnt(TopExp::LastVertex(Edge, Standard_True));
    return {gp_Vec(P1, P2)};
}

double cFace::angle_between_planar_adjacent_faces(gp_Dir n1, gp_Dir v1, gp_Dir n2, gp_Dir v2, gp_Pnt c) {

    /*
     * gp_Dir n1 - normal unit vector of face 1
     * gp_Dir v1 - half edge vector of face 1
     * gp_Dir n2 - normal unit vector of face 2
     * gp_Dir v2 - half edge vector of face 2
     * gp_Pnt c - center of edge (half edge)
     */

    gp_Dir r1 = vector_on_plane(n1, v1); // unit vector ON face 1, starting from mid of edge (90° to normal and edge vector)
    gp_Dir r2 = vector_on_plane(n2, v2); // unit vector ON face 2, starting from mid of edge (90° to normal and edge vector)
    gp_Pnt b1 = point_on_plane(r1, c);
    gp_Pnt b2 = point_on_plane(r2, c);

    if (r1.Angle(r2) < 1.0e-9) { std::cerr << "[Error] vector on planes are the same!" << std::endl; }

    auto arc = GC_MakeArcOfCircle(b1, gp_Vec(n1.XYZ()), b2); // erzeugt einen kreisbogen, genauer einen einheitskreisbogen (R=1), da b1 und b2 einen senkrechten abstand von 1 zum mittelpunkt c haben

    if (!arc.IsDone()) {
        std::cerr << "[Error] Failed creating Circle!" << std::endl;
        return 360;
    }

    TopoDS_Edge E = BRepBuilderAPI_MakeEdge(arc.Value()).Edge();
    return degree_angle_from_unit_arc_edge(E);  // length of arc of unit circle = angle

}

bool cFace::IsEnclosed(const TopTools_IndexedDataMapOfShapeListOfShape &edgeFaceMap, std::unordered_map<unsigned int, std::list<cFace *>> &id2cFaces, std::list<cFace> &cFaces) {
    return IsOffset() ? is_enclosed_offset(edgeFaceMap, id2cFaces, cFaces) : is_enclosed_shell(edgeFaceMap, id2cFaces, cFaces);
}

bool cFace::is_enclosed_offset(const TopTools_IndexedDataMapOfShapeListOfShape &edgeFaceMap, std::unordered_map<unsigned int, std::list<cFace *>> &id2cFaces, std::list<cFace> &cFaces) {

//     Offset faces are enclosed in the following cases:
//     1) Consider all adjacent faces, that are not coplanar (same face id), are not offset faces and are not connected via a seam edge. In those adjacent faces there should be two faces from the same product.
//        There is also the case that there is only one face from a product, thats totally ok. More than two is an error, because then the adjacent product consists of two intersecting shells
//        The offset face is enclosed, if the angle between the two adjacent faces is bigger than the angle between the face and that adjacent face that is connected on a "good oriented" half-edge.
//     2) If the offset face cuts a seam edge in a non-offset face, the offset face is enclosed if its on the wrong side (negative side of halfspace formed by adjacent face and its normal)

    for (auto &it: halfedges) {

        unsigned int edge_id = it.first;
        std::list<cFace *> real_adjacent_cfaces;

        // identify adjacent cfaces that fulfill the above stated conditions. If there is a seam edge combo, check 2)
        for (auto &adjacent_face: edgeFaceMap.FindFromKey(it.second)) {
            if (Kernel::hash(adjacent_face) == FaceID()) continue; // skip self reference and coplanar faces

            for (auto &possible_cface: id2cFaces[Kernel::hash(adjacent_face)])
                if (possible_cface->IsIdInHalfEdges(edge_id) && !possible_cface->IsOffset()) // skip if connection is using seam half-edge or adjacent face is offset
                    real_adjacent_cfaces.emplace_back(possible_cface);

                else if (!possible_cface->IsIdInHalfEdges(edge_id))
                    if (is_face_on_wrong_side(possible_cface)) return true; // decision 2)
        }

        // find the pairs of faces from same product
        std::unordered_map<std::basic_string<char>, std::vector<cFace *>> pairs = find_pairs_by_guid(real_adjacent_cfaces);

        // for each pair and current face, check angle condition
        for (auto &pair: pairs) {
            if (pair.second.size() > 2) {
                std::cerr << "[Error] More than two adjacent faces from same product! " << pair.second.size() << std::endl;
                continue;
            } else if (pair.second.size() < 2)
                continue;

            if (is_enclosed_by_pair_faces(pair.second[0], pair.second[1], edge_id)) return true; // decision 1)
        }
    }

    return false;
}

bool cFace::is_enclosed_shell(const TopTools_IndexedDataMapOfShapeListOfShape &edgeFaceMap, std::unordered_map<unsigned int, std::list<cFace *>> &id2cFaces, std::list<cFace> &cFaces) {


//    Shell faces are enclosed in the following cases:
//    1) Consider all adjacent faces, that are not from the same product, are not coplanar (same face id), are not offset faces and are not connected via a seam edge (last one should be possible, so throw error). In those adjacent faces there should be two faces from the same product.
//       There is also the case that there is only one face from a product, thats totally ok (if shell face is sitting coplanar on next to adjacent product). More than two is an error, because then the adjacent product consists of two intersecting shells
//       The offset face is enclosed, if the angle between the two adjacent faces is bigger than the angle between the face and that adjacent face that is connected on a "good oriented" half-edge.

    for (auto &it: halfedges) {

        unsigned int edge_id = it.first;
        std::list<cFace *> real_adjacent_cfaces;

        for (auto &adjacent_face: edgeFaceMap.FindFromKey(it.second)) {
            if (Kernel::hash(adjacent_face) == FaceID()) continue; // skip self reference and coplanar faces

            for (auto &possible_cface: id2cFaces[Kernel::hash(adjacent_face)])
                if (possible_cface->IsIdInHalfEdges(edge_id) && !possible_cface->IsOffset()) // skip if connection is using seam half-edge, adjacent face is offset or (from the same product and same shell)
                    if (RelProduct() != possible_cface->RelProduct())
                        real_adjacent_cfaces.emplace_back(possible_cface);
                    else if (ancestor->ShellID() != possible_cface->ancestor->ShellID()) // if its from the same product it must be from a different shell
                        real_adjacent_cfaces.emplace_back(possible_cface);
        }

        // find the pairs of faces from same product
        std::unordered_map<std::basic_string<char>, std::vector<cFace *>> pairs = find_pairs_by_guid(real_adjacent_cfaces);

        for (auto &pair: pairs) {
            if (pair.second.size() > 2) {
                std::cerr << "[Error] More than two adjacent faces from same product! " << pair.second.size() << std::endl;
                continue;
            } else if (pair.second.size() < 2)
                continue;

            if (is_enclosed_by_pair_faces(pair.second[0], pair.second[1], edge_id)) return true; // decision 1)
        }
    }

    return false;
}

bool cFace::is_face_on_wrong_side(cFace *other_cface) const {

    gp_Pnt c_adj = Kernel::face_center(other_cface->face);
    gp_Dir n_adj = other_cface->FixedFaceNormal();
    gp_Pnt c_off = Kernel::face_center(face);
    gp_Vec l_off(c_adj, c_off);

    return n_adj.Dot(l_off) < 0;
}

bool cFace::is_enclosed_by_pair_faces(cFace *productFace1, cFace *productFace2, unsigned int edge_id) {

    // angle between two pair faces
    gp_Dir n1 = productFace1->FixedFaceNormal();
    gp_Dir n2 = productFace2->FixedFaceNormal();
    gp_Dir v1 = productFace1->EdgeVector(edge_id);
    gp_Dir v2 = productFace2->EdgeVector(edge_id);
    gp_Pnt c = productFace1->MidPoint(edge_id);
    double angle_p = angle_between_planar_adjacent_faces(n1, v1, n2, v2, c);

    // find face to compare
    cFace *vgl = (halfedges[edge_id].Orientation() != productFace2->halfedges[edge_id].Orientation()) ? productFace1 : productFace2;

    // angle between face and "correct" pair face
    n1 = vgl->FixedFaceNormal();
    n2 = FixedFaceNormal();
    v1 = vgl->EdgeVector(edge_id);
    v2 = EdgeVector(edge_id);
    c = vgl->MidPoint(edge_id);
    double angle_v = angle_between_planar_adjacent_faces(n1, v1, n2, v2, c);

    return angle_v > angle_p;
}

std::unordered_map<std::basic_string<char>, std::vector<cFace *>> cFace::find_pairs_by_guid(const std::list<cFace *> &adjacent_cfaces) {

    std::unordered_map<std::basic_string<char>, std::vector<cFace *>> pairs;

    for (auto &real_adjacent_cface: adjacent_cfaces) {

        // In case a product consists of multiple shells, it is possible that shells are touching in one face (e.g. slab in nemetschek 2x3 fzk haus).
        // Without considering the specific shell, there would be 3 or 4 adjacent faces pushbacked into the pair vector with counter oriented normals.
        // The pair check shall only match faces with consistent orientation, so split pairs by considering shell id.
        auto guid = real_adjacent_cface->RelProduct()->guid + "_" + std::to_string(real_adjacent_cface->ancestor->ShellID());

        if (pairs.find(guid) != pairs.end())
            pairs[guid].emplace_back(real_adjacent_cface);
        else
            pairs[guid] = {real_adjacent_cface};
    }
    return pairs;
}

void cFace::UpdateFaceAdjacencies(const TopTools_IndexedDataMapOfShapeListOfShape &edgeFaceMap, std::unordered_map<unsigned int, cFace *> &id2cFace) {

    // for every half-edge, get adjacent faces
    for (auto &itr: halfedges) {
        TopTools_ListOfShape adjacent_faces = edgeFaceMap.FindFromKey(itr.second);
        adjacent_faces.Remove(face); // remove self-reference by topods_face (an opposite oriented face will not removed)
        unsigned int edge_id = itr.first;

        // get cFace of the adjacent faces
        std::list<cFace *> adjacent_cfaces;

        for (auto &adjacent_face: adjacent_faces) {

            unsigned int adjHash = Kernel::hash(adjacent_face);

            // only add if faces are in id2cFace. if not, they are in fuse shape, but not in cfaces anymore (e.g. because of deletion of duplicate faces)
            if (id2cFace.find(adjHash) == id2cFace.end()) continue;

            cFace *adjacent_cface = id2cFace[adjHash];

            // remove self-reference by id to filter also opposite oriented faces that were linked to this cface by id2cFace
            if (id == adjacent_cface->id) continue;

            // don't add adjacent faces, if they are connected by a seam edge - seam edge or seam edge - normal edge combo. So look if edge id is in maps of the "adjacent face" first
            if (adjacent_cface->IsIdInHalfEdges(edge_id))
                adjacent_cfaces.push_back(adjacent_cface);
        }

        adjacentfaces[edge_id] = adjacent_cfaces;
    }
}

void cFace::RemoveWrongOrientedAdjacentFaces() {

    for (auto &it: adjacentfaces) {
        std::list<cFace *> del; // adjacent faces to delete
        unsigned int edge_id = it.first;

        for (auto &adjacent_cface: it.second) {

            if (!adjacent_cface->is_edge_id_in_maps(edge_id)) {
                std::cerr << "ERROR: Edge id " << edge_id << " not found in map" << std::endl;
                continue;
            }

            auto twin = adjacent_cface->halfedges[edge_id];  // get adjacent half-edge of half-edge

            if (twin.Orientation() == halfedges[edge_id].Orientation())
                del.push_back(adjacent_cface);
        }

        for (auto &it_del: del)
            adjacentfaces[edge_id].remove(it_del);

    }

}

bool cFace::is_edge_id_in_maps(unsigned int ID) { return !((halfedges.find(ID) == halfedges.end()) || adjacentfaces.find(ID) == adjacentfaces.end()); }

bool cFace::SetHanging() {

    if (isHanging)
        return false;

    bool is_joined = true; // indicates if face has at least one adjacent SHELL face on all half-edges

    for (auto &it: adjacentfaces) {

        // check for existence of adjacent faces on edge
        if (it.second.empty()) {
            is_joined = false;
            break;
        }

        // there are faces on edge, so check if there is at least a non-hanging adjacent SHELL face
        bool has_non_hanging_adjacent_faces = false;
        for (const auto &adjacent_cface: it.second) {
            if (!adjacent_cface->isHanging && !adjacent_cface->IsOffset()) {
                has_non_hanging_adjacent_faces = true;
                break;
            }
        }

        if (!has_non_hanging_adjacent_faces) {
            is_joined = false;
            break;
        }
    }

    if (!is_joined) {
        isHanging = true;
        isTrash = true;
        return true;
    } else
        return false;

}

void cFace::RemoveTrashedAdjacentFaces() {

    for (auto &e: adjacentfaces)
        e.second.remove_if([](cFace *cface) { return cface->isTrash; });
}

bool cFace::FaceContainsSeamEdges() const {
    for (auto &e: Topo(face).edges())
        if (e.Orientation() == TopAbs_INTERNAL || e.Orientation() == TopAbs_EXTERNAL)
            return true;
    return false;
}

void cFace::RemoveSeamEdgesFromFace() {

    if (!FaceContainsSeamEdges()) return;

    TopExp_Explorer Ex;
    ShapeBuild_ReShape R;

    for (auto &e: Topo(face).edges())
        if (e.Orientation() == TopAbs_INTERNAL || e.Orientation() == TopAbs_EXTERNAL)
            R.Remove(e);

    face = TopoDS::Face(R.Apply(face));

    // seam edges often belong to a wire with only that edge, so remove the wire, too
    R.Clear();

    for (auto &w: Topo(face).wires())
        if (Topo(w).edges().IsEmpty())
            R.Remove(w);

    face = TopoDS::Face(R.Apply(face));
}

std::set<cFace *> cFace::UnifyableFacesMaintainHoles() {

    // get all adjacent faces (no only next neighbour) of the face on all edges fulfilling the conditions
    // starting from "this" face, find all faces from 1) the same original face with 2) the same sb type
    // doesn't add "this" face to list.

    std::set<cFace *> s;

    std::stack<cFace *> S;
    S.push(this);

    while (!S.empty()) {

        cFace *pF = S.top();
        S.pop();

        for (auto &adjacent_cfaces: pF->adjacentfaces)
            for (auto &adjacent_cface: adjacent_cfaces.second)
                if (!adjacent_cface->isTrash && adjacent_cface != this && adjacent_cface->ancestor == ancestor && adjacent_cface->sb_type == sb_type) { // accepts faces not being trash and having same sb type and original face

                    if (!corresponding_matching_superface(adjacent_cface)) continue;

                    if (s.find(adjacent_cface) == s.end())
                        S.push(adjacent_cface);

                    s.insert(adjacent_cface);
                }
    }

    return s;
}

std::set<cFace *> cFace::UnifyableFacesOvertakeOpenings() {

    // get all adjacent faces (not only next neighbour) of the face on all edges fulfilling the conditions
    // starting from "this" face, find all faces from 1) the same original face with 2) the same sb type 3) same space behind 4) same material path 5) ...
    // doesn't add "this" face to list.

    std::set<cFace *> s; // success
    std::set<cFace *> f; // fail

    std::stack<cFace *> S;
    S.push(this);

    while (!S.empty()) {

        cFace *pF = S.top();
        S.pop();

        for (auto &adjacent_cfaces: pF->adjacentfaces)
            for (auto &a: adjacent_cfaces.second) {

                if (s.find(a) != s.end() || f.find(a) != f.end()) continue; // skip if already checked

                if (unifying_criteria(a)) {
                    if (s.find(a) == s.end()) S.push(a);
                    s.insert(a);
                } else
                    f.insert(a);
            }
    }

    return s;
}

bool cFace::unifying_criteria(const cFace *a) {
    if (a->isTrash) return false;
    if (a == this) return false;
    if (a->space != space) return false;
    if (a->sb_type != sb_type) return false;
    if (a->physicalOrVirtual != physicalOrVirtual) return false;
    if (a->internalOrExternal != internalOrExternal) return false;
    if (a->space_behind != space_behind) return false;
    //if (Kernel::round_double_three_digits(a->distance) != Kernel::round_double_three_digits(distance)) return false; // already taken into account when comparing material maps
    //if (!map_compare_keys(a->materials, materials)) return false;
    if (!map_compare(a->materials, materials)) return false;
    // only faces from walls that were cut by opening have a parent_id different from 0. parent_id is hash of ancestor.
    if (!(a->ancestor == ancestor || (a->ancestor->ParentID() == ancestor->ParentID() && ancestor->ParentID() != 0))) return false;
    if (!corresponding_matching_superface(a)) return false;
    return true;
}

bool cFace::unifying_criteria_openings(const cFace *a) {
    if (!a->IsOpening()) return false;
    if (a->isTrash) return false;
    if (a == this) return false;
    if (a->space != space) return false;
    if (a->sb_type != sb_type) return false;
    if (a->physicalOrVirtual != physicalOrVirtual) return false;
    if (a->internalOrExternal != internalOrExternal) return false;
    if (a->space_behind != space_behind) return false;
    //if (Kernel::round_double_three_digits(a->distance) != Kernel::round_double_three_digits(distance)) return false; // already taken into account when comparing material maps
    // if (!map_compare_keys(a->materials, materials)) return false;
    if (!map_compare(a->materials, materials)) return false;
    if (a->ancestor != ancestor) return false;
    if (!corresponding_matching_superface(a)) return false;
    return true;
}

std::set<cFace *> cFace::UnifyableFacesOpeningsOnly() {

    // get all adjacent faces (no only next neighbour) of the face on all edges fulfilling the conditions
    // starting from "this" face, find all faces from 1) the same original face with 2) the same sb type
    // doesn't add "this" face to list.

    std::set<cFace *> s;
    std::set<cFace *> f;

    std::stack<cFace *> S;
    S.push(this);

    while (!S.empty()) {

        cFace *pF = S.top();
        S.pop();

        for (auto &adjacent_cfaces: pF->adjacentfaces)
            for (auto &a: adjacent_cfaces.second) {

                if (s.find(a) != s.end() || f.find(a) != f.end()) continue; // skip if already checked

                if (unifying_criteria_openings(a)) {
                    if (s.find(a) == s.end()) S.push(a);
                    s.insert(a);
                } else
                    f.insert(a);
            }
    }

    return s;
}

bool cFace::corresponding_matching_superface(const cFace *a) {
    // if there is a corresponding face, only unify with those adjacent faces, that also have a corresponding face, that belongs to the same superface (1st lvl)
    if (corresponding == nullptr) return true;
    if (a->corresponding == nullptr) return false;
    if (corresponding->superface != a->corresponding->superface) return false;
    return true;
}

void cFace::UnifyFaces(std::set<cFace *> &faces_to_unify) {

    if (faces_to_unify.empty()) return;

    // make compound from face and it's adjacent faces
    BRep_Builder B;
    TopoDS_Compound C;
    B.MakeCompound(C);
    B.Add(C, face);
    for (auto &adj: faces_to_unify)
        B.Add(C, adj->face);

    // unify faces in the compound
    ShapeUpgrade_UnifySameDomain U(C, true, true, true);
    U.SetAngularTolerance(1.0e-5); // Precision::Angular() * 2. Value was raised because of parallel fusing to compensate moved faces.
    U.SetLinearTolerance(Precision::Confusion() * 2);
    U.Build();

    TopoDS_ListOfShape u = Topo(U.Shape()).faces();

    if (u.Size() == 1) { // if unification was successfully, unification shape should contain only one face
        face = TopoDS::Face(u.First());
        for (auto &c: faces_to_unify) {
            c->SetIsTrash(true); // set adjacent faces to trash, so no multiple unifications are done
            c->unifying_face = this;
        }
    } else {
        std::cerr << "[Warning] Unifying failed. " << Topo(C).faces().Size() << "\t" << Info() << "\t" << IsOpening() << "\t" << (IsOpening() ? Ancestor()->Opening()->data().getArgument(0)->toString() : "null") << std::endl;
        for (auto &c: faces_to_unify)
            std::cerr << "\t" << c->Info() << "\t" << c->IsOpening() << "\t" << (c->IsOpening() ? c->Ancestor()->Opening()->data().getArgument(0)->toString() : "null") << std::endl;
    }
}

std::string cFace::STLSolidName() {
    return (std::to_string(space->id) + "_" + std::to_string(id) + "_" + std::to_string(sb_type) + "_" + IfcClass() + "_" + RelProduct()->guid);
}

void cFace::ClearMaps() {
    halfedges.clear();
    adjacentfaces.clear();
}

void cFace::ReduceEdges() {

    ShapeUpgrade_UnifySameDomain U(face, true, false, true);
    U.SetAngularTolerance(Precision::Angular() * 2);
    U.SetLinearTolerance(Precision::Confusion() * 2);
    U.Build();

    TopoDS_ListOfShape u = Topo(U.Shape()).faces();

    if (u.Size() == 1) // if unification was successfully, unification shape should contain only one face
        face = TopoDS::Face(U.Shape());
    else
        std::cerr << "[Warning] Unifying failed. " << Info() << std::endl;

}

template<typename Map>
bool cFace::map_compare_keys(Map const &lhs, Map const &rhs) {
    auto pred = [](decltype(*lhs.begin()) a, decltype(a) b) { return a.first == b.first; };
    return lhs.size() == rhs.size() && std::equal(lhs.begin(), lhs.end(), rhs.begin(), pred);
}

template<typename Map>
bool cFace::map_compare(Map const &lhs, Map const &rhs) { return lhs == rhs; } //lhs.size() == rhs.size() && std::equal(lhs.begin(), lhs.end(), rhs.begin()); } // materials double values are not rounded. "round" before comparison

void cFace::RemoveNonManifoldAdjacencyByAngle() {

    // only considers adjacent faces with correct opposite edge orientation!

    for (auto &it: adjacentfaces) {

        if (it.second.size() < 2) continue; // already only one adjacent face at all

        unsigned int edge_id = it.first;
        std::set<cFace *> L;
        for (const auto &A: it.second) {

            if (!A->is_edge_id_in_maps(edge_id)) {
                std::cerr << "[Error]: Edge id " << edge_id << " not found in map" << std::endl;
                continue;
            }
            if (A->halfedges[edge_id].Orientation() == halfedges[edge_id].Orientation())
                continue;

            L.insert(A);
        }

        if (L.size() < 2) continue; // already only one adjacent face of considered edge orientation at all

        // find adjacent face with smallest angle to face
        gp_Dir n1 = FaceNormal();
        gp_Dir v1 = EdgeVector(edge_id);
        gp_Pnt c = MidPoint(edge_id);
        cFace *correct_cface = face_with_minimal_angle_correct_oriented_edge(edge_id, n1, v1, c, L);

        if (correct_cface == nullptr) {
            std::cerr << "[Error] No fitting adjacent face found" << std::endl;
            continue;
        }

        for (auto &A: L)
            if (A != correct_cface)
                adjacentfaces[edge_id].remove(A);
    }
}

cFace *cFace::face_with_minimal_angle_correct_oriented_edge(unsigned int edge_id, gp_Dir n1, gp_Dir v1, gp_Pnt c, const std::set<cFace *> &L) {

    // takes only faces with correct opposite edge, at this point
    // otherwise do check and complement face in this function

    cFace *correct_cface = nullptr;
    double minimum_angle = 400;

    for (const auto &A: L) {
        gp_Dir n2 = A->FaceNormal();
        gp_Dir v2 = A->EdgeVector(edge_id);
        double angle = cFace::angle_between_planar_adjacent_faces(n1, v1, n2, v2, c);

        if (angle < minimum_angle) {
            minimum_angle = angle;
            correct_cface = A;
        }
    }
    return correct_cface;
}

bool cFace::IsConnectedToShell() {

    for (auto &it: adjacentfaces)
        for (auto &a: it.second)
            if (!a->IsOffset() && RelProduct() == a->RelProduct() && ancestor->ShellID() == a->ancestor->ShellID())
                return true;

    return false;
}

double cFace::SurfaceArea() const { return Kernel::area(face); }

void cFace::SetPropertiesVirtualElement() {

    internalOrExternal = SB_IE_INTERNAL;
    physicalOrVirtual = SB_PV_VIRTUAL;
    sb_type = SB_TYPE_2A;
}

void cFace::SetFixedFaceNormal(gp_Dir v) { fixed_normal = v; }

void cFace::SetPropertiesSpaceBoundary(cFace *behind_cface, bool behind_face_is_corresponding, double min_dist, double bounds_min_Z, gp_Pnt C) {

    distance = min_dist;

    if (space->is_facade) {
        internalOrExternal = SB_IE_EXTERNAL;
        //if (C.Z() - 0.05 < bounds_min_Z) internalOrExternal = SB_IE_EXTERNAL_EARTH;
        if (internalOrExternal == SB_IE_EXTERNAL_EARTH || C.Z() - 0.05 < bounds_min_Z && FixedFaceNormal().Z() < -0.95) internalOrExternal = SB_IE_EXTERNAL_EARTH;
        if (superface != nullptr) if (superface->internalOrExternal == SB_IE_EXTERNAL_EARTH) internalOrExternal = SB_IE_EXTERNAL_EARTH;
    } else {
        if (behind_cface == nullptr) internalOrExternal = SB_IE_INTERNAL; // 2b faces on the inside of an outer wall are declared as INTERNAL, too
        else internalOrExternal = behind_cface->space->is_facade ? SB_IE_EXTERNAL : SB_IE_INTERNAL;
    }

    physicalOrVirtual = SB_PV_PHYSICAL;
    sb_type = SB_TYPE_2B;
    materials[RelProduct()].insert(0);

    if (behind_cface != nullptr) {

        space_behind = behind_cface->space;
        materials[behind_cface->RelProduct()].insert(Kernel::round_double_to_n_decimal_places(min_dist, 5));

        if (behind_face_is_corresponding) {
            corresponding = behind_cface;
            sb_type = SB_TYPE_2A;
        }
    }

    // ifcopening are used as normal obejcts to get virtual sbs (e.g. fzk haus or kit institute)
    if (RelProduct()->ifcproduct->declaration().is("IfcOpeningElement"))
        physicalOrVirtual = SB_PV_VIRTUAL;
}

bool cFace::CheckParent(double tol_dist, double tol_angle) {

    if (parent == nullptr) return true;

    if (this == parent) {
        std::cerr << "[Warning] Parent face self reference. " << Info() << std::endl;
        return false;
    }

    gp_Pnt P1 = Kernel::face_center(face);
    gp_Pnt P2 = Kernel::face_center(parent->face);
    gp_Pln pln(P1, FixedFaceNormal());
    double dist = fabs(pln.Distance(P2));
    double a = fabs(FixedFaceNormal().Angle(parent->FixedFaceNormal()));

    if (RelProduct() != parent->RelProduct() || dist > tol_dist || a > tol_angle) {
        std::cerr << "[Warning] Parent face not matching. " << Info() << " - " << parent->Info() << " (" << Kernel::round_double_to_n_decimal_places(dist, 5) << "), (" << Kernel::round_double_three_digits(a * 180 / M_PI) << ")" << std::endl;
        return false;
    }

    double A1 = SurfaceArea();
    double A2 = parent->SurfaceArea();
    double dA = fabs(A1 - A2);

    if (A1 - A2 > 1e-7) {
        std::cerr << "[Warning] Parent face smaller than inner face. " << Info() << " - " << parent->Info() << " (" << Kernel::round_double_to_n_decimal_places(A1, 5) << "), (" << Kernel::round_double_to_n_decimal_places(A2, 5) << ")" << std::endl;
        return false;
    }

    if (dA < 1e-9)
        std::cout << "[Info] Parent face has same size. " << Info() << " - " << parent->Info() << " (" << Kernel::round_double_to_n_decimal_places(A1, 5) << "), (" << Kernel::round_double_to_n_decimal_places(A2, 5) << ")\n";
    else if (dA < 0.001)
        std::cout << "[Info] Parent face has approximately same size. " << Info() << " - " << parent->Info() << " (" << Kernel::round_double_to_n_decimal_places(A1, 5) << "), (" << Kernel::round_double_to_n_decimal_places(A2, 5) << ")\n";

    return true;
}

bool cFace::CheckSpaceBehind() {

    if (corresponding != nullptr) {

        if (space_behind == nullptr) {

            std::cerr << "[Warning] Face has a corresponding face but no space behind. " << id << " - " << corresponding->id << " - " << (sb_type == SB_TYPE_2B) << " - " << (physicalOrVirtual == SB_PV_VIRTUAL) << std::endl;
            return false;

        } else if (corresponding->space == nullptr) {

            std::cerr << "[Warning] Corresponding face's space is null. " << id << " - " << corresponding->id << " (" << space_behind->id << ")" << " - " << (sb_type == SB_TYPE_2B) << " - " << (physicalOrVirtual == SB_PV_VIRTUAL) << std::endl;
            return false;

        } else if (space_behind != corresponding->space) {

            std::cerr << "[Warning] Corresponding face's space does not match space behind. " << id << " - " << corresponding->id << " (" << space_behind->id << ") - (" << corresponding->RelSpace()->id << ")" << " - " << (sb_type == SB_TYPE_2B) << " - " << (physicalOrVirtual == SB_PV_VIRTUAL)
                      << std::endl;
            return false;

        }
    } //else if (space_behind != nullptr) // its an ok behaviour
    //  std::cout << "[Info] Face has no corresponding but a space behind. " << id << " - " << space_behind->id << " - " << (sb_type == SB_TYPE_2B) << " - " << (physicalOrVirtual == SB_PV_VIRTUAL) << std::endl;

    return true;
}

void cFace::CheckPointers(bool secondlvl) {

    if (ancestor == nullptr)
        std::cerr << "[Debug Info] Face has no ancestor face. " << id << std::endl;

    if (space == nullptr)
        std::cerr << "[Debug Info] Face has no space. " << id << std::endl;

    if (sb_type == SB_TYPE_2A && space_behind == nullptr)
        std::cerr << "[Debug Info] Face has no space behind. " << id << std::endl;

    //if (unifying_face == nullptr)
    //    std::cerr << "[Warning] Face has no unifying_face. " << id << std::endl;

    if (IsOpening() && parent == nullptr)
        std::cerr << "[Debug Info] Face has no parent. " << id << std::endl;

    if (sb_type == SB_TYPE_2A && corresponding == nullptr)
        std::cerr << "[Debug Info] Face has no corresponding. " << id << std::endl;

    if (secondlvl && superface == nullptr)
        std::cerr << "[Debug Info] Face has no superface. " << id << std::endl;
}

bool cFace::CheckCorresponding(double tol_area, double tol_angle) {

    if (physicalOrVirtual == SB_PV_VIRTUAL && corresponding == nullptr) // happens if space did not touch surrounding elements
        std::cout << "[Info] Virtual face does not have a corresponding face. " << id << "\n";

    if (corresponding == nullptr) return true;

    if (this != corresponding->corresponding) {
        std::string s = corresponding->corresponding == nullptr ? "null" : std::to_string(corresponding->corresponding->ID());
        std::cerr << "[Warning] Corresponding face pair not matching. " << Info() << " - " << corresponding->Info() << " | (" << RelSpace()->is_facade << ") - (" << corresponding->RelSpace()->is_facade << ") " << s << std::endl;
        return false;
    }

    if (corresponding->IfcClass() != IfcClass() && (IsIfcClass("IfcWindow") || IsIfcClass("IfcDoor"))) {
        std::cerr << "[Warning] Corresponding face pair doesn't match in ifcclass. " << Info() << " - " << corresponding->Info() << " | (" << RelSpace()->is_facade << ") - (" << corresponding->RelSpace()->is_facade << ")" << std::endl;
        return false;
    }

    double A1 = SurfaceArea();
    double A2 = corresponding->SurfaceArea();
    double a = FixedFaceNormal().Angle(corresponding->FixedFaceNormal());

    bool b = true;

    if (fabs(A1 - A2) / std::min(A1, A2) > tol_area) {
        std::cerr << "[Warning] Corresponding face pair doesn't match in area. " << Info() << " - " << corresponding->Info() << " | (" << RelSpace()->is_facade << ") - (" << corresponding->RelSpace()->is_facade << "), (" << A1
                  << ") - (" << A2 << "-> " << fabs(A1 - A2) / std::min(A1, A2) << "), (" << Kernel::round_double_to_n_decimal_places(a * 180 / M_PI, 1) << ")" << std::endl;
        b = false;
    }

    if (a < tol_angle) {
        std::cerr << "[Warning] Corresponding face pair doesn't match in angle. " << Info() << " - " << corresponding->Info() << " | (" << RelSpace()->is_facade << ") - (" << corresponding->RelSpace()->is_facade << "), (" << A1
                  << ") - (" << A2 << "-> " << fabs(A1 - A2) / std::min(A1, A2) << "), (" << Kernel::round_double_to_n_decimal_places(a * 180 / M_PI, 1) << ")" << std::endl;
        b = false;
    }

    auto n1 = Topo(face).vertices().Size();
    auto n2 = Topo(corresponding->face).vertices().Size();
    if (n1 != n2) {
        std::cerr << "[Warning] Corresponding face pair doesn't match in number of vertices. " << Info() << " - " << corresponding->Info() << " | (" << RelSpace()->is_facade << ") - (" << corresponding->RelSpace()->is_facade << "), (" << n1
                  << ") - (" << n2 << ") - (" << fabs(A1 - A2) / std::min(A1, A2) << "), (" << Kernel::round_double_to_n_decimal_places(a * 180 / M_PI, 1) << ")" << std::endl;
        //b = false;
    }

/*
    if ((v1 + v2).Magnitude() > 0.1) {
        std::cerr << "[Warning] Corresponding face pair doesn't match in normal direction. " << id << " - " << corresponding->id << " (" << RelProduct()->guid << ") - (" << corresponding->RelProduct()->guid << "), ("
                  << v1.X() << "\t" << v1.Y() << "\t" << v1.Z() << ") - (" << v2.X() << "\t" << v2.Y() << "\t" << v2.Z() << ")" << std::endl;
        b = false;
    }*/

    gp_Pnt C1 = Center();
    gp_Pnt C2 = corresponding->Center();
    double d = C1.Distance(C2);

    if (d > 1.0e-9) { // skip test for virtual face pairs without zero distance

        double dh = 0.5 * d;

        gp_Vec v1 = gp_Vec(FixedFaceNormal());
        gp_Vec v2 = gp_Vec(corresponding->FixedFaceNormal());

        gp_Pnt P1 = (gp_Vec(C1.XYZ()) + v1.Scaled(dh)).XYZ();
        gp_Pnt P2 = (gp_Vec(C2.XYZ()) + v2.Scaled(dh)).XYZ();

        // distance of center points must be smaller than distance of points moved along the normals to other center
        if (P1.Distance(C2) < d || P2.Distance(C1) < d) {
            std::cerr << "[Warning] Corresponding face pair doesn't match in normal direction. " << id << " - " << corresponding->id << " (" << RelProduct()->guid << ") - (" << corresponding->RelProduct()->guid << ") - (" << RelSpace()->is_facade << ") - (" << corresponding->RelSpace()->is_facade
                      << "), (" << v1.X() << "\t" << v1.Y() << "\t" << v1.Z() << ") - (" << v2.X() << "\t" << v2.Y() << "\t" << v2.Z() << "), (" << C1.X() << "\t" << C1.Y() << "\t" << C1.Z() << ") - (" << C2.X() << "\t" << C2.Y() << "\t" << C2.Z() << ") - (" << P1.Distance(C2) << "\t"
                      << P2.Distance(C1) << "\t" << d << ")" << std::endl;
            b = false;
        }
    }

    return b;
}

bool cFace::IsOpening() const { return Ancestor() != nullptr && ancestor->IsOpening(); }

IfcUtil::IfcBaseEntity *cFace::IfcProduct() const { return ancestor->RelProduct()->ifcproduct; }

bool cFace::IsOffset() const { return Ancestor() != nullptr && ancestor->IsOffset(); }

Product *cFace::RelProduct() const { return ancestor->RelProduct(); }

std::vector<std::pair<Product *, double>> cFace::MaterialLayers(std::set<std::pair<Product *, Product *>> &collisions) {

    std::vector<std::pair<Product *, double>> R;

    // create layers from materials map
    std::list<Layer> layers = get_layers();

    if (layers.empty()) {

        return R;

    } else if (layers.size() == 1) {

        R.push_back({std::make_pair(layers.front().RelProduct(), layers.front().Dmax() - layers.front().Dmin())});
        return R;

    } else {

        // split layers to have smallest segments
        std::list<Layer> split = get_splitlayers(layers);

        // chose one split layer for each interval
        std::vector<Layer *> dom = get_dominant_splitlayers(split, collisions);

        // recombine split layers if same material next to each other
        std::vector<Layer *> comb = recombine_splitlayers(dom);

        // export info. Product and product thickness
        for (const auto &l: comb)
            R.push_back({std::make_pair(l->RelProduct(), l->Dmax() - l->Dmin())});

        return R;
    }

}

std::list<Layer> cFace::get_layers() {

    std::list<Layer> L;

    // transform two subsequent material entries (depth) to one layer
    for (const auto &e: materials) {
        if (e.second.size() % 2 != 0) // odd
            std::cerr << "[Warning] For face " << Info() << ", " << IsOffset() << ") the product " << e.first->guid << " appears only in odd number (" << e.second.size() << ") in material list (layer cannot be closed)." << std::endl;
        else {
            std::vector<double> ds;
            for (const auto &d: e.second)
                ds.push_back(Kernel::round_double_three_digits(d));
            //std::copy(e.second.begin(), e.second.end(), std::back_inserter(ds));

            for (unsigned int i = 0; i < ds.size(); i += 2) {
                double start = ds[i];
                double end = ds[i + 1];
                double thickness = end - start;
                if (thickness < 0.001) continue;
                L.emplace_back(e.first, start, end, -(thickness)); // the thinner, the higher the rank
            }
        }
    }

    return L;
}

std::list<Layer> cFace::get_splitlayers(std::list<Layer> &layers) {

    std::set<double> S; // split points (all layer borders)
    for (auto &l: layers) {
        S.insert(l.Dmin());
        S.insert(l.Dmax());
    }

    // save (relevant) split points in layer
    std::map<Layer *, std::vector<double >> m;
    for (auto &layer: layers) {
        for (const auto &d: S)
            if (d > layer.Dmin() && d < layer.Dmax()) {
                //std::cerr << "[Warning] Collision of products in original IFC file at " << layer.P->guid << "!" << std::endl;
                m[&layer].push_back(d);
            }
        std::sort(m[&layer].begin(), m[&layer].end());
    }

    // split layer at split points
    std::list<Layer> L;
    for (auto &layer: layers) {

        if (m[&layer].empty()) // no split necessary
            L.emplace_back(layer.RelProduct(), layer.Dmin(), layer.Dmax(), layer.Rank());
        else { // split into (n_splitpoints + 1) parts
            L.emplace_back(layer.RelProduct(), layer.Dmin(), m[&layer][0], layer.Rank()); // start to first split point
            L.emplace_back(layer.RelProduct(), m[&layer].back(), layer.Dmax(), layer.Rank()); // last split point to end
            for (unsigned int i = 0; i < m[&layer].size() - 1; i++)
                L.emplace_back(layer.RelProduct(), m[&layer][i], m[&layer][i + 1], layer.Rank()); // in between
        }
    }

    return L;
}

std::vector<Layer *> cFace::get_dominant_splitlayers(std::list<Layer> &layers, std::set<std::pair<Product *, Product *>> &collisions) {

    std::vector<Layer *> L;

    std::map<std::pair<double, double>, std::vector<Layer *>> M;
    for (auto &l: layers)
        M[std::make_pair(l.Dmin(), l.Dmax())].push_back(&l);

    for (auto &m: M) {
        Layer *c = m.second[0];
        for (unsigned i = 1; i < m.second.size(); i++)
            if (m.second[i]->Rank() > c->Rank())
                c = m.second[i];
        L.push_back(c);
    }

    // collision detector
    for (auto &m: M) {
        if (m.second.size() > 1) {
/*            std::cout << "[Info] Collision of products in original IFC file between ";
            for (const auto &l: m.second)
                std::cout << "\t" << l->RelProduct()->guid;
            std::cout << " .\n";*/
            for (unsigned int i = 0; i < m.second.size(); i++)
                for (unsigned int j = 0; j < m.second.size(); j++)
                    if (i < j)
                        collisions.insert(std::make_pair(m.second[i]->RelProduct(), m.second[j]->RelProduct()));
        }
    }

    return L;
}

std::vector<Layer *> cFace::recombine_splitlayers(std::vector<Layer *> &layers) {

    std::vector<Layer *> L;

    L.push_back(layers[0]);
    for (unsigned int i = 1; i < layers.size(); i++)
        if (layers[i]->RelProduct() == L.back()->RelProduct())
            L.back()->SetDmax(layers[i]->Dmax());
        else
            L.push_back(layers[i]);

    return L;
}

IntersectionRay cFace::RayBehind(double tol) {

    // create point
    gp_Pnt C = Kernel::move_point_along_scaled_unit_vector(P, FixedFaceNormal(), tol); // move point in normal direction to ensure coplanar hits
    std::array<double, 3> Q = {C.X(), C.Y(), C.Z()};

    // create direction
    gp_Dir v = FixedFaceNormal().Reversed();
    std::array<double, 3> dir = {v.X(), v.Y(), v.Z()};

    // create length by evaluating the maximum distance in material list
    double l = 0;
    if (!materials.empty())
        for (const auto &e: materials)
            if (*e.second.rbegin() > l) l = *e.second.rbegin(); // set is ordered, so take last element

    if (l < 1.0e-9) {
        if (IsVirtual()) std::cout << "[Info] Face " << Info() << ". Ray length is zero (virtual face)." << "\n";
        else if (IsOffset()) std::cout << "[Info] Face " << Info() << ". Ray length is zero (offset face)." << "\n";
        else std::cerr << "[Warning] Face " << Info() << ". Ray length is zero." << "\n";
    }

    l += 2 * tol; // extend length so that ray start tol before bound and ends tol after transmission length

    return {Q, dir, l};
}

IntersectionRay cFace::RayBehind(double tol, double l) {

    // create point
    gp_Pnt C = Kernel::move_point_along_scaled_unit_vector(P, FixedFaceNormal(), tol); // move point in normal direction to ensure coplanar hits
    std::array<double, 3> Q = {C.X(), C.Y(), C.Z()};

    // create direction
    gp_Dir v = FixedFaceNormal().Reversed();
    std::array<double, 3> dir = {v.X(), v.Y(), v.Z()};

    l += 2 * tol; // extend length so that ray start tol before bound and ends tol after transmission length

    return {Q, dir, l};
}

bool cFace::IsPolygon() const { return Kernel::face_is_polygon(face); }

TopoDS_Face cFace::Face() const { return face; }

std::string cFace::IfcClass() const {

    if (ancestor == nullptr) {
        std::cerr << "[Warning] No ancestor. " << id << std::endl;
        return "No ancestor";
    } else return ancestor->IfcClass();
}

std::string cFace::IfcGuid() const {

    if (ancestor == nullptr) {
        std::cerr << "[Warning] No ancestor. " << id << std::endl;
        return "No ancestor";
    } else return ancestor->IfcGuid();
}

bool cFace::IsIfcClass(std::string s) const {

    if (ancestor == nullptr) {
        std::cerr << "[Warning] No ancestor. " << id << std::endl;
        return false;
    } else return ancestor->IsIfcClass(std::move(s));
}

std::list<cFace *> cFace::Subfaces() const { return subfaces; }

std::map<Product *, std::set<double>>

cFace::Materials() const { return materials; }

void cFace::SetMaterials(const std::map<Product *, std::set<double>> &M) { materials = M; }

void cFace::ClearMaterials() { materials.clear(); }

SB_internal_or_external_boundary cFace::InternalOrExternal() const { return internalOrExternal; }

SB_physical_or_virtual_boundary cFace::PhysicalOrVirtual() const { return physicalOrVirtual; }

gp_Dir cFace::AncestorNormal() const { return ancestor->Normal(); }

gp_Dir cFace::FaceNormal() const { return Kernel::face_normal(face); }

gp_Dir cFace::FixedFaceNormal() const {

    if (NormalStatus() == FACE_NORMAL_KNOWN)
        return fixed_normal;
    else {
        std::cerr << "[Error] Face has no fixed face normal, (1,0,0) is returned by default " << Info() << std::endl;
        return fixed_normal;
    }

}

Space *cFace::RelSpace() const { return space; }

SB_Type cFace::SBType() const { return sb_type; }

unsigned int cFace::ID() const { return id; }

bool cFace::IsTrash() const { return isTrash; }

bool cFace::IsInner() const { return isInner; }

bool cFace::IsHanging() const { return isHanging; }

bool cFace::IsCoplanar() const { return isCoplanar; }

bool cFace::IsVirtual() const { return Ancestor() != nullptr && ancestor->IsVirtual(); }

bool cFace::WasVisited() const { return wasVisited; }

long cFace::ProjectionAxis() const {
    gp_Dir n = FixedFaceNormal();
    std::vector<double> temp = {fabs(n.X()), fabs(n.Y()), fabs(n.Z())}; // determine projection axis (0,1 or 2), such that v[proj_axis] != 0
    return std::distance(temp.begin(), std::max_element(temp.begin(), temp.end()));
}

TopoDS_Wire cFace::OuterWire() const {
    return BRepTools::OuterWire(face);
}

std::vector<gp_Pnt2d> cFace::ProjectOuterWireTo2D() const {
    return Kernel::project_wire_to_2D(ProjectionAxis(), OuterWire());
}

bool cFace::IsConvex() const {
    std::vector<gp_Pnt2d> l = ProjectOuterWireTo2D();
    return Kernel::polygon_is_convex(l);
}

oFace *cFace::Ancestor() const { return ancestor; }

Space *cFace::SpaceBehind() const { return space_behind; }

cFace *cFace::Parent() const { return parent; }

cFace *cFace::Superface() const { return superface; }

cFace *cFace::UnifyingFace() const { return unifying_face; }

cFace *cFace::Corresponding() const { return corresponding; }

gp_Pnt cFace::PointOnFace() const { return P; }

IfcUtil::IfcBaseEntity *cFace::IfcRelSpaceBoundary() const { return ifcrelspaceboundary; }

void cFace::SetSpace(Space *s) { space = s; }

void cFace::SetSpaceBehind(Space *s) { space_behind = s; }

void cFace::SetInternalOrExternal(SB_internal_or_external_boundary e) { internalOrExternal = e; }

void cFace::SetPhysicalOrVirtual(SB_physical_or_virtual_boundary e) { physicalOrVirtual = e; }

void cFace::SetID(unsigned int n) { id = n; }

void cFace::SetSBType(SB_Type t) { sb_type = t; }

void cFace::SetIsTrash(bool b) {
    isTrash = b;
    if (b) RemoveFromSubfacesOfSuperFace();
}

void cFace::SetParent(cFace *p) { parent = p; }

void cFace::SetIsCoplanar(bool b) { isCoplanar = b; }

void cFace::SetIsInner(bool b) { isInner = b; }

void cFace::SetWasVisited(bool b) { wasVisited = b; }

void cFace::SetAncestor(oFace *o) { ancestor = o; }

void cFace::SetCorresponding(cFace *corr) {
    corresponding = corr;
    if (corr != nullptr)
        space_behind = corr->RelSpace();
}

void cFace::SetUnifyingFace(cFace *f) { unifying_face = f; }

void cFace::SetSuperface(cFace *f) { superface = f; }

void cFace::SetIfcRelSpaceBoundary(IfcUtil::IfcBaseEntity *sb) { ifcrelspaceboundary = sb; }

void cFace::SetPointOnFace(gp_Pnt Q) { P = Q; }

void cFace::AppendToSubfaces(cFace *f) { subfaces.push_back(f); }

void cFace::RemoveFromSubfaces(cFace *f) { subfaces.remove(f); }

void cFace::RemoveFromSubfacesOfSuperFace() {
    if (superface != nullptr)
        superface->RemoveFromSubfaces(this);
}

void cFace::AppendToMaterial(Product *p, double d) { materials[p].insert(d); }

bool cFace::HasHoles() const { return Topo(face).wires().Size() > 1; }

void cFace::Log(unsigned int i) const {
    bool has_ancestor = (Ancestor() != nullptr);
    std::string s = !has_ancestor ? "No Ancestor" : RelProduct() == nullptr ? "No Product" : IfcProduct()->data().toString();
    s.resize(55);
    std::cout << std::setfill('0') << std::setw(6) << i;
    std::cout << "\t" << s;
    std::cout << "\t" << std::setfill('0') << std::setw(5) << ID();
    std::cout << "\t" << std::setfill('0') << std::setw(12) << Kernel::hash(face);
    std::cout << "\tOff: " << IsOffset();
    std::cout << "\t" << face.IsNull();
    std::cout << "\t" << IsTrash();
    std::cout << "\t" << IsHanging();
    std::cout << "\t| Vis: \t" << WasVisited();
    std::cout << "\t" << IsInner();
    std::cout << "\t" << IsCoplanar();
    std::cout << "\t" << IsVirtual();
    std::cout << "\t" << SBType();
    std::cout << "\t| Int: \t" << InternalOrExternal();
    std::cout << "\t" << PhysicalOrVirtual();
    if (RelSpace() == nullptr) std::cout << "\t" << "No Space";
    else std::cout << "\t" << std::setfill('0') << std::setw(4) << RelSpace()->id;
    if (IfcRelSpaceBoundary() == nullptr) std::cout << "\t" << "No IFCSB";
    else std::cout << "\t" << IfcRelSpaceBoundary()->data().getArgument(0);
    if (Corresponding() == nullptr) std::cout << "\tCorr: " << "No Corr";
    else std::cout << "\tCorr: " << std::setfill('0') << std::setw(5) << Corresponding()->ID();
    if (Parent() == nullptr) std::cout << "\t" << "No Par";
    else std::cout << "\t" << std::setfill('0') << std::setw(5) << Parent()->ID();
    if (SpaceBehind() == nullptr) std::cout << "\t" << "NoBSp";
    else std::cout << "\t" << std::setfill('0') << std::setw(5) << SpaceBehind()->id;
    std::cout << "\t| Mat: " << std::setfill('0') << std::setw(3) << Materials().size();
    std::cout << "\t" << std::setfill('0') << std::setw(2) << halfedges.size();
    std::cout << "\t" << std::setfill('0') << std::setw(2) << adjacentfaces.size();
    unsigned int n_adj = 0;
    for (auto &as: adjacentfaces) for (auto &a: as.second) n_adj++;
    std::cout << "\t" << std::setfill('0') << std::setw(2) << n_adj;
    std::cout << "\t" << std::setfill('0') << std::setw(3) << (has_ancestor ? Ancestor()->ShellID() : -999);
    std::cout << "\t" << IsOpening();
    std::cout << "\tA: " << std::setfill('0') << std::setw(5) << SurfaceArea();
    auto n1 = FaceNormal();
    int places = 4;
    std::cout << "\tn: " << std::setfill('0')
              << std::setw(5) << Kernel::round_double_to_n_decimal_places(n1.X(), places) << "|" << std::setfill('0')
              << std::setw(5) << Kernel::round_double_to_n_decimal_places(n1.Y(), places) << "|" << std::setfill('0')
              << std::setw(5) << Kernel::round_double_to_n_decimal_places(n1.Z(), places);
    auto n = fixed_normal;
    std::cout << "\tn: " << std::setfill('0')
              << std::setw(5) << Kernel::round_double_to_n_decimal_places(n.X(), places) << "|" << std::setfill('0')
              << std::setw(5) << Kernel::round_double_to_n_decimal_places(n.Y(), places) << "|" << std::setfill('0')
              << std::setw(5) << Kernel::round_double_to_n_decimal_places(n.Z(), places);
    std::cout << "\t" << NormalStatus();
    auto PC = Kernel::face_center(face);
    std::cout << "\tC: " << std::setfill('0')
              << std::setw(5) << Kernel::round_double_to_n_decimal_places(PC.X(), places) << "|" << std::setfill('0')
              << std::setw(5) << Kernel::round_double_to_n_decimal_places(PC.Y(), places) << "|" << std::setfill('0')
              << std::setw(5) << Kernel::round_double_to_n_decimal_places(PC.Z(), places);
    std::cout << "\t" << std::setfill('0') << std::setw(3) << Topo(face).seam_edges().Size();
    std::cout << "\t" << std::setfill('0') << std::setw(3) << distance;
    std::cout << "\n";
}

std::set<cFace *> cFace::AdjacentFacesList() const {

    std::set<cFace *> L;

    for (auto &a: adjacentfaces)
        for (auto &n: a.second)
            L.insert(n);

    return L;
}

bool cFace::IsPotentiallyShading(const std::set<Product *> &products_bounding_ifcspaces) const {

    if (IsOpening()) return false;
    if (RelSpace() == nullptr) return false;
    if (!RelSpace()->is_facade) return false;
    if (SBType() == SB_TYPE_2B && DistanceToSpaceBoundaryBehind() > 1e-9) return true;

    bool has_opening_nb = false;
    if (superface != nullptr)
        for (auto &nb: Superface()->AdjacentFacesList()) {
            if (nb->IsIfcClass("IfcCurtainWall") || nb->IsOpening()) {
                has_opening_nb = true;
                break;
            }
        }

    if (SBType() == SB_TYPE_2B && (products_bounding_ifcspaces.find(RelProduct()) == products_bounding_ifcspaces.end() && has_opening_nb)) return true;
    if (Corresponding() == nullptr) return false;
    if (SBType() == SB_TYPE_2A && Corresponding()->SBType() == SB_TYPE_2A && Corresponding()->RelSpace()->is_facade) return true;
    return false;
}

void cFace::AddSuperfaceSubfaceRelationship(cFace *superface, cFace *subface) {
    superface->AppendToSubfaces(subface);
    subface->SetSuperface(superface);
}

face_normal_status cFace::NormalStatus() const { return status_normal; }

void cFace::SetPropertiesSpaceApproach(cFace *behind_cface, bool behind_face_is_corresponding, double min_dist, double bounds_min_Z, gp_Pnt C) {

    distance = min_dist;
    internalOrExternal = SB_IE_INTERNAL; // no decision possible now
    physicalOrVirtual = SB_PV_PHYSICAL;
    sb_type = SB_TYPE_2B;
    materials[RelProduct()].insert(0);

    if (behind_cface != nullptr) {
        space_behind = behind_cface->space;
        materials[behind_cface->RelProduct()].insert(Kernel::round_double_to_n_decimal_places(min_dist, 5));

        if (behind_face_is_corresponding) {
            corresponding = behind_cface;
            sb_type = SB_TYPE_2A;
            if (distance < 1e-6) physicalOrVirtual = SB_PV_VIRTUAL;
        }
    }
}

IntersectionRay cFace::RayBehindSpaceApproach(double tol, double transmission_length) {

    // create point
    gp_Pnt C = Kernel::move_point_along_scaled_unit_vector(P, FixedFaceNormal(), tol); // move point in normal direction to ensure coplanar hits
    std::array<double, 3> Q = {C.X(), C.Y(), C.Z()};

    // create direction
    gp_Dir v = FixedFaceNormal().Reversed();
    std::array<double, 3> dir = {v.X(), v.Y(), v.Z()};

    // create length by evaluating the maximum distance in material list
    double l = 0;
    if (!materials.empty())
        for (const auto &e: materials)
            if (*e.second.rbegin() > l) l = *e.second.rbegin(); // set is ordered, so take last element

    if (l < 1.0e-9) { // this happens if faces are virtual or if face is external and does not have any material info yet
        if (physicalOrVirtual == SB_PV_VIRTUAL) std::cout << "[Info] " << id << " Ray length is zero (virtual face)." << "\n";
        else l = transmission_length;
    }

    l += 2 * tol; // extend length so that ray start tol before bound and ends tol after transmission length
    return {Q, dir, l};
}

oriFace *cFace::FaceWithMinimalAngleConsiderBothOrientations(unsigned int edge_id, used_orientation o1, std::list<oriFace> &L) {

    oriFace *correct_cface = nullptr;
    double minimum_angle = 400;

    gp_Dir n1 = FaceNormal();
    gp_Dir v1 = EdgeVector(edge_id);

    if (o1 == USEDREVERSED) {
        n1.Reverse();
        v1.Reverse();
    }

    for (auto &A: L) {
        gp_Dir n2 = A.iface->cface->FaceNormal();
        gp_Dir v2 = A.iface->cface->EdgeVector(edge_id);

        if (A.orient == USEDREVERSED) {
            n2.Reverse();
            v2.Reverse();
        }

        double angle = cFace::angle_between_planar_adjacent_faces(n1, v1, n2, v2, MidPoint(edge_id));

        if (angle < minimum_angle) {
            minimum_angle = angle;
            correct_cface = &A;
        }
    }
    return correct_cface;
}

std::string cFace::Info() const { return std::to_string(ID()) + " (" + IfcClass() + ", " + IfcGuid() + ", " + std::to_string(FaceID()) + ", " + std::to_string(IsOffset()) + ", " + std::to_string(IsOpening()) + ", " + std::to_string(SurfaceArea()) + ")"; }