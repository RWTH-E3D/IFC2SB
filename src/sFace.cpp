// Copyright 2022 Eric Fichter
#include "sFace.h"

sFace::sFace(TopoDS_Face _face) : face(_face) {

    n = Kernel::face_normal(_face);
    c = Kernel::face_center(_face);

    edgeId2halfEdge.clear();
    edgeId2adjacentsFaces.clear();
    isFlagged = false;
    isTrash = false;
    isOuter = false;
    component_id = -1;

}

unsigned int sFace::face_id() { return Kernel::hash(face); }

void sFace::update_half_edges() {

    TopExp_Explorer Ex;

    for (Ex.Init(face, TopAbs_EDGE); Ex.More(); Ex.Next()) { // skip seam edges
        if (Ex.Current().Orientation() == TopAbs_INTERNAL || Ex.Current().Orientation() == TopAbs_EXTERNAL)
            continue;
        edgeId2halfEdge[Kernel::hash(Ex.Current())] = Ex.Current();
    }
}

void sFace::update_face_adjacencies(const TopTools_IndexedDataMapOfShapeListOfShape &edgeFaceMap, std::unordered_map<unsigned int, sFace *> &id2sFace) {

    // for every half-edge, get adjacent faces
    for (auto &itr: edgeId2halfEdge) {
        TopTools_ListOfShape adjacent_faces = edgeFaceMap.FindFromKey(itr.second);
        adjacent_faces.Remove(face); // remove self-reference
        unsigned int edge_id = itr.first;

        // get cFace of the adjacent faces
        std::list<sFace *> adjacent_sfaces;

        for (auto &adjacent_face: adjacent_faces) {

            Standard_Integer adjHash = adjacent_face.HashCode(INT_MAX);

            // only add if faces are in id2cFace. if not, they are in fuse shape, but not in cfaces anymore (e.g. because of deletion of duplicate faces)
            if (id2sFace.find(adjHash) != id2sFace.end()) {
                sFace *adjacent_sface = id2sFace[adjHash];

                // don't add adjacent faces, if they are connected by a seam edge - seam edge or seam edge - normal edge combo. So look if edge id is in maps of the "adjacent face" first
                if (adjacent_sface->is_id_in_edgeId2halfEdge(edge_id))
                    adjacent_sfaces.push_back(adjacent_sface);
            }
        }
        edgeId2adjacentsFaces[edge_id] = adjacent_sfaces;
    }
}

bool sFace::is_id_in_edgeId2halfEdge(unsigned int id) { return !(edgeId2halfEdge.find(id) == edgeId2halfEdge.end()); }

bool sFace::set_to_hanging() {

    if (isFlagged)
        return false;

    bool is_joined = true; // indicates if face has at least one adjacent SHELL face on all half-edges

    for (auto &it: edgeId2adjacentsFaces) {

        // check for existence of adjacent faces on edge
        if (it.second.empty()) {
            is_joined = false;
            break;
        }

        // there are faces on edge, so check if there is at least a non hanging adjacent SHELL face
        bool has_non_hanging_adjacent_faces = false;
        for (const auto &adjacent_sface: it.second) {
            if (!adjacent_sface->isFlagged) {
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
        isFlagged = true;
        isTrash = true;
        return true;
    } else
        return false;

}

void sFace::remove_deleted_adjacent_faces(const std::list<sFace *> &trash_faces) {

    for (const auto &trash_face: trash_faces)
        for (auto &e: edgeId2adjacentsFaces)
            e.second.remove(trash_face);
}

void sFace::complement() {
    face.Complement();
    n.Reverse();
    update_half_edges();
}

gp_Dir sFace::calculate_edge_vector(unsigned int edgeID) {
    TopoDS_Edge Edge = TopoDS::Edge(edgeId2halfEdge[edgeID]);
    gp_Pnt P1 = BRep_Tool::Pnt(TopExp::FirstVertex(Edge, Standard_True));
    gp_Pnt P2 = BRep_Tool::Pnt(TopExp::LastVertex(Edge, Standard_True));
    return {gp_Vec(P1, P2)};
}

gp_Pnt sFace::calculate_mid_point(unsigned int edgeID) {
    // computes the point that lies in the middle between P1 and P2
    TopoDS_Edge Edge = TopoDS::Edge(edgeId2halfEdge[edgeID]);
    gp_Pnt P1 = BRep_Tool::Pnt(TopExp::FirstVertex(Edge, Standard_True));
    gp_Pnt P2 = BRep_Tool::Pnt(TopExp::LastVertex(Edge, Standard_True));
    gp_Vec vec = 0.5 * (gp_Vec(P1.XYZ()) + gp_Vec(P2.XYZ()));
    return {vec.XYZ()};
}