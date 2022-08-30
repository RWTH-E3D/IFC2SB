// Copyright 2022 Eric Fichter
#ifndef SFACE_H
#define SFACE_H

#include "headers.h"

class Kernel;

class sFace {

public:
    sFace(TopoDS_Face _face);

    bool operator==(const sFace &F) const { return face.HashCode(INT_MAX) == F.face.HashCode(INT_MAX); }

private:

    TopoDS_Face face;
    std::unordered_map<unsigned int, TopoDS_Shape> edgeId2halfEdge;
    std::unordered_map<unsigned int, std::list<sFace *>> edgeId2adjacentsFaces;
    bool isFlagged; // marks hanging faces or visited faces
    bool isTrash; // marks trash and bad faces
    bool isOuter;
    int component_id;
    gp_Dir n; // surface normal
    gp_Pnt c; // center of face

    void update_half_edges();

    unsigned int face_id();

    void update_face_adjacencies(const TopTools_IndexedDataMapOfShapeListOfShape &edgeFaceMap, std::unordered_map<unsigned int, sFace *> &id2sFace);

    void remove_deleted_adjacent_faces(const std::list<sFace *> &trash_faces);

    bool is_id_in_edgeId2halfEdge(unsigned int id);

    bool set_to_hanging();

    gp_Dir calculate_edge_vector(unsigned int edgeID);

    gp_Pnt calculate_mid_point(unsigned int edgeID);

    void complement();

    friend class ShapeHealing;
    friend class Viewer;
};

#endif //SFACE_H