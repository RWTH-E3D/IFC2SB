// Copyright 2022 Eric Fichter
#ifndef SHAPEHEALING_H
#define SHAPEHEALING_H

#include "headers.h"

class sFace;

class ShapeHealing {

public:
    ShapeHealing(TopoDS_Shape _Shp);

    TopoDS_Shape heal_shape(bool ignore_curved_edges = false);

    void free_and_bad_edges();

private:
    TopoDS_Shape Shp;
    TopTools_IndexedDataMapOfShapeListOfShape edgeFaceMap;
    std::list<sFace> sFaces;
    std::unordered_map<unsigned int, sFace *> id2sFace;
    int component_id_counter;
    double fuse_tol;

    void remove_duplicated_hashes_from_shape();

    void curvatures(bool &has_curved_edges, bool &has_curved_faces);

    bool find_duplicated_hashes_in_shape(std::list<unsigned int> &bad_hashes);

    bool find_duplicated_hashes(std::list<unsigned int> &bad_hashes);

    void fuse();

    void fill_face_list();

    void setup_sFace_map();

    void identify_hanging_faces();

    bool identify_hanging_faces_while();

    void remove_trash_and_face_adjacency();

    void find_components();

    void find_start_face_for_component_search(sFace *&start_sface);

    TopoDS_Solid get_outer_hull_of_component(int i, bool &isGood);

    TopoDS_ListOfShape components_to_solids();

    static bool isOuterFace(const TopoDS_Face &face, gp_Pnt P, const gp_Dir &n, IntCurvesFace_ShapeIntersector &Intersector);

    static void topoHull(sFace *sface, std::stack<sFace *> &stack);

    void remove_bad_faces();
};

#endif //SHAPEHEALING_H
