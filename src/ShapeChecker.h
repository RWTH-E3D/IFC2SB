// Copyright 2022 Eric Fichter
#ifndef SHAPECHECKER_H
#define SHAPECHECKER_H

#include "headers.h"

class ShapeChecker {

public:
    ShapeChecker(const TopoDS_Shape &_shape);

    ~ShapeChecker();

    void dump_topology(int level = 0);

    void dump_topology(const TopoDS_Shape &shape, int level);

    void full_check();

    TopoDS_ListOfShape get_vertices();

    TopoDS_ListOfShape get_edges();

    TopoDS_ListOfShape get_wires();

    TopoDS_ListOfShape get_faces();

    TopoDS_ListOfShape get_shells();

    TopoDS_ListOfShape get_solids();

    TopoDS_ListOfShape get_compsolids();

    TopoDS_ListOfShape get_compounds();

    TopoDS_ListOfShape faces_from_wire(const TopoDS_Wire &W);

    bool face_is_planar(const TopoDS_Face &face) const;

    TopoDS_ListOfShape free_edges();

    //void getSubshapesOfShape(bool ignore_orientation=false);

    TopoDS_ListOfShape vertices;
    TopoDS_ListOfShape edges;
    TopoDS_ListOfShape wires;
    TopoDS_ListOfShape faces;
    TopoDS_ListOfShape shells;
    TopoDS_ListOfShape solids;
    TopoDS_ListOfShape compsolids;
    TopoDS_ListOfShape compounds;
    Standard_Boolean isValid;

private:

    static std::string shapeTypeString(const TopAbs_ShapeEnum &st);

    TopoDS_ListOfShape loop_topo(const TopAbs_ShapeEnum &topologyType, bool ignore_same_hashes = false);

    void get_results_from_analyzer(const TopoDS_Shape &shape);

    static std::string int_to_status(int i);

    TopoDS_ListOfShape map_shapes_and_ancestors(TopAbs_ShapeEnum topoTypeA, TopAbs_ShapeEnum topoTypeB, TopoDS_Shape topologicalEntity);

    void check_vertices();

    void check_edges();

    void check_wires();

    void check_faces();

    void check_shells();

    void check_vertex(const TopoDS_Vertex &V);

    void check_edge(const TopoDS_Edge &E);

    void check_wire(const TopoDS_Wire &W);

    void check_face(const TopoDS_Face &F);

    void check_shell(const TopoDS_Shell &S);

    void add_to_errors(const TopoDS_Shape &shape, BRepCheck_Status status);

    const TopoDS_Shape &shape;

    const double tol;

    BRepCheck_Analyzer *analyzer;

    std::unordered_map<const TopoDS_Shape *, std::list<BRepCheck_Status>> errors;

    void print_errors();

};


#endif //SHAPECHECKER_H
