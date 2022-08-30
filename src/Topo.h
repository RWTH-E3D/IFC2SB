// Copyright 2022 Eric Fichter
#ifndef TOPO_H
#define TOPO_H

#include "headers.h"

class Topo {

public:
    Topo(const TopoDS_Shape &_shape);

    TopoDS_ListOfShape solids();

    TopoDS_ListOfShape shells();

    TopoDS_ListOfShape faces();

    TopoDS_ListOfShape wires();

    TopoDS_ListOfShape edges();

    TopoDS_ListOfShape vertices();

    TopoDS_ListOfShape ordered_vertices_of_wire();

    TopoDS_ListOfShape ordered_edges_from_wire();

    bool all_edges_non_curved();

    TopoDS_ListOfShape seam_edges();

private:
    const TopoDS_Shape &shape;

    TopoDS_ListOfShape loop_topo(const TopAbs_ShapeEnum &topologyType);

};


#endif //TOPO_H
