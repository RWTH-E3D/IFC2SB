// Copyright 2022 Eric Fichter
#include "Topo.h"

Topo::Topo(const TopoDS_Shape &_shape) : shape(_shape) {}

TopoDS_ListOfShape Topo::solids() { return loop_topo(TopAbs_SOLID); }

TopoDS_ListOfShape Topo::shells() { return loop_topo(TopAbs_SHELL); }

TopoDS_ListOfShape Topo::faces() { return loop_topo(TopAbs_FACE); }

TopoDS_ListOfShape Topo::wires() { return loop_topo(TopAbs_WIRE); }

TopoDS_ListOfShape Topo::edges() { return loop_topo(TopAbs_EDGE); }

TopoDS_ListOfShape Topo::vertices() { return loop_topo(TopAbs_VERTEX); }

TopoDS_ListOfShape Topo::loop_topo(const TopAbs_ShapeEnum &topologyType) {

    TopExp_Explorer Ex;
    TopoDS_ListOfShape elements;

    for (Ex.Init(shape, topologyType); Ex.More(); Ex.Next())
        elements.Append(Ex.Current());

    return elements;
}

TopoDS_ListOfShape Topo::ordered_vertices_of_wire() {

    BRepTools_WireExplorer Exp(TopoDS::Wire(shape));
    TopoDS_ListOfShape vertices;

    while (Exp.More()) {
        vertices.Append(Exp.CurrentVertex());
        Exp.Next();
    }
    return vertices;
}

TopoDS_ListOfShape Topo::ordered_edges_from_wire() {

    BRepTools_WireExplorer Exp(TopoDS::Wire(shape));
    TopoDS_ListOfShape edges;

    while (Exp.More()) {
        edges.Append(Exp.Current());
        Exp.Next();
    }
    return edges;
}

bool Topo::all_edges_non_curved() {

    TopoDS_ListOfShape E = edges();

    for (const auto &edge: E)
        if (!Kernel::is_edge_line(TopoDS::Edge(edge))) return false;

    return true;
}

TopoDS_ListOfShape Topo::seam_edges() {

    TopoDS_ListOfShape L;
    for (const auto &e: edges())
        if (e.Orientation() == TopAbs_INTERNAL || e.Orientation() == TopAbs_EXTERNAL)
            L.Append(e);

    return L;
}
