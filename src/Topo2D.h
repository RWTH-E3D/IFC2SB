// Copyright 2022 Eric Fichter
#ifndef TOPO2D_H
#define TOPO2D_H

#include "headers.h"

class OrientedEdge2D;

class Kernel;

class Point2D;

class Edge2D;

class Topo2D {

public:

    Topo2D(double _tol_vertex, unsigned int _proj_axis, gp_Dir tn1, double _plane_a);

    std::list<Point2D> vertices;
    std::list<Edge2D> edges;

    bool Fuse(TopoDS_Shape &shp);

    void Add(const std::vector<gp_Pnt2d> &wire, unsigned int wire_id, typeEdge2d type);

    static bool IsZeroIsct(const double v);

private:

    double tol_isct;
    double TOL_isct;
    double tol_vertex;
    const unsigned int proj_axis;
    const gp_Dir n1;          // n1 ... normal of original face
    const double plane_a;     // plane equation constant, n[0] * x + n[1] * y + n[2] * z = a

    typedef std::vector<OrientedEdge2D> SolutionChain;

    static TopoDS_ListOfShape occ_faces(std::list<SolutionChain> &wires, std::map<SolutionChain *, TopoDS_Wire> &occ_wires, std::map<Topo2D::SolutionChain *, std::list<Topo2D::SolutionChain *>> &WireHoles);

    bool InsidePolygon(const SolutionChain &polygon, const Point2D &p);

    bool InsidePolygon(const SolutionChain &polygon, const SolutionChain &hole);

    static double Angle2D(double x1, double y1, double x2, double y2);

    gp_Pnt reproject(const Point2D &P_2D);

    static bool closed(const SolutionChain &L);

    static bool isHole(const SolutionChain &L);

    static bool orientation(const SolutionChain &p);

    static bool orientation(const std::vector<Point2D *> &l);

    std::list<std::set<Point2D *>> find_vertex_clusters(double s);

    std::map<SolutionChain *, std::list<SolutionChain *>> detect_holes_in_polygons(std::list<SolutionChain> &wires, std::list<SolutionChain> &holes);

    void add(const gp_Pnt2d &P1, const gp_Pnt2d &P2, unsigned int wire_id, typeEdge2d type);

    void calc_ee_intersections();

    void calc_nbs();

    void find_vertex_neighbours(double s);

    std::list<std::set<Point2D *>> find_chains();

    void calc_ve_intersections(double s);

    void create_new_vertices(const std::list<std::set<Point2D *>> &chains);

    void find_polygons(std::list<SolutionChain> &wires, std::list<SolutionChain> &holes);

    void inherit_type_for_vertices();

    void intersections(Edge2D &E1, Edge2D &E2) const;

    void merge_vertices(double tol);

    void occ_instances(std::list<SolutionChain> &wires, std::list<SolutionChain> &holes, std::map<Point2D *, TopoDS_Vertex> &occ_vertices, std::map<Edge2D *, TopoDS_Edge> &occ_edges, std::map<SolutionChain *, TopoDS_Wire> &occ_wires);

    void remove_degenerated_edges();

    void remove_duplicate_and_single_edges();

    void remove_hanging_edges();

    void remove_trash_vertices();

    void remove_unused_vertices();

    void reset_vertices();

    void save_add_edge(Point2D *v1, Point2D *v2, unsigned int id_wire, typeEdge2d type);

    void save_remove_trash_edges();

    void split_Edge2D(Edge2D &E);

    void split_edges(double tol);

    static bool shared_edges(std::set<Edge2D *> &edges1, std::set<Edge2D *> &edges2, std::set<Edge2D *> &free);

    std::list<Topo2D::SolutionChain> merge_adjacent_holes(std::list<SolutionChain *> &holes);

    static void remove_stacked_holes(SolutionChain &p, std::map<Topo2D::SolutionChain *, std::list<Topo2D::SolutionChain *>> &WireHoles);

    void merge_adjacent_holes(SolutionChain &p, std::map<Topo2D::SolutionChain *, std::list<Topo2D::SolutionChain *>> &WireHoles, std::list<SolutionChain> &merged_holes, std::map<Edge2D *, TopoDS_Edge> &occ_edges, std::map<SolutionChain *, TopoDS_Wire> &occ_wires);

    static std::vector<Point2D *> oriented_chain_to_oriented_vertices(const SolutionChain &p);

    static bool pointwise_self_intersection(const SolutionChain &p);

    static bool pointwise_self_intersection(const SolutionChain &p, Point2D *&v);

    void remove_pointwise_self_intersections(std::list<SolutionChain> &L);

    static bool is_point_in_chain(const SolutionChain &polygon, Point2D *v);
};

#endif //TOPO2D_H