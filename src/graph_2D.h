// Copyright 2022 Eric Fichter
#ifndef GRAPH2D_H
#define GRAPH2D_H

#include "headers.h"

class Edge2D;

class Point2D {

public:
    double X;
    double Y;
    std::set<Edge2D *> edges;
    std::set<Point2D *> nbs;
    bool wasTested;
    bool isTrash;
    typeEdge2d type;

    Point2D(double _x, double _y) : X(_x), Y(_y) {
        edges.clear();
        nbs.clear();
        wasTested = false;
        isTrash = false;
        type = TYPE2D_INNER;
    }

    Point2D() : Point2D(NAN, NAN) {}

    Point2D operator-(const Point2D &v) const { return {X - v.X, Y - v.Y}; }

    Point2D operator+(const Point2D &v) const { return {X + v.X, Y + v.Y}; }

    double operator*(const Point2D &v) const { return X * v.X + Y * v.Y; }

    Point2D operator*(double m) const { return {X * m, Y * m}; }

    double Determinant(const Point2D &v) const { return X * v.Y - Y * v.X; }

    double SquaredDistance(const Point2D &v) const {
        double dx = X - v.X;
        double dy = Y - v.Y;
        return dx * dx + dy * dy;
    }

    double Magnitude() const { return sqrt(X * X + Y * Y); }

    void Normalize() {
        double m = Magnitude();
        X = X / m;
        Y = Y / m;
    }

    Point2D UnitVector(const Point2D &end) {
        Point2D start = *this;
        Point2D v = end - start;
        v.Normalize();
        return v;
    }

    Point2D RotatedBy90DegreesCounterClockwise() const { return {-Y, X}; }

    gp_Pnt gpPnt() const { return {X, Y, 0}; }

    gp_Vec gpVec() const { return {X, Y, 0}; }

    static bool IsZeroIsct(const double v) { return fabs(v) < 1.0E-10; }

    bool Equals(const Point2D &v) const { return IsZeroIsct(X - v.X) && IsZeroIsct(Y - v.Y); }

};

struct OrientedEdge2D {
    Edge2D *edge;
    used_orientation orientation;

    OrientedEdge2D(Edge2D *_edge, used_orientation _orientation) : edge(_edge), orientation(_orientation) {}
};

class Edge2D {

public:
    Point2D *v1;
    Point2D *v2;
    const unsigned int id_wire;
    const typeEdge2d type;
    std::set<double> iscts; // list of points defined as t (n = v1 + t * (v2-v1))
    bool isTrash;
    bool isUsedForward;
    bool isUsedReversed;
    bool wasTestedForwardAsStartFace;
    bool wasTestedReversedAsStartFace;
    OrientedEdge2D nb_forward;
    OrientedEdge2D nb_reversed;

    Edge2D(Point2D *_v1, Point2D *_v2, unsigned int _id_wire, typeEdge2d _type) : v1(_v1), v2(_v2), id_wire(_id_wire), type(_type), nb_forward(OrientedEdge2D(nullptr, USEDFORWARD)), nb_reversed(OrientedEdge2D(nullptr, USEDFORWARD)) {
        iscts.clear();
        isTrash = false;
        isUsedForward = false;
        isUsedReversed = false;
        wasTestedForwardAsStartFace = false;
        wasTestedReversedAsStartFace = false;
    }

    Point2D edge_vector() const { return {*v2 - *v1}; }

    Point2D point_from_t(double t) const { return {*v1 + edge_vector() * t}; } // only valid

    double t_from_point(const Point2D &P) const {
        Point2D vec = edge_vector();
        return (vec * (P - *v1)) / (vec * vec);
    }

    bool is_degenerated_naive() const { return v1 == v2; }

    bool is_same_naive(Edge2D *e) const { return (v1 == e->v1 && v2 == e->v2) || (v1 == e->v2 && v2 == e->v1); }

    bool is_adjacent(Edge2D *e) const { return v1 == e->v1 || v1 == e->v2 || v2 == e->v1 || v2 == e->v2; }

    bool is_hanging() { return adjacent_edges_v1().empty() || adjacent_edges_v2().empty(); }

    std::set<Edge2D *> adjacent_edges() {
        std::set<Edge2D *> adj;
        adj.insert(v1->edges.begin(), v1->edges.end());
        adj.insert(v2->edges.begin(), v2->edges.end());
        adj.erase(this);
        return adj;
    }

    std::set<Edge2D *> adjacent_edges_v1() {
        std::set<Edge2D *> adj;
        adj.insert(v1->edges.begin(), v1->edges.end());
        adj.erase(this);
        return adj;
    }

    std::set<Edge2D *> adjacent_edges_v2() {
        std::set<Edge2D *> adj;
        adj.insert(v2->edges.begin(), v2->edges.end());
        adj.erase(this);
        return adj;
    }

    std::set<Edge2D *> adjacent_edges(Point2D *p) {
        std::set<Edge2D *> adj;
        adj.insert(p->edges.begin(), p->edges.end());
        adj.erase(this);
        return adj;
    }

    Point2D *first_vertex_by_orientation(used_orientation u) { return u == USEDFORWARD ? v1 : v2; }

    Point2D *last_vertex_by_orientation(used_orientation u) { return u == USEDFORWARD ? v2 : v1; }

    Point2D unit_vector_by_orientation(used_orientation u) {
        Point2D v = *last_vertex_by_orientation(u) - *first_vertex_by_orientation(u);
        v.Normalize();
        return v;
    }

    static double angle_by_arc(const Point2D &P1, const Point2D &P2, const Point2D &r1) {

        auto arc = GC_MakeArcOfCircle(P1.gpPnt(), r1.gpVec(), P2.gpPnt());

        if (!arc.IsDone()) {
            std::cerr << "[Error] Failed creating Circle!!" << std::endl;
            std::cerr << P1.X << "\t" << P1.Y << "\t" << P2.X << "\t" << P2.Y << "\t" << r1.X << "\t" << r1.Y << std::endl;
            return 360;
        }

        TopoDS_Edge E = BRepBuilderAPI_MakeEdge(arc.Value()).Edge();

        GProp_GProps props;
        BRepGProp::LinearProperties(E, props); // angle in radians
        return props.Mass() * 180.0 / M_PI; // radian to degree
    }

    void calc_nb() {
        nb_forward = nb_smallest_angle(USEDFORWARD);
        nb_reversed = nb_smallest_angle(USEDREVERSED);
        if (nb_forward.edge == nullptr) std::cerr << "HANGING1" << std::endl;
        if (nb_reversed.edge == nullptr) std::cerr << "HANGING2" << std::endl;
    }

    OrientedEdge2D nb_smallest_angle(used_orientation m) {

        Point2D *end = (m == USEDFORWARD) ? v2 : v1;

        const auto &adjacentEdges = adjacent_edges(end);

        if (adjacentEdges.empty()) return {nullptr, USEDFORWARD};

        if (adjacentEdges.size() == 1) {
            Edge2D *adj = *adjacentEdges.begin();
            return (adj->v1 == end) ? OrientedEdge2D(adj, USEDFORWARD) : OrientedEdge2D(adj, USEDREVERSED);
        }

        double amin = 400;
        OrientedEdge2D adjmin(nullptr, USEDFORWARD);

        Point2D n = unit_vector_by_orientation(m);
        Point2D vec = n.RotatedBy90DegreesCounterClockwise();
        Point2D P1 = *end - n; // distance of 1 to end

        for (const auto &adj: adjacentEdges) {
            used_orientation u = (adj->v1 == end) ? USEDFORWARD : USEDREVERSED;
            Point2D P2 = *adj->first_vertex_by_orientation(u) + adj->unit_vector_by_orientation(u); // distance of 1 to start point (common vertex)
            double a = angle_by_arc(P1, P2, vec); // create arc by two points and rotation vector for first point

            if (a >= amin) continue;

            amin = a;
            adjmin = OrientedEdge2D(adj, u);
        }
        return adjmin;
    }

    bool point_on_line(const Point2D &pt, Point2D &res) const {

        if (v1->X == v2->X && v1->Y == v2->Y) return false;

        Point2D vec = edge_vector();
        res = *v1 + vec * ((vec * (pt - *v1)) / (vec * vec));

        if (v1->Equals(res) || v2->Equals(res)) return false;

        double minx = std::min(v1->X, v2->X);
        double maxx = std::max(v1->X, v2->X);
        double miny = std::min(v1->Y, v2->Y);
        double maxy = std::max(v1->Y, v2->Y);

        return (res.X >= minx && res.X <= maxx) && (res.Y >= miny && res.Y <= maxy);
    }

    std::array<double, 4> aabb() const {
        double s = 1.0e-7;
        double xmin = std::min(v1->X, v2->X) - s;
        double xmax = std::max(v1->X, v2->X) + s;
        double ymin = std::min(v1->Y, v2->Y) - s;
        double ymax = std::max(v1->Y, v2->Y) + s;
        return {xmin, ymin, xmax, ymax};
    }
};

#endif //GRAPH2D_H
