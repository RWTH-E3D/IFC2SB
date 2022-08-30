// Copyright 2022 Eric Fichter
#include "Topo2D.h"

Topo2D::Topo2D(double _tol_vertex, unsigned int _proj_axis, gp_Dir tn1, double _plane_a) : proj_axis(_proj_axis), n1(tn1), plane_a(_plane_a) {
    tol_vertex = _tol_vertex;
    vertices.clear();
    edges.clear();
    tol_isct = 1.0e-12;
    TOL_isct = 1 - tol_isct;
}

bool Topo2D::IsZeroIsct(const double v) { return fabs(v) < 1.0E-10; }

void Topo2D::remove_unused_vertices() {
    vertices.remove_if([](Point2D &v) { return v.edges.empty(); });
}

void Topo2D::save_add_edge(Point2D *v1, Point2D *v2, unsigned int id_wire, typeEdge2d type) {
    edges.emplace_back(v1, v2, id_wire, type);
    auto e = &edges.back();

    // store edge info in vertices
    v1->edges.insert(e);
    v2->edges.insert(e);
}

void Topo2D::save_remove_trash_edges() {
    edges.remove_if([](Edge2D &e) {
        if (e.isTrash) {
            e.v1->edges.erase(&e);
            e.v2->edges.erase(&e);
            return true;
        } else return false;
    });

    remove_unused_vertices();
}

void Topo2D::reset_vertices() {
    for (auto &v: vertices) {
        v.nbs.clear();
        v.wasTested = false;
        v.isTrash = false;
    }
}

void Topo2D::inherit_type_for_vertices() {
    for (auto &e: edges) {
        e.v1->type = e.type;
        e.v2->type = e.type;
    }
}

void Topo2D::intersections(Edge2D &E1, Edge2D &E2) const {

    // https://www.codeproject.com/Tips/862988/Find-the-Intersection-Point-of-Two-Line-Segments
    const Point2D &P1 = *E1.v1;
    const Point2D &P2 = *E1.v2;
    const Point2D &Q1 = *E2.v1;
    const Point2D &Q2 = *E2.v2;

    Point2D p1q1 = Q1 - P1;
    Point2D rE1 = E1.edge_vector();
    Point2D rE2 = E2.edge_vector();
    double mE1 = rE1 * rE1;
    double mE2 = rE2 * rE2;
    double rxs = rE1.Determinant(rE2);
    double qpxr = p1q1.Determinant(rE1);

    // If r x s = 0 and (pq) x r != 0, then the two lines are parallel and non-intersecting.
    if (IsZeroIsct(rxs) && !IsZeroIsct(qpxr))
        return;

    // If r x s = 0 and (pq) x r = 0, then the two lines are collinear.
    if (IsZeroIsct(rxs) && IsZeroIsct(qpxr)) {

        // calculate t values and check if other vertices are between 0 and 1. http://stackoverflow.com/a/565282/292237
        double t1, t2;

        // segment 1
        Point2D p1q2 = Q2 - P1;
        t1 = (p1q1 * rE1) / mE1;
        t2 = (p1q2 * rE1) / mE1;
        if (t1 > tol_isct && t1 < TOL_isct) E1.iscts.insert(t1);
        if (t2 > tol_isct && t2 < TOL_isct) E1.iscts.insert(t2);

        // segment 2
        Point2D q1p1 = P1 - Q1;
        Point2D q1p2 = P2 - Q1;
        t1 = (q1p1 * rE2) / mE2;
        t2 = (q1p2 * rE2) / mE2;
        if (t1 > tol_isct && t1 < TOL_isct) E2.iscts.insert(t1);
        if (t2 > tol_isct && t2 < TOL_isct) E2.iscts.insert(t2);

        return;
    }

    // t = (pq) x s / (r x s)
    // u = (pq) x r / (r x s)
    double t1 = p1q1.Determinant(rE2) / rxs;
    double t2 = p1q1.Determinant(rE1) / rxs;

    // If r x s != 0 and 0 <= t <= 1 and 0 <= u <= 1
    // the two line segments meet at the point p + t r = q + u s.
    if (!IsZeroIsct(rxs) && (0 <= t1 && t1 <= 1) && (0 <= t2 && t2 <= 1)) {
        //Point2D intersection = P1 + rE1 * t;
        if (t1 > tol_isct && t1 < TOL_isct) E1.iscts.insert(t1);
        if (t2 > tol_isct && t2 < TOL_isct) E2.iscts.insert(t2);
        return;
    }
}

void Topo2D::calc_ee_intersections() {

    unsigned int i = 0;

    for (auto &e1: edges) {
        unsigned int j = 0;
        auto aabb1 = e1.aabb();

        for (auto &e2: edges) {
            if (i < j && e1.id_wire != e2.id_wire) {
                auto aabb2 = e2.aabb();
                if (Kernel::do_aabbs_intersect(aabb1[0], aabb1[2], aabb1[1], aabb1[3], aabb2[0], aabb2[2], aabb2[1], aabb2[3]))
                    intersections(e1, e2);
            }
            j++;
        }
        i++;
    }
}

void Topo2D::split_edges(double tol) {

    // calc intersections. they are stored in edge objects
    calc_ee_intersections(); // edge edge
    calc_ve_intersections(tol); // vertex edge

    for (auto &e: edges)
        split_Edge2D(e);

    save_remove_trash_edges(); // remove pointers and old edges
}

void Topo2D::calc_ve_intersections(double s) {

    // distance to merge vertices
    const double s2 = s * s;

    for (auto &e: edges) {

        auto aabb1 = e.aabb();

        for (const auto &v: vertices) {
            if (&v == e.v1 || &v == e.v2) continue;
            if (!Kernel::do_aabbs_intersect(aabb1[0], aabb1[2], aabb1[1], aabb1[3], v.X - s, v.X + s, v.Y - s, v.Y + s)) continue;

            Point2D P_proj(-1e9, -1e9);
            if (e.point_on_line(v, P_proj) && v.SquaredDistance(P_proj) <= s2)
                e.iscts.insert(e.t_from_point(P_proj));
        }
    }
}

void Topo2D::split_Edge2D(Edge2D &E) {

    if (E.iscts.empty()) return;

    E.isTrash = true;

    // new vertices
    std::vector<Point2D *> L;
    L.reserve(E.iscts.size() + 2);

    L.push_back(E.v1);

    for (const auto &isct: E.iscts) {
        Point2D P = E.point_from_t(isct);
        vertices.emplace_back(P.X, P.Y);
        vertices.back().type = E.type;
        if (vertices.back().type == TYPE2D_OUTER_ORIG) // new created points have lower hierarchy then original vertices. for vertex merge they must move to original vertices
            vertices.back().type = TYPE2D_OUTER;
        L.push_back(&vertices.back());
    }

    L.push_back(E.v2);

    for (unsigned int i = 0; i < L.size() - 1; i++) {
        const auto &v1 = L[i];
        const auto &v2 = L[i + 1];
        save_add_edge(v1, v2, E.id_wire, E.type);
    }
}

void Topo2D::find_vertex_neighbours(double s) {

    // distance to merge vertices, merge tolerance
    const double s2 = s * s;

    rtree_lib::RTree<Point2D *, double, 2, double> tree;

    // populate rtree
    for (auto &v: vertices) {
        double min[2] = {v.X - s, v.Y - s};
        double max[2] = {v.X + s, v.Y + s};
        tree.Insert(min, max, &v);
    }

    // find intersections
    for (auto &v: vertices) {

        double min[2] = {v.X - s, v.Y - s};
        double max[2] = {v.X + s, v.Y + s};
        tree.Search(min, max, [&v, &s2](Point2D *f) {

            // skip self-reference and already known pairs
            if (&v == f || v.nbs.find(f) != v.nbs.end())
                return true;

            if (v.SquaredDistance(*f) <= s2) {
                v.nbs.insert(f);
                f->nbs.insert(&v);
            }

            return true;
        });
    }
}

std::list<std::set<Point2D *>> Topo2D::find_chains() {

    std::list<std::set<Point2D *>> chains;

    for (auto &v: vertices) {
        if (v.wasTested || v.nbs.empty()) continue;

        std::set<Point2D *> chain;
        std::stack<Point2D *> stack;
        stack.push(&v);

        while (!stack.empty()) {

            Point2D *p = stack.top();
            stack.pop();

            if (p->wasTested) continue;

            for (auto &nb: p->nbs)
                stack.push(nb);

            chain.insert(p);
            p->wasTested = true;
        }

        chains.push_back(chain);
    }

    return chains;
}

std::list<std::set<Point2D *>> Topo2D::find_vertex_clusters(double s) {

    find_vertex_neighbours(s);
    auto chains = find_chains();
    reset_vertices();

    return chains;
}

void Topo2D::merge_vertices(double tol) {

    std::list<std::set<Point2D *>> chains = find_vertex_clusters(tol); // find chains (neighbouring vertices)
    create_new_vertices(chains);
    remove_trash_vertices();
}

void Topo2D::create_new_vertices(const std::list<std::set<Point2D *>> &chains) {

    for (auto &chain: chains) {

        double x(0), y(0);
        typeEdge2d type;

        // check if there are "blocked" entities
        auto c1 = std::count_if(chain.begin(), chain.end(), [](const Point2D *p) { return p->type == TYPE2D_OUTER; });
        auto c2 = std::count_if(chain.begin(), chain.end(), [](const Point2D *p) { return p->type == TYPE2D_HOLE; });
        auto c3 = std::count_if(chain.begin(), chain.end(), [](const Point2D *p) { return p->type == TYPE2D_OUTER_ORIG; });
        auto c_out = c1 + c3;
        auto c = c1 + c2 + c3;

        if (c == 0) { // just merge to a new inner vertex

            type = TYPE2D_INNER;
            x = std::accumulate(chain.begin(), chain.end(), 0.0, [](double sum, Point2D *p) { return sum + p->X; }) / chain.size();
            y = std::accumulate(chain.begin(), chain.end(), 0.0, [](double sum, Point2D *p) { return sum + p->Y; }) / chain.size();

        } else if (c_out > 0) { // move to center of "locked" outer vertices

            type = (c3 > 0) ? TYPE2D_OUTER_ORIG : TYPE2D_OUTER;

            unsigned int n = 0;
            for (auto &p: chain)
                if (p->type == type) {
                    x += p->X;
                    y += p->Y;
                    n++;
                }

            x = x / n;
            y = y / n;

        } else { // move to center of "locked" hole vertices

            type = TYPE2D_HOLE;
            unsigned int n = 0;
            for (auto &p: chain)
                if (p->type == type) {
                    x += p->X;
                    y += p->Y;
                    n++;
                }
            x = x / n;
            y = y / n;

        }

        x = Kernel::round_double_to_n_decimal_places(x, 6);
        y = Kernel::round_double_to_n_decimal_places(y, 6);

        // create new vertex
        vertices.emplace_back(x, y);
        Point2D *v = &vertices.back();
        v->type = type;

        // delete old vertices
        for (auto &p: chain) {
            p->isTrash = true;
            for (auto &e: p->edges) {
                v->edges.insert(e);
                if (e->v1 == p) e->v1 = v;
                if (e->v2 == p) e->v2 = v;
            }
        }
    }

}

void Topo2D::remove_trash_vertices() { vertices.remove_if([](Point2D &p) { return p.isTrash; }); }

void Topo2D::remove_degenerated_edges() {
    for (auto &e: edges)
        if (e.is_degenerated_naive())
            e.isTrash = true;
    save_remove_trash_edges();
}

void Topo2D::remove_duplicate_and_single_edges() {
    for (auto &E: edges) {

        if (E.isTrash) continue;

        std::set<Edge2D *> adj = E.adjacent_edges(); // adjacent edges

        if (adj.size() < 2) {
            E.isTrash = true;
            continue;
        }

        for (auto &a: adj) {
            if (a->isTrash) continue;
            if (E.is_same_naive(a)) {
/*                if (((E.type == TYPE2D_OUTER_ORIG || E.type == TYPE2D_OUTER) && a->type == TYPE2D_HOLE) ||
                    (a->type == TYPE2D_OUTER_ORIG || a->type == TYPE2D_OUTER) && E.type == TYPE2D_HOLE)
                    std::cerr << "Hole edge and outer edge will be merged!" << std::endl;*/
                if (E.type == TYPE2D_OUTER_ORIG) a->isTrash = true;
                else if (a->type == TYPE2D_OUTER_ORIG) E.isTrash = true;
                else if (E.type == TYPE2D_OUTER) a->isTrash = true;
                else if (a->type == TYPE2D_OUTER) E.isTrash = true;
                else if (E.type == TYPE2D_HOLE) a->isTrash = true;
                else if (a->type == TYPE2D_HOLE) E.isTrash = true;
                else E.isTrash = true;
            }
        }
    }

    save_remove_trash_edges();
}

void Topo2D::remove_hanging_edges() {

    while (true) {
        bool found = false;

        for (auto &e: edges)
            if (e.is_hanging()) {
                found = true;
                e.isTrash;
            }
        if (!found) break;

        save_remove_trash_edges();
    }
}

bool Topo2D::orientation(const std::vector<Point2D *> &l) {
    // true for counterclockwise winding
    // https://stackoverflow.com/questions/1165647/how-to-determine-if-a-list-of-polygon-points-are-in-clockwise-order
    double v = 0;

    for (unsigned int i = 0; i < l.size() - 1; i++)
        v += (l[i + 1]->X - l[i]->X) * (l[i + 1]->Y + l[i]->Y);
    v += (l[0]->X - l.back()->X) * (l[0]->Y + l.back()->Y);

    return v < 0;
}

bool Topo2D::orientation(const SolutionChain &p) {
    return orientation(oriented_chain_to_oriented_vertices(p));
}

std::vector<Point2D *> Topo2D::oriented_chain_to_oriented_vertices(const SolutionChain &p) {
    std::vector<Point2D *> pts;
    pts.reserve(p.size());
    for (auto &l: p)
        pts.push_back(l.edge->first_vertex_by_orientation(l.orientation));
    return pts;
}

bool Topo2D::closed(const SolutionChain &L) {
    Edge2D *e1 = L.back().edge;
    Edge2D *e2 = L.front().edge;

    used_orientation m1 = L.back().orientation;
    used_orientation m2 = L.front().orientation;

    Point2D *v1 = e1->last_vertex_by_orientation(m1);
    Point2D *v2 = e2->first_vertex_by_orientation(m2);

    return v1 == v2;
}

void Topo2D::calc_nbs() { for (auto &e: edges) e.calc_nb(); }

void Topo2D::find_polygons(std::list<SolutionChain> &wires, std::list<SolutionChain> &holes) {

    // returns all counter clockwise polygons (orientation is checked before solution is accepted. also in Point2D vec = n.RotatedBy90DegreesCounterClockwise(); help vector assuems counter clockwis

    while (true) {

        // get start face
        std::stack<OrientedEdge2D> stack;

        for (auto &e: edges)
            if (!e.isUsedForward && !e.wasTestedForwardAsStartFace) { // a face can only be used once as start face (in each direction)
                stack.emplace(&e, USEDFORWARD);
                e.wasTestedForwardAsStartFace = true;
                break;
            } else if (!e.isUsedReversed && !e.wasTestedReversedAsStartFace) {
                stack.emplace(&e, USEDREVERSED);
                e.wasTestedReversedAsStartFace = true;
                break;
            }

        if (stack.empty())
            break;

        SolutionChain L;
        L.reserve(20);
        std::set<Edge2D *> checked;

        while (!stack.empty()) {

            auto p = stack.top();
            stack.pop();

            auto e = p.edge;
            auto m = p.orientation;

            if (m == USEDFORWARD && e->isUsedForward || m == USEDREVERSED && e->isUsedReversed) continue; // skip if edge already part of solution
            if (checked.find(e) != checked.end()) continue; // skip if edge was already checked within this loop

            auto nb = (m == USEDFORWARD) ? e->nb_forward : e->nb_reversed;
            stack.push(nb);

            checked.insert(e);
            L.push_back(p);
        }

        // refuse solution when certain conditions meet
        if (L.size() < 2) {
            std::cerr << "[Warning] Empty or one edge only" << std::endl;
            continue;
        }
        if (!orientation(L)) {
            //std::cerr << "Wrong orientation" << std::endl;
            continue;
        }
        if (!closed(L)) {
            std::cerr << "[Warning] Not closed" << std::endl;
            continue;
        }
        if (pointwise_self_intersection(L)) {
            std::cerr << "[Warning] Pointwise self-intersection" << std::endl;
        }

        if (isHole(L))
            holes.push_back(L);
        else
            wires.push_back(L);
        for (auto &p: L)
            if (p.orientation == USEDFORWARD) p.edge->isUsedForward = true;
            else p.edge->isUsedReversed = true;

    }

}

bool Topo2D::isHole(const SolutionChain &L) {

    for (auto &p: L)
        if (p.edge->type != TYPE2D_HOLE)
            return false;
    return true;
}

double Topo2D::Angle2D(double x1, double y1, double x2, double y2) {
    /*
   Return the angle between two vectors on a plane
   The angle is from vector 1 to vector 2, positive anticlockwise
   The result is between -pi -> pi
*/
    double dtheta, theta1, theta2;

    theta1 = atan2(y1, x1);
    theta2 = atan2(y2, x2);
    dtheta = theta2 - theta1;
    while (dtheta > M_PI)
        dtheta -= 2 * M_PI;
    while (dtheta < -M_PI)
        dtheta += 2 * M_PI;

    return (dtheta);
}

bool Topo2D::InsidePolygon(const SolutionChain &polygon, const Point2D &p) {
    // https://www.eecs.umich.edu/courses/eecs380/HANDOUTS/PROJ2/InsidePoly.html

    double angle = 0;
    Point2D p1, p2;
    auto n = polygon.size();

    for (unsigned int i = 0; i < n; i++) {
        auto j = (i + 1) % n;
        auto v1 = polygon[i].edge->first_vertex_by_orientation(polygon[i].orientation);
        auto v2 = polygon[j].edge->first_vertex_by_orientation(polygon[j].orientation);
        p1.X = v1->X - p.X;
        p1.Y = v1->Y - p.Y;
        p2.X = v2->X - p.X;
        p2.Y = v2->Y - p.Y;
        angle += Angle2D(p1.X, p1.Y, p2.X, p2.Y);
    }

    return fabs(angle) > M_PI;
}

bool Topo2D::is_point_in_chain(const SolutionChain &polygon, Point2D *v) {
    for (const auto e: polygon)
        if (e.edge->first_vertex_by_orientation(e.orientation) == v) return true;
    return false;
}

bool Topo2D::InsidePolygon(const SolutionChain &polygon, const SolutionChain &hole) {

    unsigned int c = 0;
    unsigned int s = 0;

    for (const auto elem: hole) {
        Point2D *test_point = elem.edge->first_vertex_by_orientation(elem.orientation);
        if (is_point_in_chain(polygon, test_point)) continue;

        s++;
        if (InsidePolygon(polygon, *test_point))
            c++;
    }

    return c > 0.8 * s;
}

void Topo2D::occ_instances(std::list<SolutionChain> &wires, std::list<SolutionChain> &holes, std::map<Point2D *, TopoDS_Vertex> &occ_vertices, std::map<Edge2D *, TopoDS_Edge> &occ_edges, std::map<SolutionChain *, TopoDS_Wire> &occ_wires) {

    for (auto &v: vertices)
        occ_vertices[&v] = BRepBuilderAPI_MakeVertex(reproject(v)).Vertex();

    for (auto &e: edges)
        occ_edges[&e] = BRepBuilderAPI_MakeEdge(occ_vertices[e.v1], occ_vertices[e.v2]).Edge(); // we lose information here because orientation is not considered. that's ok, because faces will oriented afterwards

    for (auto &w: wires) {
        BRepBuilderAPI_MakeWire M;
        for (auto &elem: w)
            M.Add(occ_edges[elem.edge]);
        if (!M.IsDone()) std::cerr << "[Error] Wire could not be build!" << std::endl;
        occ_wires[&w] = M.Wire();
    }

    for (auto &w: holes) {
        BRepBuilderAPI_MakeWire M;
        for (auto &elem: w)
            M.Add(occ_edges[elem.edge]);
        if (!M.IsDone()) std::cerr << "[Error] Wire could not be build!" << std::endl;
        occ_wires[&w] = M.Wire();
    }

}

std::map<Topo2D::SolutionChain *, std::list<Topo2D::SolutionChain *>> Topo2D::detect_holes_in_polygons(std::list<SolutionChain> &wires, std::list<SolutionChain> &holes) {

    std::map<SolutionChain *, std::list<SolutionChain *>> WireHoles; // polygon to its "holes"

    for (auto &wire: wires) {
        for (auto &hole: holes)
            if (InsidePolygon(wire, hole))
                WireHoles[&wire].push_back(&hole);
        for (auto &wire2: wires) {
            if (&wire == &wire2) continue;
            if (InsidePolygon(wire, wire2))
                WireHoles[&wire].push_back(&wire2);
        }
    }

    return WireHoles;
}

TopoDS_ListOfShape Topo2D::occ_faces(std::list<SolutionChain> &wires, std::map<SolutionChain *, TopoDS_Wire> &occ_wires, std::map<Topo2D::SolutionChain *, std::list<Topo2D::SolutionChain *>> &WireHoles) {

    TopoDS_ListOfShape faces;

    for (auto &wire: wires) {
        TopoDS_Face Face = BRepBuilderAPI_MakeFace(occ_wires[&wire]).Face();
        if (WireHoles[&wire].empty())
            faces.Append(Face);
        else {
            BRepBuilderAPI_MakeFace MF(Face);
            for (auto &hole: WireHoles[&wire]) {

                const TopoDS_Wire &holeWire = TopoDS::Wire(occ_wires[hole]);
                TopoDS_Face temp = TopoDS::Face(BRepBuilderAPI_MakeFace(occ_wires[hole]).Face());

                if (Kernel::face_normal(Face).IsEqual(Kernel::face_normal(temp), 1.0))
                    MF.Add(TopoDS::Wire(occ_wires[hole].Complemented()));
                else
                    MF.Add(TopoDS::Wire(occ_wires[hole]));
            }
            if (!MF.IsDone()) std::cerr << "[Error] Face not created!" << std::endl;
            faces.Append(MF.Face());
        }
    }

    return faces;
}

bool Topo2D::shared_edges(std::set<Edge2D *> &edges1, std::set<Edge2D *> &edges2, std::set<Edge2D *> &free) {

    bool b = false;

    std::map<Edge2D *, std::set<unsigned int>> M;
    for (auto &e: edges1) M[e].insert(0);
    for (auto &e: edges2) M[e].insert(1);

    for (auto &m: M)
        if (m.second.size() > 1) b = true;
        else free.insert(m.first);

    return b;
}

std::list<Topo2D::SolutionChain> Topo2D::merge_adjacent_holes(std::list<SolutionChain *> &holes) {

    std::list<std::pair<std::set<Topo2D::SolutionChain *>, std::set<Edge2D *>>> cluster, temp_cluster;

    for (auto &hole: holes) { // each hole forms a cluster. a cluster knows its edges
        std::set<Topo2D::SolutionChain *> l = {hole};
        std::set<Edge2D *> r;
        for (auto &e: *hole)
            r.insert(e.edge);
        cluster.emplace_back(l, r);
    }

    unsigned int i, j;

    while (true) {

        temp_cluster = cluster;
        cluster.clear();
        bool anything_done = false;

        i = -1;
        for (auto &c1: temp_cluster) {

            i++;
            bool was_merged = false;
            j = -1;

            for (auto &c2: temp_cluster) {

                j++;

                if (i >= j) continue;

                std::set<Edge2D *> free_edges;

                if (shared_edges(c1.second, c2.second, free_edges)) {
                    std::set<Topo2D::SolutionChain *> new_wires;
                    new_wires.insert(c1.first.begin(), c1.first.end());
                    new_wires.insert(c2.first.begin(), c2.first.end());

                    if (free_edges.empty()) {
                        std::cerr << "[Warning] Both clusters are the same" << std::endl;
                        cluster.emplace_back(new_wires, c1.second);
                    } else cluster.emplace_back(new_wires, free_edges);

                    was_merged = true;
                    c2.second.clear();
                    break;
                }
            }

            if (!was_merged && !c1.second.empty())
                cluster.push_back(c1); // cluster is finished
            else anything_done = true;

        }

        if (!anything_done)
            break;
    }

    std::list<Topo2D::SolutionChain> new_holes;

    for (auto &c: cluster) {
        Topo2D::SolutionChain chain;

        // sort to have a proper wire
        std::list<Edge2D *> L = {*c.second.begin()};

        while (true) {

            Edge2D *current = L.back();
            auto c1 = L.size();

            for (auto &e: c.second)
                if (e != current && std::find(L.begin(), L.end(), e) == L.end() && e->is_adjacent(current)) {
                    L.push_back(e);
                    break;
                }

            if (c1 == L.size()) break;
        }

        // add to new solution
        for (auto &e: L)
            chain.emplace_back(e, USEDFORWARD);

        new_holes.push_back(chain);
    }

    return new_holes;
}

void Topo2D::remove_stacked_holes(SolutionChain &p, std::map<Topo2D::SolutionChain *, std::list<Topo2D::SolutionChain *>> &WireHoles) {

    if (WireHoles[&p].size() < 2) return;

    std::set<SolutionChain *> skip;
    for (auto &hole: WireHoles[&p]) {
        for (auto &comp: WireHoles[&p]) {
            if (hole == comp) continue;
            const auto &compare_holes = WireHoles[comp]; // holes of compare
            if (std::find(compare_holes.begin(), compare_holes.end(), hole) != compare_holes.end())
                skip.insert(hole);
        }
    }
    for (auto &d: skip)
        WireHoles[&p].remove(d);
}

void Topo2D::merge_adjacent_holes(SolutionChain &p, std::map<Topo2D::SolutionChain *, std::list<Topo2D::SolutionChain *>> &WireHoles, std::list<SolutionChain> &merged_holes, std::map<Edge2D *, TopoDS_Edge> &occ_edges, std::map<SolutionChain *, TopoDS_Wire> &occ_wires) {

    if (WireHoles[&p].size() < 2) return;

    std::list<Topo2D::SolutionChain> new_holes = merge_adjacent_holes(WireHoles[&p]);

    if (new_holes.size() == WireHoles[&p].size()) return;

    WireHoles[&p].clear();

    for (auto &new_hole: new_holes) {
        merged_holes.push_back(new_hole); // add to a global list
        auto &merged_hole = merged_holes.back();
        WireHoles[&p].push_back(&merged_hole); // add to wires hole list

        BRepBuilderAPI_MakeWire M; // add to occ wire map

        for (auto &elem: merged_hole)
            M.Add(occ_edges[elem.edge]);

        if (!M.IsDone()) std::cerr << "[Error] Wire could not be build! " << M.Error() << std::endl;
        occ_wires[&merged_hole] = M.Wire();
    }
}

bool Topo2D::pointwise_self_intersection(const SolutionChain &p) {

    std::set<Point2D *> pts;
    for (auto &l: p) {
        pts.insert(l.edge->v1);
        pts.insert(l.edge->v2);
    }
    return pts.size() != p.size();
}

bool Topo2D::pointwise_self_intersection(const SolutionChain &p, Point2D *&v) {

    std::set<Point2D *> pts;
    for (auto &l: p) {
        v = l.edge->first_vertex_by_orientation(l.orientation);
        if (pts.find(v) != pts.end()) return true;
        else pts.insert(v);

    }
    return false;
}

void Topo2D::remove_pointwise_self_intersections(std::list<SolutionChain> &L) {

    for (auto &chain: L) {
        Point2D *V; // intersection/split point
        if (pointwise_self_intersection(chain, V)) {

            std::cerr << "[Warning] Solve pointwise self-intersection" << std::endl;

            int s1 = -1;
            int s2 = -1;

            for (int i = 0; i < chain.size(); i++) {
                auto &E = chain[i];
                Point2D *v = E.edge->first_vertex_by_orientation(E.orientation);
                if (v == V) {
                    if (s1 == -1) s1 = i;
                    else {
                        s2 = i - 1;
                        break;
                    }
                }
            }

            SolutionChain S1, S2;

            for (int i = 0; i < chain.size(); i++)
                if (i < s1 || i > s2) S1.push_back(chain[i]);
                else S2.push_back(chain[i]);

            // only keep the counter-clockwise oriented solution. the clock wise solution is already in the chain list as separate polygon with counter-clockwise orientation
            if (orientation(S1)) L.push_back(S1);
            if (orientation(S2)) L.push_back(S2);
            chain.clear();

        }
    }

    L.remove_if([](SolutionChain &chain) { return chain.empty(); });
}

bool Topo2D::Fuse(TopoDS_Shape &shp) {

    // assumption: edges do not lie within hole

    if (edges.empty() || vertices.empty()) {
        std::cerr << "[Warning] No geometry 1!" << std::endl;
        return false;
    }

    // preprocess
    inherit_type_for_vertices();
    merge_vertices(tol_vertex);
    remove_degenerated_edges();
    remove_duplicate_and_single_edges();

    if (edges.empty() || vertices.empty()) {
        std::cerr << "[Warning] No geometry 2!" << std::endl;
        return false;
    }

    // split edges at intersections with vertices or edges
    split_edges(0.99 * tol_vertex);
    merge_vertices(tol_vertex);
    remove_degenerated_edges();

    if (edges.empty() || vertices.empty()) {
        std::cerr << "[Warning] No geometry 3!" << std::endl;
        return false;
    }

    // twice because after moving vertices, they can land on an edge. This edge must be split again
    split_edges(0.99 * tol_vertex); // split edges at vertices
    merge_vertices(tol_vertex);
    remove_degenerated_edges();
    remove_duplicate_and_single_edges();
    remove_hanging_edges();

    if (edges.empty() || vertices.empty()) {
        std::cerr << "[Warning] No geometry 4!" << std::endl;
        return false;
    }

    // find the polygons
    calc_nbs();
    std::list<SolutionChain> wires, holes;
    find_polygons(wires, holes);

    if (wires.empty()) {
        std::cerr << "[Error] No polygons built!" << std::endl;
        return false;
    }

    remove_pointwise_self_intersections(wires);
    remove_pointwise_self_intersections(holes);

    // create occ geometries
    std::map<Point2D *, TopoDS_Vertex> occ_vertices;
    std::map<Edge2D *, TopoDS_Edge> occ_edges;
    std::map<SolutionChain *, TopoDS_Wire> occ_wires;
    occ_instances(wires, holes, occ_vertices, occ_edges, occ_wires);

    if (occ_wires.empty()) {
        std::cerr << "[Error] No faces built!" << std::endl;
        return false;
    }

    auto WireHoles = detect_holes_in_polygons(wires, holes); // find holes in wires

    // remove "holes in holes" of a wire and merge adjacent holes of a wire
    std::list<SolutionChain> merged_holes;
    for (auto &p: wires) {
        remove_stacked_holes(p, WireHoles);
        merge_adjacent_holes(p, WireHoles, merged_holes, occ_edges, occ_wires);
    }

    TopoDS_ListOfShape faces = occ_faces(wires, occ_wires, WireHoles); // build faces

    if (faces.IsEmpty()) {
        std::cerr << "[Error] No faces built!" << std::endl;
        return false;
    }

    if (faces.Size() == 1) {
        bool b = Kernel::face_normal(TopoDS::Face(faces.First())).Angle(n1) > 1.5;
        shp = b ? faces.First().Reversed() : faces.First();
    } else {
        BRepBuilderAPI_Sewing sew;
        for (auto &F: faces) {
            if (Kernel::face_normal(TopoDS::Face(F)).Angle(n1) > 1.5)
                F.Reverse();
            sew.Add(F);
        }
        sew.Perform();
        shp = sew.SewedShape();

        if (Topo(shp).faces().Size() == 1)
            std::cerr << "[Info] Only one face built. " << faces.Size() << std::endl;

        if (Topo(shp).shells().Size() != 1 && Topo(shp).faces().Size() > 1)
            std::cerr << "[Error] More or less than 1 shell built!" << std::endl;
    }

    // check normals
    gp_Dir n = Kernel::face_normal(TopoDS::Face(Topo(shp).faces().First()));
    if (n.Angle(n1) > 1.5)
        shp.Reverse();

    return true;
}

gp_Pnt Topo2D::reproject(const Point2D &P_2D) {

    std::vector<double> P_3D(3);
    double t;

    if (proj_axis == 0) {
        t = n1.X();
        P_3D[0] = 0.0;
        P_3D[1] = P_2D.X;
        P_3D[2] = P_2D.Y;
    } else if (proj_axis == 1) {
        t = n1.Y();
        P_3D[0] = P_2D.X;
        P_3D[1] = 0.0;
        P_3D[2] = P_2D.Y;
    } else {
        t = n1.Z();
        P_3D[0] = P_2D.X;
        P_3D[1] = P_2D.Y;
        P_3D[2] = 0.0;
    }

    double z = plane_a;
    z -= P_3D[0] * n1.X();
    z -= P_3D[1] * n1.Y();
    z -= P_3D[2] * n1.Z();
    z /= t;
    P_3D[proj_axis] = z;
    return {P_3D[0], P_3D[1], P_3D[2]};
}

void Topo2D::Add(const std::vector<gp_Pnt2d> &wire, unsigned int wire_id, typeEdge2d type) {

    for (unsigned int i = 0; i < wire.size() - 1; i++) {
        const auto &P1 = wire[i];
        const auto &P2 = wire[i + 1];
        add(P1, P2, wire_id, type);
    }
    add(wire.back(), wire.front(), wire_id, type);
}

void Topo2D::add(const gp_Pnt2d &P1, const gp_Pnt2d &P2, unsigned int wire_id, typeEdge2d type) {

    vertices.emplace_back(P1.X(), P1.Y());
    auto v1 = &vertices.back();

    vertices.emplace_back(P2.X(), P2.Y());
    auto v2 = &vertices.back();

    edges.emplace_back(v1, v2, wire_id, type);
    auto e = &edges.back();

    v1->edges.insert(e);
    v2->edges.insert(e);
}