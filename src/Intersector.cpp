// Copyright 2022 Eric Fichter
#include "Intersector.h"

Intersector::Intersector(std::list<cFace> &cFaces, unsigned int _num_threads) : num_threads(_num_threads) {

    pairs.clear();

    for (auto &cface: cFaces) {
        faces.emplace_back();
        faces.back().cface = &cface;
    }

    // triangulate cfaces
    triangulate();

    // generate cgal triangles -> face list and triangle list
    setup_containers();

    // find neighbours of triangles
    setup_pairs();

    // intersect triangles
    intersect();
}

void Intersector::triangulate() {
    std::vector<IsctFace *> V;
    V.reserve(faces.size());
    for (auto &f: faces)
        V.push_back(&f);

#pragma omp parallel for default(none) shared(V) num_threads(num_threads)
    for (unsigned int i = 0; i < V.size(); i++) {
        cFace *cface = V[i]->cface;
        if (!cface->IsPolygon()) continue;
        BRepMesh_IncrementalMesh(V[i]->cface->face, 5.0, false, 5.0, false);
    }
}

std::list<CGAL_Triangle> Intersector::face_to_triangles(cFace *cface) {

    double tol = 1.0e-9;
    std::list<CGAL_Triangle> tris;
    TopoDS_Face Face = cface->face;
    TopLoc_Location L;
    Handle_Poly_Triangulation Poly_Triangulation = BRep_Tool::Triangulation(Face, L);

    if (Poly_Triangulation.IsNull()) {
        std::cerr << "[Warning] Null triangulation." << std::endl;
        return tris;
    }

    const Poly_Array1OfTriangle &Poly_Array1OfTriangle = Poly_Triangulation->Triangles();

    for (int i = 1; i <= Poly_Triangulation->NbTriangles(); ++i) {

        auto Poly_Triangle = Poly_Array1OfTriangle.Value(i);

        int i1, i2, i3;
        Poly_Triangle.Get(i1, i2, i3);

        gp_Pnt p1 = Poly_Triangulation->Node(i1).Transformed(L.Transformation());
        gp_Pnt p2 = Poly_Triangulation->Node(i2).Transformed(L.Transformation());
        gp_Pnt p3 = Poly_Triangulation->Node(i3).Transformed(L.Transformation());

        if (p1.IsEqual(p2, tol) || p2.IsEqual(p3, tol) || p3.IsEqual(p1, tol)) {
            std::cerr << "[Warning] Overlapping points." << std::endl;
            continue;
        }

        if (Kernel::are_points_colinear(p1, p2, p3, 1.0e-5)) {
            std::cerr << "[Warning] Colinear points." << std::endl;
            continue;
        }

        std::vector<CGAL_Point> Pnts(3);

        Pnts[0] = CGAL_Point(p1.X(), p1.Y(), p1.Z());
        Pnts[1] = CGAL_Point(p2.X(), p2.Y(), p2.Z());
        Pnts[2] = CGAL_Point(p3.X(), p3.Y(), p3.Z());

        CGAL_Triangle triangle(Pnts[0], Pnts[1], Pnts[2]);

        if (!triangle.is_degenerate())
            tris.push_back(triangle);
    }

    return tris;
}

void Intersector::setup_containers() {
    for (auto &face: faces) {

        auto tris = face_to_triangles(face.cface);
        if (tris.empty()) continue;

        for (auto &tri: tris) {
            triangles.emplace_back();
            triangles.back().T = tri;

            triangles.back().face = &faces.back();
            faces.back().triangles.push_back(&triangles.back());
        }
    }
}

void Intersector::setup_pairs() {

    rtree_lib::RTree<IsctTriangle *, double, 3, double> tree;
    double o = 1.0e-6;

    // populate rtree
    for (auto &F: triangles) {
        auto bnd = F.T.bbox();
        double min[3] = {bnd.xmin() - o, bnd.ymin() - o, bnd.zmin() - o};
        double max[3] = {bnd.xmax() + o, bnd.ymax() + o, bnd.zmax() + o};
        tree.Insert(min, max, &F);
    }

    // get neighbours
    std::set<std::set<IsctTriangle *>> M;

    for (auto &F: triangles) {
        auto bnd = F.T.bbox();
        double min[3] = {bnd.xmin() - o, bnd.ymin() - o, bnd.zmin() - o};
        double max[3] = {bnd.xmax() + o, bnd.ymax() + o, bnd.zmax() + o};

        std::set<IsctTriangle *> N;

        tree.Search(min, max, [&N, &F](IsctTriangle *found) {
            if (F.face != found->face) N.insert(found);
            return true;
        });

        for (auto &c: N) {
            std::set<IsctTriangle *> s = {&F, c};
            M.insert(s);
        }
    }

    for (auto &m: M) {
        IsctPair p;
        p.T1 = *m.begin();
        p.T2 = *m.end();
        pairs.push_back(p);
    }
}

void Intersector::intersect() {

#pragma omp parallel for default(none) shared(pairs, std::cerr) num_threads(num_threads)
    for (unsigned int i = 0; i < pairs.size(); i++) {

        CGAL_Triangle T1 = pairs[i].T1->T;
        CGAL_Triangle T2 = pairs[i].T2->T;

        if (!CGAL::do_intersect(T1, T2)) continue;

        auto result = CGAL::intersection(T1, T2);

        if (!result) continue;

        if (CGAL_Point *p = boost::get<CGAL_Point>(&*result)) {
            // CGAL_Point P = *p;

            std::cerr << "[Warning] Point not handled." << std::endl;

        } else if (CGAL_Segment *s = boost::get<CGAL_Segment>(&*result)) {
            CGAL_Segment S = *s;

            CGAL_Point p1 = S.point(0);
            CGAL_Point p2 = S.point(1);

            pairs[i].iscts.push_back(std::make_pair(p1, p2));

        } else if (CGAL_Triangle *t = boost::get<CGAL_Triangle>(&*result)) {
            CGAL_Triangle T = *t;

            CGAL_Point p1 = T.vertex(0);
            CGAL_Point p2 = T.vertex(1);
            CGAL_Point p3 = T.vertex(2);

            pairs[i].iscts.push_back(std::make_pair(p1, p2));
            pairs[i].iscts.push_back(std::make_pair(p2, p3));
            pairs[i].iscts.push_back(std::make_pair(p3, p1));

        } else if (std::vector<CGAL_Point> *v = boost::get<std::vector<CGAL_Point>>(&*result)) {

            std::cerr << "[Warning] Vector of points not handled." << std::endl;

        } else {

            std::cerr << "[Warning] Unknown intersection." << std::endl;
        }
    }
}