// Copyright 2022 Eric Fichter
#ifndef INTERSECTOR_H
#define INTERSECTOR_H

#include "headers.h"

class cFace;

class Kernel;

class Intersector {

public:
    Intersector(std::list<cFace> &cFaces, unsigned int _num_threads);

private:

    struct IsctTriangle;

    struct IsctFace {
        cFace *cface;
        std::list<IsctTriangle *> triangles;
    };

    struct IsctTriangle {
        IsctFace *face;
        CGAL_Triangle T;
    };

    struct IsctPair {
        IsctTriangle *T1;
        IsctTriangle *T2;
        std::list<std::pair<CGAL_Point, CGAL_Point>> iscts;
    };

    const unsigned int num_threads;
    std::list<IsctFace> faces;
    std::list<IsctTriangle> triangles;
    std::vector<IsctPair> pairs;

    void triangulate();

    std::list<CGAL_Triangle> face_to_triangles(cFace *cface);

    void setup_containers();

    void setup_pairs();

    void intersect();


};


#endif //INTERSECTOR_H
