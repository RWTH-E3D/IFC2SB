// Copyright 2022 Eric Fichter
#ifndef IntersectorInterface_H
#define IntersectorInterface_H

#include "intersector_includes.h"
#include "IntersectorModel.h"
#include "RayIntersector.h"

class RayIntersector;

class IntersectionRay;

class IntersectorInterface {

public:

    IntersectorInterface(const std::string &filename, const int _oct_depth, const unsigned int _min_triangles, bool _write_vtk);

    IntersectorInterface(const std::list<std::array<double, 3>> &vert, const std::list<std::array<unsigned int, 3>> &tris, const std::list<unsigned int> &attrs, const int _oct_depth, const unsigned int _min_triangles, bool _write_vtk);

    ~IntersectorInterface();

    std::multimap<double, unsigned int> Perform(IntersectionRay ray);

private:
    IntersectorModel *model;    // geometric model
    RayIntersector *intersector;    // voxel model
    const int oct_depth; // max depth
    const unsigned int min_triangles; // minimum triangles for node splitting
    unsigned int num_threads; // number of threads for parallelization
    bool write_vtk;

    void Init(const char *infile);

    void Init(const std::list<std::array<double, 3>> &vert, const std::list<std::array<unsigned int, 3>> &tris, const std::list<unsigned int> &attrs);

    static char *vtk_filename(const char *infile);

};

#endif //IntersectorInterface_H