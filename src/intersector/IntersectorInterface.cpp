// Copyright 2022 Eric Fichter
#include "IntersectorInterface.h"

IntersectorInterface::IntersectorInterface(const std::string &filename, const int _oct_depth, const unsigned int _min_triangles, bool _write_vtk) : oct_depth(_oct_depth), min_triangles(_min_triangles), write_vtk(_write_vtk) {

    intersector = nullptr;
    model = nullptr;
    Init(filename.c_str());
}

IntersectorInterface::IntersectorInterface(const std::list<std::array<double, 3>> &vert, const std::list<std::array<unsigned int, 3>> &tris, const std::list<unsigned int> &attrs, const int _oct_depth, const unsigned int _min_triangles, bool _write_vtk) : oct_depth(_oct_depth),
                                                                                                                                                                                                                                                              min_triangles(_min_triangles),
                                                                                                                                                                                                                                                              write_vtk(_write_vtk) {

    intersector = nullptr;
    model = nullptr;
    Init(vert, tris, attrs);
}

IntersectorInterface::~IntersectorInterface() {
    delete intersector;
    delete model;
};

char *IntersectorInterface::vtk_filename(const char *infile) {
    char *vtkfile = new char[strlen(infile) + 1];
    strcpy(vtkfile, infile);
    vtkfile[strlen(vtkfile) - 3] = 'v';
    vtkfile[strlen(vtkfile) - 2] = 't';
    vtkfile[strlen(vtkfile) - 1] = 'k';
    return vtkfile;
}

void IntersectorInterface::Init(const char *infile) {

    std::ifstream f(infile);
    if (!f.good()) {
        std::cerr << "[Error] File not found!";
        return;
    }

#ifdef PRINT
    std::cout << "Info: Generating octree for file " << infile << " at level depth " << oct_depth << ".\n";
    std::cout << "Info: Reading in model..." << std::flush;
#endif

    model = new IntersectorModel();
    model->load_model(infile, oct_depth);

    intersector = new RayIntersector();
    intersector->set_model(model);
    if (write_vtk) intersector->set_filename(vtk_filename(infile));
    if (intersector->gen_octree(oct_depth, min_triangles) != OCTREE_DONE) { std::cerr << " ERROR >> couldn't generate octree" << std::endl; }

#ifdef PRINT
    long white, grey;
    intersector->count(&white, &grey);
    std::cout << "Info: Number of empty voxels: " << white << ", filled voxels: " << grey << "." << std::endl;
    std::cout << "Info: Generation finished. \n";
#endif
}

void IntersectorInterface::Init(const std::list<std::array<double, 3>> &vert, const std::list<std::array<unsigned int, 3>> &tris, const std::list<unsigned int> &attrs) {

#ifdef PRINT
    std::cout << "Info: Generating octree for data set at level depth " << oct_depth << ".\n";
    std::cout << "Info: Reading in model..." << std::flush;
#endif

    model = new IntersectorModel();
    model->load_model(vert, tris, attrs, oct_depth);

    intersector = new RayIntersector();
    intersector->set_model(model);
    if (write_vtk) intersector->set_filename(vtk_filename("output.vtk"));
    if (intersector->gen_octree(oct_depth, min_triangles) != OCTREE_DONE) { std::cerr << " ERROR >> couldn't generate octree" << std::endl; }

#ifdef PRINT
    long white, grey;
    intersector->count(&white, &grey);
    std::cout << "Info: Number of empty voxels: " << white << ", filled voxels: " << grey << "." << std::endl;
    std::cout << "Info: Generation finished. \n";
#endif
}

std::multimap<double, unsigned int> IntersectorInterface::Perform(IntersectionRay ray) {
    ray.convert_lengths_to_octree_grid(model->bbox_norm, model->bbox_diffx, model->bbox_diffy, model->bbox_diffz);
    return intersector->perform(ray);
}