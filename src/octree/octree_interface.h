#ifndef OCTREE_INTERFACE_H
#define OCTREE_INTERFACE_H

#include "octree_includes.h"

class octree_interface {

public:

    static int process_stl_file(const std::string &filename, int oct_depth, bool write_vtk);

    static void process_mesh(std::vector<std::pair<std::array<double, 3>, std::set<std::string>>> &spaces_guids, std::vector<std::pair<std::array<double, 3>, std::set<unsigned int>>> &spaces_triangles, std::vector<std::set<unsigned int>> &zones,
                             const std::vector<std::tuple<double, double, double>> &vertices, const std::vector<std::tuple<int, int, int>> &faces, const std::vector<std::string> &attrs, const std::set<std::array<double, 3>> &fluid_points, int oct_depth, bool write_vtk);

};

#endif //OCTREE_INTERFACE_H