// Copyright 2022 Eric Fichter
#ifndef VIEWER_H
#define VIEWER_H

#include "headers.h"

class Viewer {

public:

    static void visualize_products(const std::list<Product> &products);

    static void visualize_products_openings(const std::list<Product> &products);

    static void visualize_orig_faces(const std::list<oFace> &orig_faces);

    static void visualize_orig_faces_fixed_normals(const std::list<oFace> &orig_faces);

    static void visualize_orig_faces_virtual_product(const std::list<oFace> &orig_faces);

    static void visualize_orig_faces_openings(const std::list<oFace> &orig_faces);

    static void visualize_cFaces(const std::list<cFace> &cFaces);

    static void visualize_face_graph(std::list<cFace> &cFaces);

    static void visualize_cFaces_fixed_normals(const std::list<cFace> &cFaces);

    static void visualize_cFaces_coplanar_and_inner(const std::list<cFace> &cFaces);

    static void visualize_cFaces_hanging(const std::list<cFace> &cFaces);

    static void visualize_cFaces_materials(const std::list<cFace> &cFaces);

    static void visualize_cFaces_openings(const std::list<cFace> &cFaces);

    static void visualize_cFaces_non_openings(const std::list<cFace> &cFaces);

    static void visualize_cFaces_holes(const std::list<cFace> &cFaces);

    static void visualize_cFaces_physical_virtual(const std::list<cFace> &cFaces);

    static void visualize_cFaces_internal_external(const std::list<cFace> &cFaces);

    static void visualize_cFaces_type(const std::list<cFace> &cFaces);

    static void visualize_cFaces_and_inners(const std::list<cFace> &cFaces, const std::list<cFace> &inners);

    static void visualize_cFace_and_shape(const cFace &cface, const TopoDS_Shape &S);

    static void visualize_cFaces_bad_corresponding(std::list<cFace> &cFaces);

    static void visualize_shadings(const std::list<cFace> &cFaces, const std::unordered_map<std::string, std::list<TopoDS_Face>>& shadings);

    static void visualize_cFace(const cFace &cface);

    static void visualize_cFace(const cFace &cface, const std::list<cFace> &cFaces);

    static void visualize_collisions(const std::set<std::pair<Product *, Product *>> &collisions);

    static void visualize_cFaces_shadings(const std::list<cFace> &cFaces);

    static void visualize_external_walls(const std::list<cFace> &cFaces);

    static void visualize_spaces(const std::list<Space> &spaces, bool show_shells, bool firstLvl);

    static void visualize_spaces_facade(const std::list<Space> &spaces, bool show_shells, bool firstLvl);

    static void visualize_cFaces_as_space_boundaries(const std::list<cFace> &cFaces);

    static void visualize_spaces_to_delete(std::list<Space> &spaces, std::list<cFace> &cFaces, bool show_shells, const std::set<Space *> &del);

    static void visualize_zones(const std::vector<std::vector<std::pair<gp_Pnt, std::set<oFace *>>>> &zones);

    static void visualize_shape(const TopoDS_Shape &S);

    static void visualize_shapelist(const TopoDS_ListOfShape &L);

    static void visualize_shape_triangulation(const TopoDS_Shape &S);

    static void visualize_component_search(const std::list<cFace> &cFaces, const std::set<std::pair<iFace *, used_orientation>> &L, const std::pair<iFace *, used_orientation> &oface);

    static void
    visualize_octree_products(std::list<oFace> &orig_faces, const std::vector<std::pair<std::array<double, 3>, std::set<std::string>>> &spaces_guids, const std::vector<std::pair<std::array<double, 3>, std::set<unsigned int>>> &spaces_triangles, const std::vector<std::set<unsigned int>> &zones,
                              const std::vector<oFace *> &tri2orig);

    static void visualize_polygon(const std::vector<std::vector<gp_Pnt2d>> &ClipWires, const std::vector<std::vector<gp_Pnt2d>> &Holes2D);

    static void visualize_clip(const ClipperLib::Paths &poly1, const ClipperLib::Paths &poly2, const ClipperLib::Paths &sol);

    static TopoDS_Face visualize_clip(const ClipperLib::Paths &poly, double z);

    static TopoDS_Face visualize_clip_sol(const ClipperLib::Paths &poly, double z);

    static void visualize_sFaces(const std::list<sFace> &sFaces);
};

#endif //VIEWER_H
