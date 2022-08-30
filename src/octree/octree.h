/*

Header for octree class.

- Ralf-Peter Mundani, Ralf-Peter.Mundani@fhgr.ch

Generate voxel model from geometric surface model.

-----------------------------------------------------------------------

Copyright (c) 2010 Technische Universit�t M�nchen.  All rights reserved.   
  
Permission to use, copy, modify and distribute this software and its   
documentation for any purpose is hereby granted without fee, provided   
that the above copyright notice and this permission notice appear in   
all copies of this software and that you do not sell the software.   
  
THE SOFTWARE IS PROVIDED "AS IS" AND WITHOUT WARRANTY OF ANY KIND,   
EXPRESS, IMPLIED OR OTHERWISE, INCLUDING WITHOUT LIMITATION, ANY   
WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.

*/

#ifndef    __OCTREE_H__
#define    __OCTREE_H__

#define    OCTREE_COLOR_WHITE           0x01
#define    OCTREE_COLOR_BLACK           0x02
#define    OCTREE_COLOR_GREY            0x04
#define    OCTREE_PREDICTED_WHITE       0x10
#define    OCTREE_PREDICTED_BLACK       0x20

#define    OCTREE_SIDE_LEFT             0x00000001
#define    OCTREE_SIDE_RIGHT            0x00000002
#define    OCTREE_SIDE_BOTTOM           0x00000010
#define    OCTREE_SIDE_TOP              0x00000020
#define    OCTREE_SIDE_BACK             0x00000100
#define    OCTREE_SIDE_FRONT            0x00000200

#define    OCTREE_DONE                    0
#define    OCTREE_ERROR                -1

#define    OCTREE_WEST            0
#define    OCTREE_EAST            1
#define    OCTREE_SOUTH           2
#define    OCTREE_NORTH           3
#define    OCTREE_BOTTOM          4
#define    OCTREE_TOP             5

/* voxel element */
typedef struct Voxel {
    double x, y, z;                // voxel midpoint
    double minx, miny, minz;    // minimal coordinates of voxel
    double maxx, maxy, maxz;    // maximal coordinates of voxel
    char color;                    // voxel attribute
    int nface;                    // amount of faces intersecting voxel
    int *findex;                // face indices
    std::bitset<64> mcode;      // voxel morton code

} Voxel;

/* voxel list */
typedef struct llist {
    Voxel *data;                // voxel data
    llist *next;                // next list element
} llist;


class octree {
public:
    /* default constructor/destructor */
    octree();

    ~octree();

    /* set geometric model */
    void set_model(model *Mod);

    /* set filename to write octree to */
    void set_filenames(char *voxfile, char *vtkfile);

    /* generate octree from geometric model */
    int gen_octree(int maxdepth);

    /* count amount of black, white and grey voxels */
    void count(long *black, long *white, long *grey);

    std::vector<std::pair<std::array<double, 3>, std::set<std::string>>> get_spaces_guids();

    std::vector<std::pair<std::array<double, 3>, std::set<unsigned int>>> get_spaces_triangles();

    std::vector<std::set<unsigned int>> get_zones();

    void set_maximum_triangles_in_a_zone(unsigned int n);

    void set_fluid_points(const std::set<std::array<double, 3>> &_fluid_points);

private:
    model *Mod;        // geometric model
    llist *in;        // stack containing voxels to be processed
    llist *out;        // stack containing voxels already finished
    double point_in_polygon_treshold;

    char *voxfile; // filename to write octree in voxel format
    char *vtkfile; // filename to write octree in vtk format

    std::vector<std::pair<std::array<double, 3>, std::set<std::string>>> spaces; // guids of products in spaces and point in air
    std::vector<std::pair<std::array<double, 3>, std::set<unsigned int>>> spaces_tri; // triangle indices of products in spaces and point in air
    std::vector<std::set<unsigned int>> zones;
    std::vector<double> spaces_minima;
    unsigned int max_triangles_in_zone;
    std::list<std::array<double, 3>> fluid_points;

    /* remove first element from stack */
    Voxel *__pop(llist *list);

    /* append element at head of stack */
    void __push(llist *list, Voxel *voxel);

    /* subdivide given voxel into eight children */
    void __subdivide_voxel(Voxel *voxel);

    /* generate child voxel for given midpoint and radius */
    Voxel *__generate_child_voxel(const Voxel *parent, const double *xyz, double radius, std::bitset<64> mcode_append);

    /* test face's bounding box for intersection with voxel */
    int __bbox_test(const Voxel *voxel, int idx);

    /* test face's vertices for being included in voxel */
    int __vertex_in_voxel_test(const Voxel *voxel, int idx);

    /* test voxel's vertices for orientation according to face's normal */
    int __voxel_orientation_test(const Voxel *voxel, int idx);

    /* test triangle for intersection with voxel */
    int __voxel_intersection_test(const Voxel *voxel, int idx);

    /* check triangle's edges for intersection with voxel's sides */
    int __line_intersects_YZ_plane(const Voxel *voxel, int idx, int side);

    int __line_intersects_XZ_plane(const Voxel *voxel, int idx, int side);

    int __line_intersects_XY_plane(const Voxel *voxel, int idx, int side);

    /* check voxel's diagonals for intersection with triangle */
    int __diagonal_intersects_polygon(const Voxel *voxel, int idx);

    /* check point S for lying inside triangle */
    int __point_in_polygon(double SX, double SY, double SZ, int idx);

    /* compute polygon's area via 1/2*|v1v2 x v1v3| */
    double __comp_polygon_area(double v1x, double v1y, double v1z,
                               double v2x, double v2y, double v2z,
                               double v3x, double v3y, double v3z);

    /* compute face's bounding box */
    void __comp_bbox();

    /* compute face's distance to origin */
    void __comp_dist();

    /* check given vertex for face's bounding box */
    void __check_vertex(int idx, double x, double y, double z);

    /* setup and initialise stack */
    void __setup_stack();

    /* cleanup (i.e. full deletion) stack */
    void __cleanup_stack(llist *list);

    /* generate uniform grid representation from octree model */
    void __generate_uniform_grid(int maxdepth);

    /* write out VTK file of voxel model */
    void __write_vox_file(const char *filename, char ***grid, const int &size);

    void __write_vtk_file(const char *filename, char ***grid, const int &size);

    void __write_adaptive_vtk_file(const char *filename);

    /* Eric neighbor search and flood fill */
    void reset_node_colors(char clr);

    void find_facade_triangles();

    int __flood_fill(std::unordered_map<std::bitset<64>, std::vector<std::bitset<64 >>

    > &hash_nb_leafs, std::unordered_map<std::bitset<64>, Voxel *> &hash_octree);

    int __get_depth_by_id(std::bitset<64> s);

    int __get_length_of_id(std::bitset<64> s);

    bool __isBoundary(std::bitset<64> s, unsigned int k);

    bool isLeaf(std::bitset<64> &mcode, std::unordered_map<std::bitset<64>, std::vector<std::bitset<64>>

    > &hash_nb_leafs);

    std::bitset<64> find_sibling_upwards(std::bitset<64> mcode, unsigned int k, unsigned int n);

    std::bitset<64> find_flood_fill_mcode_start(const std::unordered_map<std::bitset<64>, Voxel *

    > &hash_octree, std::unordered_map<std::bitset<64>, std::vector<std::bitset<64 >>> &hash_nb_leafs);

    std::unordered_map<std::bitset<64>, std::vector<std::bitset<64 >>>

    __neighbour_search();

    std::unordered_map<std::bitset<64>, Voxel *>

    create_octree_hash_map();

    std::unordered_map<std::bitset<64>, std::vector<std::bitset<64 >>>

    create_nb_hash_map();

    std::string bitset_to_clean_string(std::bitset<64> b);

    std::vector<std::bitset<64>> generate_mcode_children(std::bitset<64> mcode);

    std::vector<std::bitset<64>> generate_mcode_children_with_desired_value(const std::bitset<64> &mcode, const unsigned int &k, const unsigned int &b);

    std::vector<std::bitset<64>> get_all_boundary__leaf_nodes();

    /* test triangle for intersection with voxel */
    int __isinNode(const Voxel *voxel, int idx);

    /* exctract outer faces of facade voxels */
    void __get_facade_brep(std::unordered_map<std::bitset<64>, std::vector<std::bitset<64 >>

    > &hash_nb_leafs, std::unordered_map<std::bitset<64>, Voxel *> &hash_octree);

    /* compute direction of neighbour */
    int __get_neighbour_type(std::bitset<64> mcode, std::bitset<64> mcode_nb);

    /* append boundary faces to vectors */
    void __switch_neighbour_type(const int &nb_type, std::vector<std::array<double, 3>> &vertices, std::vector<std::array<size_t, 4>> &rectangles, llist *tmp);

    /* find duplicate vertex entries */
    // std::vector <std::array<double, 3>> create_unique_and_inverse_list_of_vertices(const std::vector <std::array<double, 3>> &vertices, std::vector<unsigned int>& inverse, const double &squared_tol);
    // void update_faces_from_inverse(std::vector <std::array<unsigned int, 4>> &faces, const std::vector<unsigned int> &inverse);

    void get_spaces(std::unordered_map<std::bitset<64>, std::vector<std::bitset<64 >>> &hash_nb_leafs, std::unordered_map<std::bitset<64>, Voxel *> &hash_octree);

    std::vector<std::list<unsigned int>> find_space_neighbours();

    void find_zones(const std::vector<std::list<unsigned int>> &nbs, unsigned int max_tri);

    bool is_point_in_voxel(const Voxel *voxel, const std::array<double, 3> &P);

    void set_non_empty_white_voxel_to_grey();

    void save_voxel_triangle_info_in_sets(const Voxel *voxel, std::set<unsigned int> &triangles, std::set<std::string> &guids);

    void save_voxel_neighbors_triangle_info_in_sets(const Voxel *voxel, std::set<unsigned int> &triangles, std::set<std::string> &guids, std::unordered_map<std::bitset<64>, std::vector<std::bitset<64>>> &hash_nb_leafs, std::unordered_map<std::bitset<64>, Voxel *> &hash_octree);

    bool voxel_contains_fluid_point(const Voxel *voxel);

    bool isZero(const double v);
};

#endif    /* !__OCTREE_H__ */