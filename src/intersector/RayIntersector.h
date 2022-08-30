// Copyright 2022 Eric Fichter
#ifndef RAY_INTERSECTOR_H
#define RAY_INTERSECTOR_H

#include "intersector_includes.h"
#include "IntersectorModel.h"

#define    OCTREE_COLOR_WHITE           0x01
//#define    OCTREE_COLOR_BLACK           0x02
#define    OCTREE_COLOR_GREY            0x04
//#define    OCTREE_PREDICTED_WHITE       0x10
//#define    OCTREE_PREDICTED_BLACK       0x20

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

//#define    MOELLER_INTERSECTION_TEST

/* voxel element */
typedef struct Voxel2 {
    double x, y, z;                // voxel midpoint
    unsigned int d;                 // voxel's tree depth
    double minx, miny, minz;    // minimal coordinates of voxel
    double maxx, maxy, maxz;    // maximal coordinates of voxel
    char color;                    // voxel attribute
    int nface;                    // amount of faces intersecting voxel
    int *findex;                // face indices
    std::array<Voxel2 *, 8> children; //children of voxel
} Voxel2;

/* voxel list */
typedef struct llist2 {
    Voxel2 *data;                // voxel data
    llist2 *next;                // next list element
} llist2;

#include "IntersectionRay.h"

class IntersectionRay;

class RayIntersector {
public:
    /* default constructor/destructor */
    RayIntersector();

    ~RayIntersector();

    /* set geometric model */
    void set_model(IntersectorModel *Mod);

    /* set filename to write octree to */
    void set_filename(char *vtkfile0);

    /* generate octree from geometric model */
    int gen_octree(int maxdepth, unsigned int min_triangles_for_division);

    /* count amount of black, white and grey voxels */
    void count(long *white, long *grey);

    std::multimap<double, unsigned int> perform(const IntersectionRay &ray);

private:
    IntersectorModel *Mod;        // geometric model
    llist2 *in;        // stack containing voxels to be processed
    llist2 *out;        // stack containing voxels already finished
    Voxel2 *root_voxel; // voxel at depth 0
    double point_in_polygon_treshold;

    char *vtkfile; // filename to write octree in vtk format

    /* remove first element from stack */
    Voxel2 *pop(llist2 *list);

    /* append element at head of stack */
    void push(llist2 *list, Voxel2 *voxel);

    /* subdivide given voxel into eight children */
    void subdivide_voxel(Voxel2 *voxel);

    /* generate child voxel for given midpoint and radius */
    Voxel2 *generate_child_voxel(const Voxel2 *parent, const double *xyz, double radius);

    /* test face's bounding box for intersection with Voxel2 */
    int bbox_test(const Voxel2 *voxel, int idx);

    /* test face's vertices for being included in Voxel2 */
    int vertex_in_voxel_test(const Voxel2 *voxel, int idx);

#ifndef MOELLER_INTERSECTION_TEST
    /* test triangle for intersection with Voxel2 */
    int voxel_intersection_test(const Voxel2 *voxel, int idx);

    /* check triangle's edges for intersection with voxel's sides */
    int line_intersects_YZ_plane(const Voxel2 *voxel, int idx, int side);

    int line_intersects_XZ_plane(const Voxel2 *voxel, int idx, int side);

    int line_intersects_XY_plane(const Voxel2 *voxel, int idx, int side);

    /* check voxel's diagonals for intersection with triangle */
    int diagonal_intersects_polygon(const Voxel2 *voxel, int idx);

    /* check point S for lying inside triangle */
    int point_in_polygon(double SX, double SY, double SZ, int idx);

    /* compute polygon's area via 1/2*|v1v2 x v1v3| */
    double comp_polygon_area(double v1x, double v1y, double v1z,
                             double v2x, double v2y, double v2z,
                             double v3x, double v3y, double v3z);
#endif

    /* compute face's bounding box */
    void comp_bbox();

    /* compute face's distance to origin */
#ifndef MOELLER_INTERSECTION_TEST
    void comp_dist();
#endif

    /* check given vertex for face's bounding box */
    void check_vertex(int idx, double x, double y, double z);

    /* setup and initialise stack */
    void setup_stack();

    /* cleanup (i.e. full deletion) stack */
    void cleanup_stack(llist2 *list);

    /* write out VTK file of voxel model */
    void write_adaptive_vtk_file();

    /* write out VTK file of voxel model */
    void write_adaptive_vtk_file_only_leafs();

    bool isLeaf(Voxel2 *voxel);

    std::list<Voxel2 *> get_intersected_voxels(const IntersectionRay &ray);

    static std::set<unsigned int> get_potential_intersected_triangles(const std::list<Voxel2 *> &voxels);

    void get_intersected_triangles(const IntersectionRay &ray, const std::set<unsigned int> &triangles_to_test, std::multimap<double, unsigned int> &hits);

    static bool isZero(double v);

    /* test triangle for intersection with voxel */
    int triangle_aabb_moeller(const Voxel2 *voxel, int idx);
};

#endif //RAY_INTERSECTOR_H