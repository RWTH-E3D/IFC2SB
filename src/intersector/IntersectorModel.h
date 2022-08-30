// Copyright 2022 Eric Fichter
#ifndef IntersectorModel_H
#define IntersectorModel_H

#include "intersector_includes.h"

#define    MODEL_ERROR                    -1

/* vertex list */
typedef struct IsctVertex {
    double x, y, z;                // vertex coordinates
} IsctVertex;

/* edge list */
typedef struct IsctEdge {
    int v1, v2;                    // vertex index list
} IsctEdge;

/* face list */
typedef struct IsctFace {
    unsigned int v1, v2, v3;                // vertex index list
    double nx, ny, nz;            // face normal
    double minx, miny, minz;    // minimal coordinates of face's bounding box
    double maxx, maxy, maxz;    // maximal coordinates of face's bounding box
    unsigned int id;
#ifndef MOELLER_INTERSECTION_TEST
    double A, dist;
#endif
} IsctFace;


class IntersectorModel {
public:
    /* default constructor/destructor */
    IntersectorModel();

    ~IntersectorModel();

    /* load object from stl file */
    void load_model(const char *filename, int oct_depth);

    void load_model(const std::list<std::array<double, 3>> &vert, const std::list<std::array<unsigned int, 3>> &tris, const std::list<unsigned int> &attrs, int oct_depth);

    int nface;            // amount of faces
    IsctVertex *vlist;        // vertex list
    IsctFace *flist;        // face list
    double minx, maxx;    // model's bounding box minimal/maximal x values
    double miny, maxy;    // model's bounding box minimal/maximal y values
    double minz, maxz;    // model's bounding box minimal/maximal z values
    double bbox_diffx;  // translation difference of geometry in x-direction
    double bbox_diffy;  // translation difference of geometry in y-direction
    double bbox_diffz;  // translation difference of geometry in z-direction
    double bbox_norm;   // normalise factor of geometry

private:
    int nvertex;        // amount of vertices
    int nedge;            // amount of edges
    IsctEdge *elist;        // edge list

    /* calculate face normals */
    void __calc_normals();

    /* normalise model to unit cube [-1,-1,-1:1,1,1] */
    void __normalise_model();

    /* compute vector cross product */
    void __vec_cross(double vec1x, double vec1y, double vec1z,
                     double vec2x, double vec2y, double vec2z,
                     double *nx, double *ny, double *nz);

    /* compute model's edges for wireframe display */
    void __comp_edges();

    /* check if given edge is already included in edge list */
    int __check_edge(int v1, int v2);

    /* compute model's bounding box */
    void __comp_BBOX();

    /* compute face's bounding box */
    void __comp_bbox();

    /* check given vertex for face's bounding box */
    void __check_vertex(int idx, double x, double y, double z);
};

#endif    //IntersectorModel_H
