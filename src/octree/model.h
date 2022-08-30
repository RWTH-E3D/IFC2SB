/*

Header for model class.

- Ralf-Peter Mundani, mundani@mytum.de

Load geometric model data as vertices and faces.

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

#ifndef __MODEL_H__
#define __MODEL_H__

#include "octree_includes.h"

#define    MODEL_ERROR                    -1

/* vertex list */
typedef struct Vertex {
    double x, y, z;                // vertex coordinates
} Vertex;

/* edge list */
typedef struct Edge {
    int v1, v2;                    // vertex index list
} Edge;

/* face list */
typedef struct Face {
    int v1, v2, v3;                // vertex index list
    double nx, ny, nz;            // face normal
    double minx, miny, minz;    // minimal coordinates of face's bounding box
    double maxx, maxy, maxz;    // maximal coordinates of face's bounding box
    double A, dist;
    std::string guid;
} Face;


class model {
public:
    /* default constructor/destructor */
    model();

    ~model();

    /* load object from TRI file */
    void load_model(const char *filename, int oct_depth, bool extend_bbox_for_floodfill);

    void load_model(const std::vector<std::tuple<double, double, double>> &vertices, const std::vector<std::tuple<int, int, int>> &faces, const std::vector<std::string> &attrs, int oct_depth, bool extend_bbox_for_floodfill);

    int nface;            // amount of faces
    Vertex *vlist;        // vertex list
    Face *flist;        // face list
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
    Edge *elist;        // edge list

    /* calculate face normals */
    void __calc_normals();

    /* normalise model to unit cube [-1,-1,-1:1,1,1] */
    void __normalise_model();

    /* compute vector cross product */
    void __vec_cross(const double vec1x, const double vec1y, const double vec1z,
                     const double vec2x, const double vec2y, const double vec2z,
                     double *nx, double *ny, double *nz);

    /* compute model's edges for wireframe display */
    void __comp_edges();

    /* check if given edge is already included in edge list */
    int __check_edge(const int v1, const int v2);

    /* compute model's bounding box */
    void __comp_BBOX();

    /* compute face's bounding box */
    void __comp_bbox();

    /* check given vertex for face's bounding box */
    void __check_vertex(const int idx, const double x, const double y, const double z);
};

#endif    /* !__MODEL_H__ */
