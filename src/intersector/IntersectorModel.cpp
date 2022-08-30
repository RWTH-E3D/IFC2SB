// Copyright 2022 Eric Fichter
#include "IntersectorModel.h"
#include "../extern/stl_reader.h"

//=============================================================================
// @ARG: none
//=============================================================================
IntersectorModel::IntersectorModel() {
    nvertex = nface = nedge = 0;
    vlist = nullptr;
    elist = nullptr;
    flist = nullptr;
    minx = maxx = miny = maxy = minz = maxz = 0;
    bbox_diffx = bbox_diffy = bbox_diffz = bbox_norm = 0;
}

//=============================================================================
// @ARG: none
//=============================================================================
IntersectorModel::~IntersectorModel() {
    if (vlist != nullptr) {
        delete[] vlist;
        vlist = nullptr;
    }
    if (elist != nullptr) {
        delete[] elist;
        elist = nullptr;
    }
    if (flist != nullptr) {
        delete[] flist;
        flist = nullptr;
    }
}

//=============================================================================
// @ARG: name of TRI file to load
//=============================================================================
void IntersectorModel::load_model(const char *filename, int oct_depth) {

    std::string s = filename;
    s = s.substr(s.length() - 3);
    if (s != "stl" && s != "STL") {
        std::cerr << "Filetype " << s << " is unknown. Abort." << std::endl;
        abort();
    }

    // read in geometry
    stl_reader::StlMesh<float, unsigned int> mesh(filename);

    /* allocate memory for vertex and face list */
    nvertex = mesh.num_vrts();
    nface = mesh.num_tris();
    vlist = new IsctVertex[nvertex];
    flist = new IsctFace[nface];

#ifdef PRINT
    std::cout << "\n... Number of vertices: " << nvertex << ", number of faces: " << nface << std::endl;
#endif

    // store vertices
    for (size_t idx_v = 0; idx_v < nvertex; ++idx_v) {
        const float *pnt = mesh.vrt_coords(idx_v);
        vlist[idx_v].x = pnt[0];
        vlist[idx_v].y = pnt[1];
        vlist[idx_v].z = pnt[2];
    }

    // store faces
    for (size_t idx_f = 0; idx_f < nface; ++idx_f) {
        const unsigned int *face = mesh.tri_corner_inds(idx_f);
        flist[idx_f].v1 = face[0];
        flist[idx_f].v2 = face[1];
        flist[idx_f].v3 = face[2];

        flist[idx_f].nx = flist[idx_f].ny = flist[idx_f].nz = 0;
        flist[idx_f].minx = flist[idx_f].miny = flist[idx_f].minz = flist[idx_f].maxx = flist[idx_f].maxy = flist[idx_f].maxz = 0;
#ifndef MOELLER_INTERSECTION_TEST
        flist[idx_f].A = 0;
        flist[idx_f].dist = 0;
#endif

        // add attribute(s)
        int val = 0;
        flist[idx_f].id = val;
    }

    /* compute face normals */
    __calc_normals();

    /* compute model's bounding box and normalise model to unit cube */
    __comp_BBOX();

    // minor margin, otherwise triangles at border are ignored
    double d = 0.01;
    minx -= d;
    maxx += d;
    miny -= d;
    maxy += d;
    minz -= d;
    maxz += d;

#ifdef PRINT
    double Lx = maxx - minx;
    double Ly = maxy - miny;
    double Lz = maxz - minz;
    double dmin = std::max(std::max(Lx, Ly), Lz) / pow(2, oct_depth); // Minimal voxel size
    std::cout << "... Original geometry data bounding box\t\t: " << Lx << " (" << minx << " -> " << maxx << ") | " << Ly << " (" << miny << " -> " << maxy << ") | " << Lz << " (" << minz << " -> " << maxz << ")" << std::endl;
    std::cout << "... Minimal voxel size\t: " << dmin << std::endl;
#endif

    __normalise_model();
    __comp_BBOX();

#ifdef PRINT
    std::cout << "... Normalized geometry data bounding box\t: " << minx << " " << maxx << " | " << miny << " " << maxy << " | " << minz << " " << maxz << "\tNorm: " << bbox_norm << std::endl;
#endif

    /* compute edges for wireframe display */
    //__comp_edges();

    /* compute face bounding boxes */
    __comp_bbox();
}

void IntersectorModel::load_model(const std::list<std::array<double, 3>> &vert, const std::list<std::array<unsigned int, 3>> &tris, const std::list<unsigned int> &attrs, int oct_depth) {

    /* allocate memory for vertex and face list */
    nvertex = vert.size();
    nface = tris.size();
    vlist = new IsctVertex[nvertex];
    flist = new IsctFace[nface];

#ifdef PRINT
    std::cout << "\n... Number of vertices: " << nvertex << ", number of faces: " << nface << std::endl;
#endif

    // store vertices
    unsigned int idx_v = 0;
    for (const auto &v: vert) {
        vlist[idx_v].x = v[0];
        vlist[idx_v].y = v[1];
        vlist[idx_v].z = v[2];
        idx_v++;
    }

    // store faces
    unsigned int idx_f = 0;
    for (const auto &t: tris) {
        flist[idx_f].v1 = t[0];
        flist[idx_f].v2 = t[1];
        flist[idx_f].v3 = t[2];

        flist[idx_f].nx = flist[idx_f].ny = flist[idx_f].nz = 0;
        flist[idx_f].minx = flist[idx_f].miny = flist[idx_f].minz = flist[idx_f].maxx = flist[idx_f].maxy = flist[idx_f].maxz = 0;
#ifndef MOELLER_INTERSECTION_TEST
        flist[idx_f].A = 0;
        flist[idx_f].dist = 0;
#endif

        idx_f++;
    }

    // add attribute
    idx_f = 0;
    for (const auto &a: attrs) {
        flist[idx_f].id = a;
        idx_f++;
    }

    // for (unsigned int i = 0; i< nface; i++)
    //    std::cerr << "MOD\t" <<flist[i].v1 << "\t" << flist[i].v2 << "\t" << flist[i].v3 << "\t" << flist[i].id << "\t" <<std::endl;

    /* compute face normals */
    __calc_normals();

    /* compute model's bounding box and normalise model to unit cube */
    __comp_BBOX();

    // minor margin, otherwise triangles at border are ignored
    double d = 0.01;
    minx -= d;
    maxx += d;
    miny -= d;
    maxy += d;
    minz -= d;
    maxz += d;

#ifdef PRINT
    double Lx = maxx - minx;
    double Ly = maxy - miny;
    double Lz = maxz - minz;
    double dmin = std::max(std::max(Lx, Ly), Lz) / pow(2, oct_depth); // Minimal voxel size
    std::cout << "... Original geometry data bounding box\t\t: " << Lx << " (" << minx << " -> " << maxx << ") | " << Ly << " (" << miny << " -> " << maxy << ") | " << Lz << " (" << minz << " -> " << maxz << ")" << std::endl;
    std::cout << "... Minimal voxel size\t: " << dmin << std::endl;
#endif

    __normalise_model();
    __comp_BBOX();

#ifdef PRINT
    std::cout << "... Normalized geometry data bounding box\t: " << minx << " " << maxx << " | " << miny << " " << maxy << " | " << minz << " " << maxz << "\tNorm: " << bbox_norm << std::endl;
#endif

    /* compute edges for wireframe display */
    //__comp_edges();

    /* compute face bounding boxes */
    __comp_bbox();
}

//=============================================================================
// @ARG: none
//=============================================================================
void IntersectorModel::__calc_normals() {
    int i;                // loop counter
    double nx, ny, nz;    // face normal's components
    double norm;        // normalisation factor

    for (i = 0; i < nface; ++i) {
        __vec_cross(vlist[flist[i].v2].x - vlist[flist[i].v1].x,
                    vlist[flist[i].v2].y - vlist[flist[i].v1].y,
                    vlist[flist[i].v2].z - vlist[flist[i].v1].z,
                    vlist[flist[i].v3].x - vlist[flist[i].v2].x,
                    vlist[flist[i].v3].y - vlist[flist[i].v2].y,
                    vlist[flist[i].v3].z - vlist[flist[i].v2].z,
                    &nx, &ny, &nz);
        norm = sqrt(nx * nx + ny * ny + nz * nz);
        flist[i].nx = nx / norm;
        flist[i].ny = ny / norm;
        flist[i].nz = nz / norm;
    }
}

//=============================================================================
// @ARG: none
//=============================================================================
void IntersectorModel::__normalise_model() {
    double diffx = (minx + maxx) / 2;    // translation difference in x-direction
    double diffy = (miny + maxy) / 2;    // translation difference in y-direction
    double diffz = (minz + maxz) / 2;    // translation difference in z-direction
    double norm = (maxx - minx) / 2;    // normalise factor
    int i;                            // loop counter

    if ((maxy - miny) / 2 > norm) { norm = (maxy - miny) / 2; }
    if ((maxz - minz) / 2 > norm) { norm = (maxz - minz) / 2; }
    for (i = 0; i < nvertex; ++i) {
        vlist[i].x = (vlist[i].x - diffx) / norm;
        vlist[i].y = (vlist[i].y - diffy) / norm;
        vlist[i].z = (vlist[i].z - diffz) / norm;
    }

    bbox_diffx = diffx;
    bbox_diffy = diffy;
    bbox_diffz = diffz;
    bbox_norm = norm;
}

//=============================================================================
// @ARG: vector vec1's components
// @ARG: vector vec2's components
// @ARG: face normal's components
//=============================================================================
void IntersectorModel::__vec_cross(const double vec1x, const double vec1y, const double vec1z,
                                   const double vec2x, const double vec2y, const double vec2z,
                                   double *nx, double *ny, double *nz) {
    *nx = vec1y * vec2z - vec1z * vec2y;
    *ny = vec1z * vec2x - vec1x * vec2z;
    *nz = vec1x * vec2y - vec1y * vec2x;
}

//=============================================================================
// @ARG: none
//=============================================================================
void IntersectorModel::__comp_edges() {
    int i;            // loop counter
    int idx(0);        // index of edge list for inserting next element

    /* allocate memory for edge list and initialise it*/
    nedge = nface * 3 / 2;
    elist = new IsctEdge[nedge];
    for (i = 0; i < nedge; ++i) {
        elist[i].v1 = MODEL_ERROR;
        elist[i].v2 = MODEL_ERROR;
    }

    /* insert edges in edge list */
    for (i = 0; i < nface; ++i) {
        if (__check_edge(flist[i].v1, flist[i].v2) != MODEL_ERROR) {
            elist[idx].v1 = flist[i].v1;
            elist[idx].v2 = flist[i].v2;
            ++idx;
        }
        if (__check_edge(flist[i].v2, flist[i].v3) != MODEL_ERROR) {
            elist[idx].v1 = flist[i].v2;
            elist[idx].v2 = flist[i].v3;
            ++idx;
        }
        if (__check_edge(flist[i].v3, flist[i].v1) != MODEL_ERROR) {
            elist[idx].v1 = flist[i].v3;
            elist[idx].v2 = flist[i].v1;
            ++idx;
        }
    }
}

//=============================================================================
// @ARG: vertex indices delimiting edge
//=============================================================================
int IntersectorModel::__check_edge(const int v1, const int v2) {
    int idx(0);        // index of edge list for checking given vertices

    while (idx < nedge) {
        if (elist[idx].v1 == MODEL_ERROR) { return idx; }
        else if (elist[idx].v1 == v1 || elist[idx].v2 == v1) {
            if (elist[idx].v1 == v2 || elist[idx].v2 == v2) { return MODEL_ERROR; }
        }
        ++idx;
    }
    return MODEL_ERROR;
}

//=============================================================================
// @ARG: none
//=============================================================================
void IntersectorModel::__comp_BBOX() {
    int i;        // loop counter

    minx = vlist[0].x;
    miny = vlist[0].y;
    minz = vlist[0].z;
    maxx = vlist[0].x;
    maxy = vlist[0].y;
    maxz = vlist[0].z;
    for (i = 1; i < nvertex; ++i) {
        if (vlist[i].x < minx) { minx = vlist[i].x; } else if (vlist[i].x > maxx) { maxx = vlist[i].x; }
        if (vlist[i].y < miny) { miny = vlist[i].y; } else if (vlist[i].y > maxy) { maxy = vlist[i].y; }
        if (vlist[i].z < minz) { minz = vlist[i].z; } else if (vlist[i].z > maxz) { maxz = vlist[i].z; }
    }
}

//=============================================================================
// @ARG: none
//=============================================================================
void IntersectorModel::__comp_bbox() {
    int i;        // loop counter

    for (i = 0; i < nface; ++i) {
        flist[i].minx = vlist[flist[i].v1].x;
        flist[i].miny = vlist[flist[i].v1].y;
        flist[i].minz = vlist[flist[i].v1].z;
        flist[i].maxx = vlist[flist[i].v1].x;
        flist[i].maxy = vlist[flist[i].v1].y;
        flist[i].maxz = vlist[flist[i].v1].z;
        __check_vertex(i, vlist[flist[i].v2].x, vlist[flist[i].v2].y, vlist[flist[i].v2].z);
        __check_vertex(i, vlist[flist[i].v3].x, vlist[flist[i].v3].y, vlist[flist[i].v3].z);
    }
}

//=============================================================================
// @ARG: face index
// @ARG: vertex coordinates to be processed
//=============================================================================
void IntersectorModel::__check_vertex(const int idx, const double x, const double y, const double z) {
    if (x < flist[idx].minx) { flist[idx].minx = x; } else if (x > flist[idx].maxx) { flist[idx].maxx = x; }
    if (y < flist[idx].miny) { flist[idx].miny = y; } else if (y > flist[idx].maxy) { flist[idx].maxy = y; }
    if (z < flist[idx].minz) { flist[idx].minz = z; } else if (z > flist[idx].maxz) { flist[idx].maxz = z; }
}