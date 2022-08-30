/*

Generate voxel model from geometric surface model.

- Ralf-Peter Mundani, Ralf-Peter.Mundani@fhgr.ch

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
#include "octree_includes.h"

#define RESIZE_VTK_OUTPUT
//#define VERIFY_NEIGHBOURSEARCH
#define MOELLER_INTERSECTION_TEST
//#define GET_FACADE_TRIANGLES
//#define GET_FACADE_MESH
#define GET_SPACES

//=============================================================================
// @ARG: none
//=============================================================================
octree::octree() {
    in = nullptr;
    out = nullptr;
    Mod = nullptr;
    voxfile = nullptr;
    vtkfile = nullptr;
    point_in_polygon_treshold = 1.0e-9;

    spaces.clear();
    spaces_tri.clear();
    zones.clear();
    spaces_minima.clear();
    max_triangles_in_zone = 0;
    fluid_points.clear();
}

//=============================================================================
// @ARG: none
//=============================================================================
octree::~octree() {
    __cleanup_stack(in);
    __cleanup_stack(out);
    Mod = nullptr;
}

//=============================================================================
// @ARG: none
//=============================================================================
void octree::set_model(model *Mod) { if (Mod != nullptr) { octree::Mod = Mod; }}

//=============================================================================
// @ARG: filename to write octree to (if not set, nothing will be written out)
//=============================================================================
void octree::set_filenames(char *voxfile, char *vtkfile) {
    octree::voxfile = voxfile;
    octree::vtkfile = vtkfile;
}

//=============================================================================
// @ARG: maximum depth of recursion
//=============================================================================
int octree::gen_octree(const int maxdepth) {
    Voxel *voxel;            // voxel to be processed
    int count;                // counter for triangles to be excluded
    int i, j;                // loop counters
    int *findex;            // indices of faces left
    unsigned char color;    // voxel's color from orientation test

#ifdef    __OCTREE_TIMER__
    std::cout << "\tGenerating octree..." << std::flush;
    MTimer timer_octree;
    timer_octree.start();
#endif    /* !__OCTREE_TIMER__ */

    /* check for model and initialise IN stack */
    if (Mod == nullptr) { return OCTREE_ERROR; }
    if (in == nullptr) { __setup_stack(); }

    /* compute for all faces its bbox and distance to origin */
    __comp_bbox();
    __comp_dist();

    while ((voxel = __pop(in)) != nullptr) {
        /* quick bounding box test */
        count = 0;
        for (i = 0; i < voxel->nface; ++i) {
            if (__bbox_test(voxel, voxel->findex[i]) == OCTREE_DONE) {
                voxel->findex[i] |= 0x80000000;
                ++count;
            }
        }

        /* check resulting triangles (if so) for intersection */
        if (count != voxel->nface) {
            for (i = 0; i < voxel->nface; ++i) {
                /* only if voxel has not been marked by __bbox_test() */
                if (!(voxel->findex[i] & 0x80000000)) {
                    /* quick vertex in voxel test */
                    if (__vertex_in_voxel_test(voxel, voxel->findex[i]) == OCTREE_COLOR_GREY) {
                        voxel->color = OCTREE_COLOR_GREY;
                    } else {
                        /* expensive intersection test */
#ifdef MOELLER_INTERSECTION_TEST
                        if (__isinNode(voxel, voxel->findex[i]) == OCTREE_COLOR_GREY) {
#else
                            if (__voxel_intersection_test(voxel, voxel->findex[i]) == OCTREE_COLOR_GREY) {
#endif
                            voxel->color = OCTREE_COLOR_GREY;
                        } else {
                            voxel->findex[i] |= 0x80000000;
                            ++count;
                        }
                    } /* END if(__vertex_in_voxel_test()) */
                } /* END if(!(voxel->findex[i] & 0x80000000)) */
            } /* END for */
        } /* END if(count != voxel->nface) */

        /* determine voxel's color in case all triangles have been discarded */
        if (count == voxel->nface) {
            i = 0;
            while (i < voxel->nface) {
                color = __voxel_orientation_test(voxel, voxel->findex[i] & 0x7fffffff);
                if (color == OCTREE_COLOR_WHITE) {
                    voxel->color = OCTREE_COLOR_WHITE;
                    break;
                } else if (color == OCTREE_COLOR_BLACK) { voxel->color = OCTREE_COLOR_BLACK; }
                ++i;
            }
            if (voxel->color == OCTREE_ERROR) {
                std::cout << "WARNING GENOCT [" << voxel->x << ":" << voxel->y << ":" << voxel->z;
                std::cout << "] (depth: " << __get_depth_by_id(voxel->mcode) << ") >> no voxel color set" << std::endl;
            }
        }

        /* push voxel back to IN or OUT stack */
        if (voxel->color != OCTREE_COLOR_GREY) {
            delete[] voxel->findex;
            voxel->nface = 0;
            voxel->findex = nullptr;
            __push(out, voxel);
        } else {
            findex = new int[voxel->nface - count];
            j = 0;
            for (i = 0; i < voxel->nface; ++i) {
                if (!(voxel->findex[i] & 0x80000000)) {
                    findex[j] = voxel->findex[i];
                    ++j;
                }
            }
            delete[] voxel->findex;
            voxel->nface = j;
            voxel->findex = nullptr;
            voxel->findex = new int[voxel->nface];
            for (i = 0; i < voxel->nface; ++i) { voxel->findex[i] = findex[i]; }
            delete[] findex;
            findex = nullptr;
            if (__get_depth_by_id(voxel->mcode) == maxdepth) { __push(out, voxel); }
            else {
                __subdivide_voxel(voxel);
                //__push(out, voxel);
            }
        }
    } /* END while */

#ifdef    __OCTREE_TIMER__
    timer_octree.stop();
    std::cout << "done. Took " << timer_octree.elapsed_cpu_time() << " s." << std::endl;
#endif    /* !__OCTREE_TIMER__ */

    //*********************************************************************************************************
    std::unordered_map<std::bitset<64>, std::vector<std::bitset<64 >>> hash_nb_leafs = __neighbour_search();
    std::unordered_map<std::bitset<64>, Voxel *> hash_octree = create_octree_hash_map(); // create hash map id->voxel
    //for (auto itr = hash_octree.begin(); itr != hash_octree.end(); itr++) std::cout << itr->first << '\t' << hash_id_leaf[(itr->second)->mcode] << '\n';
    __flood_fill(hash_nb_leafs, hash_octree);

#ifdef    GET_FACADE_TRIANGLES
    find_facade_triangles();     // for all facade nodes, find triangles lying within
#endif

#ifdef    GET_FACADE_MESH
    __get_facade_brep(hash_nb_leafs, hash_octree); // extract surface mesh from facade info
#endif

#ifdef    GET_SPACES
    get_spaces(hash_nb_leafs, hash_octree); // find spaces and return trianlges/objects linked to them
#endif
    // *********************************************************************************************************


    /* write octree as uniform grid to file */
    //__generate_uniform_grid(maxdepth);
    if (vtkfile != nullptr)
        __write_adaptive_vtk_file(vtkfile);

    return OCTREE_DONE;
}

//=============================================================================
// @ARG: linked list for fetching (i.e. deleting) voxel
//=============================================================================
Voxel *octree::__pop(llist *list) {
    Voxel *voxel;    // information of next voxel
    llist *tmp;        // next stack element

    if (list->next == nullptr) { return nullptr; }
    tmp = list->next;
    voxel = tmp->data;
    list->next = tmp->next;
    tmp->data = nullptr;
    tmp->next = nullptr;
    delete tmp;

    return voxel;
}

//=============================================================================
// @ARG: linked list for adding voxel
// @ARG: voxel data
//=============================================================================
void octree::__push(llist *list, Voxel *voxel) {
    llist *tmp;    // new stack element

    tmp = new llist;
    tmp->data = voxel;
    tmp->next = list->next;
    list->next = tmp;
}

//=============================================================================
// @ARG: voxel data
//=============================================================================
void octree::__subdivide_voxel(Voxel *voxel) {
    //double radius = pow(2.0, -(voxel->d));	// children's radius, i.e. 2^(-d)
    double radius = (voxel->maxx - voxel->minx) * 0.25;    // children's radius, i.e. 2^(-d)
    double xyz[3];                            // children's midpoint

    /* compute children's radius */

    /* compute children's midpoint and generate voxel */
    xyz[0] = voxel->x - radius;
    xyz[1] = voxel->y - radius;
    xyz[2] = voxel->z - radius;    /* BNW */
//	if(voxel->d != 1) { __push(in, __generate_child_voxel(voxel, xyz, radius)); }
    __push(in, __generate_child_voxel(voxel, xyz, radius, 0b000));
    xyz[0] = voxel->x + radius;
    xyz[1] = voxel->y - radius;
    xyz[2] = voxel->z - radius;    /* BNE */
    __push(in, __generate_child_voxel(voxel, xyz, radius, 0b100));
    xyz[0] = voxel->x - radius;
    xyz[1] = voxel->y - radius;
    xyz[2] = voxel->z + radius;    /* BSW */
//	if(voxel->d != 1) { __push(in, __generate_child_voxel(voxel, xyz, radius)); }
    __push(in, __generate_child_voxel(voxel, xyz, radius, 0b001));
    xyz[0] = voxel->x + radius;
    xyz[1] = voxel->y - radius;
    xyz[2] = voxel->z + radius;    /* BSE */
//	if(voxel->d != 1) { __push(in, __generate_child_voxel(voxel, xyz, radius)); }
    __push(in, __generate_child_voxel(voxel, xyz, radius, 0b101));
    xyz[0] = voxel->x - radius;
    xyz[1] = voxel->y + radius;
    xyz[2] = voxel->z - radius;    /* TNW */
//	if(voxel->d != 1) { __push(in, __generate_child_voxel(voxel, xyz, radius)); }
    __push(in, __generate_child_voxel(voxel, xyz, radius, 0b010));
    xyz[0] = voxel->x + radius;
    xyz[1] = voxel->y + radius;
    xyz[2] = voxel->z - radius;    /* TNE */
    __push(in, __generate_child_voxel(voxel, xyz, radius, 0b110));
    xyz[0] = voxel->x - radius;
    xyz[1] = voxel->y + radius;
    xyz[2] = voxel->z + radius;    /* TSW */
//	if(voxel->d != 1) { __push(in, __generate_child_voxel(voxel, xyz, radius)); }
    __push(in, __generate_child_voxel(voxel, xyz, radius, 0b011));
    xyz[0] = voxel->x + radius;
    xyz[1] = voxel->y + radius;
    xyz[2] = voxel->z + radius;    /* TSE */
//	if(voxel->d != 1) { __push(in, __generate_child_voxel(voxel, xyz, radius)); }
    __push(in, __generate_child_voxel(voxel, xyz, radius, 0b111));

    /* delete parent voxel */
    delete[] voxel->findex;
    voxel->findex = nullptr;
    delete voxel;
    voxel = nullptr;
}

//=============================================================================
// @ARG: parent data
// @ARG: children's midpoint
// @ARG: children's radius
//=============================================================================
Voxel *octree::__generate_child_voxel(const Voxel *parent, const double *xyz, const double radius, const std::bitset<64> mcode_append) {

    Voxel *child;    // new voxel
    int i;            // loop counter

    child = new Voxel;

    child->mcode = parent->mcode;
    child->mcode = child->mcode << 3; // shift all numbers by three bits to the left
    child->mcode |= mcode_append; // set the last three bits from 000 to child code
    child->x = xyz[0];
    child->y = xyz[1];
    child->z = xyz[2];
    child->minx = xyz[0] - radius;
    child->miny = xyz[1] - radius;
    child->minz = xyz[2] - radius;
    child->maxx = xyz[0] + radius;
    child->maxy = xyz[1] + radius;
    child->maxz = xyz[2] + radius;
    child->color = OCTREE_ERROR;
    child->nface = parent->nface;
    child->findex = new int[child->nface];
    for (i = 0; i < child->nface; ++i) { child->findex[i] = parent->findex[i]; }

    return child;
}

//=============================================================================
// @ARG: voxel data
// @ARG: face index
//=============================================================================
int octree::__bbox_test(const Voxel *voxel, const int idx) {
    if (voxel->minx >= Mod->flist[idx].maxx) { return OCTREE_DONE; }
    if (Mod->flist[idx].minx >= voxel->maxx) { return OCTREE_DONE; }
    if (voxel->miny >= Mod->flist[idx].maxy) { return OCTREE_DONE; }
    if (Mod->flist[idx].miny >= voxel->maxy) { return OCTREE_DONE; }
    if (voxel->minz >= Mod->flist[idx].maxz) { return OCTREE_DONE; }
    if (Mod->flist[idx].minz >= voxel->maxz) { return OCTREE_DONE; }

    return OCTREE_ERROR;
}

//=============================================================================
// @ARG: voxel data
// @ARG: face index
//=============================================================================
int octree::__vertex_in_voxel_test(const Voxel *voxel, const int idx) {
    int vec[3];        // index of vertices to be checked
    int i;            // loop counter

    vec[0] = Mod->flist[idx].v1;
    vec[1] = Mod->flist[idx].v2;
    vec[2] = Mod->flist[idx].v3;

    for (i = 0; i < 3; ++i) {
        if (voxel->minx <= Mod->vlist[vec[i]].x && Mod->vlist[vec[i]].x <= voxel->maxx) {
            if (voxel->miny <= Mod->vlist[vec[i]].y && Mod->vlist[vec[i]].y <= voxel->maxy) {
                if (voxel->minz <= Mod->vlist[vec[i]].z && Mod->vlist[vec[i]].z <= voxel->maxz) {
                    return OCTREE_COLOR_GREY;
                }
            }
        }
    }

    return OCTREE_ERROR;
}

//=============================================================================
// @ARG: voxel data
// @ARG: face index
//=============================================================================
int octree::__voxel_orientation_test(const Voxel *voxel, const int idx) {
    double nx = Mod->flist[idx].nx;            // face normal's x-coordinate
    double ny = Mod->flist[idx].ny;            // face normal's y-coordinate
    double nz = Mod->flist[idx].nz;            // face normal's z-coordinate
    double dist = Mod->flist[idx].dist;        // face's distance from origin
    double dvox;                            // distance of voxel's centre point from face

    /* compute distance of voxel's centre point from face
     * F: nx*x + ny*y + nz*z - dist = 0
     */
    dvox = nx * voxel->x + ny * voxel->y + nz * voxel->z - dist;

    /* check voxel's orientation (above|on|beneath face F) */
    if (dvox > 0) { return OCTREE_COLOR_WHITE; }
    else if (dvox < 0) { return OCTREE_COLOR_BLACK; }

    return OCTREE_ERROR;
}

//=============================================================================
// @ARG: voxel data
// @ARG: face index
//=============================================================================
int octree::__voxel_intersection_test(const Voxel *voxel, const int idx) {
    /* check triangle's edges for intersection with voxel's faces */
    if (__line_intersects_YZ_plane(voxel, idx, OCTREE_SIDE_LEFT) == OCTREE_DONE) { return OCTREE_COLOR_GREY; }
    if (__line_intersects_YZ_plane(voxel, idx, OCTREE_SIDE_RIGHT) == OCTREE_DONE) { return OCTREE_COLOR_GREY; }
    if (__line_intersects_XZ_plane(voxel, idx, OCTREE_SIDE_BOTTOM) == OCTREE_DONE) { return OCTREE_COLOR_GREY; }
    if (__line_intersects_XZ_plane(voxel, idx, OCTREE_SIDE_TOP) == OCTREE_DONE) { return OCTREE_COLOR_GREY; }
    if (__line_intersects_XY_plane(voxel, idx, OCTREE_SIDE_BACK) == OCTREE_DONE) { return OCTREE_COLOR_GREY; }
    if (__line_intersects_XY_plane(voxel, idx, OCTREE_SIDE_FRONT) == OCTREE_DONE) { return OCTREE_COLOR_GREY; }

    /* check voxel's diagonals for intersection with triangle */
    if (__diagonal_intersects_polygon(voxel, idx) == OCTREE_DONE) { return OCTREE_COLOR_GREY; }

    return OCTREE_ERROR;
}

//=============================================================================
// @ARG: voxel data
// @ARG: face index
// @ARG: voxel's left or right side
//=============================================================================
int octree::__line_intersects_YZ_plane(const Voxel *voxel, const int idx, const int side) {
    int vec[3];            // indices of face's vertices
    double dist;        // plane's distance from origin
    double u;            // component of line's direction vector u
    double lambda;        // line's parameter (g: x = a + lambda*u)
    int i;                // loop counter
    double SY, SZ;        // intersection point's coordinates

    vec[0] = Mod->flist[idx].v1;
    vec[1] = Mod->flist[idx].v2;
    vec[2] = Mod->flist[idx].v3;
    if (side == OCTREE_SIDE_LEFT) { dist = voxel->minx; } else { dist = voxel->maxx; }

    for (i = 0; i < 3; ++i) {
        u = Mod->vlist[vec[(i + 1) % 3]].x - Mod->vlist[vec[i]].x;
        if (u != 0) {
            lambda = (dist - Mod->vlist[vec[i]].x) / u;
            if (0 <= lambda && lambda <= 1) {
                SY = Mod->vlist[vec[i]].y + lambda * (Mod->vlist[vec[(i + 1) % 3]].y - Mod->vlist[vec[i]].y);
                SZ = Mod->vlist[vec[i]].z + lambda * (Mod->vlist[vec[(i + 1) % 3]].z - Mod->vlist[vec[i]].z);
                /* S lies truly inside plane */
                if (voxel->miny < SY && SY < voxel->maxy) {
                    if (voxel->minz < SZ && SZ < voxel->maxz) { return OCTREE_DONE; }
                }
            }
        }
    }

    return OCTREE_ERROR;
}

//=============================================================================
// @ARG: voxel data
// @ARG: face index
// @ARG: voxel's top or bottom side
//=============================================================================
int octree::__line_intersects_XZ_plane(const Voxel *voxel, const int idx, const int side) {
    int vec[3];            // indices of face's vertices
    double dist;        // plane's distance from origin
    double u;            // component of line's direction vector u
    double lambda;        // line's parameter (g: x = a + lambda*u)
    int i;                // loop counter
    double SX, SZ;        // intersection point's coordinates

    vec[0] = Mod->flist[idx].v1;
    vec[1] = Mod->flist[idx].v2;
    vec[2] = Mod->flist[idx].v3;
    if (side == OCTREE_SIDE_BOTTOM) { dist = voxel->miny; } else { dist = voxel->maxy; }

    for (i = 0; i < 3; ++i) {
        u = Mod->vlist[vec[(i + 1) % 3]].y - Mod->vlist[vec[i]].y;
        if (u != 0) {
            lambda = (dist - Mod->vlist[vec[i]].y) / u;
            if (0 <= lambda && lambda <= 1) {
                SX = Mod->vlist[vec[i]].x + lambda * (Mod->vlist[vec[(i + 1) % 3]].x - Mod->vlist[vec[i]].x);
                SZ = Mod->vlist[vec[i]].z + lambda * (Mod->vlist[vec[(i + 1) % 3]].z - Mod->vlist[vec[i]].z);
                /* S lies truly inside plane */
                if (voxel->minx < SX && SX < voxel->maxx) {
                    if (voxel->minz < SZ && SZ < voxel->maxz) { return OCTREE_DONE; }
                }
            }
        }
    }

    return OCTREE_ERROR;
}

//=============================================================================
// @ARG: voxel data
// @ARG: face index
// @ARG: voxel's front or back side
//=============================================================================
int octree::__line_intersects_XY_plane(const Voxel *voxel, const int idx, const int side) {
    int vec[3];            // indices of face's vertices
    double dist;        // plane's distance from origin
    double u;            // component of line's direction vector u
    double lambda;        // line's parameter (g: x = a + lambda*u)
    int i;                // loop counter
    double SX, SY;        // intersection point's coordinates

    vec[0] = Mod->flist[idx].v1;
    vec[1] = Mod->flist[idx].v2;
    vec[2] = Mod->flist[idx].v3;
    if (side == OCTREE_SIDE_BACK) { dist = voxel->minz; } else { dist = voxel->maxz; }

    for (i = 0; i < 3; ++i) {
        u = Mod->vlist[vec[(i + 1) % 3]].z - Mod->vlist[vec[i]].z;
        if (u != 0) {
            lambda = (dist - Mod->vlist[vec[i]].z) / u;
            if (0 <= lambda && lambda <= 1) {
                SX = Mod->vlist[vec[i]].x + lambda * (Mod->vlist[vec[(i + 1) % 3]].x - Mod->vlist[vec[i]].x);
                SY = Mod->vlist[vec[i]].y + lambda * (Mod->vlist[vec[(i + 1) % 3]].y - Mod->vlist[vec[i]].y);
                /* S lies truly inside plane */
                if (voxel->minx < SX && SX < voxel->maxx) {
                    if (voxel->miny < SY && SY < voxel->maxy) { return OCTREE_DONE; }
                }
            }
        }
    }

    return OCTREE_ERROR;
}

//=============================================================================
// @ARG: voxel data
// @ARG: face index
//=============================================================================
int octree::__diagonal_intersects_polygon(const Voxel *voxel, const int idx) {
    double lambda;            // diagonal's parameter (g: x = a + lambda*u)
    double dx, dy, dz;        // diagonal's foot point a
    double ux, uy, uz;        // diagonal's direction vector u
    double nu;                // face's normal n * diagonal's direction vector u
    int i;                    // loop counter
    double SX, SY, SZ;        // intersection point's coordinates

    dx = voxel->minx;
    ux = voxel->maxx - dx;
    for (i = 0; i < 4; ++i) {
        switch (i) {
            case 0:    /* BNW -> TSE */
                dy = voxel->miny;
                dz = voxel->minz;
                uy = voxel->maxy - dy;
                uz = voxel->maxz - dz;
                break;
            case 1: /* TNW -> BSE */
                dy = voxel->maxy;
                dz = voxel->minz;
                uy = voxel->miny - dy;
                uz = voxel->maxz - dz;
                break;
            case 2: /* BSW -> TNE */
                dy = voxel->miny;
                dz = voxel->maxz;
                uy = voxel->maxy - dy;
                uz = voxel->minz - dz;
                break;
            case 3: /* TSW -> BNE */
                dy = voxel->maxy;
                dz = voxel->maxz;
                uy = voxel->miny - dy;
                uz = voxel->minz - dz;
                break;
            default:    /* do nothing */
                break;
        }

        nu = Mod->flist[idx].nx * ux + Mod->flist[idx].ny * uy + Mod->flist[idx].nz * uz;
        if (nu != 0) {
            lambda = (Mod->flist[idx].dist - Mod->flist[idx].nx * dx - Mod->flist[idx].ny * dy - Mod->flist[idx].nz * dz) / nu;
            if (0 < lambda && lambda < 1) {
                SX = dx + lambda * ux;
                SY = dy + lambda * uy;
                SZ = dz + lambda * uz;
                if (__point_in_polygon(SX, SY, SZ, idx) == OCTREE_DONE) { return OCTREE_DONE; }
            }
        }
    }

    return OCTREE_ERROR;
}

//=============================================================================
// @ARG: point's coordinates
// @ARG: face index
//=============================================================================
int octree::__point_in_polygon(const double SX, const double SY, const double SZ, const int idx) {
    double A_diff;                                    // difference of polygons' areas
    double v1x = Mod->vlist[Mod->flist[idx].v1].x;    // node v1's x-coordinate
    double v1y = Mod->vlist[Mod->flist[idx].v1].y;    // node v1's y-coordinate
    double v1z = Mod->vlist[Mod->flist[idx].v1].z;    // node v1's z-coordinate
    double v2x = Mod->vlist[Mod->flist[idx].v2].x;    // node v2's x-coordinate
    double v2y = Mod->vlist[Mod->flist[idx].v2].y;    // node v2's y-coordinate
    double v2z = Mod->vlist[Mod->flist[idx].v2].z;    // node v2's z-coordinate
    double v3x = Mod->vlist[Mod->flist[idx].v3].x;    // node v3's x-coordinate
    double v3y = Mod->vlist[Mod->flist[idx].v3].y;    // node v3's y-coordinate
    double v3z = Mod->vlist[Mod->flist[idx].v3].z;    // node v3's z-coordinate

    /* compute area of polygon (v1,v2,v3) */
    if (isZero(Mod->flist[idx].A)) { Mod->flist[idx].A = __comp_polygon_area(v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z); }
    A_diff = Mod->flist[idx].A;

    /* compute areas of smaller polygons (v1,v2,S), (v2,v3,S), (v3,v1,S) */
    A_diff -= __comp_polygon_area(v1x, v1y, v1z, v2x, v2y, v2z, SX, SY, SZ);
    A_diff -= __comp_polygon_area(v2x, v2y, v2z, v3x, v3y, v3z, SX, SY, SZ);
    A_diff -= __comp_polygon_area(v3x, v3y, v3z, v1x, v1y, v1z, SX, SY, SZ);

    /* check if areas have 'equal' size */
    if (fabs(A_diff) < point_in_polygon_treshold) { return OCTREE_DONE; }

    return OCTREE_ERROR;
}

//=============================================================================
// @ARG: vertex v1's coordinates
// @ARG: vertex v2's coordinates
// @ARG: vertex v3's coordinates
//=============================================================================
double octree::__comp_polygon_area(const double v1x, const double v1y, const double v1z,
                                   const double v2x, const double v2y, const double v2z,
                                   const double v3x, const double v3y, const double v3z) {
    double vec1x = v2x - v1x;    // vector vec1's x-coordinate
    double vec1y = v2y - v1y;    // vector vec1's y-coordinate
    double vec1z = v2z - v1z;    // vector vec1's z-coordinate
    double vec2x = v3x - v1x;    // vector vec2's x-coordinate
    double vec2y = v3y - v1y;    // vector vec2's y-coordinate
    double vec2z = v3z - v1z;    // vector vec2's z-coordinate
    double vx, vy, vz;            // cross product vector's coordinates
    double A_poly;                // polygon's area

    vx = vec1y * vec2z - vec1z * vec2y;
    vy = vec1z * vec2x - vec1x * vec2z;
    vz = vec1x * vec2y - vec1y * vec2x;
    A_poly = sqrt(vx * vx + vy * vy + vz * vz);

    return A_poly;
}

//=============================================================================
// @ARG: none
//=============================================================================
void octree::__comp_bbox(void) {
    int i;        // loop counter

    for (i = 0; i < Mod->nface; ++i) {
        Mod->flist[i].minx = Mod->vlist[Mod->flist[i].v1].x;
        Mod->flist[i].maxx = Mod->vlist[Mod->flist[i].v1].x;
        Mod->flist[i].miny = Mod->vlist[Mod->flist[i].v1].y;
        Mod->flist[i].maxy = Mod->vlist[Mod->flist[i].v1].y;
        Mod->flist[i].minz = Mod->vlist[Mod->flist[i].v1].z;
        Mod->flist[i].maxz = Mod->vlist[Mod->flist[i].v1].z;
        __check_vertex(i, Mod->vlist[Mod->flist[i].v2].x, Mod->vlist[Mod->flist[i].v2].y, Mod->vlist[Mod->flist[i].v2].z);
        __check_vertex(i, Mod->vlist[Mod->flist[i].v3].x, Mod->vlist[Mod->flist[i].v3].y, Mod->vlist[Mod->flist[i].v3].z);
    }
}

//=============================================================================
// @ARG: none
//=============================================================================
void octree::__comp_dist(void) {
    int i;        // loop counter

    for (i = 0; i < Mod->nface; ++i) {
        Mod->flist[i].dist = Mod->flist[i].nx * Mod->vlist[Mod->flist[i].v1].x +
                             Mod->flist[i].ny * Mod->vlist[Mod->flist[i].v1].y +
                             Mod->flist[i].nz * Mod->vlist[Mod->flist[i].v1].z;
    }
}

//=============================================================================
// @ARG: face index
// @ARG: vertex coordinates to be processed
//=============================================================================
void octree::__check_vertex(const int idx, const double x, const double y, const double z) {
    if (x < Mod->flist[idx].minx) { Mod->flist[idx].minx = x; }
    else if (x > Mod->flist[idx].maxx) { Mod->flist[idx].maxx = x; }
    if (y < Mod->flist[idx].miny) { Mod->flist[idx].miny = y; }
    else if (y > Mod->flist[idx].maxy) { Mod->flist[idx].maxy = y; }
    if (z < Mod->flist[idx].minz) { Mod->flist[idx].minz = z; }
    else if (z > Mod->flist[idx].maxz) { Mod->flist[idx].maxz = z; }
}

//=============================================================================
// @ARG: none
//=============================================================================
void octree::__setup_stack(void) {
    Voxel *voxel;    // voxel data
    int i;            // loop counter

    voxel = new Voxel;
    voxel->mcode = 0b111;
    voxel->x = 0;
    voxel->y = 0;
    voxel->z = 0;
    voxel->minx = -1;
    voxel->miny = -1;
    voxel->minz = -1;
    voxel->maxx = 1;
    voxel->maxy = 1;
    voxel->maxz = 1;
    voxel->color = OCTREE_ERROR;
    voxel->nface = Mod->nface;
    voxel->findex = new int[Mod->nface];
    for (i = 0; i < Mod->nface; ++i) { voxel->findex[i] = i; }

    in = new llist;
    in->data = nullptr;
    in->next = nullptr;
    out = new llist;
    out->data = nullptr;
    out->next = nullptr;
    __push(in, voxel);
    voxel = nullptr;
}

//=============================================================================
// @ARG: linked list
//=============================================================================
void octree::__cleanup_stack(llist *list) {
    Voxel *voxel;    // voxel data

    if (list != nullptr) {
        /* delete stack elements */
        while ((voxel = __pop(list)) != nullptr) {
            if (voxel->findex != nullptr) {
                delete[] voxel->findex;
                voxel->findex = nullptr;
            }
            delete voxel;
        }
        /* delete stack head */
        delete list;
        list = nullptr;
    }
}


//=============================================================================
// @ARG: vtk filename of the adative grid to write
//=============================================================================
void octree::__write_adaptive_vtk_file(const char *filename) {

    std::ofstream os;

    /* write grid to file */
    os.open(filename, std::ofstream::out | std::ofstream::trunc);
    if (!os.is_open()) {
        std::cerr << "ERROR >> couldn't write file " << filename << std::endl;
        exit(0);
    }

    int ct_pts = 0, ct_vxls = 0, ct = 0; // counters
    llist *tmp;         // next stack element
    double v[8][3];     // voxel's node coordinates

    // count points and voxels
    if (out != nullptr && out->next != nullptr) {
        tmp = out->next;
        while (tmp != nullptr) {
            ++ct_vxls;
            tmp = tmp->next;
        }
    }
    ct_pts = 8 * ct_vxls;

    printf("\tWriting out %d points and %d voxels as unstructured VTK file.\n", ct_pts, ct_vxls);

    // write geom (condensed format)
    os << "# vtk DataFile Version 2.0\n";
    os << "MESH written by voxelisator\n";
    os << "ASCII\n";
    os << "DATASET UNSTRUCTURED_GRID\n";

    // write points
    os << "POINTS " << ct_pts << " float\n";
    if (out != nullptr && out->next != nullptr) {
        tmp = out->next;
        while (tmp != nullptr) {
            /* compute node coordinates */
            v[0][0] = tmp->data->minx;
            v[0][1] = tmp->data->miny;
            v[0][2] = tmp->data->minz;    /* BNW */
            v[1][0] = tmp->data->maxx;
            v[1][1] = tmp->data->miny;
            v[1][2] = tmp->data->minz;    /* BNE */
            v[2][0] = tmp->data->maxx;
            v[2][1] = tmp->data->maxy;
            v[2][2] = tmp->data->minz;    /* TNE */
            v[3][0] = tmp->data->minx;
            v[3][1] = tmp->data->maxy;
            v[3][2] = tmp->data->minz;    /* TNW */
            v[4][0] = tmp->data->minx;
            v[4][1] = tmp->data->miny;
            v[4][2] = tmp->data->maxz;    /* BSW */
            v[5][0] = tmp->data->maxx;
            v[5][1] = tmp->data->miny;
            v[5][2] = tmp->data->maxz;    /* BSE */
            v[6][0] = tmp->data->maxx;
            v[6][1] = tmp->data->maxy;
            v[6][2] = tmp->data->maxz;    /* TSE */
            v[7][0] = tmp->data->minx;
            v[7][1] = tmp->data->maxy;
            v[7][2] = tmp->data->maxz;    /* TSW */

#ifdef RESIZE_VTK_OUTPUT
            for (int i = 0; i < 8; ++i) os << v[i][0] * Mod->bbox_norm + Mod->bbox_diffx << " " << v[i][1] * Mod->bbox_norm + Mod->bbox_diffy << " " << v[i][2] * Mod->bbox_norm + Mod->bbox_diffz << "\n";
#else
            for (int i = 0; i < 8; ++i) os << v[i][0] << " " << v[i][1] << " " << v[i][2] << "\n";
#endif

            tmp = tmp->next;
        }
    }

    // write cells
    os << "CELLS " << ct_vxls << " " << 9 * ct_vxls << "\n";
    if (out != nullptr && out->next != nullptr) {
        tmp = out->next;
        ct = 0;
        while (tmp != nullptr) {

            os << 8 << " ";
            for (int i = 0; i < 8; ++i) os << " " << ct++;
            os << "\n";

            tmp = tmp->next;
        }
    }

    // write cell type
    os << "CELL_TYPES " << ct_vxls << "\n";
    for (int i = 0; i < ct_vxls; ++i) os << "12\n";

    // write data
    os << "CELL_DATA " << ct_vxls << "\n";
    os << "SCALARS voxel_type int\nLOOKUP_TABLE default\n";
    if (out != nullptr && out->next != nullptr) {
        tmp = out->next;
        while (tmp != nullptr) {
            os << (int) tmp->data->color << "\n";

            tmp = tmp->next;
        }
    }

    printf("\tWriteout done!\n");

    os.close();
}

//=============================================================================
// @ARG: depth of octree
//=============================================================================
void octree::__generate_uniform_grid(const int maxdepth) {
    int size = int(pow(2.0, maxdepth));    // grid size
    char ***grid;                            // uniform grid
    int i, j, k;                            // loop counters
    Voxel *voxel;                            // voxel to be processed
    double diam = pow(2.0, -(maxdepth - 1));    // diameter of smallest voxel
    int *nelem;                                // amount of 'unit' elements for different voxel sizes
    int xidx, yidx, zidx;                    // grid's start indices of voxel

    /* if no filename is given, don't write out anything */
    if (voxfile == nullptr && vtkfile == nullptr) return;

    /* compute amount of elements for different voxel sizes */
    nelem = new int[maxdepth + 1];
    nelem[maxdepth] = 1;
    for (i = maxdepth - 1; i >= 0; --i) { nelem[i] = 2 * nelem[i + 1]; }

    /* allocate memory for grid [size]x[size]x[size] */
    grid = new char **[size];
    for (i = 0; i < size; ++i) {
        grid[i] = new char *[size];
        for (j = 0; j < size; ++j) {
            grid[i][j] = new char[size];
            for (k = 0; k < size; ++k) { grid[i][j][k] = '0'; }
        }
    }

    /* fill grid with voxels */
    while ((voxel = __pop(out)) != nullptr) {
        xidx = int((1.0 + voxel->minx) / diam);
        yidx = int((1.0 + voxel->miny) / diam);
        zidx = int((1.0 + voxel->minz) / diam);
        for (i = xidx; i < xidx + nelem[__get_depth_by_id(voxel->mcode)]; ++i) {
            for (j = yidx; j < yidx + nelem[__get_depth_by_id(voxel->mcode)]; ++j) {
                for (k = zidx; k < zidx + nelem[__get_depth_by_id(voxel->mcode)]; ++k) {
                    if (voxel->color == OCTREE_COLOR_GREY) { grid[i][j][k] = '1'; }
                    else if (voxel->color == OCTREE_COLOR_BLACK) { grid[i][j][k] = '2'; }
                }
            }
        }
        if (voxel->nface != 0) {
            delete[] voxel->findex;
            voxel->findex = nullptr;
        }
        delete voxel;
    }

    /* fill OCTREE_COLOR_WHITE holes inside grid representation */
    for (i = 1; i < size - 1; ++i) {
        for (j = 1; j < size - 1; ++j) {
            for (k = 1; k < size - 1; ++k) {
                if (grid[i][j][k] == '0') {
                    if (grid[i - 1][j][k] == '2' || grid[i][j - 1][k] == '2' || grid[i][j][k - 1] == '2' ||
                        grid[i + 1][j][k] == '2' || grid[i][j + 1][k] == '2' || grid[i][j][k + 1] == '2') { grid[i][j][k] = '2'; }
                }
            }
        }
    }
    for (i = 0; i < size; ++i) {
        for (j = 0; j < size; ++j) {
            for (k = 0; k < size; ++k) { if (grid[i][j][k] == '2') { grid[i][j][k] = '1'; }}
        }
    }

    // write vox file if asked
    if (voxfile) {
        std::cout << "\tWriting VOX file to disk..." << std::flush;
        __write_vox_file(voxfile, grid, size);
#ifdef PRINT
        std::cout << "done" << std::endl;
#endif
    }

    // write vtk file if asked
    if (vtkfile) {
        std::cout << "\tWriting VTK file to disk..." << std::flush;
        __write_vtk_file(vtkfile, grid, size);
#ifdef PRINT
        std::cout << "done" << std::endl;
#endif
    }

    /* free allocated memory */
    for (i = 0; i < size; ++i) {
        for (j = 0; j < size; ++j) { delete[] grid[i][j]; }
        delete[] grid[i];
    }
    delete[] grid;
    delete[] nelem;
}


//=============================================================================
// @ARG: filename of VOX file
// @ARG: grid containing data
// @ARG: size of grid
//=============================================================================
void octree::__write_vox_file(const char *filename, char ***grid, const int &size) {
    std::ofstream os;
    int i, j, k;

    /* write vox file */
    os.open(filename, std::ofstream::out | std::ofstream::trunc);
    if (!os.is_open()) {
        std::cerr << "ERROR >> couldn't write file " << filename << std::endl;
        exit(0);
    }

    // write header
    os << size << std::endl << size << std::endl << size << std::endl;

    // write data
    for (i = 0; i < size; ++i) {
        for (j = 0; j < size; ++j) {
            for (k = 0; k < size; ++k) { os << grid[i][j][k] << "\n"; }
        }
    }

    // write data
    for (i = 0; i < size; ++i) {
        for (j = 0; j < size; ++j) {
            for (k = 0; k < size; ++k) { os << grid[i][j][k] << "\n"; }
        }
    }

    // closing file handler
    os.close();
}

//=============================================================================
// @ARG: filename of VTK file
// @ARG: grid containing data
// @ARG: size of grid
//=============================================================================
void octree::__write_vtk_file(const char *filename, char ***grid, const int &size) {
    std::ofstream os;
    int i, j, k;

    /* write grid to file */
    os.open(filename, std::ofstream::out | std::ofstream::trunc);
    if (!os.is_open()) {
        std::cerr << "ERROR >> couldn't write file " << filename << std::endl;
        exit(0);
    }

    // header
    os << "# vtk DataFile Version 2.0\nVoxelisation Visualization\nASCII\nDATASET RECTILINEAR_GRID\n";
    os << "DIMENSIONS " << size + 1 << " " << size + 1 << " " << size + 1 << "\n";
    os << "X_COORDINATES " << size + 1 << " float\n";
    for (i = 0; i < size + 1; ++i) os << i << " ";
    os << "\n";
    os << "Y_COORDINATES " << size + 1 << " float\n";
    for (i = 0; i < size + 1; ++i) os << i << " ";
    os << "\n";
    os << "Z_COORDINATES " << size + 1 << " float\n";
    for (i = 0; i < size + 1; ++i) os << i << " ";
    os << "\n";
    os << "CELL_DATA " << size * size * size << "\nSCALARS bc int 1\nLOOKUP_TABLE default\n";

    // write point data (VTK order style)
    for (k = 0; k < size; ++k) {
        for (j = 0; j < size; ++j) {
            for (i = 0; i < size; ++i) {
                os << grid[i][j][k] << "\n";
            }
        }
    }

    // closing file handler
    os.close();
}

//=============================================================================
// @ARG: amount of black voxels
// @ARG: amount of white voxels
// @ARG: amount of grey voxels
//=============================================================================
void octree::count(long *black, long *white, long *grey) {
    *black = *white = *grey = 0;

    llist *it(out);

    while (it != nullptr) {

        if (it->data != nullptr) {
            switch (it->data->color) {
                case (OCTREE_COLOR_WHITE) :
                    (*white)++;
                    break;
                case (OCTREE_COLOR_BLACK) :
                    (*black)++;
                    break;
                case (OCTREE_COLOR_GREY) :
                    (*grey)++;
                    break;
                default:
                    std::cerr << "ERROR: Undefined color found while executing octree::count(...) !" << std::endl;
                    break;
            }
        }

        it = it->next;
    }

}

/*

Flood fill octree.

- Eric Fichter, fichter@e3d.rwth-aachen.de

-----------------------------------------------------------------------

Copyright (c) 2019 RWTH Aachen University.  All rights reserved.

Permission to use, copy, modify and distribute this software and its
documentation for any purpose is hereby granted without fee, provided
that the above copyright notice and this permission notice appear in
all copies of this software and that you do not sell the software.

THE SOFTWARE IS PROVIDED "AS IS" AND WITHOUT WARRANTY OF ANY KIND,
EXPRESS, IMPLIED OR OTHERWISE, INCLUDING WITHOUT LIMITATION, ANY
WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.

*/

//=============================================================================
// @ARG: none
//=============================================================================
bool octree::isLeaf(std::bitset<64> &mcode, std::unordered_map<std::bitset<64>, std::vector<std::bitset<64>>

> &hash_nb_leafs) {
// if mcode is present in list it's a leaf node. Otherwise it is not or mcode does not exist in octree
    if (hash_nb_leafs.find(mcode) == hash_nb_leafs.end())
        return false;
    else
        return true;
}

//=============================================================================
// @ARG: none
//=============================================================================
std::unordered_map<std::bitset<64>, std::vector<std::bitset<64 >>>

octree::__neighbour_search() {
#ifdef PRINT
    std::cout << "\tNeighbour search...";
#endif
    // create hash map id->neighbour ids (already initialized)
#ifdef PRINT
    std::cout << "initialize nb hashmap...";
#endif
    std::unordered_map<std::bitset<64>, std::vector<std::bitset<64 >>> hash_nb_leafs = create_nb_hash_map(); // holds all leaf ids as key with a vector of their leaf neighbour ids
#ifdef PRINT
    std::cout << "find neighbours...";
#endif
    // iterate through all leaf voxels, find neighbour and save in second hashmap (id->vector of ids)
    llist *tmp;        // stack containing voxels to be processed
    if (out != nullptr && out->next != nullptr) {
        tmp = out->next;
        while (tmp != nullptr) {
            //if (hash_id_leaf[tmp->data->mcode]) { // only leaf nodes have to be checked
            std::bitset<64> MCODE = tmp->data->mcode; // morton code of node to be investigated for neighbours

            // for each direction find neighbours
            for (unsigned int k = 0; k < 3; k++) //0=x, 1=y, 2=z
                for (unsigned int n = 0; n < 2; n++) { // 0=west, south, bottom, 1=east, north, top

                    unsigned int bit_value_of_interest = MCODE[2 - k]; // it's the value of the bit at position of direction x/y/z

                    // ########################## 1 - Find sibling #######################
                    if (bit_value_of_interest == n) { // bit in last set in respective spatial direction compared with n, no sibling can be found

                        // ########################## 4 - Check if node is boundary node #######################
                        if (!__isBoundary(MCODE, k)) {

                            // ########################## 5 - Find off-branch nb #######################
                            // go up until sibling is found
                            std::bitset<64> mcode_sibl = find_sibling_upwards(MCODE, k, n);

                            // Step 1: Go downwards until depth of original node was reached or sibling or its children is leaf. Flip bits in respective direction for each level
                            bool neighbour_found = false;
                            unsigned int depth_mcode = __get_depth_by_id(MCODE);
                            unsigned int depth_siblg = __get_depth_by_id(mcode_sibl);

                            while (depth_siblg <= depth_mcode && !neighbour_found) {

                                if (isLeaf(mcode_sibl, hash_nb_leafs)) { // if node is a leaf node, it's the neighbour
                                    hash_nb_leafs[MCODE].push_back(mcode_sibl);
                                    neighbour_found = true;
                                } else { // else add level to morton code and flip bit in searched direction

                                    if (depth_siblg < depth_mcode) {
                                        mcode_sibl = mcode_sibl << 3; // add three zeros at the end
                                        int pos = (depth_mcode - depth_siblg - 1) * 3;
                                        mcode_sibl[2] = MCODE[2 + pos];
                                        mcode_sibl[1] = MCODE[1 + pos];
                                        mcode_sibl[0] = MCODE[0 + pos];
                                        mcode_sibl.flip(2 - k); // flip bit at position of respective direction

                                    }
                                    depth_siblg++;

                                }
                            }

                            // Step 2: If level of original node was reached, go downwards recursively and look for all nodes that have opposite of n at position k
                            if (!neighbour_found) {
                                unsigned int n_flipped = !n;

                                std::stack<std::bitset<64 >> branches_loop;
                                std::vector<std::bitset<64>> mcodes_children = generate_mcode_children_with_desired_value(mcode_sibl, k, n_flipped);

                                for (auto &i: mcodes_children)
                                    branches_loop.push(i); // push children to stack

                                while (!branches_loop.empty()) { // iterate through stack
                                    std::bitset<64> current_mcode = branches_loop.top();
                                    branches_loop.pop();

                                    if (isLeaf(current_mcode, hash_nb_leafs)) // check if leaf
                                        hash_nb_leafs[MCODE].push_back(current_mcode);

                                    else {
                                        std::vector<std::bitset<64>> mcodes_children = generate_mcode_children_with_desired_value(current_mcode, k, n_flipped);
                                        for (auto &i: mcodes_children)
                                            branches_loop.push(i); // push children to stack
                                    }
                                }

                            }

                        }
                    } else {

                        std::bitset<64> mcode_sibl = MCODE;
                        mcode_sibl.flip(2 - k);

                        // ########################## 2 - Check if sibling is leaf #######################
                        if (isLeaf(mcode_sibl, hash_nb_leafs))
                            hash_nb_leafs[MCODE].push_back(mcode_sibl);
                        else {

                            // ########################## 3 - Find all neighbours by going recursively deeper in all children and childrens children of sibling #######################
                            // example: looking for eastern nbs (k=0, n=1): node = 000, sibling = 100. All children of sibling with n=0 (bit_value_of_interest) at k=0 are neighbours
                            std::stack<std::bitset<64 >> branches_loop;
                            std::vector<std::bitset<64>> mcodes_children = generate_mcode_children_with_desired_value(mcode_sibl, k, bit_value_of_interest);
                            for (auto &i: mcodes_children)
                                branches_loop.push(i); // push children to stack
//
                            while (!branches_loop.empty()) { // iterate through stack
                                std::bitset<64> current_mcode = branches_loop.top();
                                branches_loop.pop();

                                if (isLeaf(current_mcode, hash_nb_leafs))  // check if leaf
                                    hash_nb_leafs[MCODE].push_back(current_mcode);
                                else {
                                    std::vector<std::bitset<64>> mcodes_children = generate_mcode_children_with_desired_value(current_mcode, k, bit_value_of_interest);
                                    for (auto &i: mcodes_children)
                                        branches_loop.push(i); // push children to stack

                                    /* std::vector <std::bitset<64>> mcodes_children = generate_mcode_children(current_mcode);
                                     for (unsigned int i = 0; i < mcodes_children.size(); i++)
                                         if (mcodes_children[i][2 - k] == bit_value_of_interest)
                                             branches_loop.push(mcodes_children[i]); // push children to stack*/
                                }
                            }
                        }
                    }
                }
            //}
            tmp = tmp->next;
        }
    }
#ifdef PRINT
    std::cout << "done." << std::endl;
#endif

    return hash_nb_leafs;

}

//=============================================================================
// @ARG: map of neighbours
// @ARG: map of octree
//=============================================================================
int octree::__flood_fill(std::unordered_map<std::bitset<64>, std::vector<std::bitset<64>>> &hash_nb_leafs, std::unordered_map<std::bitset<64>, Voxel *> &hash_octree) {
#ifdef PRINT
    std::cout << "\tFlood Fill Octree...";
#endif
    std::stack<std::bitset<64 >> branches_loop;

    bool ff_start_with_boundaries = true; // if set to true, bounding box of model must not be extended!

    if (ff_start_with_boundaries) {
#ifdef PRINT
        std::cout << "get start boundary nodes...";
#endif
        std::vector<std::bitset<64>> boundary_mcodes = get_all_boundary__leaf_nodes();
        for (auto &boundary_mcode: boundary_mcodes)
            branches_loop.push(boundary_mcode);
    } else {
        // get start node (000er)
#ifdef PRINT
        std::cout << "get start node...";
#endif
        std::bitset<64> mcode_start = find_flood_fill_mcode_start(hash_octree, hash_nb_leafs);
        branches_loop.push(mcode_start);
    }

    // reset node colors
#ifdef PRINT
    std::cout << "reset node colors...";
#endif
    reset_node_colors(OCTREE_COLOR_WHITE);
#ifdef PRINT
    std::cout << "start flooding...";
#endif
    std::unordered_map<std::bitset<64>, Voxel *>::const_iterator ffit;

    // OCTREE_COLOR_WHITE - not visited yet
    // OCTREE_COLOR_BLACK - already visited and empty node
    // OCTREE_COLOR_GREY - already visited and non-empty node

    while (!branches_loop.empty()) {
        ffit = hash_octree.find(branches_loop.top());
        branches_loop.pop();

        if (ffit->second->color == OCTREE_COLOR_WHITE) { // check if the node was already visited
            if (ffit->second->nface == 0) {
                ffit->second->color = OCTREE_COLOR_BLACK; // visited and empty
                std::vector<std::bitset<64 >> nbs = hash_nb_leafs[ffit->second->mcode];
                for (auto &nb: nbs) { // append neighbours to stack
                    branches_loop.push(nb);

#ifdef VERIFY_NEIGHBOURSEARCH
                    std::vector <std::bitset<64 >> nbs_test = hash_nb_leafs[*itr];
                    std::vector < std::bitset < 64 >> ::iterator
                    it_test = find(nbs_test.begin(), nbs_test.end(), ffit->second->mcode);
                    if (it_test == nbs_test.end()) {
                        std::cout << "Wrong NB " << bitset_to_clean_string(ffit->second->mcode) << " " << bitset_to_clean_string(*itr) << std::endl;
                        for (auto itr_lol = nbs_test.begin(); itr_lol != nbs_test.end(); itr_lol++)
                            std::cout << bitset_to_clean_string(*itr_lol) << std::endl;
                    }
#endif
                }

            } else
                ffit->second->color = OCTREE_COLOR_GREY; // visited and non-empty
        }
    }

#ifdef PRINT
    std::cout << "done." << std::endl;
#endif

#ifdef VERIFY_NEIGHBOURSEARCH
    std::cout << "\n############## Brute Force Test if saved neighbours are real neighbours: ###############";
    for (auto itr = hash_nb_leafs.begin(); itr != hash_nb_leafs.end(); itr++) {
        double minx = hash_octree[itr->first]->minx;
        double miny = hash_octree[itr->first]->miny;
        double minz = hash_octree[itr->first]->minz;
        double maxx = hash_octree[itr->first]->maxx;
        double maxy = hash_octree[itr->first]->maxy;
        double maxz = hash_octree[itr->first]->maxz;
        // std::cout << bitset_to_clean_string(itr->first) << '\n';
        for (unsigned int i = 0; i < (itr->second).size(); i++) {
            double cminx = hash_octree[itr->second[i]]->minx;
            double cminy = hash_octree[itr->second[i]]->miny;
            double cminz = hash_octree[itr->second[i]]->minz;
            double cmaxx = hash_octree[itr->second[i]]->maxx;
            double cmaxy = hash_octree[itr->second[i]]->maxy;
            double cmaxz = hash_octree[itr->second[i]]->maxz;
            //   std::cout << "\t" << bitset_to_clean_string(itr->second[i]) << std::endl;
            if ((fabs(minx - cmaxx) > 1e-9) && (fabs(maxx - cminx) > 1e-9) && (fabs(miny - cmaxy) > 1e-9) && (fabs(maxy - cminy) > 1e-9) && (fabs(minz - cmaxz) > 1e-9) && (fabs(maxz - cminz) > 1e-9)) {
                std::cout << minx << " " << miny << " " << minz << " " << maxx << " " << maxy << " " << maxz << std::endl;
                std::cout << cmaxx << " " << cmaxy << " " << cmaxz << " " << cminx << " " << cminy << " " << cminz << std::endl;
                abort();
            }
        }
    }
    std::cout << "...finished.\n" << std::endl;

#endif


/*
std::cout << "\n\n############## Those are all leafs and their neighbours: ###############" << std::endl;
for (auto itr = hash_nb_leafs.begin(); itr != hash_nb_leafs.end(); itr++) {
    std::cout << "A" << bitset_to_clean_string(itr->first) << '\t' << hash_octree[itr->first]->isLeaf << '\n';
    for (unsigned int i = 0; i < (itr->second).size(); i++)
        std::cout << "\t" << bitset_to_clean_string(itr->second[i]) << '\t' << hash_octree[itr->second[i]]->isLeaf << std::endl;
}
*/

    return 0;

}

//=============================================================================
// @ARG: none
//=============================================================================
void octree::find_facade_triangles() {

    std::cout << "\tLocate triangles...";

//    if (Mod->use_bof) {
//
//        std::vector<std::string> guids;
//
//        llist *tmp;
//        if (out != nullptr && out->next != nullptr) {
//            tmp = out->next;
//            while (tmp != nullptr) {
//                if (tmp->data->color == OCTREE_COLOR_GREY) {
//                    for (int i = 0; i < tmp->data->nface; ++i) {
//                        unsigned int index_tri = tmp->data->findex[i];
//                        std::string guid = Mod->flist[index_tri].guid;
//                        guids.push_back(guid);
//                    }
//                }
//                tmp = tmp->next;
//            }
//        }
//
//        sort(guids.begin(), guids.end());
//        guids.erase(std::unique(guids.begin(), guids.end()), guids.end());
//
//        std::ofstream f;
//        f.open("Faces_Guid.txt");
//        for (auto &guid : guids)
//            f << guid << std::endl;
//        f.close();
//    } else {


    std::vector<unsigned int> idx_faces;

    llist *tmp;
    if (out != nullptr && out->next != nullptr) {
        tmp = out->next;
        while (tmp != nullptr) {
            if (tmp->data->color == OCTREE_COLOR_GREY) {
                for (int i = 0; i < tmp->data->nface; ++i) {
                    unsigned int index_tri = tmp->data->findex[i];
                    idx_faces.push_back(index_tri);
                }
            }
            tmp = tmp->next;
        }
    }

    sort(idx_faces.begin(), idx_faces.end());
    idx_faces.erase(std::unique(idx_faces.begin(), idx_faces.end()), idx_faces.end());

    std::ofstream f;
    f.open("Faces_Index.txt");
    for (unsigned int idx_face: idx_faces)
        f << idx_face << std::endl;
    f.close();
    // }

    std::cout << "done" << std::endl;

}

//=============================================================================
// @ARG: color
//=============================================================================
void octree::reset_node_colors(char clr) {

    llist *tmp;
    if (out != nullptr && out->next != nullptr) {
        tmp = out->next;
        while (tmp != nullptr) {
            tmp->data->color = clr;
            tmp = tmp->next;
        }
    }
}

//=============================================================================
// @ARG: octree hash map
//=============================================================================
std::bitset<64> octree::find_flood_fill_mcode_start(const std::unordered_map<std::bitset<64>, Voxel *

> &hash_octree, std::unordered_map<std::bitset<64>, std::vector<std::bitset<64 >>> &hash_nb_leafs) {

// looks for deepest node in octree with only zeros as mcode
    std::bitset<64> mcode = 0b111; // start with root

    for (
            unsigned int i = 0;
            i < 99; i++) {
        mcode = mcode << 3;
        if (
                isLeaf(mcode, hash_nb_leafs
                ))
            break;
    }

    std::unordered_map<std::bitset<64>, Voxel *>
    ::const_iterator it;
    it = hash_octree.find(mcode);
    if ((it->second)->nface > 0) {
        std::cerr << "Start node " <<
                  bitset_to_clean_string(mcode)
                  << " is not empty: " << (it->second)->nface <<
                  std::endl;

        abort();

    }

    return
            mcode;
}

//=============================================================================
// @ARG: morton code of node
// @ARG: position of bit in searched direction
// @ARG: min or max in direction
//=============================================================================
std::bitset<64> octree::find_sibling_upwards(std::bitset<64> mcode, unsigned int k, unsigned int n) {

    std::bitset<64> mcode_sibl = mcode;
    for (int i = 0; i < __get_depth_by_id(mcode); i++) {
        mcode_sibl = mcode_sibl >> 3; // bitshift to the right
        if (__get_depth_by_id(mcode_sibl) == 0) { // root node is reached. Should not be possible because node would be a boundary node and boundary nodes are filtered before.
            std::cout << "Root was reached. " << bitset_to_clean_string(mcode_sibl) << " Should not be possible because node would be a boundary node." << std::endl;
            abort();
        }

        if (mcode_sibl[2 - k] != n)
            break;
    }

    mcode_sibl.flip(2 - k);

    return mcode_sibl;
}

//=============================================================================
// @ARG: morton code of node
//=============================================================================
std::vector<std::bitset<64>> octree::generate_mcode_children(std::bitset<64> mcode) {

    std::vector<std::bitset<64>> children;

    std::bitset<64> m = mcode << 3;
    std::bitset<64> mcode_temp_1 = m;
    std::bitset<64> mcode_temp_2 = m;
    std::bitset<64> mcode_temp_3 = m;
    std::bitset<64> mcode_temp_4 = m;
    std::bitset<64> mcode_temp_5 = m;
    std::bitset<64> mcode_temp_6 = m;
    std::bitset<64> mcode_temp_7 = m;
    std::bitset<64> mcode_temp_8 = m;
    mcode_temp_1 |= 0b000;
    mcode_temp_2 |= 0b001;
    mcode_temp_3 |= 0b010;
    mcode_temp_4 |= 0b011;
    mcode_temp_5 |= 0b100;
    mcode_temp_6 |= 0b101;
    mcode_temp_7 |= 0b110;
    mcode_temp_8 |= 0b111;
    children.push_back(mcode_temp_1);
    children.push_back(mcode_temp_2);
    children.push_back(mcode_temp_3);
    children.push_back(mcode_temp_4);
    children.push_back(mcode_temp_5);
    children.push_back(mcode_temp_6);
    children.push_back(mcode_temp_7);
    children.push_back(mcode_temp_8);

    return children;
}

//=============================================================================
// @ARG: morton code of node
//=============================================================================
std::vector<std::bitset<64>> octree::generate_mcode_children_with_desired_value(const std::bitset<64> &mcode, const unsigned int &k, const unsigned int &b) {

    // can only be 4 children
    // mcodes_children[i][2 - k] == bit_value_of_interest
    std::vector<std::bitset<64>> children;
    std::bitset<64> m = mcode << 3;
    std::bitset<64> mcode_temp_1 = m;
    std::bitset<64> mcode_temp_2 = m;
    std::bitset<64> mcode_temp_3 = m;
    std::bitset<64> mcode_temp_4 = m;

    switch (k) {
        case 0: // x-position
            if (b == 0) { // should be 0
                mcode_temp_1 |= 0b000;
                mcode_temp_2 |= 0b001;
                mcode_temp_3 |= 0b010;
                mcode_temp_4 |= 0b011;
            } else { // should be 0
                mcode_temp_1 |= 0b100;
                mcode_temp_2 |= 0b101;
                mcode_temp_3 |= 0b110;
                mcode_temp_4 |= 0b111;
            }
            break;
        case 1: // y-position
            if (b == 0) { // should be 0
                mcode_temp_1 |= 0b000;
                mcode_temp_2 |= 0b001;
                mcode_temp_3 |= 0b100;
                mcode_temp_4 |= 0b101;
            } else { // should be 1
                mcode_temp_1 |= 0b010;
                mcode_temp_2 |= 0b011;
                mcode_temp_3 |= 0b110;
                mcode_temp_4 |= 0b111;
            }
            break;
        case 2: // z-position
            if (b == 0) { // should be 0
                mcode_temp_1 |= 0b000;
                mcode_temp_2 |= 0b010;
                mcode_temp_3 |= 0b100;
                mcode_temp_4 |= 0b110;
            } else { // should be 1
                mcode_temp_1 |= 0b001;
                mcode_temp_2 |= 0b011;
                mcode_temp_3 |= 0b101;
                mcode_temp_4 |= 0b111;
            }
            break;
        default:
            abort();
    }

    children.push_back(mcode_temp_1);
    children.push_back(mcode_temp_2);
    children.push_back(mcode_temp_3);
    children.push_back(mcode_temp_4);

    return children;
}

//=============================================================================
// @ARG: morton code of node
//=============================================================================
std::string octree::bitset_to_clean_string(std::bitset<64> b) {

    std::string mcode = b.to_string();
    // get the substring
    unsigned int index = 0;
    std::string strippedstr;
    for (unsigned int i = 0; i < mcode.length(); ++i) {
        if (mcode[i] == '1') {
            index = i;
            break;
        }
    }
    return mcode.substr(index);

}

//=============================================================================
// @ARG: none
//=============================================================================
std::unordered_map<std::bitset<64>, std::vector<std::bitset<64 >>>

octree::create_nb_hash_map() {

    // initialize a hash map leaf id->neighbour leaf ids
    std::unordered_map<std::bitset<64>, std::vector<std::bitset<64 >>> hash_nb_leafs;
    llist *tmp;

    if (out != nullptr && out->next != nullptr) {
        tmp = out->next;
        while (tmp != nullptr) {
            hash_nb_leafs[tmp->data->mcode] = std::vector<std::bitset<64 >>();
            tmp = tmp->next;
        }
    }

    return hash_nb_leafs;

}

//=============================================================================
// @ARG: none
//=============================================================================
std::unordered_map<std::bitset<64>, Voxel *>

octree::create_octree_hash_map() {

#ifdef PRINT
    std::cout << "\tCreate hashmap of octree...";
#endif
    // create hash map id->voxel
    std::unordered_map<std::bitset<64>, Voxel *> hash_octree;
    llist *tmp;        // stack containing voxels to be processed

    if (out != nullptr && out->next != nullptr) {
        tmp = out->next;
        while (tmp != nullptr) {
            hash_octree.insert({tmp->data->mcode, tmp->data});
            tmp = tmp->next;
        }
    }
#ifdef PRINT
    std::cout << "done" << std::endl;
#endif
    return hash_octree;
}

//=============================================================================
// @ARG: morton code of node
// @ARG: position of bit in searched direction
//=============================================================================
bool octree::__isBoundary(std::bitset<64> s, unsigned int k) {
// id is boundary id if all bits in specific direction are all the same (except leading "111" of root)
    unsigned int ns = __get_length_of_id(s) / 3; // number of sets

    for (unsigned int i = 0; i < (2 - k); i++)
        s = s >> 1; // shift until bit with desired direction is the most right

    unsigned int last = s[0];
    for (unsigned int i = 0; i < ns - 1; i++) { // compare all bits of desired direction (without 111 of root)
        unsigned int current = s[3 * i];
        if (current != last)
            return false;
        else
            last = current;
    }

    return true;

}

//=============================================================================
// @ARG: none
//=============================================================================
std::vector<std::bitset<64>> octree::get_all_boundary__leaf_nodes() {

    std::vector<std::bitset<64>> b_nodes;

    llist *tmp;
    if (out != nullptr && out->next != nullptr) {
        tmp = out->next;
        while (tmp != nullptr) {

            if (__isBoundary(tmp->data->mcode, 0))
                b_nodes.push_back(tmp->data->mcode);
            else if (__isBoundary(tmp->data->mcode, 1))
                b_nodes.push_back(tmp->data->mcode);
            else if (__isBoundary(tmp->data->mcode, 2))
                b_nodes.push_back(tmp->data->mcode);

            tmp = tmp->next;
        }
    }

    return b_nodes;
}

//=============================================================================
// @ARG: morton code of node
//=============================================================================
int octree::__get_length_of_id(std::bitset<64> s) {

    bool bit = false;
    unsigned int count = 0;
    unsigned int last = s.size() - 1;

    while (!bit) {
        s = s << 1;
        count++;
        if (s[last] == 1) bit = true;
    }
    return s.size() - count;
}

//=============================================================================
// @ARG: morton code of node
//=============================================================================
int octree::__get_depth_by_id(std::bitset<64> s) {

    return (__get_length_of_id(s)) / 3 - 1; // root is "000" but depth 0

}

//=============================================================================
// @ARG: voxel data
// @ARG: face index
//=============================================================================
int octree::__isinNode(const Voxel *voxel, const int idx) {

    int vec[3]; // indices of face's vertices
    vec[0] = Mod->flist[idx].v1;
    vec[1] = Mod->flist[idx].v2;
    vec[2] = Mod->flist[idx].v3;

    Vertex v1 = Mod->vlist[vec[0]];
    Vertex v2 = Mod->vlist[vec[1]];
    Vertex v3 = Mod->vlist[vec[2]];

    //############### AABB-AABB Test ###########
    double xmin = std::min({v1.x, v2.x, v3.x});
    double ymin = std::min({v1.y, v2.y, v3.y});
    double zmin = std::min({v1.z, v2.z, v3.z});
    double xmax = std::max({v1.x, v2.x, v3.x});
    double ymax = std::max({v1.y, v2.y, v3.y});
    double zmax = std::max({v1.z, v2.z, v3.z});

    if (xmin > voxel->maxx) return OCTREE_ERROR;
    if (ymin > voxel->maxy) return OCTREE_ERROR;
    if (zmin > voxel->maxz) return OCTREE_ERROR;
    if (xmax < voxel->minx) return OCTREE_ERROR;
    if (ymax < voxel->miny) return OCTREE_ERROR;
    if (zmax < voxel->minz) return OCTREE_ERROR;
    //##########################################

    float triverts[3][3] = {{static_cast<float>(v1.x), static_cast<float>(v1.y), static_cast<float>(v1.z)},
                            {static_cast<float>(v2.x), static_cast<float>(v2.y), static_cast<float>(v2.z)},
                            {static_cast<float>(v3.x), static_cast<float>(v3.y), static_cast<float>(v3.z)}};

    float boxhalfsize[3] = {static_cast<float>(voxel->maxx - voxel->x), static_cast<float>(voxel->maxy - voxel->y), static_cast<float>(voxel->maxz - voxel->z)};
    float boxcenter[3] = {static_cast<float>(voxel->x), static_cast<float>(voxel->y), static_cast<float>(voxel->z)};

    if (akine_moeller_isct::triBoxOverlap(boxcenter, boxhalfsize, triverts))
        return OCTREE_COLOR_GREY;
    else
        return OCTREE_ERROR;
}

//=============================================================================
// @ARG: map of neighbours
// @ARG: map of octree
//=============================================================================
void octree::__get_facade_brep(std::unordered_map<std::bitset<64>, std::vector<std::bitset<64>>> &hash_nb_leafs, std::unordered_map<std::bitset<64>, Voxel *> &hash_octree) {

    // iterate through nodes, if nb is OCTREE_COLOR_BLACK, face between them is an outer face
    std::cout << "\tGet facade brep...";

    llist *tmp;
    std::vector<std::array<double, 3>> vertices;
    std::vector<std::array<size_t, 4>> rectangles;
    std::unordered_map<std::bitset<64>, Voxel *>::const_iterator ffit;

    if (out != nullptr && out->next != nullptr) {
        tmp = out->next;
        while (tmp != nullptr) {

            if (tmp->data->color == OCTREE_COLOR_GREY) { // only for facade boxes

                // append boundary faces if node is a boundary node
                for (unsigned int k = 0; k < 3; k++) // for x, y, and z-direction
                    if (__isBoundary(tmp->data->mcode, k)) { // check if boundary
                        unsigned int n = tmp->data->mcode[2 - k]; // find whether min or max direction

                        int nb_type = OCTREE_ERROR;
                        if (k == 0 && n == 0)
                            nb_type = OCTREE_WEST;
                        else if (k == 0 && n == 1)
                            nb_type = OCTREE_EAST;
                        else if (k == 1 && n == 0)
                            nb_type = OCTREE_SOUTH;
                        else if (k == 1 && n == 1)
                            nb_type = OCTREE_NORTH;
                        else if (k == 2 && n == 0)
                            nb_type = OCTREE_BOTTOM;
                        else if (k == 2 && n == 1)
                            nb_type = OCTREE_TOP;

                        __switch_neighbour_type(nb_type, vertices, rectangles, tmp);
                    }

                // append boundary faces if node has a neighbour node that lies outside of geometry
                std::vector<std::bitset<64 >> nbs = hash_nb_leafs[tmp->data->mcode]; // get neighbours
                for (auto mcode_nb: nbs) {
                    // morton code of neighbour
                    ffit = hash_octree.find(mcode_nb); // pointer on neighbour node
                    if (ffit->second->color == OCTREE_COLOR_BLACK) {// only for neighbours outside of geometry
                        const int nb_type = __get_neighbour_type(tmp->data->mcode, mcode_nb); // get nb type (direction)
                        __switch_neighbour_type(nb_type, vertices, rectangles, tmp
                        );
                    }
                }
            }
            tmp = tmp->next;
        }
    }


    /* std::vector<unsigned int> inverse;
     create_unique_and_inverse_list_of_vertices(vertices, inverse, 1.0e-12);
     update_faces_from_inverse(rectangles, inverse);*/

    std::ofstream v;
    v.open("Surface_Mesh_v.txt");
    for (auto &vertice: vertices)
        v << vertice[0] << "\t" << vertice[1] << "\t" << vertice[2] << "\t" << std::endl;
    v.close();


    // converting to trimesh
    std::vector<std::array<size_t, 3>> triangles;
    for (auto &rectangle: rectangles) {
        std::array<size_t, 3> t1{{rectangle[0], rectangle[1], rectangle[2]}};
        std::array<size_t, 3> t2{{rectangle[0], rectangle[2], rectangle[3]}};
        triangles.push_back(t1);
        triangles.push_back(t2);
    }

    std::ofstream f;
    f.open("Surface_Mesh_f.txt");
    for (auto &triangle: triangles)
        f << triangle[0] << "\t" << triangle[1] << "\t" << triangle[2] << "\t" << std::endl;
    f.close();

    /*
    for (unsigned int i = 0; i < vertices.size(); i++)
        std::cout << vertices[i][0] << "\t" << vertices[i][1] << "\t" << vertices[i][2] << "\t" << std::endl;

    for (unsigned int i = 0; i < rectangles.size(); i++)
        std::cout << rectangles[i][0] << "\t" << rectangles[i][1] << "\t" << rectangles[i][2] << "\t" << rectangles[i][3] << "\t" << std::endl;

    for (unsigned int i = 0; i < triangles.size(); i++)
        std::cout << triangles[i][0] << "\t" << triangles[i][1] << "\t" << triangles[i][2] << "\t" << std::endl;
    */

    std::cout << "done" << std::endl;

}
/*
//=============================================================================
// @ARG: reference on vertex list
// @ARG: reference on empty inverse list
// @ARG: tolerance for vertex merging (squared point distance)
//=============================================================================
std::vector <std::array<double, 3>> octree::create_unique_and_inverse_list_of_vertices(const std::vector <std::array<double, 3>> &vertices, std::vector<unsigned int> &inverse, const double &squared_tol) {

    std::vector<unsigned int> unique_indices;  // länge = anzahl vertex die einzigartig sind, liste gibt die indizes an in original liste
    //std::vector<unsigned int> inverse;  // länge = vertex liste, hält indizes aus unique liste
    std::vector <std::array<double, 3>> vrtcs;

    for (unsigned int idx = 0; idx < vertices.size(); idx++) {

        std::array<double, 3> v = vertices[idx]; // current vertex
        bool found = false;

        for (unsigned int idx_u = 0; idx_u < unique_indices.size(); idx_u++) { // search in unique vector

            unsigned int u = unique_indices[idx_u];
            double dist_sqrd = (v[0] - vertices[u][0]) * (v[0] - vertices[u][0]) + (v[1] - vertices[u][1]) * (v[1] - vertices[u][1]) + (v[2] - vertices[u][2]) * (v[2] - vertices[u][2]);

            if (fabs(dist_sqrd) < squared_tol) {
                found = true;
                inverse.push_back(idx_u);
                break;
            }
        }

        if (found==false) { // if no vertex was found create new entries in unique, inverse and vertices list
            unique_indices.push_back(idx);  // add index of the vertex to unique list
            inverse.push_back(unique_indices.size() - 1); // add index of current unique vertex to inverse   unique_indices.index(idx)
            vrtcs.push_back(v);
        }
    }

    return vrtcs;
}

//=============================================================================
// @ARG: reference on face list
// @ARG: reference on inverse list
//=============================================================================
void octree::update_faces_from_inverse(std::vector <std::array<unsigned int, 4>> &faces, const std::vector<unsigned int> &inverse) {

    for (unsigned int i = 0; i < faces.size(); i++) {
        faces[i][0] = inverse[faces[i][0]];
        faces[i][1] = inverse[faces[i][1]];
        faces[i][2] = inverse[faces[i][2]];
        faces[i][3] = inverse[faces[i][3]];
    }
}
*/
//=============================================================================
// @ARG: morton code of node
// @ARG: morton code of neighbour node
//=============================================================================
int octree::__get_neighbour_type(std::bitset<64> mcode, std::bitset<64> mcode_nb) {

    // shift until set starts with root (111)
    unsigned int last = mcode.size() - 1;

    for (unsigned int i = 0; i < mcode.size(); i++) {
        if (mcode[last] == 1)
            break;
        else
            mcode = mcode << 1;
    }
    for (unsigned int i = 0; i < mcode_nb.size(); i++) {
        if (mcode_nb[last] == 1)
            break;
        else
            mcode_nb = mcode_nb << 1;
    }

    for (unsigned int i = 0; i < 22; i++) {
        if (mcode[last] == 1 && mcode_nb[last] == 0)// western nb
            return OCTREE_WEST;
        else if (mcode[last] == 0 && mcode_nb[last] == 1)// eastern nb
            return OCTREE_EAST;
        else if (mcode[last - 1] == 1 && mcode_nb[last - 1] == 0)// south nb
            return OCTREE_SOUTH;
        else if (mcode[last - 1] == 0 && mcode_nb[last - 1] == 1)// north nb
            return OCTREE_NORTH;
        else if (mcode[last - 2] == 1 && mcode_nb[last - 2] == 0)// bottom nb
            return OCTREE_BOTTOM;
        else if (mcode[last - 2] == 0 && mcode_nb[last - 2] == 1)// top nb
            return OCTREE_TOP;
        else {
            mcode = mcode << 3;
            mcode_nb = mcode_nb << 3;
        }
    }

    return OCTREE_ERROR;

}

//=============================================================================
// @ARG: direction of neighbour (type)
// @ARG: global vertex list by reference to append on
// @ARG: global rectangle list by reference to append on
// @ARG: pointer on voxel list entry
//=============================================================================
void octree::__switch_neighbour_type(const int &nb_type, std::vector<std::array<double, 3>> &vertices, std::vector<std::array<size_t, 4>> &rectangles, llist *tmp) {

    switch (nb_type) {
        case OCTREE_ERROR:
            abort();
            break;
        case OCTREE_WEST: {
            // append the four vertices
            std::array<double, 3> v1{{tmp->data->minx * Mod->bbox_norm + Mod->bbox_diffx, tmp->data->miny * Mod->bbox_norm + Mod->bbox_diffy, tmp->data->minz * Mod->bbox_norm + Mod->bbox_diffz}};
            std::array<double, 3> v2{{tmp->data->minx * Mod->bbox_norm + Mod->bbox_diffx, tmp->data->miny * Mod->bbox_norm + Mod->bbox_diffy, tmp->data->maxz * Mod->bbox_norm + Mod->bbox_diffz}};
            std::array<double, 3> v3{{tmp->data->minx * Mod->bbox_norm + Mod->bbox_diffx, tmp->data->maxy * Mod->bbox_norm + Mod->bbox_diffy, tmp->data->maxz * Mod->bbox_norm + Mod->bbox_diffz}};
            std::array<double, 3> v4{{tmp->data->minx * Mod->bbox_norm + Mod->bbox_diffx, tmp->data->maxy * Mod->bbox_norm + Mod->bbox_diffy, tmp->data->minz * Mod->bbox_norm + Mod->bbox_diffz}};
            vertices.push_back(v1);
            vertices.push_back(v2);
            vertices.push_back(v3);
            vertices.push_back(v4);
            // append the face
            std::array<size_t, 4> rec{{vertices.size() - 4, vertices.size() - 3, vertices.size() - 2, vertices.size() - 1}};
            rectangles.push_back(rec);
            break;
        }
        case OCTREE_EAST: {
            // append the four vertices
            std::array<double, 3> v4{{tmp->data->maxx * Mod->bbox_norm + Mod->bbox_diffx, tmp->data->miny * Mod->bbox_norm + Mod->bbox_diffy, tmp->data->minz * Mod->bbox_norm + Mod->bbox_diffz}};
            std::array<double, 3> v3{{tmp->data->maxx * Mod->bbox_norm + Mod->bbox_diffx, tmp->data->miny * Mod->bbox_norm + Mod->bbox_diffy, tmp->data->maxz * Mod->bbox_norm + Mod->bbox_diffz}};
            std::array<double, 3> v2{{tmp->data->maxx * Mod->bbox_norm + Mod->bbox_diffx, tmp->data->maxy * Mod->bbox_norm + Mod->bbox_diffy, tmp->data->maxz * Mod->bbox_norm + Mod->bbox_diffz}};
            std::array<double, 3> v1{{tmp->data->maxx * Mod->bbox_norm + Mod->bbox_diffx, tmp->data->maxy * Mod->bbox_norm + Mod->bbox_diffy, tmp->data->minz * Mod->bbox_norm + Mod->bbox_diffz}};
            vertices.push_back(v1);
            vertices.push_back(v2);
            vertices.push_back(v3);
            vertices.push_back(v4);
            // append the face
            std::array<size_t, 4> rec{{vertices.size() - 4, vertices.size() - 3, vertices.size() - 2, vertices.size() - 1}};
            rectangles.push_back(rec);
            break;
        }
        case OCTREE_SOUTH: {
            // append the four vertices
            std::array<double, 3> v1{{tmp->data->minx * Mod->bbox_norm + Mod->bbox_diffx, tmp->data->miny * Mod->bbox_norm + Mod->bbox_diffy, tmp->data->minz * Mod->bbox_norm + Mod->bbox_diffz}};
            std::array<double, 3> v2{{tmp->data->maxx * Mod->bbox_norm + Mod->bbox_diffx, tmp->data->miny * Mod->bbox_norm + Mod->bbox_diffy, tmp->data->minz * Mod->bbox_norm + Mod->bbox_diffz}};
            std::array<double, 3> v3{{tmp->data->maxx * Mod->bbox_norm + Mod->bbox_diffx, tmp->data->miny * Mod->bbox_norm + Mod->bbox_diffy, tmp->data->maxz * Mod->bbox_norm + Mod->bbox_diffz}};
            std::array<double, 3> v4{{tmp->data->minx * Mod->bbox_norm + Mod->bbox_diffx, tmp->data->miny * Mod->bbox_norm + Mod->bbox_diffy, tmp->data->maxz * Mod->bbox_norm + Mod->bbox_diffz}};
            vertices.push_back(v1);
            vertices.push_back(v2);
            vertices.push_back(v3);
            vertices.push_back(v4);
            // append the face
            std::array<size_t, 4> rec{{vertices.size() - 4, vertices.size() - 3, vertices.size() - 2, vertices.size() - 1}};
            rectangles.push_back(rec);
            break;
        }
        case OCTREE_NORTH: {
            // append the four vertices
            std::array<double, 3> v4{{tmp->data->minx * Mod->bbox_norm + Mod->bbox_diffx, tmp->data->maxy * Mod->bbox_norm + Mod->bbox_diffy, tmp->data->minz * Mod->bbox_norm + Mod->bbox_diffz}};
            std::array<double, 3> v3{{tmp->data->maxx * Mod->bbox_norm + Mod->bbox_diffx, tmp->data->maxy * Mod->bbox_norm + Mod->bbox_diffy, tmp->data->minz * Mod->bbox_norm + Mod->bbox_diffz}};
            std::array<double, 3> v2{{tmp->data->maxx * Mod->bbox_norm + Mod->bbox_diffx, tmp->data->maxy * Mod->bbox_norm + Mod->bbox_diffy, tmp->data->maxz * Mod->bbox_norm + Mod->bbox_diffz}};
            std::array<double, 3> v1{{tmp->data->minx * Mod->bbox_norm + Mod->bbox_diffx, tmp->data->maxy * Mod->bbox_norm + Mod->bbox_diffy, tmp->data->maxz * Mod->bbox_norm + Mod->bbox_diffz}};
            vertices.push_back(v1);
            vertices.push_back(v2);
            vertices.push_back(v3);
            vertices.push_back(v4);
            // append the face
            std::array<size_t, 4> rec{{vertices.size() - 4, vertices.size() - 3, vertices.size() - 2, vertices.size() - 1}};
            rectangles.push_back(rec);
            break;
        }
        case OCTREE_BOTTOM: {
            // append the four vertices
            std::array<double, 3> v1{{tmp->data->minx * Mod->bbox_norm + Mod->bbox_diffx, tmp->data->miny * Mod->bbox_norm + Mod->bbox_diffy, tmp->data->minz * Mod->bbox_norm + Mod->bbox_diffz}};
            std::array<double, 3> v2{{tmp->data->minx * Mod->bbox_norm + Mod->bbox_diffx, tmp->data->maxy * Mod->bbox_norm + Mod->bbox_diffy, tmp->data->minz * Mod->bbox_norm + Mod->bbox_diffz}};
            std::array<double, 3> v3{{tmp->data->maxx * Mod->bbox_norm + Mod->bbox_diffx, tmp->data->maxy * Mod->bbox_norm + Mod->bbox_diffy, tmp->data->minz * Mod->bbox_norm + Mod->bbox_diffz}};
            std::array<double, 3> v4{{tmp->data->maxx * Mod->bbox_norm + Mod->bbox_diffx, tmp->data->miny * Mod->bbox_norm + Mod->bbox_diffy, tmp->data->minz * Mod->bbox_norm + Mod->bbox_diffz}};
            vertices.push_back(v1);
            vertices.push_back(v2);
            vertices.push_back(v3);
            vertices.push_back(v4);
            // append the face
            std::array<size_t, 4> rec{{vertices.size() - 4, vertices.size() - 3, vertices.size() - 2, vertices.size() - 1}};
            rectangles.push_back(rec);
            break;
        }
        case OCTREE_TOP: {
            // append the four vertices
            std::array<double, 3> v4{{tmp->data->minx * Mod->bbox_norm + Mod->bbox_diffx, tmp->data->miny * Mod->bbox_norm + Mod->bbox_diffy, tmp->data->maxz * Mod->bbox_norm + Mod->bbox_diffz}};
            std::array<double, 3> v3{{tmp->data->minx * Mod->bbox_norm + Mod->bbox_diffx, tmp->data->maxy * Mod->bbox_norm + Mod->bbox_diffy, tmp->data->maxz * Mod->bbox_norm + Mod->bbox_diffz}};
            std::array<double, 3> v2{{tmp->data->maxx * Mod->bbox_norm + Mod->bbox_diffx, tmp->data->maxy * Mod->bbox_norm + Mod->bbox_diffy, tmp->data->maxz * Mod->bbox_norm + Mod->bbox_diffz}};
            std::array<double, 3> v1{{tmp->data->maxx * Mod->bbox_norm + Mod->bbox_diffx, tmp->data->miny * Mod->bbox_norm + Mod->bbox_diffy, tmp->data->maxz * Mod->bbox_norm + Mod->bbox_diffz}};
            vertices.push_back(v1);
            vertices.push_back(v2);
            vertices.push_back(v3);
            vertices.push_back(v4);
            // append the face
            std::array<size_t, 4> rec{{vertices.size() - 4, vertices.size() - 3, vertices.size() - 2, vertices.size() - 1}};
            rectangles.push_back(rec);
            break;
        }
        default:
            abort();
    }

}

//=============================================================================
// @ARG: map of neighbours
// @ARG: map of octree
//=============================================================================
void octree::get_spaces(std::unordered_map<std::bitset<64>, std::vector<std::bitset<64>>> &hash_nb_leafs, std::unordered_map<std::bitset<64>, Voxel *> &hash_octree) {

    bool add_next_neighbours = true;

#ifdef PRINT
    std::cout << "\tGet spaces via flood fill...";
#endif

    // OCTREE_COLOR_WHITE       - not visited yet
    // OCTREE_COLOR_BLACK       - visited and empty
    // OCTREE_COLOR_GREY        - visited and non-empty
    // OCTREE_PREDICTED_WHITE   - visited, empty node of a space

    // After flood fill ...
    // OCTREE_COLOR_WHITE - not visited yet
    // OCTREE_COLOR_BLACK - already visited and empty node
    // OCTREE_COLOR_GREY - already visited and non-empty node
    // ... the environment (air) is OCTREE_COLOR_BLACK and facade nodes are OCTREE_COLOR_GREY

    //*******************************************************************************************
    // So first extract facade space.
    llist *tmp;
    std::set<std::string> facade_guids;
    std::set<unsigned int> facade_triangles;
    double min_dist = 1000;

    if (fluid_points.empty()) { // if fluid points were defined, only consider spaces points lie in
        if (out != nullptr && out->next != nullptr) {
            tmp = out->next;
            while (tmp != nullptr) {
                if (tmp->data->color == OCTREE_COLOR_GREY) {
                    save_voxel_triangle_info_in_sets(tmp->data, facade_triangles, facade_guids);

                    if (add_next_neighbours) // add neighbours too fill gaps
                        save_voxel_neighbors_triangle_info_in_sets(tmp->data, facade_triangles, facade_guids, hash_nb_leafs, hash_octree);

                } else if (tmp->data->color == OCTREE_COLOR_BLACK) { // calculate distance to -1| -1| -1
                    double dist = (tmp->data->minx + 1) * (tmp->data->minx + 1) + (tmp->data->miny + 1) * (tmp->data->miny + 1) + (tmp->data->minz + 1) * (tmp->data->minz + 1);
                    if (dist < min_dist) min_dist = dist;
                }
                tmp = tmp->next;
            }
        }
        std::array<double, 3> facade_mid{Mod->maxx * Mod->bbox_norm + Mod->bbox_diffx, Mod->maxy * Mod->bbox_norm + Mod->bbox_diffy, Mod->maxz * Mod->bbox_norm + Mod->bbox_diffz};
        spaces.emplace_back(facade_mid, facade_guids);
        spaces_tri.emplace_back(facade_mid, facade_triangles);
        spaces_minima.push_back(min_dist);
    } else {
        bool contains_fluid_point = false;

        if (out != nullptr && out->next != nullptr) {
            tmp = out->next;
            while (tmp != nullptr) {
                if (tmp->data->color == OCTREE_COLOR_GREY) {
                    save_voxel_triangle_info_in_sets(tmp->data, facade_triangles, facade_guids);

                    if (add_next_neighbours) // add neighbours too fill gaps
                        save_voxel_neighbors_triangle_info_in_sets(tmp->data, facade_triangles, facade_guids, hash_nb_leafs, hash_octree);

                } else if (tmp->data->color == OCTREE_COLOR_BLACK) { // calculate distance to -1| -1| -1
                    double dist = (tmp->data->minx + 1) * (tmp->data->minx + 1) + (tmp->data->miny + 1) * (tmp->data->miny + 1) + (tmp->data->minz + 1) * (tmp->data->minz + 1);
                    if (dist < min_dist) min_dist = dist;

                    if (!contains_fluid_point) // check if air voxel contains any fluid point
                        if (voxel_contains_fluid_point(tmp->data))
                            contains_fluid_point = true;
                }
                tmp = tmp->next;
            }
        }
        if (contains_fluid_point) {
            std::array<double, 3> facade_mid{Mod->maxx * Mod->bbox_norm + Mod->bbox_diffx, Mod->maxy * Mod->bbox_norm + Mod->bbox_diffy, Mod->maxz * Mod->bbox_norm + Mod->bbox_diffz};
            spaces.emplace_back(facade_mid, facade_guids);
            spaces_tri.emplace_back(facade_mid, facade_triangles);
            spaces_minima.push_back(min_dist);
        }
    }
    //*******************************************************************************************

    //*******************************************************************************************
    // Fill all remaining triangle nodes OCTREE_COLOR_GREY. All interior air nodes will remain as OCTREE_COLOR_WHITE
    set_non_empty_white_voxel_to_grey();
    //*******************************************************************************************

    //*******************************************************************************************
    // for all OCTREE_COLOR_WHITE nodes, find biggest node ...
    if (fluid_points.empty()) { // if fluid points were defined, only consider spaces points lie in
        if (out != nullptr && out->next != nullptr) {
            tmp = out->next;
            while (tmp != nullptr) {
                if (tmp->data->color == OCTREE_COLOR_WHITE) {
                    double dx = (tmp->data->maxx - tmp->data->minx) * Mod->bbox_norm; // scale back for criteria in metres
                    if (dx > 0.5) { // criteria a air node is considered as part of a space and not as a gap or something
                        // if voxel is big enough, start flood fill until we reach the surrounding building elements (OCTREE_COLOR_GREY)

                        // Flood Fill for Space starting from current Node
                        std::unordered_map<std::bitset<64>, Voxel *>::const_iterator ffit;
                        std::stack<std::bitset<64 >> branches_loop;
                        branches_loop.push(tmp->data->mcode);

                        double min_dist = 1000;
                        std::set<std::string> space_guids;
                        std::set<unsigned int> space_triangles;

                        // Flood Fill ...
                        while (!branches_loop.empty()) {
                            ffit = hash_octree.find(branches_loop.top());
                            branches_loop.pop();

                            if (ffit->second->color == OCTREE_COLOR_WHITE) { // its "air"
                                ffit->second->color = OCTREE_PREDICTED_WHITE; // set space color
                                auto nbs = hash_nb_leafs[ffit->second->mcode];
                                for (auto &nb: nbs)
                                    branches_loop.push(nb);

                                // calculate distance to -1| -1| -1
                                double dist = (ffit->second->minx + 1) * (ffit->second->minx + 1) + (ffit->second->miny + 1) * (ffit->second->miny + 1) + (ffit->second->minz + 1) * (ffit->second->minz + 1);
                                if (dist < min_dist) min_dist = dist;

                            } else if (ffit->second->color == OCTREE_COLOR_GREY) {
                                save_voxel_triangle_info_in_sets(ffit->second, space_triangles, space_guids);

                                if (add_next_neighbours) // add neighbours too fill gaps
                                    save_voxel_neighbors_triangle_info_in_sets(ffit->second, space_triangles, space_guids, hash_nb_leafs, hash_octree);
                            }
                        }

                        // save space defined by guids of surrounding products in octree
                        std::array<double, 3> mid{tmp->data->x * Mod->bbox_norm + Mod->bbox_diffx, tmp->data->y * Mod->bbox_norm + Mod->bbox_diffy, tmp->data->z * Mod->bbox_norm + Mod->bbox_diffz};
                        spaces.emplace_back(mid, space_guids);
                        spaces_tri.emplace_back(mid, space_triangles);
                        spaces_minima.push_back(min_dist);
                    }
                }
                tmp = tmp->next;
            }
        }
    } else {
        if (out != nullptr && out->next != nullptr) {
            tmp = out->next;
            while (tmp != nullptr) {
                if (tmp->data->color == OCTREE_COLOR_WHITE) {
                    double dx = (tmp->data->maxx - tmp->data->minx) * Mod->bbox_norm; // scale back for criteria in metres
                    if (dx > 0.5) { // criteria a air node is considered as part of a space and not as a gap or something
                        // if voxel is big enough, start flood fill until we reach the surrounding building elements (OCTREE_COLOR_GREY)

                        // Flood Fill for Space starting from current Node
                        std::unordered_map<std::bitset<64>, Voxel *>::const_iterator ffit;
                        std::stack<std::bitset<64 >> branches_loop;
                        branches_loop.push(tmp->data->mcode);

                        double min_dist = 1000;
                        std::set<std::string> space_guids;
                        std::set<unsigned int> space_triangles;

                        bool contains_fluid_point = false;

                        // Flood Fill ...
                        while (!branches_loop.empty()) {
                            ffit = hash_octree.find(branches_loop.top());
                            branches_loop.pop();

                            if (ffit->second->color == OCTREE_COLOR_WHITE) { // its "air"
                                ffit->second->color = OCTREE_PREDICTED_WHITE; // set space color
                                auto nbs = hash_nb_leafs[ffit->second->mcode];
                                for (auto &nb: nbs)
                                    branches_loop.push(nb);

                                // calculate distance to -1| -1| -1
                                double dist = (ffit->second->minx + 1) * (ffit->second->minx + 1) + (ffit->second->miny + 1) * (ffit->second->miny + 1) + (ffit->second->minz + 1) * (ffit->second->minz + 1);
                                if (dist < min_dist) min_dist = dist;

                                if (!contains_fluid_point) // check if air voxel contains any fluid point
                                    if (voxel_contains_fluid_point(ffit->second))
                                        contains_fluid_point = true;

                            } else if (ffit->second->color == OCTREE_COLOR_GREY) {
                                save_voxel_triangle_info_in_sets(ffit->second, space_triangles, space_guids);

                                if (add_next_neighbours) // add neighbours too fill gaps
                                    save_voxel_neighbors_triangle_info_in_sets(ffit->second, space_triangles, space_guids, hash_nb_leafs, hash_octree);
                            }
                        }

                        // save space defined by guids of surrounding products in octree
                        if (contains_fluid_point) {
                            std::array<double, 3> mid{tmp->data->x * Mod->bbox_norm + Mod->bbox_diffx, tmp->data->y * Mod->bbox_norm + Mod->bbox_diffy, tmp->data->z * Mod->bbox_norm + Mod->bbox_diffz};
                            spaces.emplace_back(mid, space_guids);
                            spaces_tri.emplace_back(mid, space_triangles);
                            spaces_minima.push_back(min_dist);
                        }
                    }
                }
                tmp = tmp->next;
            }
        }
    }
    //*******************************************************************************************

    //*******************************************************************************************
    // sort spaces by distance to origin, minima vector is not sorted!
    std::multimap<double, unsigned int> M; // store distance to origin

    for (unsigned int i = 0; i < spaces.size(); i++)
        M.insert(std::pair<double, unsigned int>(spaces_minima[i], i));

    std::vector<std::pair<std::array<double, 3>, std::set<std::string>>> new_spaces;
    std::vector<std::pair<std::array<double, 3>, std::set<unsigned int>>> new_s_triangles;

    new_spaces.reserve(spaces.size());
    new_s_triangles.reserve(spaces.size());

    for (const auto &m: M) {
        new_spaces.emplace_back(spaces[m.second]);
        new_s_triangles.emplace_back(spaces_tri[m.second]);
    }

    spaces = new_spaces;
    spaces_tri = new_s_triangles;
    //*******************************************************************************************

    //*******************************************************************************************
    auto nbs = find_space_neighbours();
    find_zones(nbs, max_triangles_in_zone);
    //*******************************************************************************************

#ifdef PRINT
    std::cout << "done" << std::endl;
#endif
}

void octree::save_voxel_neighbors_triangle_info_in_sets(const Voxel *voxel, std::set<unsigned int> &triangles, std::set<std::string> &guids, std::unordered_map<std::bitset<64>, std::vector<std::bitset<64>>> &hash_nb_leafs, std::unordered_map<std::bitset<64>, Voxel *> &hash_octree) {

    auto mocode_nbs = hash_nb_leafs[voxel->mcode];
    for (const auto &mcode_nb: mocode_nbs) {
        auto nb = hash_octree.find(mcode_nb)->second;
        if (nb->color == OCTREE_COLOR_GREY)
            save_voxel_triangle_info_in_sets(nb, triangles, guids);
    }
}

void octree::save_voxel_triangle_info_in_sets(const Voxel *voxel, std::set<unsigned int> &triangles, std::set<std::string> &guids) {

    for (int i = 0; i < voxel->nface; ++i) {
        unsigned int index_tri = voxel->findex[i];
        guids.insert(Mod->flist[index_tri].guid);
        triangles.insert(index_tri);
    }
}

void octree::set_non_empty_white_voxel_to_grey() {

    llist *tmp;
    if (out != nullptr && out->next != nullptr) {
        tmp = out->next;
        while (tmp != nullptr) {
            if (tmp->data->color == OCTREE_COLOR_WHITE)
                if (tmp->data->nface > 0)
                    tmp->data->color = OCTREE_COLOR_GREY;
            tmp = tmp->next;
        }
    }
}

bool octree::voxel_contains_fluid_point(const Voxel *voxel) {

    for (const auto P: fluid_points)
        if (is_point_in_voxel(voxel, P))
            return true;
    return false;
}

std::vector<std::pair<std::array<double, 3>, std::set<std::string>>> octree::get_spaces_guids() { return spaces; }

std::vector<std::pair<std::array<double, 3>, std::set<unsigned int>>> octree::get_spaces_triangles() { return spaces_tri; }

std::vector<std::set<unsigned int>> octree::get_zones() { return zones; }

std::vector<std::list<unsigned int>> octree::find_space_neighbours() {

    std::vector<std::set<unsigned int>> nbs(spaces.size());

    // lookup table GUID to spaces, the GUID is used in
    std::map<std::string, std::vector<unsigned int>> temp;
    for (unsigned int i = 0; i < spaces.size(); i++)
        for (const auto &guid: spaces[i].second)
            temp[guid].push_back(i);

    // get neighbours of each space
    for (const auto &t: temp)
        for (const auto &space: t.second) {
            for (const auto &nb: t.second)
                nbs[space].insert(nb);
            nbs[space].erase(space); // remove self reference
        }

    std::vector<std::list<unsigned int>> nbs_sorted(spaces.size());
    for (int i = 0; i < nbs.size(); i++)
        for (const auto &n: nbs[i])
            nbs_sorted[i].push_back(n);

    return nbs_sorted;
}

void octree::find_zones(const std::vector<std::list<unsigned int>> &nbs, unsigned int max_tri) {

    std::vector<int> space2zone(spaces.size(), -1);
    std::vector<std::set<unsigned int>> zones_content; // count triangles

    for (unsigned int i = 0; i < spaces.size(); i++) {

        // skip if already belonging to zone
        if (space2zone[i] != -1) continue;

        // add start face to new zone
        zones.emplace_back();
        zones.back().insert(i);
        space2zone[i] = zones.size() - 1;
        std::set<unsigned int> visited = {i};

        // create guid/triangle content of zone
        zones_content.emplace_back();
        zones_content.back().insert(spaces_tri[i].second.begin(), spaces_tri[i].second.end());

        // add nbs to deque
        std::deque<unsigned int> face_deque;
        for (const auto &n: nbs[i])
            face_deque.push_back(n);

        while (!face_deque.empty()) {

            unsigned int space = face_deque.front();
            face_deque.pop_front();

            // skip space if already belonging to a zone or visited in this cluster round
            if (space2zone[space] != -1 || visited.find(space) != visited.end()) continue;

            // check content of zone, if we add the space
            auto t = zones_content.back();
            t.insert(spaces_tri[space].second.begin(), spaces_tri[space].second.end());

            // add space to zone, if max content is not reached
            if (t.size() < max_tri) {
                zones.back().insert(space);
                space2zone[space] = zones.size() - 1;
                visited.insert(space);
                zones_content.back().insert(spaces_tri[space].second.begin(), spaces_tri[space].second.end());

                for (unsigned int j = 0; j < nbs[space].size(); j++)
                    face_deque.push_back(j);
            }
        }
    }
}

void octree::set_maximum_triangles_in_a_zone(unsigned int n) { max_triangles_in_zone = n; }

bool octree::is_point_in_voxel(const Voxel *voxel, const std::array<double, 3> &P) {

    if (P[0] > voxel->maxx) return false;
    if (P[1] > voxel->maxy) return false;
    if (P[2] > voxel->maxz) return false;
    if (P[0] < voxel->minx) return false;
    if (P[1] < voxel->miny) return false;
    if (P[2] < voxel->minz) return false;
    return true;
}

void octree::set_fluid_points(const std::set<std::array<double, 3>> &_fluid_points) {
    for (const auto &P: _fluid_points) {
        auto x = (P[0] - Mod->bbox_diffx) / Mod->bbox_norm;
        auto y = (P[1] - Mod->bbox_diffy) / Mod->bbox_norm;
        auto z = (P[2] - Mod->bbox_diffz) / Mod->bbox_norm;

        // in case fluid point is outside of the model, shift it to the border
        // this can happen if fluidpoint of coarse octree is used in a finer octree
        // (voxel width at the border is smaller)
        if (x < -1) x = -1;
        if (x > 1) x = 1;
        if (y < -1) y = -1;
        if (y > 1) y = 1;
        if (z < -1) z = -1;
        if (z > 1) z = 1;

        fluid_points.push_back({x, y, z});
    }
}

bool octree::isZero(const double v) { return v >= -1.0e-9 && v <= 1.0e-9; }