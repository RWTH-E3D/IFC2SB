// Copyright 2022 Eric Fichter
#include "RayIntersector.h"
#include "../extern/aabb_triang_akenine-moeller.h"

//=============================================================================
// @ARG: none
//=============================================================================
RayIntersector::RayIntersector() {
    in = nullptr;
    out = nullptr;
    Mod = nullptr;
    vtkfile = nullptr;
    root_voxel = nullptr;
    point_in_polygon_treshold = 1.0e-9;
}

//=============================================================================
// @ARG: none
//=============================================================================
RayIntersector::~RayIntersector() {
    cleanup_stack(in);
    cleanup_stack(out);
    Mod = nullptr;
    vtkfile = nullptr;
    root_voxel = nullptr;
}

//=============================================================================
// @ARG: none
//=============================================================================
void RayIntersector::set_model(IntersectorModel *Mod0) { if (Mod0 != nullptr) { Mod = Mod0; }}

//=============================================================================
// @ARG: filename to write octree to (if not set, nothing will be written out)
//=============================================================================
void RayIntersector::set_filename(char *vtkfile0) { vtkfile = vtkfile0; }

//=============================================================================
// @ARG: maximum depth of recursion
//=============================================================================
int RayIntersector::gen_octree(const int maxdepth, const unsigned int min_triangles_for_division) {
    Voxel2 *voxel;            // voxel to be processed
    int count;                // counter for triangles to be excluded
    int i, j;                // loop counters
    int *findex;            // indices of faces left
    unsigned char color;    // voxel's color from orientation test

    /* check for model and initialise IN stack */
    if (Mod == nullptr) { return OCTREE_ERROR; }
    if (in == nullptr) { setup_stack(); }

    /* compute for all faces its bbox and distance to origin */
    comp_bbox();

#ifndef MOELLER_INTERSECTION_TEST
    comp_dist();
#endif

    while ((voxel = pop(in)) != nullptr) {

        count = 0;

        /* quick bounding box test */
        //std::cerr << voxel->nface << std::endl;
        for (i = 0; i < voxel->nface; ++i) {
            if (bbox_test(voxel, voxel->findex[i]) == OCTREE_DONE) {
                voxel->findex[i] |= 0x80000000;
                ++count;
            }
        }

        /* check resulting triangles (if so) for intersection. if at least one triangle in node -> set to grey */
        if (count != voxel->nface) {
            for (i = 0; i < voxel->nface; ++i) {
                /* only if voxel has not been marked by __bbox_test() */
                if (!(voxel->findex[i] & 0x80000000)) {
                    /* quick vertex in voxel test */
                    if (vertex_in_voxel_test(voxel, voxel->findex[i]) == OCTREE_COLOR_GREY) {
                        voxel->color = OCTREE_COLOR_GREY;
                    } else {
                        /* expensive intersection test */
#ifdef MOELLER_INTERSECTION_TEST
                        if (triangle_aabb_moeller(voxel, voxel->findex[i]) == OCTREE_COLOR_GREY) {
#else
                            if (voxel_intersection_test(voxel, voxel->findex[i]) == OCTREE_COLOR_GREY) {
#endif
                            voxel->color = OCTREE_COLOR_GREY;
                        } else {
                            voxel->findex[i] |= 0x80000000;
                            ++count;
                        }
                    }
                }
            }
        }

        /* no triangle in voxel -> set to white */
        if (count == voxel->nface)
            voxel->color = OCTREE_COLOR_WHITE;

        /* push voxel back to IN or OUT stack */
        if (voxel->color != OCTREE_COLOR_GREY) {
            delete[] voxel->findex;
            voxel->nface = 0;
            voxel->findex = nullptr;
            push(out, voxel);
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
            if (voxel->d == maxdepth || voxel->nface < min_triangles_for_division) { push(out, voxel); }
            else {
                subdivide_voxel(voxel);
                push(out, voxel);
            }
        }
    } /* END while */

    /* write octree as uniform grid to file */
    if (vtkfile != nullptr)
        write_adaptive_vtk_file_only_leafs();

    return OCTREE_DONE;
}

//=============================================================================
// @ARG: linked list for fetching (i.e. deleting) voxel
//=============================================================================
Voxel2 *RayIntersector::pop(llist2 *list) {
    Voxel2 *voxel;    // information of next voxel
    llist2 *tmp;        // next stack element

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
void RayIntersector::push(llist2 *list, Voxel2 *voxel) {
    llist2 *tmp;    // new stack element

    tmp = new llist2;
    tmp->data = voxel;
    tmp->next = list->next;
    list->next = tmp;
}

//=============================================================================
// @ARG: voxel data
//=============================================================================
void RayIntersector::subdivide_voxel(Voxel2 *voxel) {
    //double radius = pow(2.0, -(voxel->d));	// children's radius, i.e. 2^(-d)
    double radius = (voxel->maxx - voxel->minx) * 0.25;    // children's radius, i.e. 2^(-d)
    double xyz[3];                            // children's midpoint

    /* compute children's radius */

    /* compute children's midpoint and generate Voxel2 */
    xyz[0] = voxel->x - radius;
    xyz[1] = voxel->y - radius;
    xyz[2] = voxel->z - radius;    /* BNW */
    auto child1 = generate_child_voxel(voxel, xyz, radius);

    xyz[0] = voxel->x + radius;
    xyz[1] = voxel->y - radius;
    xyz[2] = voxel->z - radius;    /* BNE */
    auto child2 = generate_child_voxel(voxel, xyz, radius);

    xyz[0] = voxel->x - radius;
    xyz[1] = voxel->y - radius;
    xyz[2] = voxel->z + radius;    /* BSW */
    auto child3 = generate_child_voxel(voxel, xyz, radius);

    xyz[0] = voxel->x + radius;
    xyz[1] = voxel->y - radius;
    xyz[2] = voxel->z + radius;    /* BSE */
    auto child4 = generate_child_voxel(voxel, xyz, radius);

    xyz[0] = voxel->x - radius;
    xyz[1] = voxel->y + radius;
    xyz[2] = voxel->z - radius;    /* TNW */
    auto child5 = generate_child_voxel(voxel, xyz, radius);

    xyz[0] = voxel->x + radius;
    xyz[1] = voxel->y + radius;
    xyz[2] = voxel->z - radius;    /* TNE */
    auto child6 = generate_child_voxel(voxel, xyz, radius);

    xyz[0] = voxel->x - radius;
    xyz[1] = voxel->y + radius;
    xyz[2] = voxel->z + radius;    /* TSW */
    auto child7 = generate_child_voxel(voxel, xyz, radius);

    xyz[0] = voxel->x + radius;
    xyz[1] = voxel->y + radius;
    xyz[2] = voxel->z + radius;    /* TSE */
    auto child8 = generate_child_voxel(voxel, xyz, radius);

    voxel->children[0] = child1;
    voxel->children[1] = child2;
    voxel->children[2] = child3;
    voxel->children[3] = child4;
    voxel->children[4] = child5;
    voxel->children[5] = child6;
    voxel->children[6] = child7;
    voxel->children[7] = child8;

    push(in, child1);
    push(in, child2);
    push(in, child3);
    push(in, child4);
    push(in, child5);
    push(in, child6);
    push(in, child7);
    push(in, child8);

    /* delete parent voxel */
//    delete[] voxel->findex;
//    voxel->findex = nullptr;
//    delete voxel;
//    voxel = nullptr;

    /* delete info from parent voxel */
    delete[] voxel->findex;
    voxel->findex = nullptr;

}

//=============================================================================
// @ARG: parent data
// @ARG: children's midpoint
// @ARG: children's radius
//=============================================================================
Voxel2 *RayIntersector::generate_child_voxel(const Voxel2 *parent, const double *xyz, const double radius) {

    Voxel2 *child;    // new voxel
    int i;            // loop counter

    child = new Voxel2;

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
    child->d = parent->d + 1;
    child->findex = new int[child->nface];
    for (i = 0; i < child->nface; ++i) { child->findex[i] = parent->findex[i]; }
    for (i = 0; i < 8; i++) { child->children[i] = nullptr; }

    return child;
}

//=============================================================================
// @ARG: voxel data
// @ARG: face index
//=============================================================================
int RayIntersector::bbox_test(const Voxel2 *voxel, const int idx) {
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
int RayIntersector::vertex_in_voxel_test(const Voxel2 *voxel, const int idx) {
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


#ifndef MOELLER_INTERSECTION_TEST
//=============================================================================
// @ARG: voxel data
// @ARG: face index
//=============================================================================
int RayIntersector::voxel_intersection_test(const Voxel2 *voxel, const int idx) {
    /* check triangle's edges for intersection with voxel's faces */
    if (line_intersects_YZ_plane(voxel, idx, OCTREE_SIDE_LEFT) == OCTREE_DONE) { return OCTREE_COLOR_GREY; }
    if (line_intersects_YZ_plane(voxel, idx, OCTREE_SIDE_RIGHT) == OCTREE_DONE) { return OCTREE_COLOR_GREY; }
    if (line_intersects_XZ_plane(voxel, idx, OCTREE_SIDE_BOTTOM) == OCTREE_DONE) { return OCTREE_COLOR_GREY; }
    if (line_intersects_XZ_plane(voxel, idx, OCTREE_SIDE_TOP) == OCTREE_DONE) { return OCTREE_COLOR_GREY; }
    if (line_intersects_XY_plane(voxel, idx, OCTREE_SIDE_BACK) == OCTREE_DONE) { return OCTREE_COLOR_GREY; }
    if (line_intersects_XY_plane(voxel, idx, OCTREE_SIDE_FRONT) == OCTREE_DONE) { return OCTREE_COLOR_GREY; }

    /* check voxel's diagonals for intersection with triangle */
    if (diagonal_intersects_polygon(voxel, idx) == OCTREE_DONE) { return OCTREE_COLOR_GREY; }

    return OCTREE_ERROR;
}

//=============================================================================
// @ARG: voxel data
// @ARG: face index
// @ARG: voxel's left or right side
//=============================================================================
int RayIntersector::line_intersects_YZ_plane(const Voxel2 *voxel, const int idx, const int side) {
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
int RayIntersector::line_intersects_XZ_plane(const Voxel2 *voxel, const int idx, const int side) {
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
int RayIntersector::line_intersects_XY_plane(const Voxel2 *voxel, const int idx, const int side) {
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
int RayIntersector::diagonal_intersects_polygon(const Voxel2 *voxel, const int idx) {
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
                if (point_in_polygon(SX, SY, SZ, idx) == OCTREE_DONE) { return OCTREE_DONE; }
            }
        }
    }

    return OCTREE_ERROR;
}

//=============================================================================
// @ARG: point's coordinates
// @ARG: face index
//=============================================================================
int RayIntersector::point_in_polygon(const double SX, const double SY, const double SZ, const int idx) {
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
    if (isZero(Mod->flist[idx].A)) { Mod->flist[idx].A = comp_polygon_area(v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z); }
    A_diff = Mod->flist[idx].A;

    /* compute areas of smaller polygons (v1,v2,S), (v2,v3,S), (v3,v1,S) */
    A_diff -= comp_polygon_area(v1x, v1y, v1z, v2x, v2y, v2z, SX, SY, SZ);
    A_diff -= comp_polygon_area(v2x, v2y, v2z, v3x, v3y, v3z, SX, SY, SZ);
    A_diff -= comp_polygon_area(v3x, v3y, v3z, v1x, v1y, v1z, SX, SY, SZ);

    /* check if areas have 'equal' size */
    if (fabs(A_diff) < point_in_polygon_treshold) { return OCTREE_DONE; }

    return OCTREE_ERROR;
}

//=============================================================================
// @ARG: vertex v1's coordinates
// @ARG: vertex v2's coordinates
// @ARG: vertex v3's coordinates
//=============================================================================
double RayIntersector::comp_polygon_area(const double v1x, const double v1y, const double v1z,
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
#endif

//=============================================================================
// @ARG: none
//=============================================================================
void RayIntersector::comp_bbox() {
    int i;        // loop counter

    for (i = 0; i < Mod->nface; ++i) {
        Mod->flist[i].minx = Mod->vlist[Mod->flist[i].v1].x;
        Mod->flist[i].maxx = Mod->vlist[Mod->flist[i].v1].x;
        Mod->flist[i].miny = Mod->vlist[Mod->flist[i].v1].y;
        Mod->flist[i].maxy = Mod->vlist[Mod->flist[i].v1].y;
        Mod->flist[i].minz = Mod->vlist[Mod->flist[i].v1].z;
        Mod->flist[i].maxz = Mod->vlist[Mod->flist[i].v1].z;
        check_vertex(i, Mod->vlist[Mod->flist[i].v2].x, Mod->vlist[Mod->flist[i].v2].y, Mod->vlist[Mod->flist[i].v2].z);
        check_vertex(i, Mod->vlist[Mod->flist[i].v3].x, Mod->vlist[Mod->flist[i].v3].y, Mod->vlist[Mod->flist[i].v3].z);
    }
}

//=============================================================================
// @ARG: none
//=============================================================================
#ifndef MOELLER_INTERSECTION_TEST
void RayIntersector::comp_dist() {
    int i;        // loop counter

    for (i = 0; i < Mod->nface; ++i) {
        Mod->flist[i].dist = Mod->flist[i].nx * Mod->vlist[Mod->flist[i].v1].x +
                             Mod->flist[i].ny * Mod->vlist[Mod->flist[i].v1].y +
                             Mod->flist[i].nz * Mod->vlist[Mod->flist[i].v1].z;
    }
}
#endif

//=============================================================================
// @ARG: face index
// @ARG: vertex coordinates to be processed
//=============================================================================
void RayIntersector::check_vertex(const int idx, const double x, const double y, const double z) {
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
void RayIntersector::setup_stack(void) {
    Voxel2 *voxel;    // voxel data
    int i;            // loop counter

    voxel = new Voxel2;
    voxel->d = 0;
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
    for (i = 0; i < 8; i++) { voxel->children[i] = nullptr; }

    root_voxel = voxel;

    in = new llist2;
    in->data = nullptr;
    in->next = nullptr;
    out = new llist2;
    out->data = nullptr;
    out->next = nullptr;
    push(in, voxel);
    voxel = nullptr;
}

//=============================================================================
// @ARG: linked list
//=============================================================================
void RayIntersector::cleanup_stack(llist2 *list) {
    Voxel2 *voxel;    // voxel data

    if (list != nullptr) {
        /* delete stack elements */
        while ((voxel = pop(list)) != nullptr) {
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
void RayIntersector::write_adaptive_vtk_file() {

    std::ofstream os;

    /* write grid to file */
    os.open(vtkfile, std::ofstream::out | std::ofstream::trunc);
    if (!os.is_open()) {
        std::cerr << "ERROR >> couldn't write file " << vtkfile << std::endl;
        exit(0);
    }

    int ct_pts = 0, ct_vxls = 0, ct = 0; // counters
    llist2 *tmp;         // next stack element
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

#ifdef PRINT
    std::cout << "\tWriting out " << ct_pts << " points and " << ct_vxls << " voxels as unstructured VTK file " << vtkfile << ".\n";
#endif

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

#ifdef PRINT
    printf("\tWriteout done!\n");
#endif

    os.close();
}

//=============================================================================
// @ARG: vtk filename of the adative grid to write
//=============================================================================
void RayIntersector::write_adaptive_vtk_file_only_leafs() {

    std::ofstream os;

    /* write grid to file */
    os.open(vtkfile, std::ofstream::out | std::ofstream::trunc);
    if (!os.is_open()) {
        std::cerr << "ERROR >> couldn't write file " << vtkfile << std::endl;
        exit(0);
    }

    int ct_pts = 0, ct_vxls = 0, ct = 0; // counters
    llist2 *tmp;         // next stack element
    double v[8][3];     // voxel's node coordinates

    // count points and voxels
    if (out != nullptr && out->next != nullptr) {
        tmp = out->next;
        while (tmp != nullptr) {
            if (isLeaf(tmp->data))
                ++ct_vxls;
            tmp = tmp->next;
        }
    }
    ct_pts = 8 * ct_vxls;

#ifdef PRINT
    std::cout << "\tWriting out " << ct_pts << " points and " << ct_vxls << " voxels as unstructured VTK file " << vtkfile << ".\n";
#endif

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
            if (isLeaf(tmp->data)) {
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
            }
            tmp = tmp->next;
        }
    }

    // write cells
    os << "CELLS " << ct_vxls << " " << 9 * ct_vxls << "\n";
    if (out != nullptr && out->next != nullptr) {
        tmp = out->next;
        ct = 0;
        while (tmp != nullptr) {
            if (isLeaf(tmp->data)) {
                os << 8 << " ";
                for (int i = 0; i < 8; ++i) os << " " << ct++;
                os << "\n";
            }

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
            if (isLeaf(tmp->data))
                os << (int) tmp->data->color << "\n";

            tmp = tmp->next;
        }
    }

#ifdef PRINT
    printf("\tWriteout done!\n");
#endif

    os.close();
}

//=============================================================================
// @ARG: amount of black voxels
// @ARG: amount of white voxels
// @ARG: amount of grey voxels
//=============================================================================
void RayIntersector::count(long *white, long *grey) {
    *white = *grey = 0;

    llist2 *it(out);

    while (it != nullptr) {

        if (it->data != nullptr) {
            switch (it->data->color) {
                case (OCTREE_COLOR_WHITE) :
                    (*white)++;
                    break;
                case (OCTREE_COLOR_GREY) :
                    (*grey)++;
                    break;
                default:
                    std::cerr << "ERROR: Undefined color found while executing RayIntersector::count(...) !" << std::endl;
                    break;
            }
        }

        it = it->next;
    }

}

bool RayIntersector::isLeaf(Voxel2 *voxel) { return voxel->children[0] == nullptr; }

std::multimap<double, unsigned int> RayIntersector::perform(const IntersectionRay &ray) {

    std::multimap<double, unsigned int> hits;

    std::list<Voxel2 *> voxels = get_intersected_voxels(ray);
    std::set<unsigned int> triangles_to_test = get_potential_intersected_triangles(voxels);
    if (triangles_to_test.empty()) return hits;
    get_intersected_triangles(ray, triangles_to_test, hits);

    //for (const auto &hit: hits)
    //    std::cerr << hit.first << " " << hit.second << std::endl;

    return hits;
}

std::list<Voxel2 *> RayIntersector::get_intersected_voxels(const IntersectionRay &ray) {

    std::list<Voxel2 *> voxels;
    std::stack<Voxel2 *> s;

    for (int i = 0; i < 8; i++)
        if (root_voxel->children[i] != nullptr)
            s.push(root_voxel->children[i]);

    while (!s.empty()) {

        Voxel2 *voxel = s.top();
        s.pop();

        if (!ray.intersects_voxel(voxel)) continue;

        if (isLeaf(voxel)) {
            if (voxel->nface > 0)
                voxels.push_back(voxel);
        } else
            for (int i = 0; i < 8; i++)
                s.push(voxel->children[i]);
    }

    return voxels;
}

std::set<unsigned int> RayIntersector::get_potential_intersected_triangles(const std::list<Voxel2 *> &voxels) {

    // detect triangles to test for intersection
    std::set<unsigned int> triangles_to_test;

    for (const auto &voxel: voxels)
        for (int i = 0; i < voxel->nface; ++i)
            triangles_to_test.insert(voxel->findex[i]);

    return triangles_to_test;
}

void RayIntersector::get_intersected_triangles(const IntersectionRay &ray, const std::set<unsigned int> &triangles_to_test, std::multimap<double, unsigned int> &hits) {

    std::set<unsigned int> ids; // ids (attribute indicating cFace* of triangle) of triangles that were already insert in hit map

    for (const auto &i: triangles_to_test) {

        const IsctFace &tri = Mod->flist[i];

        if (ids.find(tri.id) != ids.end()) // if triangle belongs to same cface as a triangle already in hit list, no need to check intersection
            continue;

        double dist;

        if (ray.intersects_triangle(Mod->vlist[tri.v1], Mod->vlist[tri.v2], Mod->vlist[tri.v3], dist)) {
            hits.insert(std::pair<double, unsigned int>(dist * Mod->bbox_norm, tri.id)); // rescale dist
            ids.insert(tri.id);
        }
    }
}

bool RayIntersector::isZero(const double v) { return v >= -1.0e-9 && v <= 1.0e-9; }

int RayIntersector::triangle_aabb_moeller(const Voxel2 *voxel, const int idx) {

    IsctFace f = Mod->flist[idx];
    IsctVertex v1 = Mod->vlist[f.v1];
    IsctVertex v2 = Mod->vlist[f.v2];
    IsctVertex v3 = Mod->vlist[f.v3];

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