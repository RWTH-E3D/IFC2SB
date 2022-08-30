// Copyright 2022 Eric Fichter
#ifndef INTERSECTIONRAY_H
#define INTERSECTIONRAY_H

#include "intersector_includes.h"
#include "RayIntersector.h"

//! Describes a ray as a finite line segment.
class IntersectionRay {

public:
    IntersectionRay(std::array<double, 3> P0, std::array<double, 3> v0, double l0);

    double distance_from_origin(std::array<double, 3> Q) const;

    //! Returns true if ray intersects a voxel.
    bool intersects_voxel(const Voxel2 *voxel) const;

    //! Rescales ray to fit normalized octree model.
    void convert_lengths_to_octree_grid(double norm, double dx, double dy, double dz);

    //! Returns true if ray intersects a triangle defined by its vertices.
    //! If true, sets distance t from ray origin to intersection point.
    bool intersects_triangle(const IsctVertex &v1, const IsctVertex &v2, const IsctVertex &v3, double &t) const;

    //! Get start point.
    std::array<double, 3> get_start_point();

    //! Get end point.
    std::array<double, 3> get_end_point();

private:
    std::array<double, 3> P; // origin, start point
    std::array<double, 3> Q; // end point

    std::array<double, 3> v; // direction
    std::array<double, 3> v_inv; // 1/v for performance reasons

    double l; // length

    //! tolerance for ray - triangle intersection.
    double tol;

    void normalize_vector();

    void calc_end_point();

    void calc_v_helper();

    static bool point_in_voxel(const std::array<double, 3> &p, const Voxel2 *voxel) ;

    bool intersects_aabb(const Voxel2 *voxel) const;

    static std::array<double, 3> subtract(const IsctVertex &a, const IsctVertex &b) ;

    static std::array<double, 3> subtract(const std::array<double, 3> &a, const IsctVertex &b) ;

    static std::array<double, 3> add(const std::array<double, 3> &a, const std::array<double, 3> &b) ;

    static std::array<double, 3> multiply(const std::array<double, 3> &a, double d) ;

    static std::array<double, 3> cross(const std::array<double, 3> &a, const std::array<double, 3> &b) ;

    static double dot(const std::array<double, 3> &a, const std::array<double, 3> &b) ;
};

#endif //INTERSECTIONRAY_H