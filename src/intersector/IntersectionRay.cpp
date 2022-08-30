// Copyright 2022 Eric Fichter
#include "IntersectionRay.h"

IntersectionRay::IntersectionRay(std::array<double, 3> P0, std::array<double, 3> v0, double l0) : P(P0), v(v0), l(l0) {

    tol = 1.0e-9;
    normalize_vector();
    calc_v_helper();
    calc_end_point();
}

void IntersectionRay::convert_lengths_to_octree_grid(double norm, double dx, double dy, double dz) {
    P = {(P[0] - dx) / norm, (P[1] - dy) / norm, (P[2] - dz) / norm};
    Q = {(Q[0] - dx) / norm, (Q[1] - dy) / norm, (Q[2] - dz) / norm};
    l = l / norm;
}

void IntersectionRay::normalize_vector() {
    double m = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    v = {v[0] / m, v[1] / m, v[2] / m};
}

void IntersectionRay::calc_end_point() {

    // IEEE 754 floating-point properties implicitly handles the edge case where a component of the direction is zero.
    // t values in intersection calculation will be inf and -inf, causing the test to "fail" (return false).

    double x = P[0] + l * v[0];
    double y = P[1] + l * v[1];
    double z = P[2] + l * v[2];
    Q = {x, y, z};
}

void IntersectionRay::calc_v_helper() { v_inv = {1.0f / v[0], 1.0f / v[1], 1.0f / v[2]}; }

double IntersectionRay::distance_from_origin(std::array<double, 3> p) const {
    double a = P[0] - p[0];
    double b = P[1] - p[1];
    double c = P[2] - p[2];
    return sqrt(a * a + b * b + c * c);
}

bool IntersectionRay::intersects_voxel(const Voxel2 *voxel) const {

    // check if start or end point are in node
    if (point_in_voxel(P, voxel) || point_in_voxel(Q, voxel))
        return true;

    return intersects_aabb(voxel);
}

bool IntersectionRay::point_in_voxel(const std::array<double, 3> &p, const Voxel2 *voxel) {
    if (p[0] > voxel->maxx) return false;
    if (p[1] > voxel->maxy) return false;
    if (p[2] > voxel->maxz) return false;
    if (p[0] < voxel->minx) return false;
    if (p[1] < voxel->miny) return false;
    if (p[2] < voxel->minz) return false;
    return true;
}

bool IntersectionRay::intersects_aabb(const Voxel2 *voxel) const {

    // https://gamedev.stackexchange.com/questions/18436/most-efficient-aabb-vs-ray-collision-algorithms
    // https://tavianator.com/2011/ray_box.html

    // P ---v--> tmin (BBOXMIN) ---> tmax (BBOXMAX)
    // tmin (BBOXMIN) ---> tmax (BBOXMAX) | P ---v-->       (tmax < 0)

    // calculate distances from origin to voxel vertices in specific direction (x,y,z) along ray direction
    // (length of ray until intersection)
    double t1 = (voxel->minx - P[0]) * v_inv[0]; // to xmin
    double t2 = (voxel->maxx - P[0]) * v_inv[0]; // to xmax
    double t3 = (voxel->miny - P[1]) * v_inv[1]; // to ymin
    double t4 = (voxel->maxy - P[1]) * v_inv[1]; // to ymax
    double t5 = (voxel->minz - P[2]) * v_inv[2]; // to zmin
    double t6 = (voxel->maxz - P[2]) * v_inv[2]; // to zmax

    // tmin and tmax are the distances from origin to intersection with voxel
    double tmin = std::max(std::max(std::min(t1, t2), std::min(t3, t4)), std::min(t5, t6));
    double tmax = std::min(std::min(std::max(t1, t2), std::max(t3, t4)), std::max(t5, t6));

    if (tmax < 0)  // ray is intersecting AABB, but the whole AABB is behind us. t would be tmax
        return false;

    if (tmin > tmax) // ray doesn't intersect AABB. t would be tmax
        return false;

    // t is now tmin (f origin is in box, tmin is negative)
    if (tmin > l) // ray intersects after end point of ray
        return false;

    return true;
}

bool IntersectionRay::intersects_triangle(const IsctVertex &v0, const IsctVertex &v1, const IsctVertex &v2, double &t) const {

    // https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
    // Möller–Trumbore intersection algorithm

    std::array<double, 3> edge1 = subtract(v1, v0);
    std::array<double, 3> edge2 = subtract(v2, v0);
    std::array<double, 3> h = cross(v, edge2);

    std::array<double, 3> q{};
    double a, f, u, k;

    a = dot(edge1, h);

    if (a > -tol && a < tol) // ray is parallel to triangle
        return false;

    f = 1.0 / a;
    std::array<double, 3> s = subtract(P, v0);
    u = f * dot(s, h);

    if (u < 0.0 || u > 1.0)
        return false;

    q = cross(s, edge1);
    k = f * dot(v, q);

    if (k < 0.0 || u + k > 1.0)
        return false;

    // At this stage we can compute distance t from origin to intersection point.
    t = f * dot(edge2, q);

    if (t > l || t < -tol)
        return false;

    // intersection point on the line.
    //   P_isct = add(P, multiply(v, t));

    return true;
}

std::array<double, 3> IntersectionRay::subtract(const IsctVertex &a, const IsctVertex &b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }

std::array<double, 3> IntersectionRay::subtract(const std::array<double, 3> &a, const IsctVertex &b) { return {a[0] - b.x, a[1] - b.y, a[2] - b.z}; }

std::array<double, 3> IntersectionRay::add(const std::array<double, 3> &a, const std::array<double, 3> &b) { return {a[0] + b[0], a[1] + b[1], a[2] + b[2]}; }

std::array<double, 3> IntersectionRay::multiply(const std::array<double, 3> &a, double d) { return {a[0] * d, a[1] * d, a[2] * d}; }

std::array<double, 3> IntersectionRay::cross(const std::array<double, 3> &a, const std::array<double, 3> &b) {

    std::array<double, 3> c{};
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
    return c;

}

double IntersectionRay::dot(const std::array<double, 3> &a, const std::array<double, 3> &b) { return std::inner_product(std::begin(a), std::end(a), std::begin(b), 0.0); }

std::array<double, 3> IntersectionRay::get_start_point() { return P; }

std::array<double, 3> IntersectionRay::get_end_point() { return Q; }