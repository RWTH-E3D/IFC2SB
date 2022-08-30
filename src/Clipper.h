// Copyright 2022 Eric Fichter
#ifndef CLIPPER_H
#define CLIPPER_H

#include "headers.h"

class Clipper {

public:
    //! Constructor.
    //! tF1 ... face to be clipped
    //! tF2 ... projected and clipping face
    //! tn1 ... normal of tF1
    //! tn2 ... normal of tF2
    Clipper(TopoDS_Face tF1, gp_Dir tn1, bool _dimension);

    bool success;
    bool same;
    TopoDS_ListOfShape Result;   // Resulting faces

    typedef std::vector<gp_Pnt2d> chain2D;
    typedef std::vector<chain2D> chains2D;

    chains2D Result2DWires; // no distinguishing is made between holes and normal wires, because resulting wires will always be one F1
    chain2D F1_wire2D;
    chains2D F1_holes2D;

    void clip(const TopoDS_Face &F2, gp_Dir n2);

    static void calculate_algebraic_values(gp_Pnt P, gp_Dir n, double &temp_a, unsigned int &temp_proj_axis);

private:

    const bool dimension;     // dimension of result (true -> 3d, false 2d)
    const TopoDS_Face F1;     // F1 ... face to be clipped
    const gp_Dir n1;          // n1 ... normal of face to be clipped
    const gp_Vec d1;          // Vector of n1, plane normal
    double a;                 // plane equation constant, n[0] * x + n[1] * y + n[2] * z = a
    unsigned int proj_axis;
    const gp_Pnt P_plane;     // Random point on F1
    const unsigned long long int f; // scaling_factor, resulting number can not be higher than int max 2147483647
    double f_inv;             // inverse scaling factor
    double tol_area;          // if area of intersection is same, skip result
    double tol_length;        // removal of colinear points and external edges

    struct Point2D {
        double X;
        double Y;
    };

    typedef std::vector<Point2D> PointChain2D;

    PointChain2D P1_2D;
    std::vector<PointChain2D> P1_holes_2D;

    void calculate_algebraic_values();

    void clipping(const PointChain2D &PolyPoints1, const PointChain2D &PolyPoints2, const std::vector<PointChain2D> &HolePoints1, const std::vector<PointChain2D> &HolePoints2, std::list<PointChain2D> &common_2D_polygons, std::list<PointChain2D> &common_2D_holes);

    static std::vector<gp_Pnt> vertex_to_point_list(const TopoDS_ListOfShape &V);

    static TopoDS_Wire polygon_to_wire(const std::vector<gp_Pnt> &polygon_3D);

    static TopoDS_Face polygon_to_face(const std::vector<gp_Pnt> &polygon_3D, const TopoDS_ListOfShape &Holes);

    static std::vector<double> aabb_polygon(const PointChain2D &poly);

    static bool do_boxes_intersect(const std::vector<double> &A, const std::vector<double> &B);

    static std::list<TopoDS_ListOfShape> ordered_vertices_of_hole_wires(const TopoDS_Face &F, const TopoDS_Wire &outerWire);

    static std::list<std::vector<gp_Pnt>> vertex_to_point_list_holes(const std::list<TopoDS_ListOfShape> &V);

    PointChain2D wire_to_projected_points_1(const TopoDS_Wire &W) const;

    PointChain2D wire_to_projected_points_2(const TopoDS_Wire &W) const;

    std::vector<PointChain2D> face_to_projected_holes_1(const TopoDS_Face &F, const TopoDS_Wire &W) const;

    std::vector<PointChain2D> face_to_projected_holes_2(const TopoDS_Face &F, const TopoDS_Wire &W) const;

    gp_Pnt project_point_onto_plane(gp_Vec P) const;

    std::vector<gp_Pnt> projects_points_onto_face_plane(const std::vector<gp_Pnt> &R) const;

    Point2D project(gp_Pnt P) const;

    PointChain2D project_points_to_2D(const std::vector<gp_Pnt> &R) const;

    std::list<std::vector<gp_Pnt>> projects_points_onto_face_plane_holes(const std::list<std::vector<gp_Pnt>> &R) const;

    std::vector<PointChain2D> project_points_to_2D_holes(const std::list<std::vector<gp_Pnt>> &R) const;

    gp_Pnt project_inv(Point2D P_2D) const;

    std::vector<std::vector<gp_Pnt>> reproject_2D_to_3D(const std::list<PointChain2D> &R) const;
};

#endif //CLIPPER_H
