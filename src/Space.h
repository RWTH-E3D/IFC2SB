// Copyright 2022 Eric Fichter
#ifndef SPACE_H
#define SPACE_H

#include "headers.h"

class cFace;

class Space {

public:
    Space(unsigned int _id, std::set<cFace *> &_faces, bool update_faces_after_sewing);

    bool operator==(const Space &S) const { return id == S.id; }

    unsigned int id;
    bool is_facade;

    //! Returns volume of space calculated by its shell.
    double Volume() const;

    //! Calculates some ifc space quantities based on space's SBs.
    void CalcQuantities();

    //! Enrich space with space boundary.
    void AddSecondLvlFace(cFace *cface);

    //! Returns second level space boundaries.
    std::set<cFace *> SecondLevel() const;

    //! Returns first level space boundaries.
    std::set<cFace *> FirstLevel() const;

    //! Returns space shape.
    TopoDS_Shell Shell() const;

    //! Clear first level space boundaries.
    void ClearFirstLevel();

    //! Print attributes.
    void Log() const;

    //! Returns IfcProducts present in space.
    std::set<IfcUtil::IfcBaseClass *> ProductsInSpace();

private:
    TopoDS_Shell shell;
    std::set<cFace *> SecondLvl;
    std::set<cFace *> FirstLvl;
    double aabb_xmin;
    double aabb_xmax;
    double aabb_ymin;
    double aabb_ymax;
    double aabb_zmin;
    double aabb_zmax;
    double volume;
    double floor_elevation;
    bool wasFlipped; // indicates if faces had to be flipped to get an outward oriented shell
    std::unique_ptr<BRepClass3d_SolidClassifier> classifier;
    //std::unique_ptr<IntCurvesFace_ShapeIntersector> intersector;
    std::set<Space *> nbs;

    // IFC related
    IfcUtil::IfcBaseClass *storey;
    std::set<std::string> old_space_guids;
    std::map<std::string, std::string> old_space_info;
    std::set<IfcUtil::IfcBaseClass *> IfcRelContainedInSpatialStructures;
    std::set<IfcUtil::IfcBaseClass *> IfcRelDefinesByProperties;
    std::set<IfcUtil::IfcBaseClass *> IfcRelAssigns;
    std::set<IfcUtil::IfcBaseClass *> IfcRelAssociates;
    IfcUtil::IfcBaseClass * IfcRelDefinesByType;
    double NetVolume; // Net volume enclosed by the space, excluding the volume of construction elements inside the space.
    double NetFloorArea; // Sum of all usable floor areas covered by the space. It excludes the area covered by elements inside the space (columns, inner walls, built-in's etc.), slab openings, or other protruding elements.
    double NetWallArea; // Sum of all wall (and other vertically bounding elements, like columns) areas bounded by the space. It excludes the area covered by elements inside the wall area (doors, windows, other openings, etc.).
    double NetCeilingArea; // Sum of all ceiling areas of the space. It excludes the area covered by elements inside the space (columns, inner walls, etc.). The ceiling area is the real (and not the projected) area (e.g. in case of sloped ceilings).

    TopoDS_Shell build_shell_from_cFaces(bool update_faces_after_sewing);

    bool is_point_in_aabb(gp_Pnt P) const;

    bool is_point_in_aabb(gp_Pnt P, double tol) const;

    bool is_point_in_shell(gp_Pnt P);

    gp_Pnt random_point_on_cface_shell();

    void create_solid_classifier_from_shell();

    // void create_intersector_from_cfaces(double tol);

    void get_second_level();

    void calc_floor_elevation();

    friend class Kernel;
};

#endif //SPACE_H
