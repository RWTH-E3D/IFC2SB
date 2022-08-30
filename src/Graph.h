// Copyright 2022 Eric Fichter
#ifndef GRAPH_H
#define GRAPH_H

#include "headers.h"

class Product;

class oFace;

class cFace;

class Space;

class Kernel;

//! Structures the graph-based Space Boundary generation method.
//! Kernel functions are called.
//! Stores ifc model and all other important topological and geometric objects.
class Graph {

public:
    //! Stores ifc model and all other important topological and geometric objects.
    //! _input: Input path and name of ifc file.
    //! _input: Output path and name of result file.
    //! _num_threads: Number of threads used for parallelization.
    //! _stl: If true, output will be an stl with first level space boundaries. Else IFC.
    //! _space_split: To be deprecated. If true, use old IfcSpaces to create virtual space boundaries.
    //! _max_transmission_length: Heat transmission length to define, if SB is of type 2a or 2b.
    //! _first_level_only: If true, only IfcRelSpaceBoundary1stLevel objects are created..
    Graph(std::string _input,
          std::string _output,
          unsigned int _num_threads,
          bool _stl,
          bool _space_split,
          double _max_transmission_length,
          bool _first_level_only,
          bool _remove_holes_by_face_split,
          bool _decompose_concave_polygons,
          bool _simplify_fenestration_faces,
          bool _calculate_shadings,
          bool _use_ifcopeningelelements_for_virtual_boundaries,
          bool _use_spaces_for_virtual_boundaries
    );

    //! Starts space boundary generation
    bool run();

private:

    //! IfcOpenShell IFC model. To be parsed and enriched.
    std::unique_ptr<IfcParse::IfcFile> model;

    //! Number of threads used for parallelization.
    const unsigned int num_threads;

    //! Input path and name of ifc file.
    const std::string input;

    //! Output path and name of result file.
    const std::string output;

    //! If true, output will be an stl with first level space boundaries. Else IFC.
    bool stl;

    //! To be deprecated. If true, use old IfcSpaces to create virtual space boundaries.
    bool space_split;

    //! If true, only IfcRelSpaceBoundary1stLevel objects are created.
    bool first_level_only;

    //! Product instance containing e.g. geometry of IFC products.
    std::list<Product> products;

    //! IFC classes to be included in the processing.
    std::set<std::string> include_entities;

    //! IFC classes to be included in the stl export, in addition to include_entities.
    std::set<std::string> cfd_entities;

    //! IFC classes that should be considered for face offsetting.
    std::set<std::string> offset_entities;

    //! Faces coming from IFC file or derived from them by offsetting.
    std::list<oFace> ifc_faces;

    //! First level space boundaries and other generated faces (that are deleted at some point).
    std::list<cFace> faces_1st_level;

    //! Faces not bordering an air volume (internal faces) and (redundant) faces coplanar to other faces.
    std::list<cFace> faces_non_sb;

    //! Second level space boundaries.
    std::list<cFace> faces_2nd_level;

    //! Shading faces.
    std::unordered_map<std::string, std::list<TopoDS_Face>> shadings;

    //! "Trash can" for faces.
    std::list<cFace> faces_trash;

    //! Spaces (closed volumes surrounded by faces).
    std::list<Space> spaces;

    //! Minimum coordinates of AABB of all products.
    gp_XYZ bounds_min;

    //! Maximum coordinates of AABB of all products.
    gp_XYZ bounds_max;

    //! IfcSpaces from the original ifc file and their shapes.
    ifcspaceInfoList ifcspaces;

    //! Heat transmission length to define, if SB is of type 2a or 2b.
    double max_transmission_length;

    //! Offsetting length of faces.
    double face_offset_length;

    //! Tolerance for Boolean Operations (added to standard tolerance of shape).
    double fuzzy_tol;

    //! Integrates windows and doors into wall faces.
    bool integrate_openings_into_walls;

    //! Split second level space boundaries with holes into faces without holes.
    bool remove_holes_by_face_split;

    //! Decompose concave second level space boundaries to convexity.
    bool decompose_concave_polygons;

    //! Simplify openings with more than four edges by using a recantgular with same surface area instead.
    bool simplify_fenestration_faces;

    //! Detect shading faces.
    bool calculate_shadings;

    //! Virtual boundaries are created on the faces of "virtual" IfcOpeningElements.
    bool use_ifcopening_elements_for_virtual_boundaries;

    //! Virtual boundaries are created by contact faces of ifcspaces.
    bool use_spaces_for_virtual_boundaries;

    //! Rank system for ifc classes applied in cases of coplanar faces. Lower rank is preferred as choice for relatedbuildingelement.
    std::unordered_map<std::string, unsigned int> ranks;

    unsigned int space_id_counter;
    unsigned int face_id_counter;
    IfcGeom::IteratorSettings settings;
    std::set<IfcUtil::IfcBaseEntity *> virtual_products;

    bool prepare_products(Kernel &K, std::unordered_map<IfcUtil::IfcBaseClass *, TopoDS_Shape> &opening_shapes);

    bool prepare_faces(Kernel &K, std::unordered_map<IfcUtil::IfcBaseClass *, TopoDS_Shape> &opening_shapes);

    bool calc_faces_first_level(Kernel &K);

    void calc_faces_first_level_normals_known(const TopoDS_Shape& fuse);

    void calc_faces_first_level_normals_unknown(const TopoDS_Shape& fuse, Kernel &K);

    bool calc_faces_second_level(Kernel &K);

    bool process_spaces(Kernel &K);

    bool process_model(Kernel &K);
};

#endif //GRAPH_H