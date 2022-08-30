// Copyright 2022 Eric Fichter
#ifndef KERNEL_H
#define KERNEL_H

#include "headers.h"

class Product;

class oFace;

class cFace;

class Kernel {

public:

    Kernel(unsigned int _num_threads);

    static unsigned int hash(const TopoDS_Shape &s);

    static gp_XY surface_uv_coordinates_at_point(const Handle(Geom_Surface) &surface, gp_Pnt P);

    static gp_Dir face_normal_at_point(const TopoDS_Face &face, gp_Pnt P, bool &success);

    static gp_Dir face_normal(const TopoDS_Face &face);

    static gp_Pnt face_center(const TopoDS_Face &face);

    static void move(TopoDS_Shape &S, gp_Vec v);

    static TopoDS_Shape moved(TopoDS_Shape &S, gp_Vec v);

    static void rotate(TopoDS_Shape &S, gp_Pnt P, gp_Dir n, double angle);

    static Bnd_Box aabb(const TopoDS_Shape &shape, double gap);

    static Bnd_OBB obb(const TopoDS_Shape &shape);

    static double volume(const TopoDS_Shape &shape);

    static double length(const TopoDS_Edge &e);

    static double area(const TopoDS_Shape &shape);

    static TopoDS_Shape geom_object_to_shape(IfcGeom::Element<real_t> *geom_object);

    static TopoDS_Shape replace_subshapes_in_shape(const std::list<std::pair<TopoDS_Shape, TopoDS_Shape>> &old_and_new_subshapes, TopoDS_Shape &shape);

    static TopoDS_Compound compound_from_shape_list(const TopoDS_ListOfShape &L);

    static bool are_points_colinear(gp_Pnt p1, gp_Pnt p2, gp_Pnt p3, double tol);

    static bool are_points_colinear(gp_Pnt2d p1, gp_Pnt2d p2, gp_Pnt2d p3, double tol);

    static bool face_is_polygon(const TopoDS_Face &face);

    static bool face_is_planar(const TopoDS_Face &face);

    static bool polygon_is_convex(const std::vector<gp_Pnt2d> &l);

    static TopoDS_Shape polygonize_shape_2a_curvature(const TopoDS_Shape &shp, double linear_deflection = 0.9, bool isRelative = true, double angular_deflection = 0.9, bool isInParallel = false);

    static TopoDS_Shape polygonize_shape(const TopoDS_Shape &shp, double linear_deflection = 0.9, bool isRelative = true, double angular_deflection = 0.9, bool isInParallel = false);

    static TopoDS_Shape shells_to_solids(const TopoDS_Shape &shp, bool complement = true);

    static TopoDS_Shape shape_copy(const TopoDS_Shape &S);

    static gp_Pnt shape_center(const TopoDS_Shape &S);

    static double round_double_to_n_decimal_places(double d, unsigned int n);

    static bool is_edge_line(const TopoDS_Edge &E);

    static double round_double_one_digit(double d);

    static double round_double_two_digits(double d);

    static double round_double_three_digits(double d);

    static gp_Pnt point_on_face(TopoDS_Face &inputFace, const gp_Dir &n, std::string s = "");

    static gp_Pnt move_point_along_scaled_unit_vector(gp_Pnt c, gp_Dir n, double s);

    static bool is_point_in_face(const TopoDS_Face &face, gp_Pnt P, double tol);

    static gp_Pnt mid_uv_point_on_surface(const TopoDS_Face &face);

    static gp_Pnt random_point_on_face_using_triangulation(const TopoDS_Shape &shp);

    static gp_Pnt intersection_face_line(const TopoDS_Face &face, gp_Lin line, double tol);

    static bool do_aabbs_intersect(double A_Xmin, double A_Xmax, double A_Ymin, double A_Ymax, double B_Xmin, double B_Xmax, double B_Ymin, double B_Ymax);

    static void check_product_tolerances(std::list<Product> &products, double critical);

    static void subtract_openings_from_products(std::list<Product> &products, std::unordered_map<IfcUtil::IfcBaseClass *, TopoDS_Shape> &opening_shapes);

    static void heal_products(std::list<Product *> &products_to_heal);

    static TopoDS_Shape best_fitting_bbox(const TopoDS_ListOfShape &L);

    static TopoDS_Shape best_fitting_bbox(const TopoDS_Shape &S);

    bool geom_context_info(std::unique_ptr<IfcParse::IfcFile> &model) const;

    void simplify_products(std::list<Product> &products) const;

    void cut_products(std::list<Product> &products, const std::set<std::string> &C, const double fuzzy_tol) const;

    void fuse_products(std::list<Product> &products, const double tol) const;

    static void polygonize_products(std::list<Product> &products);

    static bool read_ifc_file(const std::string &ifc_path, std::unique_ptr<IfcParse::IfcFile> &model);

    bool check_ifc_schema(std::unique_ptr<IfcParse::IfcFile> &model);

    bool get_length_conversion_factor(std::unique_ptr<IfcParse::IfcFile> &model);

    void find_relevant_products(std::unique_ptr<IfcParse::IfcFile> &model, const std::set<std::string> &include_entities, std::set<std::string> &non_void_filling_products, bool integrate_openings_into_walls, bool use_ifcopeningelelements_for_virtual_boundaries);

    bool generate_shapes_from_ifc_guids(std::unique_ptr<IfcParse::IfcFile> &model, IfcGeom::IteratorSettings settings, std::list<Product> &products, const std::set<std::string> &include_guids, gp_XYZ &bounds_min, gp_XYZ &bounds_max) const;

    void add_decomposed_entities_to_products(std::unique_ptr<IfcParse::IfcFile> &model, IfcGeom::IteratorSettings settings, std::list<Product> &products, const std::set<std::string> &include_entities);

    static std::list<Product *> check_shells(std::list<Product> &products);

    static std::list<Product *> check_shells(std::list<Product *> &products);

    std::set<std::string> find_openings_in_products(std::list<Product> &products);

    std::unordered_map<IfcUtil::IfcBaseClass *, TopoDS_Shape> generate_opening_shapes(std::unique_ptr<IfcParse::IfcFile> &model, IfcGeom::IteratorSettings settings, const std::set<std::string> &guids) const;

    static bool identify_non_planar_products(std::list<Product> &products);

    static void remove_products(std::list<Product> &products);

    static void check_products(std::list<Product> &products, double l_warn, double l_crit);

    static std::vector<gp_Pnt> vertex_to_point_list(const TopoDS_ListOfShape &V);

    static std::vector<gp_Pnt2d> project_wire_to_2D(long a, const TopoDS_Wire &W);

    static double determinant_three_points(gp_Pnt2d V1, gp_Pnt2d V2, gp_Pnt2d V3);

    template<typename SetType>
    static bool is_element_in_set(std::set<SetType> const &s, SetType const &element);

    static TopoDS_Face remove_seam_edges(const TopoDS_Face &f);

    static TopoDS_Shape mold(const std::list<Space> &spaces);

private:

    const unsigned int num_threads;
    unsigned int n_digits_round_double_ifc_write;
    unsigned int round_double_ifc_write;
    double round_ifc_write_inv;
    IFC_SCHEMA ifcSchema;
    double conv_fctr; // conversion factor from metre to ifc length unit <prefix>-metre
    std::map<TopAbs_ShapeEnum, std::string> shapeEnum_to_string;

    struct opening_info {
        TopoDS_Face face;
        Product *product;
        unsigned short int shell_id;
        unsigned int parent_id;
        IfcUtil::IfcBaseClass *opening_product;
        face_normal_status normal_status;

        opening_info(TopoDS_Face face, Product *product, unsigned short int shell_id, unsigned int parent_id, IfcUtil::IfcBaseClass *opening_product, face_normal_status normal_status) :
                face(face), product(product), shell_id(shell_id), parent_id(parent_id), opening_product(opening_product), normal_status(normal_status) {}
    };

    //! Structure holding some face data.
    struct site_face {
        TopoDS_Face face;
        Bnd_Box aabb;
        gp_Dir n;

        site_face(TopoDS_Face face, Bnd_Box aabb, gp_Dir n) : face(face), aabb(aabb), n(n) {}
    };

    struct segmentation_prism {
        TopoDS_Shape S;
        TopoDS_Face base;
        gp_Dir n;
        gp_Pnt c;

        segmentation_prism(TopoDS_Shape S, TopoDS_Face base, gp_Dir n, gp_Pnt c) : S(std::move(S)), base(std::move(base)), n(n), c(c) {}
    };

    typedef rtree_lib::RTree<cFace *, double, 3, double> cface_tree3D;

    typedef std::unordered_map<IfcUtil::IfcBaseClass *, std::pair<TopoDS_Shape, IfcUtil::IfcBaseClass * >> OWP; // key: opening product, value: a pair of opening shape and window/door product

    bool generate_shapes_from_ifc_classes(std::unique_ptr<IfcParse::IfcFile> &model, IfcGeom::IteratorSettings settings, std::list<Product> &products, const std::set<std::string> &include_entities, gp_XYZ &bounds_min, gp_XYZ &bounds_max) const;

    static void collect_original_faces(std::list<Product> &products, std::list<oFace> &orig_faces);

    static void collect_original_faces(Product &product, std::list<oFace> &orig_faces);

    void get_virtual_faces(std::unique_ptr<IfcParse::IfcFile> &model, ifcspaceInfoList &ifcspaces, std::list<oFace> &ifc_faces, std::list<Product> &products, std::set<IfcUtil::IfcBaseEntity *> &virtual_products);

    template<typename Schema>
    void get_virtual_faces_worker(ifcspaceInfoList &spaces, std::list<oFace> &orig_faces, std::list<Product> &products, std::set<IfcUtil::IfcBaseEntity *> &virtual_products);

    static void add_offset_faces(std::list<oFace> &orig_faces, double offset, const std::set<std::string> &classes_for_face_extension);

    static bool offset_faces_outer_wire_cut(const TopoDS_Face &face, TopoDS_Face &offset_face, double offset);

    static bool offset_wire_cut(const TopoDS_Wire &wire, TopoDS_Wire &offsetWire, double offset, GeomAbs_JoinType joinType);

    static bool fuse_original_faces(TopoDS_Shape &fuse, std::list<oFace> &orig_faces, std::list<cFace> &cFaces, double fuzzy_tol, unsigned int &fid);

    static bool fuse_original_faces(TopoDS_Shape &fuse, std::set<oFace *> &orig_faces, std::list<cFace> &cFaces, double fuzzy_tol);

    static void correct_face_normals(std::list<cFace> &cFaces);

    static void update_half_edges(std::list<cFace> &cFaces);

    static void identify_enclosed_faces(const TopoDS_Shape &fuse, std::list<cFace> &cFaces);

    static void remove_trash(std::list<cFace> &cFaces, std::list<cFace> &cFaces_trash);

    static void identify_duplicate_faces(std::list<cFace> &cFaces, std::unordered_map<std::string, unsigned int> ranks);

    static bool compare_normals(const std::list<cFace *> &cFaces);

    static std::multimap<int, cFace *> sort_mat_mat_faces(const std::list<cFace *> &mat_faces, std::unordered_map<std::string, unsigned int> &ranks);

    static void check_duplicate_faces(std::list<cFace> &cFaces);

    static void update_face_adjacencies(std::list<cFace> &cFaces, const TopoDS_Shape &fuse);

    static std::unordered_map<unsigned int, cFace *> setup_face_cFace_map(std::list<cFace> &cFaces);

    static void remove_adjacency_by_orientation(std::list<cFace> &cFaces);

    static void identify_hanging_faces(std::list<cFace> &cFaces);

    static bool identify_hanging_faces_while(std::list<cFace> &cFaces);

    static void remove_trash_and_face_adjacency(std::list<cFace> &cFaces);

    static void find_spaces(std::list<cFace> &cFaces, std::list<Space> &spaces, unsigned int &space_id);

    void find_spaces_normals_unknown(std::list<cFace> &cFaces, std::list<Space> &spaces, unsigned int &space_id, unsigned int &fid);

    static void identify_sb_types(std::list<cFace> &cFaces, std::list<Space> &spaces, gp_XYZ bounds_min, double depth_2b);

    static void identify_facade_space(std::list<Space> &spaces);

    static void remove_seam_edges_from_faces(std::list<cFace> &cFaces);

    static void unify_cFaces(std::list<cFace> &cFaces, bool close_filled_wall_holes);

    static void remove_trash_from_spaces(std::list<cFace> &cFaces, std::list<Space> &spaces);

    static void triangulate_cfaces(std::list<cFace> &cFaces);

    bool write_spaces_to_stl(std::list<cFace> &cFaces, std::list<Space> &spaces, const std::string &input_ifc_filepath, const std::string &output_ifc_file, const std::string &path, const std::string &mode, bool ascii, bool ignore_facade_space,
                             const std::map<std::string, TopoDS_Shape> &additional_shapes, double fuzzy_tol);

    template<typename Schema>
    void write_spaces_to_stl_worker(std::list<cFace> &cFaces, std::list<Space> &spaces, const std::string &input_ifc_filepath, const std::string &output_ifc_file, std::string path, const std::string &mode, bool ascii, bool ignore_facade_space,
                                    const std::map<std::string, TopoDS_Shape> &additional_shapes, double fuzzy_tol);

    static void clear_cface_maps(std::list<cFace> &cFaces);

    static void process_space_shells_for_export(std::list<Space> &spaces, bool skip_facade);

    static TopoDS_Shape process_shape_for_export(TopoDS_Shape S);

    static void write_ifc_file(std::unique_ptr<IfcParse::IfcFile> &model, const std::string &output_filename);

    static void add_virtual_elements_to_model(std::unique_ptr<IfcParse::IfcFile> &model, const std::set<IfcUtil::IfcBaseEntity *> &virtual_products);

    bool add_entities_to_model(std::unique_ptr<IfcParse::IfcFile> &model, std::list<Space> &spaces, const std::set<IfcUtil::IfcBaseEntity *> &virtual_products, const std::unordered_map<std::string, std::list<TopoDS_Face>> &shadings);

    template<typename Schema>
    void add_entities_to_model_worker(std::unique_ptr<IfcParse::IfcFile> &model, std::list<Space> &spaces, const std::set<IfcUtil::IfcBaseEntity *> &virtual_products, const std::unordered_map<std::string, std::list<TopoDS_Face>> &shadings, IfcUtil::IfcBaseClass *IfcGeomReprContext_Model);

    static void remove_entities_from_model(std::unique_ptr<IfcParse::IfcFile> &model, const std::set<std::string> &remove_classes);

    static void remove_relationships_from_model(std::unique_ptr<IfcParse::IfcFile> &model);

    static void remove_empty_IfcRelDefines_from_model(std::unique_ptr<IfcParse::IfcFile> &model);

    static void remove_unused_instances_from_model(std::unique_ptr<IfcParse::IfcFile> &model);

    template<typename Schema>
    IfcUtil::IfcBaseClass *space_storey(const std::set<cFace *> &cfaces, IfcUtil::IfcBaseClass *backup_storey);

    template<typename Schema>
    typename Schema::IfcConnectionGeometry *create_connection_geometry(cFace *cface, gp_Trsf trsf, bool complement_faces, const std::string &schema);

    template<typename Schema>
    typename Schema::IfcProductDefinitionShape *createRepresentation(const TopoDS_Solid &solid, const std::string &schema, typename Schema::IfcGeometricRepresentationSubContext *IfcGeometricRepresentationSubContext);

    template<typename Schema>
    typename Schema::IfcSpace *
    space_to_ifc_space(Space &space, std::unique_ptr<IfcParse::IfcFile> &model, typename Schema::IfcOwnerHistory *IfcOwnerHistory, std::unordered_map<typename Schema::IfcBuildingStorey *, gp_Trsf> &storeyTrsf, const std::string &schema, typename Schema::IfcAxis2Placement *IfcAxis2Placement,
                       typename Schema::IfcGeometricRepresentationSubContext *IfcGeometricRepresentationSubContext, typename Schema::IfcSIUnit *VolumeUnit, typename Schema::IfcSIUnit *AreaUnit);

    template<typename Schema>
    typename Schema::IfcRelSpaceBoundary *cface_to_ifc_rel_space_boundary(cFace *cface, typename Schema::IfcSpaceBoundarySelect *IfcSpaceBoundarySelect, typename Schema::IfcOwnerHistory *IfcOwnerHistory, gp_Trsf trsf, bool complement_faces, bool second_lvl, const std::string &schema);

    template<typename Schema>
    typename Schema::IfcRelSpaceBoundary *cface_to_ifc_rel_space_boundary_shading(cFace *cface, typename Schema::IfcSpaceBoundarySelect *IfcSpaceBoundarySelect, typename Schema::IfcOwnerHistory *IfcOwnerHistory, gp_Trsf trsf, const std::string &schema, std::string description);

    static void reduce_edges_cFaces(std::list<cFace> &cFaces);

    double round_double_for_ifc_write(double d) const;

    template<typename Schema>
    void find_relevant_products_worker(std::unique_ptr<IfcParse::IfcFile> &model, const std::set<std::string> &include_entities, std::set<std::string> &non_void_filling_products, bool integrate_openings_into_walls, bool use_ifcopeningelelements_for_virtual_boundaries);

    template<typename Schema>
    void find_openings_in_products_worker(std::list<Product> &products, std::set<std::string> &guids);

    template<typename T>
    inline void free_container(T &p_container);

    template<typename Schema>
    void add_entities_to_model_worker_Ifc2x3(std::unique_ptr<IfcParse::IfcFile> &model, std::list<Space> &spaces, const std::set<IfcUtil::IfcBaseEntity *> &virtual_products, IfcUtil::IfcBaseClass *IfcGeomReprContext_Model);

    template<typename Schema>
    typename Schema::IfcSpace *
    space_to_ifc_space_Ifc2x3(Space &space, std::unique_ptr<IfcParse::IfcFile> &model, typename Schema::IfcOwnerHistory *IfcOwnerHistory, std::unordered_map<typename Schema::IfcBuildingStorey *, gp_Trsf> &storeyTrsf, const std::string &schema, typename Schema::IfcAxis2Placement *IfcAxis2Placement,
                              typename Schema::IfcGeometricRepresentationSubContext *IfcGeometricRepresentationSubContext, typename Schema::IfcSIUnit *VolumeUnit, typename Schema::IfcSIUnit *AreaUnit);

    template<typename Schema>
    typename Schema::IfcRelSpaceBoundary *cface_to_ifc_rel_space_boundary_Ifc2x3(cFace *cface, typename Schema::IfcSpace *IfcSpace, typename Schema::IfcOwnerHistory *IfcOwnerHistory, gp_Trsf trsf, const std::string &schema);

    static void log_spaces(const std::list<Space> &spaces);

    static void identify_blocked_faces(std::list<cFace> &cFaces);

    void create_void_filling_cuts_in_products(std::list<Product> &products, std::list<oFace> &orig_faces, std::unordered_map<IfcUtil::IfcBaseClass *, TopoDS_Shape> &opening_shapes, bool simplifyOpeningElement, std::unique_ptr<IfcParse::IfcFile> &model, IfcGeom::IteratorSettings settings) const;

    static bool create_void_filling_cuts_in_products_worker(Product &product, std::unordered_map<IfcUtil::IfcBaseClass *, TopoDS_Shape> &opening_shapes, bool simplifyOpeningElement, std::set<oFace *> &rem, std::list<opening_info> &ins, std::unordered_map<std::string, Product *> &togenerate);

    static OWP cvfcipw_collect_opening_shapes(const Product &product, std::unordered_map<IfcUtil::IfcBaseClass *, TopoDS_Shape> &opening_shapes, bool simplifyOpeningElement);

    static void cvfcipw_check_for_common_volume(Product &product, std::list<opening_info> &ins, std::unordered_map<std::string, Product *> &togenerate, size_t n, OWP &OSW);

    static void cvfcipw_commons_between_opening_and_wall_faces(Product &product, BOPAlgo_CellsBuilder &CellBuilder, TopTools_ListOfShape &take, TopTools_ListOfShape &avoi, const TopoDS_Compound &Shell, const std::set<unsigned int> &ids, const OWP &OSW, std::list<opening_info> &ins,
                                                               const std::pair<const unsigned int, std::list<oFace *>> &m);

    static std::set<unsigned int>
    cvfcipw_cut_openings_from_wall_faces(Product &product, BOPAlgo_CellsBuilder &CellBuilder, TopTools_ListOfShape &take, TopTools_ListOfShape &avoi, const TopoDS_Compound &Shell, const OWP &OSW,
                                         std::list<opening_info> &ins, const std::pair<const unsigned int, std::list<oFace *>> &m);

    static TopoDS_Compound cvfcipw_compound_from_orig_faces(size_t n, const std::list<oFace *> &orig_faces);

    static void cvfcipw_perform_fuse(const std::string &guid, BOPAlgo_CellsBuilder &CellBuilder, const TopoDS_Compound &Shell, const OWP &OSW);

    static TopoDS_Shape unify_shape(TopoDS_Shape &shp);

    static void update_space_neighbourships(std::list<Space> &spaces, double s);

    static Bnd_OBB obb(const TopoDS_ListOfShape &L);

    static TopoDS_Shape obb_to_shape(const Bnd_OBB &obb);

    template<typename Schema>
    bool is_product_part_of_ifccurtain_wall(typename Schema::IfcProduct *IfcProduct);

    template<typename Schema>
    std::list<typename Schema::IfcProduct *> get_related_objects_of_curtain_wall(typename Schema::IfcCurtainWall *IfcCurtainWall);

    void add_curtain_walls_to_products(std::unique_ptr<IfcParse::IfcFile> &model, IfcGeom::IteratorSettings settings, std::list<Product> &products);

    template<typename Schema>
    void add_curtain_walls_to_products_worker(std::unique_ptr<IfcParse::IfcFile> &model, IfcGeom::IteratorSettings settings, std::list<Product> &products);

    static Bnd_Box create_aabb_from_shape_list(const TopoDS_ListOfShape &L, double gap);

    static TopoDS_Shape aabb_to_shape(const Bnd_Box &bnd);

    static void calc_plane_uv_vectors(gp_Pnt O, gp_Pnt P, gp_Dir n, gp_Dir &u_plane, gp_Dir &v_plane);

    static void calc_point_uv_parameters_on_plane(gp_Pnt P, gp_Pnt O, gp_Vec u_plane, gp_Vec v_plane, double &u, double &v);

    template<typename Schema>
    void get_length_conversion_factor_worker(std::unique_ptr<IfcParse::IfcFile> &model);

    static void remove_products(std::list<Product> &products, const std::list<Product *> &products_to_remove);

    static bool find_duplicate_hashes_in_orig_faces(std::list<oFace> &orig_faces, std::set<unsigned int> &bad_hashes, std::set<unsigned int> &bad_hashes_location);

    static void remove_invalid_spaces(std::list<Space> &spaces, std::list<cFace> &cFaces);

    static void initialize_space_solid_classifiers(std::list<Space> &spaces);

    std::map<std::string, TopoDS_Shape> add_additional_shapes(std::unique_ptr<IfcParse::IfcFile> &model, IfcGeom::IteratorSettings settings, const std::set<std::string> &include_entities, const std::set<std::string> &cfd_entities, bool heal_shapes, bool triangulate_shapes) const;

    void log_orig_faces(const std::list<oFace> &orig_faces);

    template<typename Schema>
    void round_ifc_cartesian_points_traversal(std::unique_ptr<IfcParse::IfcFile> &model, IfcUtil::IfcBaseClass *E);

    template<typename Schema>
    typename Schema::IfcExternalSpatialElement *space_to_ifc_external_spatial_element(Space &space, std::unique_ptr<IfcParse::IfcFile> &model, typename Schema::IfcOwnerHistory *IfcOwnerHistory, std::unordered_map<typename Schema::IfcBuilding *, gp_Trsf> &buildingTrsf, const std::string &schema,
                                                                                      typename Schema::IfcAxis2Placement *IfcAxis2Placement, typename Schema::IfcGeometricRepresentationSubContext *IfcGeometricRepresentationSubContext);

    static void log_cFaces(const std::list<cFace> &cFaces);

    void log_cFace(const cFace &cface, unsigned int i);

    static void log_products(const std::list<Product> &products);

    static bool is_point_in_aabb(double A_Xmin, double A_Xmax, double A_Ymin, double A_Ymax, double A_Zmin, double A_Zmax, double P_X, double P_Y, double P_Z);

    static double overlapping_volume_aabbs(double A_Xmin, double A_Xmax, double A_Ymin, double A_Ymax, double A_Zmin, double A_Zmax, double B_Xmin, double B_Xmax, double B_Ymin, double B_Ymax, double B_Zmin, double B_Zmax);

    void link_old_IfcSpace_data(std::unique_ptr<IfcParse::IfcFile> &model, IfcGeom::IteratorSettings settings, std::list<Space> &spaces, ifcspaceInfoList &old_spaces);

    void split_spaces_by_IfcSpaces(std::list<Space> &spaces, std::list<cFace> &cFaces, std::list<oFace> &orig_faces, std::list<Product> &products, ifcspaceInfoList &old_spaces, unsigned int &fid);

    ifcspaceInfoList create_shapes_of_IfcSpaces(std::unique_ptr<IfcParse::IfcFile> &model, IfcGeom::IteratorSettings settings) const;

    static std::unordered_map<Space *, std::list<ifcspaceInfoList::iterator>> IfcSpaces_rtree(ifcspaceInfoList &old_spaces, std::list<Space> &spaces);

    static void link_spaces_and_IfcSpaces(std::list<Space> &spaces, std::unordered_map<Space *, std::list<ifcspaceInfoList::iterator>> &comp, bool reconsider);

    template<typename Schema>
    void old_spaces_split_space(Space *space, std::list<cFace> &cFaces, std::list<oFace> &orig_faces, std::list<Space> &spaces, std::list<Product> &products, unsigned int &space_id, unsigned int &fid, std::unordered_map<Space *, std::list<ifcspaceInfoList::iterator>> &comp,
                                std::set<cFace *> &del_cface, std::set<Space *> &del_space);

    template<typename Schema>
    void old_spaces_split_and_create_new_spaces(std::list<cFace> &cFaces, std::list<oFace> &orig_faces, std::list<Product> &products, std::list<Space> &spaces, std::unordered_map<Space *, std::list<ifcspaceInfoList::iterator>> &comp, unsigned int &fid);

    static NCollection_DataMap<TopoDS_Face, TopoDS_Face> split_spaces_by_faces(const TopoDS_ListOfShape &cutFaces, TopoDS_ListOfShape &subspaces, double V_crit, double h_crit);

    static std::list<std::list<TopoDS_Face>> virtual_face_pairs(TopoDS_ListOfShape &subspaces, const NCollection_DataMap<TopoDS_Face, TopoDS_Face> &H);

    static int double_to_scaled_int(double d, int s);

    template<typename Schema>
    void link_infos_of_spaces(std::unique_ptr<IfcParse::IfcFile> &model, std::list<Space> &spaces, std::unordered_map<Space *, std::list<ifcspaceInfoList::iterator>> &comp);

    static void modify_view_definition(std::unique_ptr<IfcParse::IfcFile> &model);

    template<typename Schema>
    bool does_product_decompose_another_product(typename Schema::IfcProduct *IfcProduct, bool consider_allowed_classes, const std::set<std::string> &inc);

    template<typename Schema>
    void add_decomposed_entities_to_products_worker(std::unique_ptr<IfcParse::IfcFile> &model, IfcGeom::IteratorSettings settings, std::list<Product> &products, const std::set<std::string> &include_entities);

    template<typename Schema>
    std::set<std::string> guids_of_products_decomposing_another_product(typename Schema::IfcProduct *IfcProduct, const std::set<std::string> &inc);

    TopoDS_ListOfShape products_to_shapelist(std::unique_ptr<IfcParse::IfcFile> &model, IfcGeom::IteratorSettings settings, const std::set<std::string> &guids) const;

    static std::string print_time(double t, const std::string &s1, const std::string &s2);

    static std::vector<std::string> split_string_by_character(const std::string &s, char delimiter);

    static void check_containment_in_fuse(const TopoDS_Shape &fuse, std::list<cFace> &cFaces, bool check_edges);

    static void check_corresponding_face_pairs(std::list<cFace> &cFaces);

    static void check_fixed_normal(std::list<cFace> &cFaces);

    static void check_redundant_faces(std::list<cFace> &cFaces, std::list<cFace> &innerFaces);

    static void check_adjacency_self_reference(std::list<cFace> &cFaces);

    static void check_holes(std::list<cFace> &cFaces);

    void check_loops(std::list<cFace> &cFaces, double tol);

    static void check_duplicate_vertices(std::list<cFace> &cFaces);

    static void check_superface(std::list<cFace> &cFaces);

    static void check_super_and_subfaces(std::list<cFace> &cFaces);

    static void check_convex(std::list<cFace> &cFaces);

    static void remove_corresponding_trash(std::list<cFace> &cFaces);

    static void get_zones_from_octree(std::vector<std::pair<std::array<double, 3>, std::set<std::string>>> &spaces_guids, std::vector<std::pair<std::array<double, 3>, std::set<unsigned int>>> &spaces_triangles, std::vector<std::set<unsigned int>> &zones, std::list<oFace> &orig_faces,
                                      double maximum_voxel_size, int maximum_octree_levels, double write_vtk, std::vector<oFace *> &tri2orig, bool apply_refinement_octrees);

    static void refined_octree(const std::array<double, 3> &fluid_point, std::set<unsigned int> &old_space_triangles, std::set<std::string> &old_spaces_guids, const std::vector<std::tuple<double, double, double>> &global_vertices, const std::vector<std::tuple<int, int, int>> &global_faces,
                               const std::vector<std::string> &global_attrs, double maximum_voxel_size, int maximum_octree_levels);

    static std::vector<std::pair<std::vector<gp_Pnt>, std::set<oFace *>>> get_origfaces_from_octree(std::list<oFace> &orig_faces, double maximum_voxel_size, int maximum_octree_levels, double write_vtk, bool apply_refinement_octrees);

    std::vector<std::tuple<std::array<double, 3>, TopoDS_Shape, std::list<cFace>>> perform_fuse_on_zone_builders(std::vector<BOPAlgo_Builder> &B, std::vector<std::list<oFace *>> &O, std::vector<NCollection_DataMap<TopoDS_Face, TopoDS_Face>> &N) const;

    static void log_zones(const std::vector<std::tuple<std::array<double, 3>, TopoDS_Shape, std::list<cFace>>> &zones);

    static void find_spaces_filter_start_face_by_fluid_points(std::list<cFace> &GlobalcFaces, std::list<cFace> &cFaces, std::list<Space> &spaces, unsigned int &space_id, const std::vector<gp_Pnt> &Points_air);

    static void remove_triangulation(std::list<oFace> &orig_faces);

    std::vector<std::tuple<std::vector<gp_Pnt>, TopoDS_Shape, std::list<cFace>>> fuse_zones_in_parallel(std::vector<std::pair<std::vector<gp_Pnt>, std::set<oFace *>>> zones_orig, double fuzzy_tol) const;

    std::vector<std::tuple<std::vector<gp_Pnt>, TopoDS_Shape, std::list<cFace>>> get_cfaces_from_octree(std::list<oFace> &orig_faces, double maximum_voxel_size, int maximum_octree_levels, double write_vtk, double fuzzy_tol, bool apply_refinement_octrees);

    static void nullify_orig_faces(std::list<oFace> &orig_faces);

    static void check_manifoldness(const std::list<cFace> &cFaces);

    static void remove_non_manifold_adjacency_by_angle(std::list<cFace> &cFaces);

    static TopoDS_ListOfShape create_extension_faces(const TopoDS_Edge &edge, gp_Dir n1, gp_Dir n2, double length);

    static bool check_faces(std::list<cFace> &cFaces, bool mark_as_trash, double l_warn, double l_crit);

    static bool check_face(const TopoDS_Face &F, unsigned int n, const std::string &g, double l_warn, double l_crit);

    static void check_orig_faces(std::list<oFace> &orig_faces);

    static bool fuse_cFaces(TopoDS_Shape &fuse, const std::list<cFace> &old_cFaces, std::list<cFace> &new_cFaces, double fuzzy_tol, unsigned int &fid);

    static void analyze_orig_faces(std::list<oFace> &orig_faces, double l_warn, double l_crit);

    static void check_fuse(const TopoDS_Shape &fuse);

    static void check_redundant_cFaces(std::list<cFace> &cFaces);

    static void check_adjacency(const TopoDS_Shape &fuse, std::list<cFace> &cFaces, bool b1 = false);

    static bool check_shape_tolerances(const TopoDS_Shape &S, double critical);

    static TopoDS_Shape fuse_shape(TopoDS_Shape &S, double tol, const std::string &s);

    static bool check_ifc_context(IfcGeom::Element<real_t> *geom_object);

    static TopoDS_Shape prism_from_face(const TopoDS_Face &F, gp_Dir n, double l, bool rm_tf, bool rm_bf);

    static void add_prism_faces_by_guid(std::list<oFace> &orig_faces, double offset, const std::set<std::string> &offset_guids);

    static void add_prism_faces(std::list<oFace> &orig_faces, double offset, const std::set<std::string> &classes_for_face_extension);

    static void identify_decoupled_offset_faces(std::list<cFace> &cFaces);

    static void add_offset_face_reversed_duplicates(std::list<cFace> &cFaces, unsigned int &fid);

    static void nonoffset_inner_or_coplanar_faces(std::list<cFace> &cFaces, std::list<cFace> &L, bool consider_duplicates = false);

    static void nonoffset_hanging_faces(std::list<cFace> &cFaces, std::list<cFace> &L);

    static void add_clipping_box_to_tree(gp_Vec v, double l, cFace *f, cface_tree3D &tree_inner);

    static void find_clipping_pairs(cfaceSetMap &M_planar, cfaceSetMap &M_curved, std::list<cFace> &cFaces, std::list<cFace> &innerFaces, double depth_2b);

    static TopoDS_ListOfShape clip_faces(std::list<cFace> &cFaces, cfaceSetMap &M_planar, cfaceSetMap &M_curved, double depth_2b, double fuzzy_tol);

    static TopoDS_ListOfShape clip_faces(const TopoDS_Face &F1, const TopoDS_Face &F2, double depth_2b);

    static bool fuse_clipped_faces(TopoDS_Shape &fuse, std::list<cFace> &cFaces, const TopoDS_ListOfShape &clips, std::list<cFace> &new_cFaces, double fuzzy_tol, unsigned int &fid);

    static void find_clipping_pairs_2(cfaceSetMap &M_planar, cfaceSetMap &M_curved, std::list<cFace> &cFaces, double depth_2b, double min_angle);

    static std::vector<gp_Pnt> remove_colinear_points_from_polygon(std::vector<gp_Pnt> Pnts, double tol);

    static std::vector<gp_Pnt2d> remove_colinear_points_from_polygon(std::vector<gp_Pnt2d> Pnts, double tol);

    static std::vector<gp_Pnt> remove_coincident_points_from_polygon(std::vector<gp_Pnt> Pnts, double tol);

    static std::vector<gp_Pnt2d> remove_coincident_points_from_polygon(std::vector<gp_Pnt2d> Pnts, double tol);

    static double unsigned_angle_between_points(gp_Pnt p1, gp_Pnt mid, gp_Pnt p3);

    static double unsigned_angle_between_points(gp_Pnt2d p1, gp_Pnt2d mid, gp_Pnt2d p3);

    static void check_containment_of_fuse_in_cFaces(const TopoDS_Shape &fuse, std::list<cFace> &cFaces);

    static TopoDS_Wire rebuild_wire(const TopoDS_Wire &w);

    static TopoDS_Face rebuild_face(const TopoDS_Face &f, bool &success);

    void rebuild_origfaces(std::list<oFace> &orig_faces) const;

    static std::vector<gp_Pnt> sort_points_from_polygon_by_min(std::vector<gp_Pnt> Pnts, cFace *cface);

    static std::vector<gp_Pnt> sort_points_from_polygon_by_dist(std::vector<gp_Pnt> Pnts, cFace *cface);

    static void check_parent_faces(std::list<cFace> &cFaces);

    static void check_adjacency_null(std::list<cFace> &cFaces);

    static void check_space_behind(std::list<cFace> &cFaces);

    static void check_pointers(std::list<cFace> &cFaces, bool secondlvl);

    static void check_enum_intext(std::list<cFace> &cFaces);

    static bool are_points_coincident(gp_Pnt p1, gp_Pnt p2, double tol);

    static bool are_points_coincident(gp_Pnt2d p1, gp_Pnt2d p2, double tol);

    static void change_cface_sbtype(std::list<cFace> &cFaces, std::list<cFace> &cFaces_2nd);

    static void spaces_introduce_second_level_cfaces(std::list<Space> &spaces);

    void clip_faces(std::list<cFace> &cFaces, std::list<cFace> &sec, cfaceSetMap &M_planar, cfaceSetMap &M_curved, double depth_2b, double fuzzy_tol, unsigned int &fid) const;

    static TopoDS_Shape clip_face_occ(cFace &cface, const std::set<cFace *> &planars, const std::set<cFace *> &curved, double length, double fuzzy_tol);

    static TopoDS_Shape clip_face(cFace &cface, const std::set<cFace *> &planars, double length, double fuzzy_tol);

    static void update_face_adjacencies(std::list<cFace *> &cFaces, const TopoDS_Shape &fuse);

    static std::unordered_map<unsigned int, cFace *> setup_face_cFace_map(std::list<cFace *> &cFaces);

    static void trash_2b_opening_faces(std::list<cFace> &cFaces);

    static void create_missing_parent_faces(std::list<cFace> &cFaces, std::list<oFace> &orig_faces, unsigned int &fid);

    static double area_by_triangulation(const TopoDS_Shape &shp);

    void site(std::unique_ptr<IfcParse::IfcFile> &model, IfcGeom::IteratorSettings settings, std::list<cFace> &cFaces, std::list<Space> &spaces, unsigned int &fid) const;

    TopoDS_ListOfShape site_ifcsite(std::unique_ptr<IfcParse::IfcFile> &model, IfcGeom::IteratorSettings settings) const;

    std::map<cFace *, std::set<cFace * >> site_clip(std::list<cFace> &cFaces, Space *facade_space, const rtree_lib::RTree<site_face *, double, 3, double> &tree, double site_min, double site_max, unsigned int &fid) const;

    static void site_update_face_relationships(std::list<cFace> &cFaces, Space *facade_space, std::map<cFace *, std::set<cFace * >> &oldToNew);

    static gp_Pnt project_point_onto_plane(gp_Pnt P, gp_Dir n, gp_Pnt A);

    static TopoDS_Wire rebuild_wire_projection(const TopoDS_Wire &w, gp_Dir n, gp_Pnt A);

    static void check_ifcproduct_pointers(std::unique_ptr<IfcParse::IfcFile> &model, std::list<Space> &spaces);

    static void add_attribut_strings(Argument *arg, const std::string &name, Space &space);

    template<typename Schema>
    typename Schema::IfcTable *sb_ifctable(std::list<Space> &spaces, typename Schema::IfcSIUnit *AreaUnit);

    template<typename Schema>
    std::list<typename Schema::IfcTable *> sb_material_ifctable(std::list<Space> &spaces, typename Schema::IfcSIUnit *AreaUnit);

    void complete_material_list_with_inner_faces(std::list<cFace> &cFaces, std::list<cFace> &cFaces_2ndLvl);

    void triangulate_cfaces_for_ray_tracer(std::list<cFace> &cFaces) const;

    static std::list<std::array<gp_Pnt, 3>> triangles_from_face(const TopoDS_Face &face);

    static void create_mesh_for_ray_tracer(std::list<cFace> &cFaces, std::list<std::array<double, 3>> &vertices, std::list<std::array<unsigned int, 3>> &triangles, std::list<unsigned int> &attributes, std::unordered_map<unsigned int, cFace *> &M);

    static void create_rays_for_ray_tracer(std::list<cFace> &cFaces_2ndLvl, std::vector<cFace *> &send, std::vector<IntersectionRay> &rays, double tol);

    static void create_rays_for_ray_tracer(std::list<cFace> &cFaces_2ndLvl, std::vector<cFace *> &send, std::vector<IntersectionRay> &rays, double tol, double l);

    static void perform_ray_tracing(IntersectorInterface &intersector, const std::vector<cFace *> &send, std::vector<IntersectionRay> &rays, double tol, std::unordered_map<unsigned int, cFace *> &M);

    static void perform_ray_tracing(IntersectorInterface &intersector, const std::vector<cFace *> &send, std::vector<IntersectionRay> &rays, double tol, std::unordered_map<unsigned int, cFace *> &M, double min_angle, gp_XYZ bounds_min);

    static void perform_ray_tracing_first_level(IntersectorInterface &intersector, const std::vector<cFace *> &send, std::vector<IntersectionRay> &rays, double tol, std::unordered_map<unsigned int, cFace *> &M, gp_XYZ bounds_min);

    void calculate_points_on_faces(std::list<cFace> &cFaces_2ndLvl) const;

    void identify_sb_types_ray_tracing(std::list<cFace> &cFaces, gp_XYZ bounds_min, double transmission_length, double min_angle);

    void identify_sb_types_ray_tracing_first_level(std::list<cFace> &cFaces, gp_XYZ bounds_min, double transmission_length);

    void revise_is_external_attribute(std::unique_ptr<IfcParse::IfcFile> &model, std::list<Space> &spaces);

    static void remove_holes(std::list<cFace> &cFaces, unsigned int &fid);

    static void decompose_concave_polygons(std::list<cFace> &cFaces, unsigned int &fid);

    static bool subpolygon_check(const cxdsb::ConcavePolygon &subPolygon, std::list<std::vector<gp_Pnt2d>> &result, const std::vector<gp_Pnt2d> &L, const cFace &cface);

    static void cface_linked_naive_triangulation(std::list<cFace> &cFaces, cFace &cface, unsigned int &fid, const std::vector<cFace *> &inner_boundaries, unsigned int &c);

    static void cface_linked_hole_removal(std::list<cFace> &cFaces, cFace &cface, unsigned int &fid, const std::vector<cFace *> &inner_boundaries, unsigned int &c);

    static void cface_linked_concave_decomposing(std::list<cFace> &cFaces, cFace &cface, unsigned int &fid, const std::vector<cFace *> &inner_boundaries, unsigned int &c);

    static void decompose_concave_polygons_triangulation(std::list<cFace> &cFaces, unsigned int &fid);

    static void simplify_fenestration_faces(std::list<cFace> &cFaces, unsigned int &fid);

    static void split_corresponding_equivalent(const cFace &cface, const std::set<cFace *> &NEW, unsigned int &fid, std::list<cFace> &cFaces);

    static gp_Pnt project_inv(double X, double Y, long proj_axis, gp_Dir n, gp_Pnt R);

    template<typename Schema>
    void revise_is_external_attribute_worker(const std::set<IfcUtil::IfcBaseEntity *> &entities, bool b);

    template<typename Schema>
    void link_relationships(const Space &space, typename Schema::IfcSpace *IfcSpace);

    template<typename Schema>
    void set_old_space_info(Space &space, typename Schema::IfcSpace *IfcSpace);

    static void clip_face_add_multiple_faces(cFace &cface, const TopoDS_Shape &C, std::list<cFace> &sec, unsigned int &fid);

    static void clip_face_add_single_face(cFace &cface, std::list<cFace> &sec, unsigned int &fid);

    void calculate_shading(std::list<cFace> &cFaces, std::list<Space> &spaces, std::unordered_map<std::string, std::list<TopoDS_Face>> &shadings, std::unique_ptr<IfcParse::IfcFile> &model, IfcGeom::IteratorSettings settings);

    static std::list<TopoDS_Face> calculate_shading_facade(std::list<cFace> &cFaces, std::list<Space> &spaces);

    std::list<TopoDS_Face> calculate_shading_extern(std::list<Space> &spaces, std::unique_ptr<IfcParse::IfcFile> &model, IfcGeom::IteratorSettings settings) const;

    static TopoDS_Shape sew_cfaces(const std::set<cFace *> &faces, double tol, bool nonmanifold);

    static void update_adjacence_map_after_unification(std::list<cFace> &cFaces);

    static double minimal_distance_between_shapes(const TopoDS_Shape &b1, const TopoDS_Shape &b2);

    static bool accept_solid_split(const TopoDS_ListOfShape &solids, double V_crit, double h_crit);

    static std::list<ifcspaceInfoList::iterator> link_spaces_and_IfcSpaces_worker(const Space &space, const std::list<ifcspaceInfoList::iterator> &comp, const TopoDS_Solid &new_solid);

    bool create_connection_geometry_loop(cFace *cface, const TopoDS_ListOfShape &vertices, std::vector<gp_Pnt> &Pnts, double tol) const;

    template<typename Schema>
    typename IfcTemplatedEntityList<typename Schema::IfcCartesianPoint>::ptr create_connection_geometry_cartesian_points(const std::vector<gp_Pnt> &Pnts, gp_Pnt O, gp_Dir u_plane, gp_Dir v_plane);

    static gp_Pnt2d centroid_2D_polygon(const std::vector<gp_Pnt2d> &poly);

    static bool is_point_inside_2D_polygon(const std::vector<gp_Pnt2d> &poly, const gp_Pnt2d &p);

    static double Angle2D(double x1, double y1, double x2, double y2);

    // space clipper approach
    bool
    add_entities_to_model_space_approach(std::unique_ptr<IfcParse::IfcFile> &model, std::list<Space> &spaces, const std::set<IfcUtil::IfcBaseEntity *> &virtual_products, const std::unordered_map<std::string, std::list<TopoDS_Face>> &shadings, std::unordered_map<Space *, ifcspaceInfo *> &space_map);

    template<typename Schema>
    void add_entities_to_model_space_approach_worker(std::unique_ptr<IfcParse::IfcFile> &model, std::list<Space> &spaces, const std::set<IfcUtil::IfcBaseEntity *> &virtual_products, const std::unordered_map<std::string, std::list<TopoDS_Face>> &shadings,
                                                     IfcUtil::IfcBaseClass *IfcGeomReprContext_Model, std::unordered_map<Space *, ifcspaceInfo *> &space_map);

    static void find_clipping_pairs_all(std::map<cFace *, std::set<cFace * >> &M_planar, std::map<cFace *, std::set<cFace * >> &M_curved, std::list<cFace> &space_faces, std::list<cFace> &element_faces, double depth_2b);

    static void find_clipping_pairs_space_approach(std::map<cFace *, std::set<cFace * >> &M_planar, std::map<cFace *, std::set<cFace * >> &M_curved, std::list<cFace> &cFaces, double depth_2b, double min_angle);

    static void
    collect_spaces_and_ifcfaces(std::list<Product> &products, ifcspaceInfoList &ifcspaces, std::list<oFace> &ifcspace_ifc_faces, std::list<cFace> &ifcspace_faces, std::list<Space> &spaces, unsigned int &space_id_counter, std::unordered_map<Space *, ifcspaceInfo *> &M, unsigned int &fid);

    void identify_sb_types_ray_tracing_space_approach(std::list<cFace> &cFaces, gp_XYZ bounds_min, double transmission_length, double min_angle);

    static void perform_ray_tracing_space_approach(IntersectorInterface &intersector, const std::vector<cFace *> &send, std::vector<IntersectionRay> &rays, double tol, std::unordered_map<unsigned int, cFace *> &M, double min_angle, gp_XYZ bounds_min);

    void perform_ray_tracing_space_approach_elements(std::list<cFace> &element_faces, std::list<cFace> &cFaces_2ndLvl, double transmission_length);

    static void create_rays_for_ray_tracer_space_approach(std::list<cFace> &cFaces_2ndLvl, std::vector<cFace *> &send, std::vector<IntersectionRay> &rays, double tol, double transmission_length);

    static void perform_ray_tracing_space_approach(IntersectorInterface &intersector, const std::vector<cFace *> &send, std::vector<IntersectionRay> &rays, double tol, std::unordered_map<unsigned int, cFace *> &M);

    void create_virtual_elements_space_approach(std::list<cFace> &cFaces_2ndLvl, std::list<Product> &products, std::list<oFace> &ifcspace_ifc_faces);

    template<typename Schema>
    void create_virtual_elements_space_approach_worker(std::list<cFace> &cFaces_2ndLvl, std::list<Product> &products, std::list<oFace> &ifcspace_ifc_faces);

    static void check_corresponding_face_pairs_space_approach(std::list<cFace> &cFaces);

    static unsigned int number_of_faces_of_unknown_orientation(const std::list<cFace> &cFaces);

    static void ifaces_nbs(iFace &f1, used_orientation o1, std::unordered_map<cFace *, iFace *> &M);

    static void ifaces_nbs(iFace &iface, std::unordered_map<cFace *, iFace *> &M);

    static std::stack<oriFace> find_spaces_normals_unknown_start_face(std::queue<oriFace> &q);

    static bool find_spaces_normals_unknown_accept_solution(const std::set<oriFace> &L);

    static void identify_duplicate_faces_unknown_normals(std::list<cFace> &cFaces, std::unordered_map<std::string, unsigned int> ranks);

    static double find_spaces_normals_unknown_volume(const std::set<oriFace> &L);

    static void add_face_reversed(std::list<cFace> &cFaces, cFace *cface, unsigned int &fid);

    static std::multimap<double, std::set<oriFace>> find_spaces_normals_unknown_recursion(std::list<iFace> &ifaces, const std::list<cFace> &cFaces, std::queue<oriFace> &q);

    static void check_closed_space_edge_id(std::list<Space> &spaces);

    static void find_spaces_normals_unknown_save_spaces(std::list<cFace> &cFaces, std::list<Space> &spaces, unsigned int &space_id, unsigned int &fid, const std::multimap<double, std::set<oriFace>> &components, std::unordered_map<cFace *, iFace *> &M);

    static void find_spaces_normals_unknown_save_spaces_component(const std::pair<const double, std::set<oriFace>> &i, std::list<cFace> &cFaces, std::list<Space> &spaces, unsigned int &space_id, unsigned int &fid, const std::multimap<double, std::set<oriFace>> &components);

    static void find_spaces_normals_unknown_face_usages(const std::multimap<double, std::set<oriFace>> &components);

    static std::queue<oriFace> find_spaces_normals_unknown_queue(std::list<iFace> &ifaces);

    static std::string remove_first_and_last_char(std::string s);

    void bad_openings_to_box_faces(std::list<oFace> &orig_faces, std::unique_ptr<IfcParse::IfcFile> &model, IfcGeom::IteratorSettings settings, std::list<Product> &products, std::unordered_map<std::string, Product *> &todo) const;

    static bool clip_faces_surface_check(const TopoDS_Shape &C, const cFace &cface, bool level = false);

    static gp_Pnt mid_point(gp_Pnt P1, gp_Pnt P2);

    static void remove_inner_window_and_door_faces(std::list<cFace> &cFaces);

    static double polygonize_shape_2a_curvature_distance(const TopoDS_Shape &F_tool, const TopoDS_Shape &common, const segmentation_prism &p);

    friend class Graph;

    friend class Clip;
};

#endif //KERNEL_H
