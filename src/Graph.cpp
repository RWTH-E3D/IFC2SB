// Copyright 2022 Eric Fichter
#include "Graph.h"

Graph::Graph(std::string _input,
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
             bool _use_spaces_for_virtual_boundaries) :
        input(std::move(_input)),
        output(std::move(_output)),
        num_threads(_num_threads),
        stl(_stl),
        space_split(_space_split),
        max_transmission_length(_max_transmission_length),
        first_level_only(_first_level_only),
        remove_holes_by_face_split(_remove_holes_by_face_split),
        decompose_concave_polygons(_decompose_concave_polygons),
        simplify_fenestration_faces(_simplify_fenestration_faces),
        calculate_shadings(_calculate_shadings),
        use_ifcopening_elements_for_virtual_boundaries(_use_ifcopeningelelements_for_virtual_boundaries),
        use_spaces_for_virtual_boundaries(_use_spaces_for_virtual_boundaries) {

    products.clear();
    ifc_faces.clear();
    faces_1st_level.clear();
    faces_non_sb.clear();
    faces_2nd_level.clear();
    faces_trash.clear();
    spaces.clear();
    virtual_products.clear();
    shadings.clear();
    ifcspaces.clear();
    bounds_min.SetCoord(0, 0, 0);
    bounds_max.SetCoord(0, 0, 0);

    include_entities = {
            "IfcBeam",
            //"IfcBearing",
            // "IfcBuildingElementProxy",
            //"IfcChimney",
            "IfcColumn",
            //"IfcCovering",
            "IfcCurtainWall",
            //"IfcDeepFoundation",
            "IfcDoor",
            //"IfcFooting",
            "IfcMember",
            "IfcPlate",
            // "IfcRailing",
            //"IfcRamp",
            //"IfcRampFlight",
            "IfcRoof",
            //"IfcShadingDevice",
            "IfcSlab",
            //"IfcStair",
            //"IfcStairFlight",
            "IfcWall",
            "IfcWindow"
    };

    ranks["IfcSlab"] = 0;
    ranks["IfcRoof"] = 1;
    ranks["IfcWall"] = 2;
    ranks["IfcDoor"] = 3;
    ranks["IfcWindow"] = 4;
    ranks["IfcColumn"] = 5;
    ranks["IfcBeam"] = 6;
    ranks["IfcCovering"] = 7;
    ranks["IfcCurtainWall"] = 8;
    ranks["IfcRailing"] = 9;
    ranks["IfcPlate"] = 10;
    ranks["IfcMember"] = 11;
    ranks["IfcFooting"] = 12;
    ranks["IfcBearing"] = 13;
    ranks["IfcRamp"] = 14;
    ranks["IfcRampFlight"] = 15;
    ranks["IfcStair"] = 16;
    ranks["IfcStairFlight"] = 17;
    ranks["IfcDeepFoundation"] = 18;
    ranks["IfcChimney"] = 19;
    ranks["IfcShadingDevice"] = 20;
    ranks["IfcVirtualElement"] = 21;

    cfd_entities = {"IfcProduct"};
    offset_entities = {"IfcWall"};
    fuzzy_tol = 0.0001;
    face_offset_length = 0;
    space_id_counter = 0;
    face_id_counter = 0;
    integrate_openings_into_walls = true;

    settings.set(IfcGeom::IteratorSettings::FASTER_BOOLEANS, true);
    settings.set(IfcGeom::IteratorSettings::SEW_SHELLS, true);
    settings.set(IfcGeom::IteratorSettings::USE_WORLD_COORDS, false);
    settings.set(IfcGeom::IteratorSettings::DISABLE_TRIANGULATION, true);
    settings.set(IfcGeom::IteratorSettings::STRICT_TOLERANCE, true);
    settings.set(IfcGeom::IteratorSettings::DISABLE_OPENING_SUBTRACTIONS, integrate_openings_into_walls);
}

bool Graph::run() {

    Kernel K(num_threads);

    std::unordered_map<IfcUtil::IfcBaseClass *, TopoDS_Shape> opening_shapes;
    if (!prepare_products(K, opening_shapes))
        return false;

    if (!prepare_faces(K, opening_shapes))
        return false;

    if (!calc_faces_first_level(K))
        return false;

    if (!process_spaces(K))
        return false;

    if (stl) return true;

    if (!first_level_only)
        if (!calc_faces_second_level(K))
            return false;

    if (!process_model(K))
        return false;

    // check results written to model
    IfcCheck ifccheck(model.get(), num_threads);

#ifdef VISUALIZATION
    // Viewer::visualize_products(products);
    // Viewer::visualize_orig_faces(ifc_faces);
    // Viewer::visualize_cFaces(faces_1st_level);
    // Viewer::visualize_spaces(spaces, faces_1st_level, true);
    // Viewer::visualize_cFaces_as_space_boundaries(faces_2nd_level);
#endif

    return true;
}

// Products
bool Graph::prepare_products(Kernel &K, std::unordered_map<IfcUtil::IfcBaseClass *, TopoDS_Shape> &opening_shapes) {

    std::cout << "\n### Products" << std::endl;

    //***************************************************************
    // Parsing Ifc file to model
    if (!Kernel::read_ifc_file(input, model))
        return false;
    //***************************************************************

    //***************************************************************
    // Identify ifc schema of model and store in kernel
    if (!K.check_ifc_schema(model))
        return false;
    //***************************************************************

    //***************************************************************
    if (!K.get_length_conversion_factor(model))
        return false;
    if (!K.geom_context_info(model))
        return false;
    //***************************************************************

    //***************************************************************
    // Create container for products, that are not void filling products (e.g. walls)
    std::set<std::string> non_void_filling_products;
    //***************************************************************

    //***************************************************************
    // Find products and their guids that are from a desired class and are a) not void filling products or b) void filling products
    K.find_relevant_products(model, include_entities, non_void_filling_products, integrate_openings_into_walls, use_ifcopening_elements_for_virtual_boundaries);
    //***************************************************************

    //***************************************************************
    // Generating shapes from specified Ifc entities
    if (!K.generate_shapes_from_ifc_guids(model, settings, products, non_void_filling_products, bounds_min, bounds_max))
        return false;
    //***************************************************************

    //***************************************************************
    // Simplify products (windows and doors that lacked an IfcOpeningElement) by bounding boxes
    K.simplify_products(products);
    //***************************************************************

    //***************************************************************
    // Big tolerances of vertices, edges and faces can lead to failing fuse
    Kernel::check_product_tolerances(products, fuzzy_tol);
    //***************************************************************

    //***************************************************************
    // Add entities to products that are composed of sub-entities, e.g. IfcCurtainWall. The shape can be simplified
    K.add_decomposed_entities_to_products(model, settings, products, include_entities);
    //***************************************************************

    //***************************************************************
    // Generating shapes from specified Ifc entities
    // if (!K.generate_shapes_from_ifc_classes(model, products, include_entities, bounds_min, bounds_max))
    //     return false;
    //***************************************************************

    if (integrate_openings_into_walls) {
        //***************************************************************
        // For all created products (non void filling products) find non-filled openings that will be cut from product
        std::set<std::string> opening_guids = K.find_openings_in_products(products);
        // Create shapes of openings
        opening_shapes = K.generate_opening_shapes(model, settings, opening_guids);
        // Cut opening shapes from products
        Kernel::subtract_openings_from_products(products, opening_shapes);
        //***************************************************************

        //***************************************************************
        // Clear container, part I
        for (auto &product: products) product.IfcOpeningElements.clear();
        non_void_filling_products.clear(); // K.free_container(normal_products);
        opening_guids.clear(); // K.free_container(opening_guids);
        //***************************************************************
    }

    //***************************************************************
    // Add agglomerations of IfcCurtainWall to products, Not necessary anymore
    //if (consider_curtain_walls) K.add_curtain_walls_to_products(model, products);
    //***************************************************************

    //***************************************************************
    // Check if all faces belong to a shell and shells are closed and correctly oriented
    std::list<Product *> products_with_faulty_shape = Kernel::check_shells(products);
    if (!products_with_faulty_shape.empty()) {
        // Try to heal shapes
        Kernel::heal_products(products_with_faulty_shape);
        // Check again
        products_with_faulty_shape = Kernel::check_shells(products_with_faulty_shape);
        if (!products_with_faulty_shape.empty())
            for (auto &p: products_with_faulty_shape) {
                TopoDS_ListOfShape L;
                L.Append(p->shape);
                TopoDS_Shape box = Kernel::best_fitting_bbox(L);
                if (!box.IsNull()) {
                    p->shape = box;
                    std::cerr << "[Info] Use bounding box as shape for " << p->ifcproduct->data().toString() << std::endl;
                }
            }
            // Delete products that couldn't be healed
            //K.remove_products(products, products_with_faulty_shape);
        else
            std::cout << "\t[Info] All issues resolved." << std::endl;
    }
    products_with_faulty_shape.clear();
    //***************************************************************

    //***************************************************************
    // Checks if shapes consist of planar faces only. If not, planarize them
    if (Kernel::identify_non_planar_products(products)) {
        Kernel::polygonize_products(products);
        if (Kernel::identify_non_planar_products(products))
            return false;
    }
    //***************************************************************

    //***************************************************************
    // Cut some shapes
    std::set<std::string> cut;
    if (!cut.empty())
        K.cut_products(products, cut, fuzzy_tol);
    //***************************************************************

    //***************************************************************
    // Analyze product shapes, set flag in product (true if product shape is correct)
    double l_crit = 1.0e-6;
    Kernel::check_products(products, 5 * l_crit, l_crit);
    //***************************************************************

    //***************************************************************
    // Fuse products with false flag
    K.fuse_products(products, l_crit);
    //***************************************************************

    //***************************************************************
    // Analyze product shapes, set flag in product (true if product shape is correct)
    Kernel::check_products(products, 2 * l_crit, l_crit);
    //***************************************************************

    //***************************************************************
    // Remove products by flag
    Kernel::remove_products(products);
    //***************************************************************

    //***************************************************************
    // Create shapes of IfcSpaces
    ifcspaces = K.create_shapes_of_IfcSpaces(model, settings);
    //***************************************************************

    return !products.empty();
}

// Original Faces
bool Graph::prepare_faces(Kernel &K, std::unordered_map<IfcUtil::IfcBaseClass *, TopoDS_Shape> &opening_shapes) {

    std::cout << "\n### Faces" << std::endl;

    //***************************************************************
    // Save all faces of product shapes
    Kernel::collect_original_faces(products, ifc_faces);
    //***************************************************************

    //***************************************************************
    // Get virtual faces by space x space clipping
    if (use_spaces_for_virtual_boundaries)
        K.get_virtual_faces(model, ifcspaces, ifc_faces, products, virtual_products);
    //Viewer::visualize_orig_faces_virtual_product(ifc_faces);
    //***************************************************************

    //***************************************************************
    // Add extended faces
    if (face_offset_length > 1.0e-6) {
        // if (!stl)
        //    K.add_offset_faces(ifc_faces, face_offset_length, offset_entities);
        Kernel::add_prism_faces(ifc_faces, face_offset_length, offset_entities); // prism extrusion also extrudes faces adjacent to inner wire!
        //std::set<std::string> offset_guids = {"2ptk1k7qn8_Qk22vjh$0DE"};//, "3jjW3rL656ex34Gws22EfM"
        //Kernel::add_prism_faces_by_guid(ifc_faces, face_offset_length, offset_guids);
    }
    //***************************************************************

    //***************************************************************
    // Create opening faces in walls
    if (integrate_openings_into_walls) {
        K.create_void_filling_cuts_in_products(products, ifc_faces, opening_shapes, false, model, settings);
        //  for (auto &f: ifc_faces)
        //      if (f.IsIfcClass("IfcDoor") || f.IsIfcClass("IfcWindow"))
        //          if.SetOpening(f.RelProduct()->ifcproduct);
    }
    //***************************************************************

    //***************************************************************
    // Clear container, part II
    for (auto &product: products) {
        product.VoidFillingElements.clear();
        product.orig_faces.clear();
        product.shape = TopoDS_Shape();
    }
    opening_shapes.clear(); // K.free_container(opening_shapes);
    //***************************************************************

    //***************************************************************
    // Find duplicate hash codes in ifc_faces
    Kernel::analyze_orig_faces(ifc_faces, 3 * fuzzy_tol, 1.5 * fuzzy_tol);
    Kernel::check_orig_faces(ifc_faces);
    std::set<unsigned int> bad_hashes, bad_hashes_location;
    if (!Kernel::find_duplicate_hashes_in_orig_faces(ifc_faces, bad_hashes, bad_hashes_location)) {
        std::cerr << "[Warning] Duplicate hashes (" << bad_hashes.size() << ", " << bad_hashes_location.size() << ")." << std::endl;
        for (auto &h: bad_hashes) std::cerr << "\tBad hash: " << h << std::endl;
        for (auto &h: bad_hashes_location) std::cerr << "\tBad loc hash: " << h << std::endl;
        for (auto &f: ifc_faces) if (bad_hashes.find(Kernel::hash(f.face)) != bad_hashes.end()) std::cerr << "\tBad hash face: " << f.Info() << "\t" << f.face.Location().HashCode(INT_MAX) << std::endl;
    }
    bad_hashes.clear();
    //***************************************************************

    std::cout << std::flush;

    return !ifc_faces.empty();
}

// First level faces_1st_level
bool Graph::calc_faces_first_level(Kernel &K) {

    std::cout << "\n### Boundaries" << std::endl;

    //***************************************************************
    // Remove triangulation attached to faces, to reduce fusing time
    Kernel::remove_triangulation(ifc_faces);
    //***************************************************************

    //***************************************************************
    // Fuse all faces
//    TopoDS_Shape fuse;
//    // K.fuse_original_faces_parallel(fuse, ifc_faces, faces_1st_level, fuzzy_tol);
//    if (!K.fuse_original_faces(fuse, ifc_faces, faces_1st_level, fuzzy_tol))
//        return false;
//    K.nullify_orig_faces(ifc_faces);
    //***************************************************************

    //***************************************************************
    // Fuse all faces
    TopoDS_Shape fuse;
    // K.fuse_original_faces_parallel(fuse, ifc_faces, faces_1st_level, fuzzy_tol);
    if (!Kernel::fuse_original_faces(fuse, ifc_faces, faces_1st_level, fuzzy_tol, face_id_counter))
        return false;

    while (true) {
        if (Kernel::check_faces(faces_1st_level, true, fuzzy_tol, 0.5 * fuzzy_tol)) break;

        std::cerr << "[Error] Check found faulty faces!" << std::endl;
        Kernel::remove_trash(faces_1st_level, faces_trash);
        std::list<cFace> temp_cFaces = faces_1st_level;
        faces_1st_level.clear();
        if (!Kernel::fuse_cFaces(fuse, temp_cFaces, faces_1st_level, fuzzy_tol, face_id_counter)) return false;
    }

    Kernel::nullify_orig_faces(ifc_faces);
    //***************************************************************

    //***************************************************************
    // Check some face properties
    // if (!K.check_faces(faces_1st_level)) std::cerr << "[Error] Check found faulty faces!" << std::endl;
    //***************************************************************

    //***************************************************************
    // Check fuse
    Kernel::check_fuse(fuse);
    Kernel::check_fixed_normal(faces_1st_level);
    Kernel::check_redundant_cFaces(faces_1st_level);
    Kernel::check_containment_in_fuse(fuse, faces_1st_level, false);
    Kernel::check_containment_of_fuse_in_cFaces(fuse, faces_1st_level);
    Kernel::check_duplicate_vertices(faces_1st_level);
    //***************************************************************

    //***************************************************************
    // this mode has been used without virtual faces
    bool offsets_have_fixed_normal = false;
    if (offsets_have_fixed_normal)
        for (auto &cface: faces_1st_level)
            if (cface.IsOffset()) cface.SetNormalStatus(FACE_NORMAL_KNOWN);
    //***************************************************************

    //***************************************************************
    if (Kernel::number_of_faces_of_unknown_orientation(faces_1st_level) == 0)
        calc_faces_first_level_normals_known(fuse);
    else
        calc_faces_first_level_normals_unknown(fuse, K);
    //***************************************************************

    return !faces_1st_level.empty();
}

void Graph::calc_faces_first_level_normals_known(const TopoDS_Shape &fuse) {

    std::cout << "All face normals are known." << std::endl;

    //***************************************************************
    // If face normals were changed during fusing, correct them
    Kernel::correct_face_normals(faces_1st_level);
    //***************************************************************

    //***************************************************************
    // faces_1st_level get their non-seam half-edges from their face
    Kernel::update_half_edges(faces_1st_level);
    //***************************************************************

    //***************************************************************
    // Check some adjacence properties between fuse and faces_1st_level
    Kernel::check_adjacency(fuse, faces_1st_level);
    //***************************************************************

    //***************************************************************
    // flags enclosed faces as trash
    Kernel::identify_enclosed_faces(fuse, faces_1st_level);
    //***************************************************************

    //***************************************************************
    // get copy of inner faces
    if (!stl)
        Kernel::nonoffset_inner_or_coplanar_faces(faces_1st_level, faces_non_sb);
    //***************************************************************

    //***************************************************************
    // deletes trash faces from faces_1st_level and moves them to faces_trash
    Kernel::remove_trash(faces_1st_level, faces_trash);
    // Attention. Now there are faces in fuse shape, that are not in faces_1st_level anymore
    //***************************************************************

    //***************************************************************
    // identifies coplanar faces. Will keep one of them if necessary
    Kernel::identify_duplicate_faces(faces_1st_level, ranks);
    //***************************************************************

    //***************************************************************
    // get copy of inner faces (material-material)
    if (!stl)
        Kernel::nonoffset_inner_or_coplanar_faces(faces_1st_level, faces_non_sb);
    //***************************************************************

    //***************************************************************
    // deletes trash faces from faces_1st_level and moves them to faces_trash
    Kernel::remove_trash(faces_1st_level, faces_trash);
    //***************************************************************

    //***************************************************************
    // checks for duplicate HashCodes in faces_1st_level' TopoDS_Faces
    Kernel::check_duplicate_faces(faces_1st_level);
    // In faces_1st_level there are no faces with same HashCode anymore
    //***************************************************************

    //***************************************************************
    // create network of cface - halfedge - cface
    Kernel::update_face_adjacencies(faces_1st_level, fuse);
    //***************************************************************

    //***************************************************************
    // Add reversed duplicates of offset faces, this will lead to duplicated hashes again
    Kernel::add_offset_face_reversed_duplicates(faces_1st_level, face_id_counter);
    //***************************************************************

    //***************************************************************
    // adjacencies are removed, if faces can not build a proper shell
    Kernel::remove_adjacency_by_orientation(faces_1st_level);
    //***************************************************************

    //***************************************************************
    // identify faces that don't have adjacent faces on all edges
    Kernel::identify_hanging_faces(faces_1st_level);
    //***************************************************************

    //***************************************************************
    // get copy of hanging faces
    if (!stl)
        Kernel::nonoffset_hanging_faces(faces_1st_level, faces_non_sb);
    //***************************************************************

    //***************************************************************
    // deletes trash faces from faces_1st_level and moves them to faces_trash
    // adjacence info of a cface' edges are removed, if adjacent cface is trash
    Kernel::remove_trash_and_face_adjacency(faces_1st_level);
    //***************************************************************

    //***************************************************************
    // identifies offset faces that are not in contact with shell anymore
    Kernel::identify_decoupled_offset_faces(faces_1st_level);
    //***************************************************************

    //***************************************************************
    // identify blocked faces, lying behind non hanging offset faces
    // K.identify_blocked_faces(faces_1st_level);
    //***************************************************************

    //***************************************************************
    Kernel::remove_trash_and_face_adjacency(faces_1st_level);
    //***************************************************************

    //***************************************************************
    // check faces having a non-manifold connection to other faces on an edge
    Kernel::check_manifoldness(faces_1st_level);
    //***************************************************************

    //***************************************************************
    // faces having a non-manifold connection to other faces maintain only on adjacent face
    // (with opposite edge direction)
    Kernel::remove_non_manifold_adjacency_by_angle(faces_1st_level);
    //***************************************************************

    //***************************************************************
    // check faces having a non-manifold connection to other faces on an edge
    Kernel::check_manifoldness(faces_1st_level);
    //***************************************************************

    //***************************************************************
    // find components in building graph aka spaces
    Kernel::check_adjacency_self_reference(faces_1st_level);
    Kernel::find_spaces(faces_1st_level, spaces, space_id_counter);
    Kernel::check_redundant_faces(faces_1st_level, faces_non_sb);
    Kernel::check_closed_space_edge_id(spaces);
    Kernel::check_fixed_normal(faces_1st_level);
    //***************************************************************
}

void Graph::calc_faces_first_level_normals_unknown(const TopoDS_Shape &fuse, Kernel &K) {

    std::cout << "There are unknown face normals." << std::endl;

    //***************************************************************
    // faces_1st_level get their non-seam half-edges from their face
    Kernel::update_half_edges(faces_1st_level);
    //***************************************************************

    //***************************************************************
    // Check some adjacence properties between fuse and faces_1st_level
    Kernel::check_adjacency(fuse, faces_1st_level);
    //***************************************************************

    //***************************************************************
    // identifies coplanar faces. Will keep one of them if necessary
    Kernel::identify_duplicate_faces_unknown_normals(faces_1st_level, ranks);
    //***************************************************************

    //***************************************************************
    // get copy of inner faces (material-material)
    if (!stl)
        Kernel::nonoffset_inner_or_coplanar_faces(faces_1st_level, faces_non_sb);
    //***************************************************************

    //***************************************************************
    // deletes trash faces from faces_1st_level and moves them to faces_trash
    Kernel::remove_trash(faces_1st_level, faces_trash);
    //***************************************************************

    //***************************************************************
    // checks for duplicate HashCodes in faces_1st_level' TopoDS_Faces
    Kernel::check_duplicate_faces(faces_1st_level);
    // In faces_1st_level there are no faces with same HashCode anymore
    //***************************************************************

    //***************************************************************
    // create network of cface - halfedge - cface
    Kernel::update_face_adjacencies(faces_1st_level, fuse);
    //***************************************************************

    //***************************************************************
    // identify faces that don't have adjacent faces on all edges
    Kernel::identify_hanging_faces(faces_1st_level);
    //***************************************************************

    //***************************************************************
    // get copy of hanging faces
    if (!stl)
        Kernel::nonoffset_hanging_faces(faces_1st_level, faces_non_sb);
    //***************************************************************

    //***************************************************************
    // deletes trash faces from faces_1st_level and moves them to faces_trash
    // adjacence info of a cface' edges are removed, if adjacent cface is trash
    Kernel::remove_trash_and_face_adjacency(faces_1st_level);
    //***************************************************************

    //***************************************************************
    // identifies offset faces that are not in contact with shell anymore
    Kernel::identify_decoupled_offset_faces(faces_1st_level);
    //***************************************************************

    //***************************************************************
    Kernel::remove_trash_and_face_adjacency(faces_1st_level);
    //***************************************************************

    //***************************************************************
    // find components in building graph aka spaces
    Kernel::check_adjacency_self_reference(faces_1st_level);
    K.find_spaces_normals_unknown(faces_1st_level, spaces, space_id_counter, face_id_counter);
    //***************************************************************

    //***************************************************************
    // get copy of inner faces
    if (!stl)
        Kernel::nonoffset_inner_or_coplanar_faces(faces_1st_level, faces_non_sb);
    //***************************************************************

    //***************************************************************
    Kernel::remove_trash_and_face_adjacency(faces_1st_level);
    //***************************************************************

    //***************************************************************
    Kernel::check_redundant_faces(faces_1st_level, faces_non_sb);
    Kernel::check_fixed_normal(faces_1st_level);
    Kernel::check_closed_space_edge_id(spaces);
    //***************************************************************
}

// Spaces
bool Graph::process_spaces(Kernel &K) {

    std::cout << "\n### Spaces" << std::endl;

    //***************************************************************
    // get outer hull of building
    Kernel::identify_facade_space(spaces);
    //***************************************************************

    //***************************************************************
    // split spaces by IfcSpaces from original model
    if (space_split) {
        K.split_spaces_by_IfcSpaces(spaces, faces_1st_level, ifc_faces, products, ifcspaces, face_id_counter);
        for (const auto &p: products)
            if (p.IsIfcClass("IfcVirtualElement"))
                virtual_products.insert(p.ifcproduct);
    }
    Kernel::check_fixed_normal(faces_1st_level);
    // all normals must be known at this point
    //***************************************************************

    //***************************************************************
    // calculate solid classifiers of spaces used in point in shell operations
    Kernel::initialize_space_solid_classifiers(spaces);
    //***************************************************************

    //***************************************************************
    // remove invalid spaces, e.g. railings hanging in the air
    Kernel::remove_invalid_spaces(spaces, faces_1st_level);
    Kernel::nonoffset_inner_or_coplanar_faces(faces_1st_level, faces_non_sb, true);
    Kernel::remove_trash_and_face_adjacency(faces_1st_level);
    //***************************************************************

    //***************************************************************
    // remove window and door faces occurring as inner or coplanar faces.
    // this is done for bad situations where opening integration into wall
    // did not work and faces prevent wall faces from unifying
    Kernel::remove_inner_window_and_door_faces(faces_non_sb);
    //***************************************************************

    //***************************************************************
    // remove seam edges from faces
    Kernel::remove_seam_edges_from_faces(faces_1st_level);
    //***************************************************************

    //***************************************************************
    // unify first level faces_1st_level with adjacent neighbours
    // Attention: before calling this function, all trash faces must be removed)
    if (!stl)
        Kernel::unify_cFaces(faces_1st_level, true);
    else
        Kernel::unify_cFaces(faces_1st_level, false);
    Kernel::update_adjacence_map_after_unification(faces_1st_level);
    //***************************************************************

    //***************************************************************
    // clear cface's adjacence maps
    //K.clear_cface_maps(faces_1st_level); needed for shading
    //***************************************************************

    //***************************************************************
    // remove the trashed faces from spaces before deleting the cfaces
    // deletes trash faces from faces_1st_level and moves them to faces_trash
    Kernel::remove_trash_from_spaces(faces_1st_level, spaces);
    Kernel::remove_trash(faces_1st_level, faces_trash);
    Kernel::check_adjacency_null(faces_1st_level);
    //***************************************************************

    //***************************************************************
    // In some cases, when opening face can't be unified by wall face, it's missing a parentface
    // Introduces new origfaces and cfaces. Origface is pointing to RelProduct of the face with missing parent
    if (!stl)
        Kernel::create_missing_parent_faces(faces_1st_level, ifc_faces, face_id_counter);
    //***************************************************************

    //***************************************************************
    // merge colinear edges in faces (somehow it's not done by unify_cFaces())
    Kernel::reduce_edges_cFaces(faces_1st_level);
    Kernel::check_duplicate_vertices(faces_1st_level);
    //***************************************************************

    //***************************************************************
    // process space shell for export as IfcSpace
    Kernel::process_space_shells_for_export(spaces, false);
    //***************************************************************

    //***************************************************************
    // write stl
    if (stl) {
        std::string stl_mode = "productWiseInOne";
        std::map<std::string, TopoDS_Shape> additional_shapes = K.add_additional_shapes(model, settings, include_entities, cfd_entities, true, true);
        // K.triangulate_cfaces(faces_1st_level);
        if (!K.write_spaces_to_stl(faces_1st_level, spaces, input, output, "stl/", stl_mode, true, false, additional_shapes, fuzzy_tol))
            return false;
        return true;
    }
    //***************************************************************

    //***************************************************************
    // prepare first level space boundaries for ifc export
    if (first_level_only) {

        //***************************************************************
        // calculate ray origins on faces
        K.calculate_points_on_faces(faces_1st_level);
        //***************************************************************

        //***************************************************************
        // set enumeration types of sb
        K.identify_sb_types_ray_tracing_first_level(faces_1st_level, bounds_min, max_transmission_length);
        //***************************************************************

        //***************************************************************
        Kernel::check_parent_faces(faces_1st_level);
        //***************************************************************

        return true;
    }
    //***************************************************************

    //***************************************************************
    // split facade faces by IfcSite
    K.site(model, settings, faces_1st_level, spaces, face_id_counter);
    Kernel::remove_trash_from_spaces(faces_1st_level, spaces);
    Kernel::remove_trash(faces_1st_level, faces_trash);
    //***************************************************************

    //***************************************************************
    Kernel::check_parent_faces(faces_1st_level);
    Kernel::check_fixed_normal(faces_1st_level);
    //***************************************************************

    return !spaces.empty();
}

// Second level faces_1st_level
bool Graph::calc_faces_second_level(Kernel &K) {

    std::cout << "\n### Boundaries second level" << std::endl;

    //***************************************************************
    // project non boundary faces onto boundary faces and split them by clipping
    cfaceSetMap M_planar, M_curved; // faces with planar/linear faces or edges and faces with curved faces or edges
    Kernel::find_clipping_pairs(M_planar, M_curved, faces_1st_level, faces_non_sb, max_transmission_length); // non SB faces clipping faces_1st_level, virtual faces are excluded
    Kernel::find_clipping_pairs_2(M_planar, M_curved, faces_1st_level, max_transmission_length, 95 * M_PI / 180); //  SB faces mutual clipping, virtual faces are excluded
    K.clip_faces(faces_1st_level, faces_2nd_level, M_planar, M_curved, max_transmission_length, 1.0e-5, face_id_counter);
    Kernel::correct_face_normals(faces_2nd_level);
    M_curved.clear();
    M_planar.clear();
    //***************************************************************

    //***************************************************************
    // Check some face properties
    Kernel::check_faces(faces_2nd_level, false, 3 * fuzzy_tol, 1.5 * fuzzy_tol);
    Kernel::check_duplicate_vertices(faces_2nd_level);
    Kernel::check_fixed_normal(faces_2nd_level);
    //***************************************************************

    //***************************************************************
    // Fill second level face list in spaces
    Kernel::spaces_introduce_second_level_cfaces(spaces);
    //***************************************************************

    //***************************************************************
    // calculate ray origins on faces
    K.calculate_points_on_faces(faces_2nd_level);
    //***************************************************************

    //***************************************************************
    // set enumeration types of sb
    K.identify_sb_types_ray_tracing(faces_2nd_level, bounds_min, max_transmission_length, 179 * M_PI / 180);
    //K.initialize_space_intersectors(spaces); // intersector is based on cfaces of space, not on shell. So must be up-to-date
    //K.identify_sb_types_ray(faces_2nd_level, spaces, bounds_min, max_transmission_length, 179 * M_PI / 180, true, true);
    //***************************************************************

    //***************************************************************
    // consider internal and coplanar faces for material list
    K.complete_material_list_with_inner_faces(faces_non_sb, faces_2nd_level);
    //***************************************************************

    //***************************************************************
    // unify faces_1st_level with adjacent neighbours
    // Attention: before calling this function, all trash faces must be removed
    Kernel::unify_cFaces(faces_2nd_level, true);
    //***************************************************************

    //***************************************************************
    // flag non-2a opening faces as trash
    Kernel::trash_2b_opening_faces(faces_2nd_level);
    //***************************************************************

    //***************************************************************
    // remove the trashed faces from spaces before deleting the cfaces
    // deletes trash faces from faces_1st_level and moves them to faces_trash
    Kernel::remove_trash_from_spaces(faces_2nd_level, spaces);
    Kernel::remove_trash(faces_2nd_level, faces_trash);
    //***************************************************************

    //***************************************************************
    // merge colinear edges in faces (somehow it's not done by unify_cFaces())
    Kernel::reduce_edges_cFaces(faces_2nd_level);
    //***************************************************************

    //***************************************************************
    if (remove_holes_by_face_split) {
        Kernel::check_holes(faces_2nd_level);
        Kernel::remove_holes(faces_2nd_level, face_id_counter);
        Kernel::remove_corresponding_trash(faces_2nd_level);
        Kernel::remove_trash_from_spaces(faces_2nd_level, spaces);
        Kernel::remove_trash(faces_2nd_level, faces_trash);
    }
    //***************************************************************

    //***************************************************************
    if (decompose_concave_polygons && remove_holes_by_face_split) {
        Kernel::check_convex(faces_2nd_level);
        Kernel::decompose_concave_polygons(faces_2nd_level, face_id_counter);
        Kernel::remove_corresponding_trash(faces_2nd_level);
        Kernel::remove_trash_from_spaces(faces_2nd_level, spaces);
        Kernel::remove_trash(faces_2nd_level, faces_trash);
    }
    //***************************************************************

    //***************************************************************
    if (simplify_fenestration_faces) {
        Kernel::simplify_fenestration_faces(faces_2nd_level, face_id_counter);
        Kernel::remove_corresponding_trash(faces_2nd_level);
        Kernel::remove_trash_from_spaces(faces_2nd_level, spaces);
        Kernel::remove_trash(faces_2nd_level, faces_trash);
    }
    //***************************************************************

    //***************************************************************
    Kernel::check_pointers(faces_2nd_level, true);
    K.check_loops(faces_2nd_level, 1.0e-5);
    Kernel::check_parent_faces(faces_1st_level);
    Kernel::check_parent_faces(faces_2nd_level);
    Kernel::check_corresponding_face_pairs(faces_2nd_level);
    //Kernel::check_corresponding_face_pairs(faces_2nd_level);
    Kernel::check_space_behind(faces_2nd_level);
    Kernel::check_holes(faces_2nd_level);
    Kernel::check_convex(faces_2nd_level);
    Kernel::check_duplicate_vertices(faces_2nd_level);
    //***************************************************************

    //***************************************************************
    // set cfaces to first level SBs
    Kernel::check_super_and_subfaces(faces_1st_level);
    Kernel::check_superface(faces_2nd_level);
    Kernel::change_cface_sbtype(faces_1st_level, faces_2nd_level);
    Kernel::check_enum_intext(faces_1st_level);
    Kernel::check_enum_intext(faces_2nd_level);
    //***************************************************************

    //***************************************************************
    if (calculate_shadings)
        K.calculate_shading(faces_2nd_level, spaces, shadings, model, settings);
    //***************************************************************

    return !faces_2nd_level.empty();
}

// Model
bool Graph::process_model(Kernel &K) {

    std::cout << "\n### Model" << std::endl;

    //***************************************************************
    // link guids and psets of old spaces to new ones
    K.link_old_IfcSpace_data(model, settings, spaces, ifcspaces);
    //***************************************************************

    //***************************************************************
    // remove original IfcSpaces and IfcRelSpaceBoundaries from model
    std::set<std::string> remove_classes = {"IfcSpace", "IfcRelSpaceBoundary", "IfcExternalSpatialElement", "IfcVirtualElement"}; // "IfcSpaceType",
    Kernel::remove_entities_from_model(model, remove_classes);
    //***************************************************************

    //***************************************************************
    // add IfcVirtualElements to model
    Kernel::add_virtual_elements_to_model(model, virtual_products);
    //***************************************************************

    //***************************************************************
    // removes invalid pointer to entities in model
    Kernel::check_ifcproduct_pointers(model, spaces);
    //***************************************************************

    //***************************************************************
    // add IfcSpaces and IfcRelSpaceBoundaries to model
    if (!K.add_entities_to_model(model, spaces, virtual_products, shadings))
        return false;
    //***************************************************************

    //***************************************************************
    // remove unused relationship instances (e.g. psets)
    Kernel::remove_relationships_from_model(model);
    //***************************************************************

    //***************************************************************
    // correct isExternal attribut for IfcWalls
    K.revise_is_external_attribute(model, spaces);
    //***************************************************************

    //***************************************************************
    // add "SpaceBoundary2ndLevelAddOnView" tag to header description
    if (K.ifcSchema == IFC2X3)
        Kernel::modify_view_definition(model);
    //***************************************************************

    //***************************************************************
    // export the current model into new ifc file
    Kernel::write_ifc_file(model, output);
    //***************************************************************

    return true;
}