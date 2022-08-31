// Copyright 2022 Eric Fichter
#include "Clip.h"

Clip::Clip(std::string _input,
           std::string _output,
           unsigned int _num_threads,
           double _max_transmission_length,
           bool _remove_holes_by_face_split,
           bool _decompose_concave_polygons,
           bool _simplify_fenestration_faces,
           bool _calculate_shadings,
           bool _use_ifcopeningelelements_for_virtual_boundaries) :
        input(std::move(_input)),
        output(std::move(_output)),
        num_threads(_num_threads),
        max_transmission_length(_max_transmission_length),
        remove_holes_by_face_split(_remove_holes_by_face_split),
        decompose_concave_polygons(_decompose_concave_polygons),
        simplify_fenestration_faces(_simplify_fenestration_faces),
        calculate_shadings(_calculate_shadings),
        use_ifcopeningelelements_for_virtual_boundaries(_use_ifcopeningelelements_for_virtual_boundaries) {

    products.clear();
    element_ifc_faces.clear();
    ifcspace_faces.clear();
    element_faces.clear();
    faces_2nd_level.clear();
    faces_trash.clear();
    spaces.clear();
    virtual_products.clear();
    shadings.clear();
    ifcspaces.clear();
    ifcspace_ifc_faces.clear();
    space_map.clear();
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

    fuzzy_tol = 0.0001;
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

bool Clip::run() {

    Kernel K(num_threads);

    std::unordered_map<IfcUtil::IfcBaseClass *, TopoDS_Shape> opening_shapes;
    if (!prepare_products(K, opening_shapes))
        return false;

    if (!prepare_building_element_faces(K, opening_shapes))
        return false;

    if (!process_spaces(K))
        return false;

    if (!calc_faces_second_level(K))
        return false;

    if (!process_model(K))
        return false;

    // check results written to model
    IfcCheck ifccheck(model.get(), num_threads);

#ifdef VISUALIZATION
    // Viewer::visualize_products(products);
    // Viewer::visualize_orig_faces(element_ifc_faces);
    // Viewer::visualize_cFaces(ifcspace_faces);
    // Viewer::visualize_spaces(spaces, ifcspace_faces, true);
    // Viewer::visualize_cFaces_as_space_boundaries(faces_2nd_level);
#endif

    return true;
}

// Products
bool Clip::prepare_products(Kernel &K, std::unordered_map<IfcUtil::IfcBaseClass *, TopoDS_Shape> &opening_shapes) {

    std::cout << "\n### Products" << std::endl;

    //***************************************************************
    // Parsing Ifc file to model
    if (!K.read_ifc_file(input, model))
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
    K.find_relevant_products(model, include_entities, non_void_filling_products, integrate_openings_into_walls, use_ifcopeningelelements_for_virtual_boundaries);
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
    // Add entities to products that are composed of sub-entities, e.g. IfcCurtainWall. The shape can be simplified
    K.add_decomposed_entities_to_products(model, settings, products, include_entities);
    //***************************************************************

    if (integrate_openings_into_walls) {
        //***************************************************************
        // For all created products (non void filling products) find non-filled openings that will be cut from product
        std::set<std::string> opening_guids = K.find_openings_in_products(products);
        // Create shapes of openings
        opening_shapes = K.generate_opening_shapes(model, settings, opening_guids);
        // Cut opening shapes from products
        K.subtract_openings_from_products(products, opening_shapes);
        //***************************************************************


        //***************************************************************
        // Clear container, part I
        for (auto &product: products) product.IfcOpeningElements.clear();
        non_void_filling_products.clear(); // K.free_container(normal_products);
        opening_guids.clear(); // K.free_container(opening_guids);
        //***************************************************************
    }

    // add site for external earth
    boost::shared_ptr<IfcEntityList> IfcEntityList = nullptr;
    std::set<std::string> site_guid;
    gp_XYZ fake;
    try { IfcEntityList = model->instances_by_type("IfcSite"); }
    catch (...) {}

    if (IfcEntityList != nullptr)
        for (auto E: *IfcEntityList) {
            std::string s = E->data().getArgument(0)->toString();
            site_guid.insert(s.substr(1, s.size() - 2));
        }

    if (!site_guid.empty()) K.generate_shapes_from_ifc_guids(model, settings, products, site_guid, fake, fake);

    return !products.empty();
}

bool Clip::prepare_building_element_faces(Kernel &K, std::unordered_map<IfcUtil::IfcBaseClass *, TopoDS_Shape> &opening_shapes) {

    std::cout << "\n### Building Element Faces" << std::endl;

    //***************************************************************
    // Save all faces of product shapes
    K.collect_original_faces_clip(products, element_ifc_faces);
    //***************************************************************

    //***************************************************************
    // Get virtual faces by space x space clipping
    //K.get_virtual_faces(model, ifcspaces, element_ifc_faces, products);
    //Viewer::visualize_orig_faces_virtual_product(element_ifc_faces);
    //***************************************************************

    //***************************************************************
    // Create opening faces in walls
    if (integrate_openings_into_walls)
        K.create_void_filling_cuts_in_products(products, element_ifc_faces, opening_shapes, false, model, settings);
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
    for (auto &ifcface: element_ifc_faces)
        element_faces.emplace_back(TopoDS::Face(K.shape_copy(ifcface.face)), &ifcface, 0); // copy, otherwise occ triangulation created segfault sometimes in parallel
    //***************************************************************

    std::cout << std::flush;

    return !element_ifc_faces.empty();
}

// Spaces
bool Clip::process_spaces(Kernel &K) {

    std::cout << "\n### Spaces" << std::endl;

    //***************************************************************
    // Create shapes of IfcSpaces
    ifcspaces = K.create_shapes_of_IfcSpaces(model, settings);
    //***************************************************************

    //***************************************************************
    // Store IfcSpaces as products and collect their faces as oFaces and cFaces. Create spaces from IfcSpaces.
    K.collect_spaces_and_ifcfaces(products, ifcspaces, ifcspace_ifc_faces, ifcspace_faces, spaces, space_id_counter, space_map, face_id_counter);
    //***************************************************************

    if (ifcspace_faces.empty()) return false;

    //***************************************************************
    // calculate solid classifiers of spaces used in point in shell operations
    K.initialize_space_solid_classifiers(spaces);
    //***************************************************************

    //***************************************************************
    // split facade faces by IfcSite
    K.site(model, settings, ifcspace_faces, spaces, face_id_counter);
    K.remove_trash_from_spaces(ifcspace_faces, spaces);
    K.remove_trash(ifcspace_faces, faces_trash);
    //***************************************************************

    //***************************************************************
    K.check_parent_faces(ifcspace_faces);
    //***************************************************************

    return !spaces.empty();
}

// Second level ifcspace_faces
bool Clip::calc_faces_second_level(Kernel &K) {

    std::cout << "\n### Boundaries second level" << std::endl;

    //***************************************************************
    // project non boundary faces onto boundary faces and split them by clipping
    std::map<cFace *, std::set<cFace * >> M_planar, M_curved; // faces with planar/linear faces or edges and faces with curved faces or edges
    K.find_clipping_pairs_space_approach(M_planar, M_curved, ifcspace_faces, max_transmission_length, 95 * M_PI / 180); //  SB faces mutual clipping, virtual faces are excluded
    K.find_clipping_pairs_all(M_planar, M_curved, ifcspace_faces, element_faces, max_transmission_length);
    K.clip_faces(ifcspace_faces, faces_2nd_level, M_planar, M_curved, max_transmission_length, 1.0e-5, face_id_counter);
    K.correct_face_normals(faces_2nd_level);
    M_curved.clear();
    M_planar.clear();
    //***************************************************************

    //***************************************************************
    // Check some face properties
    K.check_faces(faces_2nd_level, false, 3 * fuzzy_tol, 1.5 * fuzzy_tol);
    K.check_duplicate_vertices(faces_2nd_level);
    K.check_fixed_normal(faces_2nd_level);
    //***************************************************************

    //***************************************************************
    // Fill second level face list in spaces
    K.spaces_introduce_second_level_cfaces(spaces);
    //***************************************************************

    //***************************************************************
    // calculate ray origins on faces
    K.calculate_points_on_faces(faces_2nd_level);
    //***************************************************************

    //***************************************************************
    // set enumeration types of sb
    K.identify_sb_types_ray_tracing_space_approach(faces_2nd_level, bounds_min, max_transmission_length, 179 * M_PI / 180);
    //***************************************************************

    //***************************************************************
    // consider internal and coplanar faces for material list
    K.perform_ray_tracing_space_approach_elements(element_faces, faces_2nd_level, max_transmission_length);
    //***************************************************************

    //***************************************************************
    // Pseudo facade space
    std::set<cFace *> spaceFaces;
    for (auto &cface: faces_2nd_level)
        if (cface.InternalOrExternal() != SB_IE_INTERNAL && cface.SBType() == SB_TYPE_2A) {
            if (cface.Corresponding() != nullptr) spaceFaces.insert(cface.Corresponding());
            else std::cerr << "[Warning] No corresponding for external 2a face!" << std::endl;
        }
    spaces.emplace_back(space_id_counter, spaceFaces, false);
    for (auto &spaceFace: spaceFaces) spaceFace->SetSpace(&(spaces.back()));
    //space_map[&spaces.back()] = &spaces.back();
    space_id_counter++;
    spaces.back().is_facade = true;
    //***************************************************************

    //***************************************************************
    for (auto &cface: faces_2nd_level) {

        // space behind
        if (cface.InternalOrExternal() != SB_IE_INTERNAL && cface.SBType() == SB_TYPE_2A) {
            if (cface.Corresponding() != nullptr) cface.SetSpaceBehind(cface.Corresponding()->RelSpace());
            else std::cerr << "[Warning] No corresponding for external 2a face!" << std::endl;
        }
        // external_earth for physical 2a
        if (cface.InternalOrExternal() == SB_IE_EXTERNAL && cface.PhysicalOrVirtual() == SB_PV_PHYSICAL && cface.Corresponding() != nullptr)
            if (cface.Corresponding()->RelSpace()->is_facade)
                if (cface.Corresponding()->PointOnFace().Z() - 0.05 < bounds_min.Z() && cface.Corresponding()->FixedFaceNormal().Z() < -0.95) cface.SetInternalOrExternal(SB_IE_EXTERNAL_EARTH);
    }
    //***************************************************************

    //***************************************************************
    // unify ifcspace_faces with adjacent neighbours
    // Attention: before calling this function, all trash faces must be removed
    K.unify_cFaces(faces_2nd_level, true);
    //***************************************************************

    //***************************************************************
    // flag non-2a opening faces as trash
    K.trash_2b_opening_faces(faces_2nd_level);
    //***************************************************************

    //***************************************************************
    // remove the trashed faces from spaces before deleting the cfaces
    // deletes trash faces from ifcspace_faces and moves them to faces_trash
    K.remove_trash_from_spaces(faces_2nd_level, spaces);
    K.remove_trash(faces_2nd_level, faces_trash);
    //***************************************************************

    //***************************************************************
    // merge colinear edges in faces (somehow it's not done by unify_cFaces())
    K.reduce_edges_cFaces(faces_2nd_level);
    //***************************************************************

    //***************************************************************
    K.create_virtual_elements_space_approach(faces_2nd_level, products, ifcspace_ifc_faces);
    //***************************************************************

    //***************************************************************
    for (const auto &p: products)
        if (p.IsIfcClass("IfcVirtualElement"))
            virtual_products.insert(p.ifcproduct);
    //***************************************************************

    //***************************************************************
    if (remove_holes_by_face_split) {
        K.check_holes(faces_2nd_level);
        K.remove_holes(faces_2nd_level, face_id_counter);
        K.remove_corresponding_trash(faces_2nd_level);
        K.remove_trash_from_spaces(faces_2nd_level, spaces);
        K.remove_trash(faces_2nd_level, faces_trash);
    }
    //***************************************************************

    //***************************************************************
    if (decompose_concave_polygons && remove_holes_by_face_split) {
        K.check_convex(faces_2nd_level);
        K.decompose_concave_polygons(faces_2nd_level, face_id_counter);
        K.remove_corresponding_trash(faces_2nd_level);
        K.remove_trash_from_spaces(faces_2nd_level, spaces);
        K.remove_trash(faces_2nd_level, faces_trash);

        K.decompose_concave_polygons_triangulation(faces_2nd_level, face_id_counter);
        K.remove_corresponding_trash(faces_2nd_level);
        K.remove_trash_from_spaces(faces_2nd_level, spaces);
        K.remove_trash(faces_2nd_level, faces_trash);
    }
    //***************************************************************

    //***************************************************************
    if (simplify_fenestration_faces) {
        K.simplify_fenestration_faces(faces_2nd_level, face_id_counter);
        K.remove_corresponding_trash(faces_2nd_level);
        K.remove_trash_from_spaces(faces_2nd_level, spaces);
        K.remove_trash(faces_2nd_level, faces_trash);
    }
    //***************************************************************

    //***************************************************************
    K.check_pointers(faces_2nd_level, true);
    K.check_loops(faces_2nd_level, 1.0e-5);
    K.check_parent_faces(ifcspace_faces);
    K.check_parent_faces(faces_2nd_level);
    K.check_corresponding_face_pairs_space_approach(faces_2nd_level);
    K.check_space_behind(faces_2nd_level);
    K.check_holes(faces_2nd_level);
    K.check_convex(faces_2nd_level);
    K.check_duplicate_vertices(faces_2nd_level);
    //***************************************************************

    //***************************************************************
    // set cfaces to first level SBs
    K.check_super_and_subfaces(ifcspace_faces);
    K.check_superface(faces_2nd_level);
    K.change_cface_sbtype(ifcspace_faces, faces_2nd_level);
    K.check_enum_intext(ifcspace_faces);
    K.check_enum_intext(faces_2nd_level);
    //***************************************************************

    //***************************************************************
    if (calculate_shadings)
        K.calculate_shading(faces_2nd_level, spaces, shadings, model, settings);
    //***************************************************************

    return !faces_2nd_level.empty();
}

// Model
bool Clip::process_model(Kernel &K) {

    std::cout << "\n### Model" << std::endl;

    //***************************************************************
    // remove original IfcSpaces and IfcRelSpaceBoundaries from model
    std::set<std::string> remove_classes = {"IfcRelSpaceBoundary", "IfcExternalSpatialElement", "IfcVirtualElement"};
    K.remove_entities_from_model(model, remove_classes);
    //***************************************************************

    //***************************************************************
    // add IfcVirtualElements to model
    K.add_virtual_elements_to_model(model, virtual_products);
    //***************************************************************

    //***************************************************************
    for (auto &space: spaces)
        space.ClearFirstLevel();
    //***************************************************************

    //***************************************************************
    // removes invalid pointer to entities in model
    K.check_ifcproduct_pointers(model, spaces);
    //***************************************************************

    //***************************************************************
    // add IfcSpaces and IfcRelSpaceBoundaries to model
    if (!K.add_entities_to_model_space_approach(model, spaces, virtual_products, shadings, space_map))
        return false;
    //***************************************************************

    //***************************************************************
    // add "SpaceBoundary2ndLevelAddOnView" tag to header description
    if (K.ifcSchema == IFC2X3)
        K.modify_view_definition(model);
    //***************************************************************

    //***************************************************************
    // export the current model into new ifc file
    K.write_ifc_file(model, output);
    //***************************************************************

    return true;
}