// Copyright 2022 Eric Fichter
#include "Kernel.h"

// Concerns mostly Product and oFace objects only. Longer functions with deeper geometric processing.

bool Kernel::check_ifc_schema(std::unique_ptr<IfcParse::IfcFile> &model) {

    auto start = std::chrono::high_resolution_clock::now();

    if (model->schema()->name() == "IFC2X3")
        ifcSchema = IFC2X3;
    else if (model->schema()->name() == "IFC4")
        ifcSchema = IFC4;
    else if (model->schema()->name() == "IFC4X1")
        ifcSchema = IFC4X1;
    else if (model->schema()->name() == "IFC4X2")
        ifcSchema = IFC4X2;
    else if (model->schema()->name() == "IFC4X3_RC1")
        ifcSchema = IFC4X3_RC1;
    else {
        std::cerr << "[Error] Schema " << model->schema()->name() << " not implemented yet! " << std::endl;
        return false;
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Check ifc schema", model->schema()->name());

    return true;
}

bool Kernel::check_ifc_context(IfcGeom::Element<real_t> *geom_object) {
    //Reference 3D representation that is not part of the Body representation. This is used, e.g., for opening geometries, if there are to be excluded from an implicit Boolean operation.
    //Body      3D Body representation, e.g. as wireframe, surface, or solid model, of an element
    std::string t = geom_object->context();
    boost::algorithm::to_lower(t);
    if (t != "body" && t != "reference") {
        std::string s = geom_object->type() + " " + geom_object->guid();
        std::cerr << "[Warning] Wrong context: " << geom_object->context() << ". " << s << std::endl;
        return false;
    } else return true;
}

bool Kernel::generate_shapes_from_ifc_guids(std::unique_ptr<IfcParse::IfcFile> &model, IfcGeom::IteratorSettings settings, std::list<Product> &products, const std::set<std::string> &include_guids, gp_XYZ &bounds_min, gp_XYZ &bounds_max) const {

    auto start = std::chrono::high_resolution_clock::now();

    IfcGeom::attribute_filter attribute_filter;
    attribute_filter.include = true;
    attribute_filter.traverse = false;
    attribute_filter.attribute_name = "GlobalId";
    attribute_filter.populate(include_guids);
    std::vector<IfcGeom::filter_t> filter_funcs;
    filter_funcs.emplace_back(boost::ref(attribute_filter));

    IfcGeom::Iterator<real_t> geom_iterator(settings, model.get(), filter_funcs, num_threads);
    std::set<std::string> contexts;

    std::string unit_name;

    if (!geom_iterator.initialize()) {
        std::cerr << "No geometrical entities found. Initialization failed." << std::endl;
        return false;
    } else {
        unit_name = geom_iterator.unit_name();
        std::for_each(unit_name.begin(), unit_name.end(), [](char &c) { c = ::tolower(c); });
    }

    do {
        IfcGeom::Element<real_t> *geom_object = geom_iterator.get();
        if (!check_ifc_context(geom_object)) continue;
        products.emplace_back(geom_object->product(), geom_object->guid(), geom_object_to_shape(geom_object));
        contexts.insert(geom_object->context());
    } while (geom_iterator.next());

    // store aabb values of building
    bounds_min = geom_iterator.bounds_min();
    bounds_max = geom_iterator.bounds_max();

    std::string context;
    for (const auto &c: contexts) context += " " + c;

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Create shapes", "unit: " + unit_name + ", " + std::to_string(geom_iterator.unit_magnitude()) + ", context:" + context + ", " + std::to_string(products.size()));

    return true;
}

void Kernel::bad_openings_to_box_faces(std::list<oFace> &orig_faces, std::unique_ptr<IfcParse::IfcFile> &model, IfcGeom::IteratorSettings settings, std::list<Product> &products, std::unordered_map<std::string, Product *> &todo) const {

    std::cout << "[Info] Create shapes of window and door instances that could not be processed, yet.\n";

    std::set<std::string> include_guids;
    std::unordered_map<std::string, TopoDS_Shape> aabbs;
    unsigned int shell_id = 99999;

    for (auto &t: todo) {
        include_guids.insert(t.first);
        Product *wall = t.second;
        aabbs.insert({t.first, aabb_to_shape(aabb(wall->shape, 0))});
    }

    IfcGeom::attribute_filter attribute_filter;
    attribute_filter.include = true;
    attribute_filter.traverse = false;
    attribute_filter.attribute_name = "GlobalId";
    attribute_filter.populate(include_guids);
    std::vector<IfcGeom::filter_t> filter_funcs;
    filter_funcs.emplace_back(boost::ref(attribute_filter));

    IfcGeom::Iterator<real_t> geom_iterator(settings, model.get(), filter_funcs, num_threads);
    std::set<std::string> contexts;

    if (!geom_iterator.initialize()) {
        std::cerr << "[Warning] No geometrical entities found. Initialization failed." << std::endl;
        return;
    }

    do {
        IfcGeom::Element<real_t> *geom_object = geom_iterator.get();
        if (!check_ifc_context(geom_object)) continue;

        TopoDS_Shape S = aabb_to_shape(aabb(geom_object_to_shape(geom_object), 0));
        //TopoDS_Shape S = best_fitting_bbox(geom_object_to_shape(geom_object));
        if (S.IsNull()) continue;

        TopoDS_Shape comm = BRepAlgoAPI_Common(aabbs[geom_object->guid()], S).Shape();
        if (!comm.IsNull()) {
            double V_c = volume(comm);
            double V_s = volume(S);
            double q = fabs(V_s - V_c) / V_s;
            if (q < 0.18) {
                S = comm; // take the common of the aabb of the wall and the box of the window, if the volume difference of window box and common does not exceed certain a value (high value means that window box would be downsized too much, maybe creating holes)
                std::cout << "[Info] Common passed volume criteria. " << geom_object->guid() << " " << round_double_three_digits(q) << "\n";
            } else
                std::cout << "[Info] Common did not pass volume criteria. " << geom_object->guid() << " " << round_double_three_digits(q) << "\n";
        } else
            std::cout << "[Info] Common is null.\n";

        //products.emplace_back(geom_object->product(), geom_object->guid(), S);
        //products.back().valid = true;
        //collect_original_faces_openings(products.back(), orig_faces);

        // add to ifcfaces
        auto shells = Topo(S).shells();
        auto wall = todo[geom_object->guid()];
        auto window = geom_object->product();

        for (auto &shell: shells) {
            for (auto &face: Topo(shell).faces()) {
                orig_faces.emplace_back(TopoDS::Face(face), wall, shell_id);
                wall->orig_faces.push_back(&orig_faces.back());
                orig_faces.back().SetParentID(1);
                orig_faces.back().SetOpening(window);
                orig_faces.back().SetNormalStatus(FACE_NORMAL_KNOWN);
            }
            shell_id++;
        }

    } while (geom_iterator.next());

}

std::list<Product *> Kernel::check_shells(std::list<Product> &products) {

    auto start = std::chrono::high_resolution_clock::now();

    std::list<Product *> faulty_products;

#ifdef PARALLEL_PROCESSING
    tbb::parallel_for_each(products.begin(), products.end(), [&](Product &product) { if (!product.CheckShells()) faulty_products.push_back(&product); });
#else
    for (auto &product: products)
        if (!product.CheckShells()) faulty_products.push_back(&product);
#endif

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Check for broken shells", std::to_string(products.size()));

    if (!faulty_products.empty()) {
        std::cerr << "\t[Warning] Products with faulty shape:" << std::endl;
        for (const auto &p: faulty_products)
            std::cerr << "\t" << p->ifcproduct->data().toString() << std::endl;
    }

    return faulty_products;
}

std::list<Product *> Kernel::check_shells(std::list<Product *> &products) {

    std::cout << "Check for broken shells ";

    auto start = std::chrono::high_resolution_clock::now();

    std::list<Product *> faulty_products;

#ifdef PARALLEL_PROCESSING
    tbb::parallel_for_each(products.begin(), products.end(), [&](Product *product) { if (!product->CheckShells()) faulty_products.push_back(product); });
#else
    for (auto &product: products)
        if (!product->CheckShells()) faulty_products.push_back(product);
#endif

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "(" << products.size() << ") \t... \tElapsed time: " << elapsed.count() << " s\n";
    std::cout << print_time(elapsed.count(), "Check for broken shells", std::to_string(products.size()));

    if (!faulty_products.empty()) {
        std::cerr << "\t[Warning] Products with faulty shape:" << std::endl;
        for (const auto &p: faulty_products)
            std::cerr << "\t" << p->ifcproduct->data().toString() << std::endl;
    }

    return faulty_products;
}

void Kernel::find_relevant_products(std::unique_ptr<IfcParse::IfcFile> &model, const std::set<std::string> &include_entities, std::set<std::string> &non_void_filling_products, bool integrate_openings_into_walls, bool use_ifcopeningelelements_for_virtual_boundaries) {

    auto start = std::chrono::high_resolution_clock::now();

    if (ifcSchema == IFC2X3)
        find_relevant_products_worker<Ifc2x3>(model, include_entities, non_void_filling_products, integrate_openings_into_walls, use_ifcopeningelelements_for_virtual_boundaries);
    else if (ifcSchema == IFC4)
        find_relevant_products_worker<Ifc4>(model, include_entities, non_void_filling_products, integrate_openings_into_walls, use_ifcopeningelelements_for_virtual_boundaries);
    else if (ifcSchema == IFC4X1)
        find_relevant_products_worker<Ifc4x1>(model, include_entities, non_void_filling_products, integrate_openings_into_walls, use_ifcopeningelelements_for_virtual_boundaries);
    else if (ifcSchema == IFC4X2)
        find_relevant_products_worker<Ifc4x2>(model, include_entities, non_void_filling_products, integrate_openings_into_walls, use_ifcopeningelelements_for_virtual_boundaries);
    else if (ifcSchema == IFC4X3_RC1)
        find_relevant_products_worker<Ifc4x3_rc1>(model, include_entities, non_void_filling_products, integrate_openings_into_walls, use_ifcopeningelelements_for_virtual_boundaries);
    else
        std::cerr << "[Error] Schema not implemented yet! " << model->schema()->name() << std::endl;

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Identify relevant products", std::to_string(non_void_filling_products.size()));

}

template<typename Schema>
void Kernel::find_relevant_products_worker(std::unique_ptr<IfcParse::IfcFile> &model, const std::set<std::string> &include_entities, std::set<std::string> &non_void_filling_products, bool integrate_openings_into_walls, bool use_ifcopeningelelements_for_virtual_boundaries) {

    // remove duplicate guids in products by replacing them
    ice::remove_duplicate_guids<Schema>(model);

    /*
     IfcElement (IfcWindow)
     IfcElement.FillsVoids -> IfcRelFillsElements
     IfcRelFillsElement.RelatingOpeningElement -> IfcOpeningElement
     IfcOpeningElement.VoidsElements -> IfcRelVoidsElements
     IfcRelVoidsElement.RelatingBuildingElement -> IfcElement (IfcWall)
     */

    for (const auto &c: include_entities) {

        boost::shared_ptr<IfcEntityList> IfcEntityList = nullptr;
        try { IfcEntityList = model->instances_by_type(c); }
        catch (...) {}

        if (IfcEntityList != nullptr)
            for (auto E: *IfcEntityList) { // for the entities, check if they are child of IfcElement class, only those can have or be an void filling product

                if (does_product_decompose_another_product<Schema>(E->as<typename Schema::IfcProduct>(), false, include_entities)) continue; // skip products decomposing parent products (e.g. windows in IfcCurtainWall or IfcMember in IfcStair)

                if (E->declaration().is("IfcElement")) {

                    if (integrate_openings_into_walls) {

                        if (E->as<typename Schema::IfcElement>()->FillsVoids()->size() == 0) {

                            if (E->declaration().is("IfcDoor") || E->declaration().is("IfcWindow"))
                                std::cerr << "[Warning] IfcElement is an IfcWindow or IfcDoor, but is not declared as void filling opening (" << E->data().toString().substr(0, 55) << " ...)." << std::endl;

                            non_void_filling_products.insert(E->as<typename Schema::IfcProduct>()->GlobalId()); // entity is not a void filling element

                        }
                    } else
                        non_void_filling_products.insert(E->as<typename Schema::IfcProduct>()->GlobalId()); // assume no entity is a void filling element
                } else
                    non_void_filling_products.insert(E->as<typename Schema::IfcProduct>()->GlobalId()); // entity is no element and therefore no void filling
            }
    }

    // this code applies if ifcopeningelements shall be used as normal objects to get virtual sbs (e.g. fzk haus or kit institute)
    if (use_ifcopeningelelements_for_virtual_boundaries) {
        boost::shared_ptr<IfcEntityList> IfcEntityList = nullptr;
        try { IfcEntityList = model->instances_by_type("IfcOpeningElement"); }
        catch (...) {}

        if (IfcEntityList != nullptr)
            for (auto E: *IfcEntityList) {

                auto IfcOpeningElement = E->as<typename Schema::IfcOpeningElement>();

                if (IfcOpeningElement->HasFillings()->size() == 0) {

                    auto IfcRelVoidsElements = IfcOpeningElement->VoidsElements();
                    for (auto tIfcRelVoidsElement: *IfcRelVoidsElements) {
                        auto IfcRelVoidsElement = tIfcRelVoidsElement->template as<typename Schema::IfcRelVoidsElement>();
                        auto RelatingBuildingElement = IfcRelVoidsElement->RelatingBuildingElement();
                        if (RelatingBuildingElement->declaration().is("IfcSlab")) {
                            non_void_filling_products.insert(E->as<typename Schema::IfcProduct>()->GlobalId()); // entity is not a void filling element
                            std::cout << "[Info] Add IfcOpeningElement as base for virtual space boundaries (" << E->data().toString().substr(0, 55) << " ...).\n";
                        }
                    }
                }
            }
    }
}

bool Kernel::identify_non_planar_products(std::list<Product> &products) {

    auto start = std::chrono::high_resolution_clock::now();

#ifdef PARALLEL_PROCESSING
    tbb::parallel_for_each(products.begin(), products.end(), [&](Product &product) {

        Topo top(product.shape);
        TopoDS_ListOfShape faces = top.faces();

        for (const auto &face: faces)
            if (!face_is_polygon(TopoDS::Face(face))) {
                product.hasOnlyPolygons = false;
                break;
            }
    });
#else
    for (auto &product: products) {
        Topo top(product.shape);
        TopoDS_ListOfShape faces = top.faces();

        for (auto & face : faces)
            if (!face_is_polygon(TopoDS::Face(face))) {
                product.hasOnlyPolygons = false;
                break;
            }
    }
#endif

    unsigned int count = 0;
    for (const auto &product: products)
        if (!product.hasOnlyPolygons) count++;

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Identify non-planar products", std::to_string(count) + "/" + std::to_string(products.size()));

    return count > 0;
}

void Kernel::polygonize_products(std::list<Product> &products) {

    auto start = std::chrono::high_resolution_clock::now();

#ifdef PARALLEL_PROCESSING
    tbb::parallel_for_each(products.begin(), products.end(), [&](Product &product) {
        if (!product.hasOnlyPolygons) {
            TopoDS_Shape copy = shape_copy(product.shape); // without copying function runs into segmentation fault
            TopoDS_Shape S = polygonize_shape_2a_curvature(copy);
            if (S.IsNull()) {
                std::cout << "[Info] Polygonize regarding corresponding faces not possible. " << product.guid << "\n";
                S = polygonize_shape(copy);
            }
            std::cout << "[Info] Polygonize " << product.guid << "\n";
            if (!S.IsNull()) {
                product.shape = S;
                product.hasOnlyPolygons = true;
            } else
                std::cerr << "[Warning] Triangulation failed!" << std::endl;
        }
    });
#else
    for (auto &product: products) {
        if (!product.hasOnlyPolygons) {
            TopoDS_Shape copy = shape_copy(product.shape); // without copying function runs into segmentation fault
            TopoDS_Shape S = polygonize_shape_2a_curvature(copy);
            if (S.IsNull()) {
                std::cout << "[Info] Polygonize regarding corresponding faces not possible. " << product.guid << "\n";
                S = polygonize_shape(copy);
            }
            std::cout << "[Info] Polygonize " << product.guid << "\n";
            if (!S.IsNull()) {
                product.shape = S;
                product.hasOnlyPolygons = true;
            } else
                std::cerr << "[Warning] Triangulation failed!" << std::endl;
        }
    }
#endif

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Polygonize products", std::to_string(products.size()));

}

void Kernel::collect_original_faces_clip(std::list<Product> &products, std::list<oFace> &orig_faces) {

    auto start = std::chrono::high_resolution_clock::now();

    for (auto &product: products) {
        for (auto &face: Topo(product.shape).faces()) {
            orig_faces.emplace_back(TopoDS::Face(face), &product, 0);
            product.orig_faces.push_back(&orig_faces.back());
            orig_faces.back().SetNormalStatus(FACE_NORMAL_KNOWN);
        }
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Collect original faces", std::to_string(orig_faces.size()));
}

void Kernel::collect_original_faces(std::list<Product> &products, std::list<oFace> &orig_faces) {

    auto start = std::chrono::high_resolution_clock::now();

    for (auto &product: products)
        collect_original_faces(product, orig_faces);

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Collect original faces", std::to_string(orig_faces.size()));
}

void Kernel::collect_original_faces(Product &product, std::list<oFace> &orig_faces) {

    unsigned int shell_id = 0;
    auto shells = Topo(product.shape).shells();

    for (auto &shell: shells) {
        for (auto &face: Topo(shell).faces()) {
            orig_faces.emplace_back(TopoDS::Face(face), &product, shell_id);
            product.orig_faces.push_back(&orig_faces.back());
            if (product.valid) orig_faces.back().SetNormalStatus(FACE_NORMAL_KNOWN);
            else std::cout << "[Info] Add face of non-valid product " << product.guid << ".\n";
        }
        shell_id++;
    }

    if (shells.Size() == 0)
        for (auto &face: Topo(product.shape).faces()) {
            orig_faces.emplace_back(TopoDS::Face(face), &product, shell_id);
            product.orig_faces.push_back(&orig_faces.back());
            std::cout << "[Info] Add non-shell face " << product.guid << ".\n";
        }
}

void Kernel::get_virtual_faces(std::unique_ptr<IfcParse::IfcFile> &model, ifcspaceInfoList &ifcspaces, std::list<oFace> &ifc_faces, std::list<Product> &products, std::set<IfcUtil::IfcBaseEntity *> &virtual_products) {

    auto start = std::chrono::high_resolution_clock::now();

    if (ifcSchema == IFC2X3)
        get_virtual_faces_worker<Ifc2x3>(ifcspaces, ifc_faces, products, virtual_products);
    else if (ifcSchema == IFC4)
        get_virtual_faces_worker<Ifc4>(ifcspaces, ifc_faces, products, virtual_products);
    else if (ifcSchema == IFC4X1)
        get_virtual_faces_worker<Ifc4x1>(ifcspaces, ifc_faces, products, virtual_products);
    else if (ifcSchema == IFC4X2)
        get_virtual_faces_worker<Ifc4x2>(ifcspaces, ifc_faces, products, virtual_products);
    else if (ifcSchema == IFC4X3_RC1)
        get_virtual_faces_worker<Ifc4x3_rc1>(ifcspaces, ifc_faces, products, virtual_products);
    else {
        std::cerr << "[Error] Schema not implemented yet! " << model->schema()->name() << std::endl;
        return;
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Calculate virtual faces by space clipping", std::to_string(ifc_faces.size()));
}

template<typename Schema>
void Kernel::get_virtual_faces_worker(ifcspaceInfoList &spaces, std::list<oFace> &orig_faces, std::list<Product> &products, std::set<IfcUtil::IfcBaseEntity *> &virtual_products) {

    rtree_lib::RTree<ifcspaceInfo *, double, 3, double> tree;

    // fill rtree
    for (auto &space: spaces) {
        Bnd_Box bnd = space.bbox;
        double min[3] = {bnd.CornerMin().X(), bnd.CornerMin().Y(), bnd.CornerMin().Z()};
        double max[3] = {bnd.CornerMax().X(), bnd.CornerMax().Y(), bnd.CornerMax().Z()};
        tree.Insert(min, max, &space);
    }

    // find intersections
    std::set<std::set<ifcspaceInfo *>> pairs;

    for (auto &space: spaces) {
        Bnd_Box bnd = space.bbox;
        double min[3] = {bnd.CornerMin().X(), bnd.CornerMin().Y(), bnd.CornerMin().Z()};
        double max[3] = {bnd.CornerMax().X(), bnd.CornerMax().Y(), bnd.CornerMax().Z()};
        tree.Search(min, max, [&space, &pairs](ifcspaceInfo *found_space) {
            if (space.guid != found_space->guid)
                pairs.insert({&space, found_space});
            return true;
        });
    }

    TopoDS_ListOfShape V; // virtual faces

    for (auto &pair: pairs) {
        ifcspaceInfo *s1 = *pair.begin();
        ifcspaceInfo *s2 = *pair.rbegin();

        for (const auto &face1: Topo(s1->shape).faces()) {

            gp_Dir n1 = face_normal(TopoDS::Face(face1));

            for (const auto &face2: Topo(s2->shape).faces()) {

                gp_Dir n2 = face_normal(TopoDS::Face(face2));

                if (n1.Angle(n2) < 3.14) continue;
                if (minimal_distance_between_shapes(face1, face2) > 1.0e-4) continue;

                Clipper clip(TopoDS::Face(face1), n1, true);
                clip.clip(TopoDS::Face(face2), n2);
                if (clip.success) V.Append(clip.Result);
                else if (clip.same) V.Append(face1);
            }
        }
    }

    for (const auto &face: V) {
        auto IfcVirtualElement = ice::IfcVirtualElement<Schema>();
        products.emplace_back(IfcVirtualElement, IfcVirtualElement->GlobalId(), TopoDS::Face(face));
        orig_faces.emplace_back(TopoDS::Face(face), &products.back(), 0);
        virtual_products.insert(IfcVirtualElement);
    }
}

void Kernel::add_offset_faces(std::list<oFace> &orig_faces, double offset, const std::set<std::string> &classes_for_face_extension) {

    auto start = std::chrono::high_resolution_clock::now();

    std::list<oFace> offsets_of_orig_faces;

    for (const auto &orig_face: orig_faces) {

        //if (area(orig_face.face) < 2) continue;

        bool to_extend = false;
        for (const auto &cl: classes_for_face_extension)
            if (orig_face.IsIfcClass(cl)) {
                to_extend = true;
                break;
            }

        if (!to_extend)
            continue;

        // skip wall faces that are product of cut with window/door
        if (orig_face.Opening() != nullptr) continue;

        TopoDS_Face offset_face;
        if (offset_faces_outer_wire_cut(orig_face.face, offset_face, offset)) {
            if (!orig_face.Normal().IsEqual(face_normal(offset_face), 0.1))
                offset_face.Reverse();

            offsets_of_orig_faces.emplace_back(offset_face, orig_face.RelProduct(), orig_face.ShellID());
            offsets_of_orig_faces.back().SetIsOffset(true);

            if (!orig_face.Normal().IsEqual(offsets_of_orig_faces.back().Normal(), 0.1))
                std::cerr << "Error: Face normal changed. " << std::endl;
        }
    }

    for (const auto &f: offsets_of_orig_faces)
        orig_faces.push_back(f);

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Add offset faces", std::to_string(offsets_of_orig_faces.size()) + "/" + std::to_string(orig_faces.size()));

}

template<typename Schema>
bool Kernel::does_product_decompose_another_product(typename Schema::IfcProduct *IfcProduct, bool consider_allowed_classes, const std::set<std::string> &inc) {

    // return false if product does not decompose another product (of a considered ifc class)
    if (consider_allowed_classes) {
        auto IfcRelAggregates = IfcProduct->Decomposes();
        for (const auto &IfcRelAggregate: *IfcRelAggregates)
            for (const auto &c: inc)
                if (IfcRelAggregate->RelatingObject()->declaration().is(c))
                    return true;
    } else {
        auto IfcRelAggregates = IfcProduct->Decomposes();
        for (const auto &IfcRelAggregate: *IfcRelAggregates)
            return true;
    }

    return false;
}

void Kernel::add_decomposed_entities_to_products(std::unique_ptr<IfcParse::IfcFile> &model, IfcGeom::IteratorSettings settings, std::list<Product> &products, const std::set<std::string> &include_entities) {

    auto start = std::chrono::high_resolution_clock::now();

    if (ifcSchema == IFC2X3)
        add_decomposed_entities_to_products_worker<Ifc2x3>(model, settings, products, include_entities);
    else if (ifcSchema == IFC4)
        add_decomposed_entities_to_products_worker<Ifc4>(model, settings, products, include_entities);
    else if (ifcSchema == IFC4X1)
        add_decomposed_entities_to_products_worker<Ifc4x1>(model, settings, products, include_entities);
    else if (ifcSchema == IFC4X2)
        add_decomposed_entities_to_products_worker<Ifc4x2>(model, settings, products, include_entities);
    else if (ifcSchema == IFC4X3_RC1)
        add_decomposed_entities_to_products_worker<Ifc4x3_rc1>(model, settings, products, include_entities);
    else
        std::cerr << "[Error] Schema not implemented yet! " << model->schema()->name() << std::endl;

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Add decomposed products", std::to_string(products.size()));
}

template<typename Schema>
void Kernel::add_decomposed_entities_to_products_worker(std::unique_ptr<IfcParse::IfcFile> &model, IfcGeom::IteratorSettings settings, std::list<Product> &products, const std::set<std::string> &include_entities) {

    // map to check if product does already exists and therefore already has a shape on its own
    std::unordered_map<IfcUtil::IfcBaseEntity *, Product *> M;
    for (auto &product: products)
        M[product.ifcproduct] = &product;

    // entities that are converted to OBB along with its decomposing entities
    const std::set<std::string> &to_be_simplified_classes = {"IfcCurtainWall", "IfcRoof"};

    for (const auto &c: include_entities) {

        boost::shared_ptr<IfcEntityList> IfcEntityList = nullptr;
        try { IfcEntityList = model->instances_by_type(c); }
        catch (...) {}

        if (IfcEntityList != nullptr)
            for (auto E: *IfcEntityList) {

                auto P = E->as<typename Schema::IfcProduct>();

                if (does_product_decompose_another_product<Schema>(P, false, include_entities)) continue;

                auto guids = guids_of_products_decomposing_another_product<Schema>(P, include_entities);
                if (guids.empty()) continue;

                TopoDS_ListOfShape L = products_to_shapelist(model, settings, guids);
                if (L.IsEmpty()) continue;

                // existing and kind of geometrical processing
                bool exists = (M.find(P) != M.end());
                bool to_be_simplified = false;
                for (const auto &tbsc: to_be_simplified_classes)
                    if (P->declaration().is(tbsc)) {
                        to_be_simplified = true;
                        break;
                    }

                if (exists)
                    L.Append(M[P]->shape);

                TopoDS_Shape S;
                if (to_be_simplified) S = best_fitting_bbox(L);
                else S = compound_from_shape_list(L);

                if (S.IsNull()) continue;

                if (exists)
                    M[P]->shape = S;
                else
                    products.emplace_back(P, P->GlobalId(), S);
            }
    }
}

TopoDS_ListOfShape Kernel::products_to_shapelist(std::unique_ptr<IfcParse::IfcFile> &model, IfcGeom::IteratorSettings settings, const std::set<std::string> &guids) const {

    IfcGeom::attribute_filter attribute_filter;
    attribute_filter.include = true;
    attribute_filter.traverse = false;
    attribute_filter.attribute_name = "GlobalId";
    attribute_filter.populate(guids);
    std::vector<IfcGeom::filter_t> filter_funcs;
    filter_funcs.emplace_back(boost::ref(attribute_filter));

    IfcGeom::Iterator<real_t> geom_iterator(settings, model.get(), filter_funcs, num_threads);
    geom_iterator.initialize();

    TopoDS_ListOfShape L;

    do {
        IfcGeom::Element<real_t> *geom_object = geom_iterator.get();
        if (!check_ifc_context(geom_object)) continue;
        auto S = geom_object_to_shape(geom_object);
        if (!S.IsNull()) L.Append(S);
    } while (geom_iterator.next());

    return L;
}

template<typename Schema>
std::set<std::string> Kernel::guids_of_products_decomposing_another_product(typename Schema::IfcProduct *IfcProduct, const std::set<std::string> &inc) {

    // returns entities that are decomposing the product
    std::set<std::string> guids;

    auto IfcRelAggregates = IfcProduct->IsDecomposedBy();

    for (auto &IfcRelAggregate: *IfcRelAggregates) {
        auto RelatedObjects = IfcRelAggregate->RelatedObjects();
        for (auto &RelatedObject: *RelatedObjects) {
            auto P = RelatedObject->template as<typename Schema::IfcProduct>();
            for (const auto &c: inc)
                if (P->declaration().is(c)) {
                    guids.insert(P->GlobalId());
                    break;
                }
        }
    }

    return guids;
}

void Kernel::remove_triangulation(std::list<oFace> &orig_faces) {

    for (const auto &f: orig_faces)
        BRepTools::Clean(f.face);
}

void Kernel::create_void_filling_cuts_in_products(std::list<Product> &products, std::list<oFace> &orig_faces, std::unordered_map<IfcUtil::IfcBaseClass *, TopoDS_Shape> &opening_shapes,
                                                  bool simplifyOpeningElement, std::unique_ptr<IfcParse::IfcFile> &model, IfcGeom::IteratorSettings settings) const {

    auto start = std::chrono::high_resolution_clock::now();

    std::vector<std::set<oFace *>> removes; // orig_faces that need to be removed
    std::vector<std::list<opening_info>> inserts; // information for orig_faces to be created (face, product, shell id, parent id, opening product)
    std::unordered_map<std::string, Product *> bad_openings; // openings that could not operated with a wall and the wall product

#ifdef PARALLEL_PROCESSING
    std::vector<Product *> V;
    V.reserve(products.size());
    for (auto &product: products)
        V.push_back(&product);

#pragma omp parallel for default(none) shared(V, removes, inserts, bad_openings, opening_shapes, simplifyOpeningElement) num_threads(num_threads)
    for (auto p: V) {
        if (p->VoidFillingElements.empty()) continue;
        std::set<oFace *> rem;
        std::list<opening_info> ins;
        std::unordered_map<std::string, Product *> bad;

        bool changes = create_void_filling_cuts_in_products_worker(*p, opening_shapes, simplifyOpeningElement, rem, ins, bad);
        if (changes) {
#pragma omp critical
            {
                removes.push_back(rem);
                inserts.push_back(ins);
                bad_openings.insert(bad.begin(), bad.end());
            }
        } else {
#pragma omp critical
            bad_openings.insert(bad.begin(), bad.end());
        }
    }
#else
    for (auto &product: products) {
        if (product.VoidFillingElements.empty()) continue;
        std::set<oFace *> rem;
        std::list<opening_info> ins;
        bool changes = create_void_filling_cuts_in_products_worker(product, opening_shapes, simplifyOpeningElement, rem, ins);
        if (changes) {
            removes.push_back(rem);
            inserts.push_back(ins);
        }
    }
#endif

    // remove
    std::set<oFace *> s;
    for (auto &rem: removes)
        s.insert(rem.begin(), rem.end());
    orig_faces.remove_if([&s](oFace &O) { return s.find(&O) != s.end(); });

    // insert
    for (auto &ins: inserts)
        for (auto &T: ins) {
            orig_faces.emplace_back(T.face, T.product, T.shell_id);
            orig_faces.back().SetParentID(T.parent_id);
            orig_faces.back().SetOpening(T.opening_product);
            orig_faces.back().SetNormalStatus(T.normal_status);
        }

    // generate window/door geometries (bad_openings) that could not be integrated into wall e.g. because wall already has hole and simplify them
    if (!bad_openings.empty())
        bad_openings_to_box_faces(orig_faces, model, settings, products, bad_openings);

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Add window/door faces to walls", std::to_string(products.size()) + "/" + std::to_string(orig_faces.size()));
}

Kernel::OWP Kernel::cvfcipw_collect_opening_shapes(const Product &product, std::unordered_map<IfcUtil::IfcBaseClass *, TopoDS_Shape> &opening_shapes, bool simplifyOpeningElement) {

    OWP OSW;

    for (auto &OpeningElement: product.VoidFillingElements)
        if (opening_shapes.find(OpeningElement.first) != opening_shapes.end()) {
            if (simplifyOpeningElement) {
                auto bnd = aabb(opening_shapes[OpeningElement.first], 1.0e-5);
                OSW[OpeningElement.first] = std::make_pair(aabb_to_shape(bnd), OpeningElement.second);
            } else {
                // Opening elements having multiple solids (e.g. by extruding a face in two opposite directions creating two solids) made some problems.
                // A cut result was e.g. a face consisting of two spatial separated faces connected by an edge (in plane of the extrusion face).
                // Unify all solids of the shape and their faces beforehand, to minimize the number of faces (6 instead of 12 for example)
                OSW[OpeningElement.first] = std::make_pair(unify_shape(opening_shapes[OpeningElement.first]), OpeningElement.second);
                //OSW[OpeningElement.first] = std::make_pair(opening_shapes[OpeningElement.first], OpeningElement.second);
            }
        }

    return OSW;
}

void Kernel::cvfcipw_check_for_common_volume(Product &product, std::list<opening_info> &ins, std::unordered_map<std::string, Product *> &togenerate, size_t n, OWP &OSW) {

    unsigned int c = -1;

    for (auto it = OSW.cbegin(); it != OSW.cend();) {
        TopoDS_Shape comm = BRepAlgoAPI_Common(product.shape, it->second.first).Shape();
        if (volume(comm) < 0.001) {
            c++;
            std::string guid = remove_first_and_last_char(it->second.second->data().getArgument(0)->toString());
            std::cout << "[Info] Parent element " << product.guid << " (" << product.ifcproduct->declaration().name()
                      << ") of " << guid << " (" << it->second.second->declaration().name()
                      << ") linked by " << remove_first_and_last_char(it->first->data().getArgument(0)->toString()) << " (" << it->first->declaration().name()
                      << ") already has a hole where opening should be cut.\n";

            // add simplified window/door box faces later
            togenerate.insert({guid, &product});

            // add it now as simplified box
/*            auto box = aabb(it->second.first, 0);
            gp_Pnt P1 = box.CornerMin();
            gp_Pnt P2 = box.CornerMax();
            double z = 0.5 * (P1.Z() + P2.Z());
            double w = 0.10;
            auto S = BRepPrimAPI_MakeBox(gp_Pnt(P1.X(), P1.Y(), z - w), gp_Pnt(P2.X(), P2.Y(), z + w)).Shape();
            for (auto &f: Topo(aabb_to_shape(box)).faces()) // simplified window
                ins.emplace_back(TopoDS::Face(f), &product, c + 2 * n, 0, it->second.second, FACE_NORMAL_KNOWN);*/

            OSW.erase(it++);
        } else ++it;
    }
}

void Kernel::cvfcipw_commons_between_opening_and_wall_faces(Product &product, BOPAlgo_CellsBuilder &CellBuilder, TopTools_ListOfShape &take, TopTools_ListOfShape &avoi, const TopoDS_Compound &Shell, const std::set<unsigned int> &ids,
                                                            const OWP &OSW, std::list<opening_info> &ins, const std::pair<const unsigned int, std::list<oFace *>> &m) {

    // Commons between openings and wall faces
    for (const auto &O: OSW) {
        take.Clear();
        avoi.Clear();
        CellBuilder.RemoveAllFromResult();
        take.Append(Shell);
        take.Append(O.second.first);
        for (const auto &O2: OSW) {
            if (O.first == O2.first) break;
            avoi.Append(O2.second.first);
        }
        CellBuilder.AddToResult(take, avoi);

        // It is possible to just add all faces of the CellBuilder.Shape() to the orig_faces ...
        //  add_cut_faces_to_orig_face(orig_faces, Topo(CellBuilder.Shape()).faces(), &product, m.first, O.second.second);
        // But I want to check the face normals
        for (auto &orig_face: m.second) {

            auto NewFaces = CellBuilder.Modified(orig_face->face);
            if (NewFaces.IsEmpty()) continue;

            auto itr = TopTools_ListIteratorOfListOfShape(NewFaces);
            while (itr.More()) {

                if (face_normal(TopoDS::Face(itr.Value())).Angle(orig_face->Normal()) > 0.1)
                    itr.Value().Complement();

                // add to orig_faces
                if (ids.find(hash(itr.Value())) == ids.end()) // if ifc in ifc the walls do already have openings before openingelements are subtracted,
                    // result of common operation are side faces of wall in the opening. Those are already taken into account by cut. So skip them, otherwise duplicate faces and ids
                    ins.emplace_back(TopoDS::Face(itr.Value()), &product, m.first, hash(orig_face->face), O.second.second, orig_face->NormalStatus());

                itr.Next();
            }
        }
    }
}

std::set<unsigned int>
Kernel::cvfcipw_cut_openings_from_wall_faces(Product &product, BOPAlgo_CellsBuilder &CellBuilder, TopTools_ListOfShape &take, TopTools_ListOfShape &avoi, const TopoDS_Compound &Shell, const OWP &OSW,
                                             std::list<opening_info> &ins, const std::pair<const unsigned int, std::list<oFace *>> &m) {

    // cut openings from wall faces
    take.Clear();
    avoi.Clear();
    take.Append(Shell);
    for (const auto &O: OSW)
        avoi.Append(O.second.first);
    CellBuilder.AddToResult(take, avoi);

    // add to orig_faces
    std::set<unsigned int> ids;

    // fast way, without tracking parentship:
    // for (auto &face : Topo(CellBuilder.Shape()).faces())
    //    orig_faces.emplace_back(TopoDS::Face(face), &product, m.first);

    // slow way, tracking parentship:
    for (auto &orig_face: m.second) {
        auto NewFaces = CellBuilder.Modified(orig_face->face);

        // add to orig_faces
        if (NewFaces.IsEmpty()) {
            ins.emplace_back(TopoDS::Face(orig_face->face), &product, m.first, orig_face->ParentID(), orig_face->Opening(), orig_face->NormalStatus());
            ids.insert(hash(orig_face->face));
        } else {
            auto itr = TopTools_ListIteratorOfListOfShape(NewFaces);
            while (itr.More()) {
                ins.emplace_back(TopoDS::Face(itr.Value()), &product, m.first, hash(orig_face->face), orig_face->Opening(), orig_face->NormalStatus());
                ids.insert(hash(itr.Value()));
                itr.Next();
            }
        }
    }
    return ids;
}

TopoDS_Compound Kernel::cvfcipw_compound_from_orig_faces(size_t n, const std::list<oFace *> &orig_faces) {

    // Create Compound from set of orig faces
    TopoDS_Compound Shell;
    BRep_Builder B;
    B.MakeCompound(Shell);
    for (auto &orig_face: orig_faces) {
        B.Add(Shell, orig_face->face);
        orig_face->SetShellID(orig_face->ShellID() + n); // change of shell_id is necessary. Otherwise the faces in cut that were not altered will be deleted within orig_faces.remove(*orig_face)
    }
    return Shell;
}

void Kernel::cvfcipw_perform_fuse(const std::string &guid, BOPAlgo_CellsBuilder &CellBuilder, const TopoDS_Compound &Shell, const OWP &OSW) {

    // Add Walls and Openings to the Fuse list and Perform Fuse
    TopTools_ListOfShape L;
    L.Append(Shell);
    for (const auto &O: OSW)
        L.Append(O.second.first);
    CellBuilder.SetArguments(L);
    CellBuilder.Perform();

    if (CellBuilder.HasWarnings()) {
        std::cerr << "Warnings for " << guid << ":" << std::endl;
        CellBuilder.DumpWarnings(std::cerr);
    }
    if (CellBuilder.HasErrors()) {
        std::cerr << "Errors for " << guid << ":" << std::endl;
        CellBuilder.DumpErrors(std::cerr);
    }
}

bool Kernel::create_void_filling_cuts_in_products_worker(Product &product, std::unordered_map<IfcUtil::IfcBaseClass *, TopoDS_Shape> &opening_shapes, bool simplifyOpeningElement, std::set<oFace *> &rem, std::list<opening_info> &ins, std::unordered_map<std::string, Product *> &togenerate) {

    // Sort orig faces of the products by shell id
    std::unordered_map<unsigned int, std::list<oFace *>> M; // key: shell_id, value: faces
    for (auto &orig_face: product.orig_faces)
        M[orig_face->ShellID()].push_back(orig_face);

    // Collect opening shapes
    auto OSW = cvfcipw_collect_opening_shapes(product, opening_shapes, simplifyOpeningElement); // key: IfcOpeningElement, value: a pair of IfcOpeningElement shape and IfcProduct (window/door)
    if (OSW.empty()) return false;

    // Check if opening element has common volume with product. if not, parent already has a hole even before opening subtraction. remove from OSW then but add openingelement faces to solution
    cvfcipw_check_for_common_volume(product, ins, togenerate, M.size(), OSW);
    if (OSW.empty()) return true;

    // Setup CellBuilder
    BOPAlgo_CellsBuilder CellBuilder;
    TopTools_ListOfShape take;
    TopTools_ListOfShape avoi;

    for (auto &m: M) { // for every set of orig faces belonging to a shell

        CellBuilder.Clear();

        // Create Compound from set of orig faces
        TopoDS_Compound Shell = cvfcipw_compound_from_orig_faces(M.size(), m.second);

        // Fuse shell and shape of IfcOpeningElement
        cvfcipw_perform_fuse(product.guid, CellBuilder, Shell, OSW);

        // Cut openings from wall faces
        std::set<unsigned int> ids = cvfcipw_cut_openings_from_wall_faces(product, CellBuilder, take, avoi, Shell, OSW, ins, m);

        // Get Commons between openings and wall faces
        cvfcipw_commons_between_opening_and_wall_faces(product, CellBuilder, take, avoi, Shell, ids, OSW, ins, m);

        // Delete all old faces belonging to current shell
        for (auto &orig_face: m.second)
            rem.insert(orig_face);
    }

    return true;
}

std::unordered_map<IfcUtil::IfcBaseClass *, TopoDS_Shape> Kernel::generate_opening_shapes(std::unique_ptr<IfcParse::IfcFile> &model, IfcGeom::IteratorSettings settings, const std::set<std::string> &guids) const {

    auto start = std::chrono::high_resolution_clock::now();

    std::unordered_map<IfcUtil::IfcBaseClass *, TopoDS_Shape> opening_shapes;

    IfcGeom::attribute_filter attribute_filter;
    attribute_filter.include = true;
    attribute_filter.traverse = false;
    attribute_filter.attribute_name = "GlobalId";
    attribute_filter.populate(guids);
    std::vector<IfcGeom::filter_t> filter_funcs;
    filter_funcs.emplace_back(boost::ref(attribute_filter));

    IfcGeom::Iterator<real_t> geom_iterator(settings, model.get(), filter_funcs, num_threads); // construct iterator
    unsigned int n = 0;

    if (!geom_iterator.initialize()) {
        std::cerr << "No geometrical entities found. Initialization failed." << std::endl;
        return opening_shapes;
    }

    do {
        IfcGeom::Element<real_t> *geom_object = geom_iterator.get();
        if (!check_ifc_context(geom_object)) continue;
        //if (geom_object->guid() == "04RcXu1mKxpGnJw_7oRBLn" || geom_object->guid() == "0LM8GvGe$G3dlW4mZ4aA9R") continue;
        opening_shapes[geom_object->product()] = geom_object_to_shape(geom_object);
        n++;
    } while (geom_iterator.next());

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Generate shapes of openings", std::to_string(n) + "/" + std::to_string(guids.size()));

    return opening_shapes;
}

TopoDS_Shape Kernel::geom_object_to_shape(IfcGeom::Element<real_t> *geom_object) {

    const auto *o = dynamic_cast<const IfcGeom::BRepElement<real_t> *>(geom_object);
    TopoDS_Shape shape = o->geometry().as_compound();
    gp_Trsf trsf = o->transformation().data();
    shape.Move(trsf);
    return shape;
}

void Kernel::remove_products(std::list<Product> &products, const std::list<Product *> &products_to_remove) {

    auto start = std::chrono::high_resolution_clock::now();

    for (auto &p: products_to_remove)
        products.remove(*p);

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Remove products", std::to_string(products.size()));
}

void Kernel::remove_products(std::list<Product> &products) {

    auto start = std::chrono::high_resolution_clock::now();

    products.remove_if([](Product &p) {
        if (!p.valid) {
            std::cerr << "[Warning] Remove product " << p.guid << std::endl;
            return true;
        } else
            return false;
    });

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Remove products", std::to_string(products.size()));
}

void Kernel::heal_products(std::list<Product *> &products_to_heal) {

    std::cout << "Heal products ";

    auto start = std::chrono::high_resolution_clock::now();

    for (const auto &p: products_to_heal) {
        TopoDS_Shape S = ShapeHealing(p->shape).heal_shape();
        if (S.IsNull()) std::cerr << "[Warning] Shape healing failed." << std::endl;
        else p->shape = S;
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "(" << products_to_heal.size() << ") \t... \tElapsed time: " << elapsed.count() << " s\n";
}

TopoDS_Shape Kernel::replace_subshapes_in_shape(const std::list<std::pair<TopoDS_Shape, TopoDS_Shape>> &old_and_new_subshapes, TopoDS_Shape &shape) {

    ShapeBuild_ReShape R;
    for (const auto &p: old_and_new_subshapes)
        R.Replace(p.first, p.second);
    return R.Apply(shape);
}

bool Kernel::find_duplicate_hashes_in_orig_faces(std::list<oFace> &orig_faces, std::set<unsigned int> &bad_hashes, std::set<unsigned int> &bad_hashes_location) {

    auto start = std::chrono::high_resolution_clock::now();

    std::set<unsigned int> a; // all hashes
    std::set<std::pair<unsigned int, unsigned int>> b; // incl. location

    for (auto &orig_face: orig_faces) {
        unsigned int h = Kernel::hash(orig_face.face);
        unsigned int l = orig_face.face.Location().HashCode(INT_MAX);

        if (a.find(h) != a.end())
            bad_hashes.insert(h);
        a.insert(h);

        auto p = std::make_pair(h, l);
        if (b.find(p) != b.end())
            bad_hashes_location.insert(h);
        b.insert(p);
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Find duplicate hashes in original faces", std::to_string(bad_hashes.size()) + "/" + std::to_string(bad_hashes_location.size()));

    return bad_hashes.empty();
}

ifcspaceInfoList Kernel::create_shapes_of_IfcSpaces(std::unique_ptr<IfcParse::IfcFile> &model, IfcGeom::IteratorSettings settings) const {

    auto start = std::chrono::high_resolution_clock::now();

    std::set<std::string> include_entities = {"IfcSpace"};
    ifcspaceInfoList old_spaces;

    IfcGeom::entity_filter entity_filter;
    entity_filter.include = true;
    entity_filter.traverse = false;
    entity_filter.entity_names = include_entities;
    std::vector<IfcGeom::filter_t> filter_funcs;
    filter_funcs.emplace_back(boost::ref(entity_filter));

    IfcGeom::Iterator<real_t> geom_iterator(settings, model.get(), filter_funcs, num_threads);

    if (!geom_iterator.initialize()) {
        std::cerr << "No geometrical entities found. Initialization failed." << std::endl;
        return old_spaces;
    }

    do {
        IfcGeom::Element<real_t> *geom_object = geom_iterator.get();
        if (!check_ifc_context(geom_object)) continue;
        TopoDS_Shape S = geom_object_to_shape(geom_object);
        auto bnd = aabb(S, 0.01);

        S = ShapeHealing(S).heal_shape();
        if (S.IsNull()) std::cerr << "[Warning] Shape healing failed." << std::endl;
        auto L = Topo(S).shells();

        TopoDS_Solid solid;

        if (!S.IsNull() && !L.IsEmpty()) {
            solid = BRepBuilderAPI_MakeSolid(TopoDS::Shell(L.First())).Solid();
            if (volume(solid) < 0) solid.Complement();
        }

        old_spaces.emplace_back(geom_object->product(), geom_object->guid(), solid, bnd, volume(solid));

    } while (geom_iterator.next());

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), " Create shapes of IfcSpaces", std::to_string(old_spaces.size()));

    return old_spaces;
}

void Kernel::check_orig_faces(std::list<oFace> &orig_faces) {

    auto start = std::chrono::high_resolution_clock::now();

    std::map<Product *, std::list<oFace * >> O;
    for (auto &origface: orig_faces)
        O[origface.RelProduct()].push_back(&origface);

    for (auto &p: O) {

        std::map<unsigned int, std::list<unsigned int>> M;
        std::list<unsigned int> L;
        std::set<unsigned int> S;

        for (auto &origface: p.second) {
            unsigned int id = hash(origface->face);
            M[origface->ShellID()].push_back(id);
            L.push_back(id);
            S.insert(id);
        }

        if (L.size() == S.size()) continue;

        for (const auto &s: S) {
            auto n = std::count(L.begin(), L.end(), s);
            if (n > 1) std::cerr << "[Warning] Face " << s << " occurs multiple times (" << n << ") in product " << p.first->guid << std::endl;
        }

        for (const auto &m: M)
            for (const auto &s: S) {
                auto n = std::count(m.second.begin(), m.second.end(), s);
                if (n > 1) std::cerr << "[Warning] Face " << s << " occurs multiple times (" << n << ") in shell " << m.first << std::endl;
            }
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Check original faces", std::to_string(orig_faces.size()));
}

void Kernel::check_products(std::list<Product> &products, double l_warn, double l_crit) {

    auto start = std::chrono::high_resolution_clock::now();

#ifdef PARALLEL_PROCESSING
    tbb::parallel_for_each(products.begin(), products.end(), [&](Product &p) {

        bool is_good = true;

        if (p.shape.IsNull()) {
            std::cerr << "[Warning] Null shape " << hash(p.shape) << "\t" << p.guid << std::endl;
            is_good = false;
        } else {
            double V = volume(p.shape);
            if (V < 1.0e-8) {
                std::cerr << "[Warning] Shape small volume " << hash(p.shape) << "\t" << p.guid << "\t" << V << std::endl;
                is_good = false;
            }

            BRepCheck_Analyzer A(p.shape);
            if (!A.IsValid()) {
                std::cerr << "[Warning] Shape not valid " << hash(p.shape) << "\t" << p.guid << std::endl;
                is_good = false;
            }

            auto shells = Topo(p.shape).shells();
            if (shells.IsEmpty()) {
                std::cerr << "[Warning] No shell " << hash(p.shape) << "\t" << p.guid << std::endl;
                is_good = false;
            } else if (!p.CheckShells()) is_good = false;

            for (auto &s: shells) {
                BRepCheck_Analyzer B(s);
                if (!B.IsValid()) {
                    std::cerr << "[Warning] Shell not valid " << hash(s) << "\t" << p.guid << std::endl;
                    is_good = false;
                }

                ShapeAnalysis_Shell a;
                a.LoadShells(TopoDS::Shell(s));
                a.CheckOrientedShells(s, true, true);
                if (a.HasBadEdges()) {
                    std::cerr << "[Warning] Shell bad edges " << hash(s) << "\t" << p.guid << std::endl;
                    is_good = false;
                }
                if (a.HasFreeEdges()) {
                    std::cerr << "[Warning] Shell free edges " << hash(s) << "\t" << p.guid << std::endl;
                    is_good = false;
                }
            }

//        TopTools_IndexedDataMapOfShapeListOfShape faceShellMap;
//        TopExp::MapShapesAndAncestors(p.shape, TopAbs_FACE, TopAbs_SHELL, faceShellMap);

            for (auto &f: Topo(p.shape).faces()) {
                BRepCheck_Analyzer B(f);
                if (!B.IsValid()) {
                    std::cerr << "[Warning] Face not valid " << hash(f) << "\t" << p.guid << std::endl;
                    is_good = false;
                }

//            const TopTools_ListOfShape &fshells = faceShellMap.FindFromKey(TopoDS::Face(f));
//            if (fshells.IsEmpty()) {
//                std::cerr << "[Warning] Face has no shell " << hash(f) << "\t" << p.guid << std::endl;
//                is_good = false;
//            }

                if (!check_face(TopoDS::Face(f), 0, p.Info(), l_warn, l_crit)) is_good = false;
            }

            for (auto &w: Topo(p.shape).wires()) {
                BRepCheck_Analyzer B(w);
                if (!B.IsValid()) {
                    std::cerr << "[Warning] Wire not valid " << hash(w) << "\t" << p.guid << std::endl;
                    is_good = false;
                }
            }

            for (auto &e: Topo(p.shape).edges()) {
                BRepCheck_Analyzer B(e);
                if (!B.IsValid()) {
                    std::cerr << "[Warning] Edge not valid " << hash(e) << "\t" << p.guid << std::endl;
                    is_good = false;
                }
            }

            for (auto &v: Topo(p.shape).vertices()) {
                BRepCheck_Analyzer B(v);
                if (!B.IsValid()) {
                    std::cerr << "[Warning] Vertex not valid " << hash(v) << "\t" << p.guid << std::endl;
                    is_good = false;
                }
            }
        }

        p.valid = is_good;
    });
#else
    for (auto &p: products) {

        bool is_good = true;

        if (p.shape.IsNull()) {
            std::cerr << "[Warning] Null shape " << hash(p.shape) << "\t" << p.guid << std::endl;
            is_good = false;
        } else {
            double V = volume(p.shape);
            if (V < 1.0e-8) {
                std::cerr << "[Warning] Shape small volume " << hash(p.shape) << "\t" << p.guid << "\t" << V << std::endl;
                is_good = false;
            }

            BRepCheck_Analyzer A(p.shape);
            if (!A.IsValid()) {
                std::cerr << "[Warning] Shape not valid " << hash(p.shape) << "\t" << p.guid << std::endl;
                is_good = false;
            }

            auto shells = Topo(p.shape).get_shells();
            if (shells.IsEmpty()) {
                std::cerr << "[Warning] No shell " << hash(p.shape) << "\t" << p.guid << std::endl;
                is_good = false;
            } else if (!p.CheckShells()) is_good = false;

            for (auto &s: shells) {
                BRepCheck_Analyzer B(s);
                if (!B.IsValid()) {
                    std::cerr << "[Warning] Shell not valid " << hash(s) << "\t" << p.guid << std::endl;
                    is_good = false;
                }

                ShapeAnalysis_Shell a;
                a.LoadShells(TopoDS::Shell(s));
                a.CheckOrientedShells(s, true, true);
                if (a.HasBadEdges()) {
                    std::cerr << "[Warning] Shell bad edges " << hash(s) << "\t" << p.guid << std::endl;
                    is_good = false;
                }
                if (a.HasFreeEdges()) {
                    std::cerr << "[Warning] Shell free edges " << hash(s) << "\t" << p.guid << std::endl;
                    is_good = false;
                }
            }

//        TopTools_IndexedDataMapOfShapeListOfShape faceShellMap;
//        TopExp::MapShapesAndAncestors(p.shape, TopAbs_FACE, TopAbs_SHELL, faceShellMap);

            for (auto &f: Topo(p.shape).faces()) {
                BRepCheck_Analyzer B(f);
                if (!B.IsValid()) {
                    std::cerr << "[Warning] Face not valid " << hash(f) << "\t" << p.guid << std::endl;
                    is_good = false;
                }

//            const TopTools_ListOfShape &fshells = faceShellMap.FindFromKey(TopoDS::Face(f));
//            if (fshells.IsEmpty()) {
//                std::cerr << "[Warning] Face has no shell " << hash(f) << "\t" << p.guid << std::endl;
//                is_good = false;
//            }

                if (!check_face(TopoDS::Face(f), 0, p.Info(), l_warn, l_crit)) is_good = false;
            }

            for (auto &w: Topo(p.shape).get_wires()) {
                BRepCheck_Analyzer B(w);
                if (!B.IsValid()) {
                    std::cerr << "[Warning] Wire not valid " << hash(w) << "\t" << p.guid << std::endl;
                    is_good = false;
                }
            }

            for (auto &e: Topo(p.shape).get_edges()) {
                BRepCheck_Analyzer B(e);
                if (!B.IsValid()) {
                    std::cerr << "[Warning] Edge not valid " << hash(e) << "\t" << p.guid << std::endl;
                    is_good = false;
                }
            }

            for (auto &v: Topo(p.shape).get_vertices()) {
                BRepCheck_Analyzer B(v);
                if (!B.IsValid()) {
                    std::cerr << "[Warning] Vertex not valid " << hash(v) << "\t" << p.guid << std::endl;
                    is_good = false;
                }
            }
        }

        p.flag = is_good;
    }
#endif

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Check products", std::to_string(products.size()));
}

void Kernel::analyze_orig_faces(std::list<oFace> &orig_faces, double l_warn, double l_crit) {

    auto start = std::chrono::high_resolution_clock::now();

#ifdef PARALLEL_PROCESSING
    tbb::parallel_for_each(orig_faces.begin(), orig_faces.end(), [&](oFace &o) {
        BRepCheck_Analyzer A(o.face);
        if (!A.IsValid()) std::cerr << "[Warning] Original face not valid " << hash(o.face) << "\t" << o.IfcGuid() << std::endl;
        check_face(o.face, o.ShellID(), o.Info(), l_warn, l_crit);
    });
#else
    for (auto &o: orig_faces) {
        BRepCheck_Analyzer A(o.face);
        if (!A.IsValid()) std::cerr << "[Warning] Original face not valid " << hash(o.face) << "\t" << o.IfcGuid() << std::endl;
        check_face(o.face, o.ShellID(), o.Info(), l_warn, l_crit);
    }
#endif

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Analyze original faces", std::to_string(orig_faces.size()));
}

void Kernel::rebuild_origfaces(std::list<oFace> &orig_faces) const {

    auto start = std::chrono::high_resolution_clock::now();

#ifdef PARALLEL_PROCESSING
    std::vector<oFace *> V;
    V.reserve(orig_faces.size());
    for (auto &o: orig_faces)
        V.push_back(&o);

#pragma omp parallel for default(none) shared(V, std::cerr) num_threads(num_threads)
    for (unsigned int i = 0; i < V.size(); i++) {
        oFace *O = V[i];
        bool success;
        TopoDS_Face f = rebuild_face(O->face, success);
        if (success) O->face = f;
        else std::cerr << "[Warning] Rebuild failed." << O->IfcGuid() << std::endl;
    }
#else
    for (auto &O: orig_faces) {
        bool success;
        TopoDS_Face f = rebuild_face(O.face, success);
        if (success) O.face = f;
        else std::cerr << "[Warning] Face could not be rebuilt." << O.IfcGuid() << std::endl;
    }
#endif

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Rebuild faces", std::to_string(orig_faces.size()));
}

void Kernel::log_products(const std::list<Product> &products) {
    unsigned int i = 0;
    for (const auto &p: products) {
        p.Log(i);
        i++;
    }
}