// Copyright 2022 Eric Fichter
#include "IfcCheck.h"

IfcCheck::IfcCheck(IfcParse::IfcFile *_model, const unsigned int _num_threads) : model(_model), num_threads(_num_threads) {

    settings.set(IfcGeom::IteratorSettings::FASTER_BOOLEANS, true);
    settings.set(IfcGeom::IteratorSettings::SEW_SHELLS, true);
    settings.set(IfcGeom::IteratorSettings::USE_WORLD_COORDS, false);
    settings.set(IfcGeom::IteratorSettings::DISABLE_TRIANGULATION, true);
    settings.set(IfcGeom::IteratorSettings::DISABLE_OPENING_SUBTRACTIONS, false);

    check_model();
}

void IfcCheck::check_model() {

    std::cout << "\n### Check" << std::endl;
    std::cout << "Check model\n";

    auto start = std::chrono::high_resolution_clock::now();

    // check IfcSpaces
    check_ifc_spaces();

    // check IfcExternalSpatialElement
    check_ifc_external_spatial_elements();

    // check IfcSpaceBoundaries
    check_ifc_space_boundaries("IfcRelSpaceBoundary");
    check_ifc_space_boundaries("IfcRelSpaceBoundary1stLevel");
    check_ifc_space_boundaries("IfcRelSpaceBoundary2ndLevel");

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "\tElapsed time: " << elapsed.count() << " s\n";
}

void IfcCheck::check_ifc_spaces() {

    std::set<std::string> include_entities = {"IfcSpace"};
    std::list<IfcSpace> spaces;

    //***************************************************************
    // Shape generation
    IfcGeom::entity_filter entity_filter;
    entity_filter.include = true;
    entity_filter.traverse = false;
    entity_filter.entity_names = include_entities;
    std::vector<IfcGeom::filter_t> filter_funcs;
    filter_funcs.emplace_back(boost::ref(entity_filter));

    IfcGeom::Iterator<real_t> geom_iterator(settings, model, filter_funcs, num_threads);

    if (!geom_iterator.initialize()) {
        std::cout << "No IfcSpace geometries found." << std::endl;
        return;
    }

    do {
        IfcGeom::Element<real_t> *geom_object = geom_iterator.get();
        spaces.emplace_back(geom_object->product(), Kernel::geom_object_to_shape(geom_object), geom_object->guid());
    } while (geom_iterator.next());
    //***************************************************************


    //***************************************************************
    // Checks
    check_ifc_spaces_shell_number(spaces); // Number of shells
    check_ifc_spaces_shell(spaces); // Shell check
    evaluation(spaces, "Space"); // Evaluation
    //***************************************************************
}

void IfcCheck::check_ifc_external_spatial_elements() {

    std::set<std::string> include_entities = {"IfcExternalSpatialElement"};
    std::list<IfcExternalSpatialElement> spatials;

    //***************************************************************
    // Shape generation
    IfcGeom::entity_filter entity_filter;
    entity_filter.include = true;
    entity_filter.traverse = false;
    entity_filter.entity_names = include_entities;
    std::vector<IfcGeom::filter_t> filter_funcs;
    filter_funcs.emplace_back(boost::ref(entity_filter));

    IfcGeom::Iterator<real_t> geom_iterator(settings, model, filter_funcs, num_threads);

    if (!geom_iterator.initialize()) {
        std::cout << "No IfcExternalSpatialElement geometries found." << std::endl;
        return;
    }

    do {
        IfcGeom::Element<real_t> *geom_object = geom_iterator.get();
        spatials.emplace_back(geom_object->product(), Kernel::geom_object_to_shape(geom_object), geom_object->guid());
    } while (geom_iterator.next());
    //***************************************************************


    //***************************************************************
    // Checks
    check_ifc_spaces_shell_number(spatials); // Number of shells
    check_ifc_spaces_shell(spatials); // Shell check
    evaluation(spatials, "ExternalSpatialElement"); // Evaluation
    //***************************************************************
}

template<typename T>
void IfcCheck::check_ifc_spaces_shell_number(std::list<T> &spaces) {

    for (auto &space: spaces) {
        unsigned int n = Topo(space.shape).shells().Size();
        if (n != 1)
            space.errors.push_back(WRONG_NUMBER_OF_SHELLS);
    }
}

template<typename T>
void IfcCheck::check_ifc_spaces_shell(std::list<T> &spaces) {

    TopExp_Explorer Ex;

    for (auto &space: spaces) {

        for (Ex.Init(space.shape, TopAbs_SHELL); Ex.More(); Ex.Next()) {

            BRepCheck_Shell c = BRepCheck_Shell(TopoDS::Shell(Ex.Current()));
            BRepCheck_Status Closed = c.Closed();
            BRepCheck_Status Orientation = c.Orientation();

            if (Closed != BRepCheck_NoError)
                space.errors.push_back(NOT_CLOSED);

            if (Orientation != BRepCheck_NoError)
                space.errors.push_back(NON_CONSISTENT_ORIENTATION);

            if (Kernel::volume(Ex.Current()) < 0)
                space.errors.push_back(NEGATIVE_VOLUME);
        }
    }
}

template<typename T>
void IfcCheck::evaluation(const std::list<T> &L, const std::string &name) {

    for (const auto &O: L) {
        if (O.errors.empty())
            std::cout << "\t" << name << " " << O.guid << ": ok\n";
        else {
            std::cout << "\t" << name << " " << O.guid << ":\n";
            for (const auto &error: O.errors)
                std::cout << "\t\t" << error_to_string[error] << "\n";
        }
    }
}

void IfcCheck::check_ifc_space_boundaries(const std::string &ifcclass) {

#define IfcSchema Ifc4 // TODO ALLGM mit tyname schema

    IfcGeom::Kernel K(model);
    std::list<IfcSpaceBoundary> boundaries;

    //***************************************************************
    // Shape generation
    boost::shared_ptr<IfcEntityList> IfcEntityList = nullptr;
    try { IfcEntityList = model->instances_by_type_excl_subtypes(ifcclass); }
    catch (...) {}

    if (IfcEntityList == nullptr) {
        std::cout << "No " << ifcclass << " found." << std::endl;
        return;
    } else
        std::cout << IfcEntityList->size() << " instances of " << ifcclass << " found." << std::endl;

    for (auto e: *IfcEntityList){
        if( e->declaration().name() == "IfcRelSpaceBoundary"){
            std::cout << "[Info] IfcRelSpaceBoundary " << e->data().getArgument(0)->toString() << " will be skipped" << std::endl;
            continue;
        }
        auto IfcRelSpaceBoundary = e->as<Ifc4::IfcRelSpaceBoundary1stLevel>();
        //if (IfcRelSpaceBoundary->hasParentBoundary()) continue;

        gp_Trsf trsf;
        std::string relating_space_type;

        try { auto rs = IfcRelSpaceBoundary->RelatingSpace(); }
        catch (...) {
            std::cout << "[Warning] There is no valid ifc entity linked via RelatingSpace. SpaceBoundary " << IfcRelSpaceBoundary->GlobalId() << " will be skipped" << std::endl;
            continue;
        }

        if (IfcRelSpaceBoundary->RelatingSpace()->declaration().is("IfcSpace"))
            K.convert_placement(IfcRelSpaceBoundary->RelatingSpace()->as<IfcSchema::IfcSpace>()->ObjectPlacement(), trsf);
        else if (IfcRelSpaceBoundary->RelatingSpace()->declaration().is("IfcExternalSpatialElement"))
            K.convert_placement(IfcRelSpaceBoundary->RelatingSpace()->as<IfcSchema::IfcExternalSpatialElement>()->ObjectPlacement(), trsf);
        else {
            std::cout << "[Warning] Relating Space is no IfcSpaceBoundarySelect. SpaceBoundary " << IfcRelSpaceBoundary->GlobalId() << " will be skipped" << std::endl;
            continue;
        }

        auto ConnectionGeometry = IfcRelSpaceBoundary->ConnectionGeometry()->as<IfcSchema::IfcConnectionGeometry>();
        auto ConnectionSurfaceGeometry = ConnectionGeometry->as<IfcSchema::IfcConnectionSurfaceGeometry>();

        if (ConnectionSurfaceGeometry->hasSurfaceOnRelatedElement()) {
            std::cout << "[Warning] ConnectionSurfaceGeometry has SurfaceOnRelatedElement. SpaceBoundary " << IfcRelSpaceBoundary->GlobalId() << " will be skipped" << std::endl;
            continue;
        }

        std::vector<IfcGeom::IfcRepresentationShapeItem> Shape_List;
        try { Shape_List = K.convert(ConnectionSurfaceGeometry->SurfaceOnRelatingElement()); }
        catch (...) {
            std::cerr << "[Warning] Error in geometry. SpaceBoundary " << IfcRelSpaceBoundary->GlobalId() << " will be skipped" << std::endl;
            continue;
        }

        TopoDS_Shape S = BRepBuilderAPI_Transform(Shape_List.begin()->Shape(), trsf, false);
        boundaries.emplace_back(IfcRelSpaceBoundary, IfcRelSpaceBoundary->RelatingSpace(), S, IfcRelSpaceBoundary->GlobalId());
    }
    //***************************************************************


    //***************************************************************
    // Checks
    check_ifc_spaceboundaries_shell(boundaries);
    evaluation(boundaries, ifcclass); // Evaluation
    //***************************************************************

}

void IfcCheck::check_ifc_spaceboundaries_shell(std::list<IfcSpaceBoundary> &boundaries) {

    std::unordered_map<IfcUtil::IfcBaseClass *, std::list<IfcSpaceBoundary *>> C;

    for (auto &boundary: boundaries)
        C[boundary.space].push_back(&boundary);

    TopExp_Explorer Ex;

    for (auto &item: C) {

        // Fuse
        BOPAlgo_Builder B;
        B.SetFuzzyValue(2.0e-5);
        B.SetGlue(BOPAlgo_GlueShift);
        B.SetRunParallel(true);
        for (auto &boundary: item.second)
            B.AddArgument(boundary->shape);
        B.Perform();
        TopoDS_Shape Shp = B.Shape();

        // Sew
        BRepBuilderAPI_Sewing sew;
        sew.SetTolerance(1.0e-5);
        for (auto &F: Topo(Shp).faces())
            sew.Add(F);
        sew.Perform();
        Shp = sew.SewedShape();

        // Check shells
        for (Ex.Init(Shp, TopAbs_SHELL); Ex.More(); Ex.Next()) {

            BRepCheck_Shell c = BRepCheck_Shell(TopoDS::Shell(Ex.Current()));
            BRepCheck_Status Closed = c.Closed();
            BRepCheck_Status Orientation = c.Orientation();

            for (auto &boundary: item.second) {
                if (Closed != BRepCheck_NoError)
                    boundary->errors.push_back(NOT_CLOSED);

                if (Orientation != BRepCheck_NoError)
                    boundary->errors.push_back(NON_CONSISTENT_ORIENTATION);

                if (Kernel::volume(Ex.Current()) < 0)
                    boundary->errors.push_back(NEGATIVE_VOLUME);

                if (Kernel::area(boundary->shape) < 0.0025)
                    boundary->errors.push_back(SMALL_AREA);
            }
        }
    }
}
