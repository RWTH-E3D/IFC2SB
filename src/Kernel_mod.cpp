// Copyright 2022 Eric Fichter
#include "Kernel.h"

// Concerns the ifc model, like removing and adding of objects.

bool Kernel::read_ifc_file(const std::string &ifc_path, std::unique_ptr<IfcParse::IfcFile> &model) {

    auto start = std::chrono::high_resolution_clock::now();

    model = std::make_unique<IfcParse::IfcFile>(ifc_path);

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Parse ifc file", "");

    if (!model->good()) {
        std::cout << "Unable to parse .ifc file. " << model->good().value() << "\n";
        return false;
    } else {
        return true;
    }
}

void Kernel::link_old_IfcSpace_data(std::unique_ptr<IfcParse::IfcFile> &model, IfcGeom::IteratorSettings settings, std::list<Space> &spaces, ifcspaceInfoList &old_spaces) {

    auto start = std::chrono::high_resolution_clock::now();

    if (old_spaces.empty()) return;

    // create rtree and save intersections
    auto comp = IfcSpaces_rtree(old_spaces, spaces);
    if (comp.empty()) return;

    // evaluate and filter according to overlap volumes
    link_spaces_and_IfcSpaces(spaces, comp, false);
    if (comp.empty()) return;

    // save guids of old spaces in space
    for (auto &space: spaces)
        if (comp.find(&space) != comp.end())
            for (const auto &subspace: comp[&space])
                space.old_space_guids.insert(subspace->guid);

    // save attributes of old spaces in space
    for (auto &space: spaces)
        if (comp.find(&space) != comp.end())
            for (const auto &subspace: comp[&space]) {
                add_attribut_strings(subspace->product->data().getArgument(2), "Name", space);
                add_attribut_strings(subspace->product->data().getArgument(3), "Description", space);
                add_attribut_strings(subspace->product->data().getArgument(4), "ObjectType", space);
                add_attribut_strings(subspace->product->data().getArgument(7), "LongName", space);
                add_attribut_strings(subspace->product->data().getArgument(9), "PredefinedType", space);
            }

    // link psets
    if (ifcSchema == IFC2X3)
        link_infos_of_spaces<Ifc2x3>(model, spaces, comp);
    else if (ifcSchema == IFC4)
        link_infos_of_spaces<Ifc4>(model, spaces, comp);
    else if (ifcSchema == IFC4X1)
        link_infos_of_spaces<Ifc4x1>(model, spaces, comp);
    else if (ifcSchema == IFC4X2)
        link_infos_of_spaces<Ifc4x2>(model, spaces, comp);
    else if (ifcSchema == IFC4X3_RC1)
        link_infos_of_spaces<Ifc4x3_rc1>(model, spaces, comp);
    else {
        std::cerr << "[Error] Schema not implemented yet! " << model->schema()->name() << std::endl;
        return;
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Link data of original IfcSpaces", std::to_string(old_spaces.size()) + "/" + std::to_string(old_spaces.size()));
}

void Kernel::remove_entities_from_model(std::unique_ptr<IfcParse::IfcFile> &model, const std::set<std::string> &remove_classes) {

    auto start = std::chrono::high_resolution_clock::now();

    for (const auto &C: remove_classes) {
        boost::shared_ptr<IfcEntityList> entities;
        try { entities = model->instances_by_type(C); }
        catch (...) { continue; }
        ice::remove_entities(model, entities);
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Remove entities from model", "");
}

void Kernel::remove_relationships_from_model(std::unique_ptr<IfcParse::IfcFile> &model) {

    auto start = std::chrono::high_resolution_clock::now();

    remove_empty_IfcRelDefines_from_model(model);
    remove_unused_instances_from_model(model);

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Remove relationships from model", "");
}

void Kernel::remove_empty_IfcRelDefines_from_model(std::unique_ptr<IfcParse::IfcFile> &model) {

    // remove empty IfcRelDefines that connected IfcSpaces and IfcSpaceType
    boost::shared_ptr<IfcEntityList> entities;
    try { entities = model->instances_by_type("IfcRelDefines"); }
    catch (...) {}

    if (entities->size() > 0)
        for (auto it = entities->end() - 1; it >= entities->begin(); --it) {
            IfcUtil::IfcBaseClass *const entity = *it;

            bool null_objs = entity->data().getArgument(4)->isNull(); // RelatedObjects
            bool null_type = entity->data().getArgument(5)->isNull(); // RelatingType

            if (null_objs || null_type)
                ice::remove_entity(model, entity);

            else if (!null_objs) {
                unsigned int len_objs = entity->data().getArgument(4)->size();
                if (len_objs == 0)
                    ice::remove_entity(model, entity);
            }
        }
}

void Kernel::remove_unused_instances_from_model(std::unique_ptr<IfcParse::IfcFile> &model) {

    // delete by inverses
    std::map<std::string, std::set<std::string>> M2 = {
            {"IfcProductDefinitionShape", {"ShapeOfProduct"}},
            {"IfcLocalPlacement",         {"PlacesObject"}},
            {"IfcShapeModel",             {"OfProductRepresentation", "OfShapeAspect", "RepresentationMap"}}
    };

    for (const auto &m: M2) {
        boost::shared_ptr<IfcEntityList> entities;
        try { entities = model->instances_by_type(m.first); }
        catch (...) {}

        if (entities->size() > 0)
            for (auto it = entities->end() - 1; it >= entities->begin(); --it) {
                auto *entity = (IfcUtil::IfcBaseEntity *) *it;

                bool keep = false;

                for (const auto &i: m.second) {
                    auto Inverses = entity->get_inverse(i);
                    if (Inverses->size() > 0) {
                        keep = true;
                        break;
                    }
                }
                if (!keep) ice::remove_entity(model, entity);
            }
    }

    // delete by direct attribute
    std::map<std::string, std::string> M3 = {
            {"IfcPresentationLayerAssignment", "AssignedItems"}
    };

    for (const auto &m: M3) {
        boost::shared_ptr<IfcEntityList> entities;
        try { entities = model->instances_by_type(m.first); }
        catch (...) {}

        if (entities->size() > 0)
            for (auto it = entities->end() - 1; it >= entities->begin(); --it) {
                auto *entity = (IfcUtil::IfcBaseEntity *) *it;
                auto Attribute = entity->get(m.second);
                if (Attribute->size() == 0) ice::remove_entity(model, entity);
            }
    }
}

void Kernel::write_ifc_file(std::unique_ptr<IfcParse::IfcFile> &model, const std::string &output_filename) {

    auto start = std::chrono::high_resolution_clock::now();

    std::ofstream f(output_filename);
    f << *model;

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Write Ifc file", "");
}

void Kernel::add_virtual_elements_to_model(std::unique_ptr<IfcParse::IfcFile> &model, const std::set<IfcUtil::IfcBaseEntity *> &V) {

    auto start = std::chrono::high_resolution_clock::now();

    for (const auto &v: V)
        model->addEntity(v);

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Add virtual elements to model", std::to_string(V.size()));
}

void Kernel::check_ifcproduct_pointers(std::unique_ptr<IfcParse::IfcFile> &model, std::list<Space> &spaces) {

    auto start = std::chrono::high_resolution_clock::now();

    // collect pointers
    std::set<Product *> P;
    for (auto &space: spaces) {
        for (const auto &cface: space.SecondLvl) P.insert(cface->RelProduct());
        for (const auto &cface: space.FirstLvl) P.insert(cface->RelProduct());
    }

    // find bad pointers
    for (auto &p: P) {
        IfcUtil::IfcBaseClass *entity;
        try {
            entity = model->instance_by_guid(p->guid);
        }
        catch (...) {
            std::cerr << "[Warning] No entity with guid " << p->guid << " found in ifc file. Remove pointer" << std::endl;
            p->ifcproduct = nullptr;
        }
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Check model pointers", std::to_string(spaces.size()));
}

void Kernel::revise_is_external_attribute(std::unique_ptr<IfcParse::IfcFile> &model, std::list<Space> &spaces) {

    auto start = std::chrono::high_resolution_clock::now();

    std::set<IfcUtil::IfcBaseEntity *> externals;
    std::set<IfcUtil::IfcBaseEntity *> all;

    // Is External: Indication whether the element is designed for use in the exterior (TRUE) or not (FALSE). If (TRUE) it is an external element and faces the outside of the building.
    for (const auto &space: spaces)
        if (space.is_facade)
            for (const auto &cface: space.FirstLvl)
                if (cface->IfcProduct() == nullptr) std::cerr << "[Warning] Face has no product. " << cface->Info() << std::endl;
                else externals.insert(cface->IfcProduct());
    //else
    //for (const auto &cface: space.FirstLvl)
    //         if (cface->SpaceBehind() != nullptr) { // also use space behind info. maybe clashes with ifc definition
    //            if (cface->SpaceBehind()->is_facade)
    //                externals.insert(cface->IfcProduct());
    //             else
    // internals.insert(cface->IfcProduct());
    //        }

    for (const auto &space: spaces)
        for (const auto &cface: space.FirstLvl)
            if (cface->IfcProduct() == nullptr) std::cerr << "[Warning] Face has no product. " << cface->Info() << std::endl;
            else all.insert(cface->IfcProduct());

/*    auto it = internals.begin();
    while (it != internals.end())
        if (externals.find(*it) != externals.end()) it = internals.erase(it);
        else it++;*/

    // TODO currently value is set in the pset not considering whether pset is referenced by e.g. a type (ifcwalltype). so changing the isexternal value may also apply to other objects. Instead create individual pset

    if (ifcSchema == IFC2X3) {
        revise_is_external_attribute_worker<Ifc2x3>(all, false);
        revise_is_external_attribute_worker<Ifc2x3>(externals, true);
    } else if (ifcSchema == IFC4) {
        revise_is_external_attribute_worker<Ifc4>(all, false);
        revise_is_external_attribute_worker<Ifc4>(externals, true);
    } else if (ifcSchema == IFC4X1) {
        revise_is_external_attribute_worker<Ifc4x1>(all, false);
        revise_is_external_attribute_worker<Ifc4x1>(externals, true);
    } else if (ifcSchema == IFC4X2) {
        revise_is_external_attribute_worker<Ifc4x2>(all, false);
        revise_is_external_attribute_worker<Ifc4x2>(externals, true);
    } else if (ifcSchema == IFC4X3_RC1) {
        revise_is_external_attribute_worker<Ifc4x3_rc1>(all, false);
        revise_is_external_attribute_worker<Ifc4x3_rc1>(externals, true);
    } else
        std::cerr << "[Error] Schema not implemented yet! " << model->schema()->name() << std::endl;

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Correct isExternal attribut in product psets", std::to_string(spaces.size()));
}

void Kernel::modify_view_definition(std::unique_ptr<IfcParse::IfcFile> &model) {

    std::vector<std::string> add = {"SpaceBoundary2ndLevelAddOnView"};
    std::vector<std::string> del = {"SpaceBoundary1stLevelAddOnView"};

    auto desc = model->header().file_description().description();

    for (auto &x: desc) {

        std::size_t found = x.find("ViewDefinition");

        if (found == std::string::npos) continue;

        //construct a vector of strings from original string, extract arguments from the main string
        size_t end_pos = x.find(']') - x.find('[') - 1;

        if (end_pos <= 0)
            break;

        std::string s = x.substr(x.find('[') + 1, end_pos);
        std::regex rgx("[,][\\s]");
        std::sregex_token_iterator iter(s.begin(), s.end(), rgx, -1);
        std::sregex_token_iterator end;
        std::vector<std::string> vd;

        for (; iter != end; ++iter) {
            if (iter->length() <= 0)
                continue;
            vd.push_back(*iter);
        }

        for (const auto &str: add) {
            if (std::find(vd.begin(), vd.end(), str) != vd.end())
                continue;
            vd.push_back(str);
        }

        // construct a new viewDefinition to replace the old one
        std::string viewDef;
        viewDef = "ViewDefinition [";
        bool b = false;

        for (const auto &str: vd) {

            if (std::find(del.begin(), del.end(), str) != del.end()) continue;

            if (!b)
                viewDef += str;
            else {
                viewDef += "," + str;
                b = true;
            }

        }

        viewDef += ']';
        x = viewDef;
        break;
    }

    model->header().file_description().description(desc);
}

bool Kernel::add_entities_to_model(std::unique_ptr<IfcParse::IfcFile> &model, std::list<Space> &spaces, const std::set<IfcUtil::IfcBaseEntity *> &virtual_products, const std::unordered_map<std::string, std::list<TopoDS_Face>> &shadings) {

    auto start = std::chrono::high_resolution_clock::now();

    // check if there is at least one IfcBuildingStorey in ifc file
    if (model->instances_by_type("IfcBuildingStorey")->size() == 0) {
        std::cerr << "[Error] No storeys found in ifc file!" << std::endl;
        return false;
    }

    // check if there is at least one IfcGeometricRepresentationContext for 3D models in ifc file
    IfcUtil::IfcBaseClass *IfcGeomReprContext_Model = nullptr;
    for (auto &c: *model->instances_by_type_excl_subtypes("IfcGeometricRepresentationContext")) {
        if (c->data().getArgument(1)->toString() == "'Model'" && c->data().getArgument(2)->toString() == "3") {
            IfcGeomReprContext_Model = c;
            break;
        }
    }
    if (IfcGeomReprContext_Model == nullptr) {
        std::cerr << "[Error] No valid representation context found in ifc file!" << std::endl;
        return false;
    }

    // start model enrichment based on schema
    if (ifcSchema == IFC2X3)
        add_entities_to_model_worker_Ifc2x3<Ifc2x3>(model, spaces, virtual_products, IfcGeomReprContext_Model);
    else if (ifcSchema == IFC4)
        add_entities_to_model_worker<Ifc4>(model, spaces, virtual_products, shadings, IfcGeomReprContext_Model);
    else if (ifcSchema == IFC4X1)
        add_entities_to_model_worker<Ifc4x1>(model, spaces, virtual_products, shadings, IfcGeomReprContext_Model);
    else if (ifcSchema == IFC4X2)
        add_entities_to_model_worker<Ifc4x2>(model, spaces, virtual_products, shadings, IfcGeomReprContext_Model);
    else if (ifcSchema == IFC4X3_RC1)
        add_entities_to_model_worker<Ifc4x3_rc1>(model, spaces, virtual_products, shadings, IfcGeomReprContext_Model);
    else {
        std::cerr << "[Error] Schema not implemented yet! " << model->schema()->name() << std::endl;
        return false;
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Add entities to model", "");

    return true;
}

bool Kernel::get_length_conversion_factor(std::unique_ptr<IfcParse::IfcFile> &model) {

    auto start = std::chrono::high_resolution_clock::now();

    if (ifcSchema == IFC2X3)
        get_length_conversion_factor_worker<Ifc2x3>(model);
    else if (ifcSchema == IFC4)
        get_length_conversion_factor_worker<Ifc4>(model);
    else if (ifcSchema == IFC4X1)
        get_length_conversion_factor_worker<Ifc4x1>(model);
    else if (ifcSchema == IFC4X2)
        get_length_conversion_factor_worker<Ifc4x2>(model);
    else if (ifcSchema == IFC4X3_RC1)
        get_length_conversion_factor_worker<Ifc4x3_rc1>(model);
    else {
        std::cerr << "[Error] Schema not implemented yet! " << model->schema()->name() << std::endl;
        return false;
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Get length unit conversion factor", "f: " + std::to_string(conv_fctr) + ", n: " + std::to_string(n_digits_round_double_ifc_write));

    return true;
}

bool Kernel::geom_context_info(std::unique_ptr<IfcParse::IfcFile> &model) const {

    boost::shared_ptr<IfcEntityList> IfcEntityList = nullptr;
    try { IfcEntityList = model->instances_by_type_excl_subtypes("IfcGeometricRepresentationContext"); }
    catch (...) {}

    if (IfcEntityList != nullptr)
        for (auto E: *IfcEntityList)
            std::cout << "\t#" << E->data().id() << "\t" << E->data().getArgument(1)->toString() << "\t" << E->data().getArgument(3)->toString() << "\t" << std::stod(E->data().getArgument(3)->toString()) / conv_fctr << "\n";
    else {
        std::cerr << "[Error] No context!" << std::endl;
        return false;
    }

    return true;
}

template<typename Schema>
void Kernel::add_entities_to_model_worker_Ifc2x3(std::unique_ptr<IfcParse::IfcFile> &model, std::list<Space> &spaces, const std::set<IfcUtil::IfcBaseEntity *> &virtual_products, IfcUtil::IfcBaseClass *IfcGeomReprContext_Model) {

    // set IfcGeometricRepresentationContext *************************
    typename Schema::IfcGeometricRepresentationContext *IfcGeomeReprCont_Copy = IfcGeomReprContext_Model->as<typename Schema::IfcGeometricRepresentationContext>();
    //*****************************************************************

    // calculate rotation matrices for each storey ********************
    // Because spaces are located relative to storey,
    // calculate the position by taken the storey component from the world coordinates
    IfcGeom::Kernel kernel(model.get());
    std::unordered_map<typename Schema::IfcBuildingStorey *, gp_Trsf> storeyTrsf;
    for (const auto &entity: *model->instances_by_type("IfcBuildingStorey")) {
        typename Schema::IfcBuildingStorey *IfcBuildingStorey = entity->as<typename Schema::IfcBuildingStorey>();
        gp_Trsf trsf;
        kernel.convert_placement(IfcBuildingStorey->ObjectPlacement(), trsf);
        trsf.Invert();
        storeyTrsf[IfcBuildingStorey] = trsf;
    }
    //*****************************************************************

    // add meta data **************************************************
    auto CreationDate = time(nullptr);
    typename Schema::IfcPerson *IfcPerson = ice::IfcPerson<Schema>("Eric", "Fichter");
    typename Schema::IfcOrganization *IfcOrganization = ice::IfcOrganization<Schema>("E3D - RWTH Aachen University");
    typename Schema::IfcPersonAndOrganization *IfcPersonAndOrganization = ice::IfcPersonAndOrganization<Schema>(IfcPerson, IfcOrganization);
    typename Schema::IfcApplication *IfcApplication = ice::IfcApplication<Schema>(IfcOrganization, IFC2SB_VERSION, IFC2SB_FULLNAME, IFC2SB_NAME);
    typename Schema::IfcOwnerHistory *IfcOwnerHistory = ice::IfcOwnerHistory_Ifc2x3<Schema>(IfcPersonAndOrganization, IfcApplication, IfcPersonAndOrganization, IfcApplication, CreationDate);
    model->addEntity(IfcOwnerHistory);
    //*****************************************************************

    // attributes for virtual elements ********************************
    for (const auto &v: virtual_products) {
        auto IfcVirtualElement = v->template as<typename Schema::IfcVirtualElement>();
        IfcVirtualElement->setOwnerHistory(IfcOwnerHistory);
        IfcVirtualElement->setDescription("created by original IfcSpaces");
    }
    //*****************************************************************

    // reusable data **************************************************
    std::string schema = Schema::get_schema().name();
    typename Schema::IfcCartesianPoint *IfcCartesianPoint = ice::IfcCartesianPoint<Schema>(0, 0, 0);
    typename Schema::IfcDirection *IfcDirection_Axis = ice::IfcDirection<Schema>(0, 0, 1);
    typename Schema::IfcDirection *IfcDirection_RefDirection = ice::IfcDirection<Schema>(1, 0, 0);
    typename Schema::IfcAxis2Placement *IfcAxis2Placement = ice::IfcAxis2Placement3D<Schema>(IfcCartesianPoint, IfcDirection_Axis, IfcDirection_RefDirection);
    //typename Schema::IfcAxis2Placement3D *IfcAxis2Placement3D = ice::IfcAxis2Placement3D<Schema>(IfcCartesianPoint, IfcDirection_Axis, IfcDirection_RefDirection);
    typename Schema::IfcGeometricRepresentationContext *IfcGeometricRepresentationContext = ice::IfcGeometricRepresentationContext<Schema>("Model", IfcGeomeReprCont_Copy->WorldCoordinateSystem(), 0);
    if (IfcGeomeReprCont_Copy->hasTrueNorth()) IfcGeometricRepresentationContext->setTrueNorth(IfcGeomeReprCont_Copy->TrueNorth());
    IfcGeometricRepresentationContext->setPrecision(round_ifc_write_inv);
    typename Schema::IfcGeometricRepresentationSubContext *IfcGeometricRepresentationSubContext = ice::IfcGeometricRepresentationSubContext<Schema>("Body", "Model", IfcGeometricRepresentationContext);
    typename Schema::IfcSIUnit *VolumeUnit = ice::IfcSIUnit<Schema>(Schema::IfcUnitEnum::IfcUnit_VOLUMEUNIT, Schema::IfcSIUnitName::IfcSIUnitName_CUBIC_METRE);
    typename Schema::IfcSIUnit *AreaUnit = ice::IfcSIUnit<Schema>(Schema::IfcUnitEnum::IfcUnit_AREAUNIT, Schema::IfcSIUnitName::IfcSIUnitName_SQUARE_METRE);
    //*****************************************************************

    // add spaces
    for (auto &space: spaces) {

        if (space.is_facade) continue;

        typename Schema::IfcSpace *IfcSpace = space_to_ifc_space_Ifc2x3<Schema>(space, model, IfcOwnerHistory, storeyTrsf, schema, IfcAxis2Placement, IfcGeometricRepresentationSubContext, VolumeUnit, AreaUnit);
        model->addEntity(IfcSpace);

        // add space boundaries
        gp_Trsf trsf;
        kernel.convert_placement(IfcSpace->ObjectPlacement(), trsf);
        trsf.Invert();
        for (auto &cface: space.SecondLvl) {
            typename Schema::IfcRelSpaceBoundary *IfcRelSpaceBoundary = cface_to_ifc_rel_space_boundary_Ifc2x3<Schema>(cface, IfcSpace, IfcOwnerHistory, trsf, schema);
            model->addEntity(IfcRelSpaceBoundary);
            cface->SetIfcRelSpaceBoundary(IfcRelSpaceBoundary);
        }
    }
}

template<typename Schema>
typename Schema::IfcRelSpaceBoundary *Kernel::cface_to_ifc_rel_space_boundary_Ifc2x3(cFace *cface, typename Schema::IfcSpace *IfcSpace, typename Schema::IfcOwnerHistory *IfcOwnerHistory, gp_Trsf trsf, const std::string &schema) {

    // Name
    std::string Name = "2ndLevel";

    // Description
    std::string Description = cface->SBType() == SB_TYPE_2B ? "2b" : "2a";

    // IfcElement
    typename Schema::IfcElement *IfcElement = nullptr;
    if (!cface->IsOpening()) {
        if (cface->IfcProduct() != nullptr)
            IfcElement = cface->IfcProduct()->as<typename Schema::IfcElement>();
        else
            std::cerr << "[Info] Space boundary is related to non-existent product with guid " << cface->RelProduct()->guid << "." << std::endl;
    } else
        IfcElement = cface->Ancestor()->Opening()->as<typename Schema::IfcElement>();

    // IfcConnectionGeometry
    typename Schema::IfcConnectionGeometry *IfcConnectionGeometry = create_connection_geometry<Schema>(cface, trsf, true, schema);

    // IfcInternalOrExternalEnumValue
    typename Schema::IfcInternalOrExternalEnum::Value IfcInternalOrExternalEnumValue;
    if (cface->InternalOrExternal() == SB_IE_INTERNAL) IfcInternalOrExternalEnumValue = Schema::IfcInternalOrExternalEnum::IfcInternalOrExternal_INTERNAL;
    else if (cface->InternalOrExternal() == SB_IE_EXTERNAL) IfcInternalOrExternalEnumValue = Schema::IfcInternalOrExternalEnum::IfcInternalOrExternal_EXTERNAL;
    else if (cface->InternalOrExternal() == SB_IE_EXTERNAL_EARTH) IfcInternalOrExternalEnumValue = Schema::IfcInternalOrExternalEnum::IfcInternalOrExternal_EXTERNAL;
    else if (cface->InternalOrExternal() == SB_IE_EXTERNAL_WATER) IfcInternalOrExternalEnumValue = Schema::IfcInternalOrExternalEnum::IfcInternalOrExternal_EXTERNAL;
    else if (cface->InternalOrExternal() == SB_IE_EXTERNAL_FIRE) IfcInternalOrExternalEnumValue = Schema::IfcInternalOrExternalEnum::IfcInternalOrExternal_EXTERNAL;
    else IfcInternalOrExternalEnumValue = Schema::IfcInternalOrExternalEnum::IfcInternalOrExternal_NOTDEFINED;

    // IfcPhysicalOrVirtualEnumValue
    typename Schema::IfcPhysicalOrVirtualEnum::Value IfcPhysicalOrVirtualEnumValue;
    if (cface->PhysicalOrVirtual() == SB_PV_PHYSICAL) IfcPhysicalOrVirtualEnumValue = Schema::IfcPhysicalOrVirtualEnum::IfcPhysicalOrVirtual_PHYSICAL;
    else if (cface->PhysicalOrVirtual() == SB_PV_VIRTUAL) IfcPhysicalOrVirtualEnumValue = Schema::IfcPhysicalOrVirtualEnum::IfcPhysicalOrVirtual_VIRTUAL;
    else IfcPhysicalOrVirtualEnumValue = Schema::IfcPhysicalOrVirtualEnum::IfcPhysicalOrVirtual_NOTDEFINED;

    // IfcRelSpaceBoundary
    typename Schema::IfcRelSpaceBoundary *IfcRelSpaceBoundary = ice::IfcRelSpaceBoundary_Ifc2x3<Schema>(IfcOwnerHistory, Name, Description, IfcSpace, IfcElement, IfcConnectionGeometry, IfcPhysicalOrVirtualEnumValue, IfcInternalOrExternalEnumValue);

    return IfcRelSpaceBoundary;
}

template<typename Schema>
void Kernel::revise_is_external_attribute_worker(const std::set<IfcUtil::IfcBaseEntity *> &entities, bool b) {

    auto start = std::chrono::high_resolution_clock::now();

    // set isExternal in psets
    for (const auto &product: entities) {
        auto rel_defines = product->get_inverse("IsDefinedBy");
        for (const auto &prop: *rel_defines) {

            // Only interested in property sets
            if (!prop->declaration().is("IfcRelDefinesByProperties")) continue;
            typename Schema::IfcRelDefinesByProperties *IfcRelDefinesByProperties = prop->template as<typename Schema::IfcRelDefinesByProperties>();
            auto p = IfcRelDefinesByProperties->RelatingPropertyDefinition();

            if (p) {
                if (!p->declaration().is("IfcPropertySet")) continue;

                typename Schema::IfcPropertySet *IfcPropertySet = p->template as<typename Schema::IfcPropertySet>();
                //if (IfcPropertySet->Name() != "Pset_WallCommon") continue;
                auto h = IfcPropertySet->HasProperties();

                for (const auto &v: *h) {
                    typename Schema::IfcProperty *IfcProperty = v->template as<typename Schema::IfcProperty>();
                    if (!IfcProperty->declaration().is("IfcPropertySingleValue")) continue;
                    typename Schema::IfcPropertySingleValue *IfcPropertySingleValue = IfcProperty->template as<typename Schema::IfcPropertySingleValue>();
                    if (IfcPropertySingleValue->Name() != "IsExternal") continue;
                    typename Schema::IfcBoolean *IfcBoolean = ice::IfcBoolean<Schema>(b);
                    IfcPropertySingleValue->setNominalValue(IfcBoolean);
                    std::cout << "[Info] Changed isExternal attribute to " << std::boolalpha << b << " for product " << product->data().getArgument(0)->toString() << " in pset " << IfcPropertySet->data().getArgument(0)->toString() << " (" << IfcPropertySingleValue->data().toString() << ").\n";
                }
            }
        }
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Revise IsExternal attribute", "");
}

template<typename Schema>
void Kernel::link_infos_of_spaces(std::unique_ptr<IfcParse::IfcFile> &model, std::list<Space> &spaces, std::unordered_map<Space *, std::list<ifcspaceInfoList::iterator>> &comp) {

    // store subtypes of IfcRelAssigns
    std::set<std::string> assigns = {"HasAssignments"};
    for (auto &space: spaces)
        if (comp.find(&space) != comp.end())
            for (const auto &s: comp[&space]) {
                auto IfcSpace = s->product->template as<typename Schema::IfcSpace>();
                for (const auto &assign: assigns) {
                    auto IfcRelAssigns = IfcSpace->get_inverse(assign);
                    for (const auto &IfcRelAssign: *IfcRelAssigns)
                        space.IfcRelAssigns.insert(IfcRelAssign);
                }
            }

    // store subtypes of IfcRelAssociates
    std::set<std::string> associates = {"HasAssociations"};
    for (auto &space: spaces)
        if (comp.find(&space) != comp.end())
            for (const auto &s: comp[&space]) {
                auto IfcSpace = s->product->template as<typename Schema::IfcSpace>();
                for (const auto &associate: associates) {
                    auto IfcRelAssociates = IfcSpace->get_inverse(associate);
                    for (const auto &IfcRelAssociate: *IfcRelAssociates)
                        space.IfcRelAssociates.insert(IfcRelAssociate);
                }
            }

    // store psets
    for (auto &space: spaces)
        if (comp.find(&space) != comp.end())
            for (const auto &s: comp[&space]) {
                auto IfcSpace = s->product->template as<typename Schema::IfcSpace>();
                auto IfcRelDefines = IfcSpace->get_inverse("IsDefinedBy");
                for (const auto &IfcRelDefine: *IfcRelDefines) {
                    if (!IfcRelDefine->declaration().is("IfcRelDefinesByProperties")) continue;
                    auto IfcRelDefinesByProperties = IfcRelDefine->template as<typename Schema::IfcRelDefinesByProperties>();
                    auto RelatingPropertyDefinition = IfcRelDefinesByProperties->RelatingPropertyDefinition();

                    // skip all IfcElementQuantity
                    if (RelatingPropertyDefinition->declaration().is("IfcElementQuantity")) continue;

                    // skip BaseQuantities quantity set
                    // std::string Name = RelatingPropertyDefinition->data().getArgument(2)->toString();
                    // if (Name == "'BaseQuantities'")
                    //     continue;
                    space.IfcRelDefinesByProperties.insert(IfcRelDefine);
                }
            }

    // store types
    for (auto &space: spaces)
        if (comp.find(&space) != comp.end())
            for (const auto &s: comp[&space]) {
                auto IfcSpace = s->product->template as<typename Schema::IfcSpace>();

                // 2x3 handling
                boost::shared_ptr<IfcEntityList> IfcRelDefines;
                try { IfcRelDefines = IfcSpace->get_inverse("IsTypedBy"); }
                catch (...) { IfcRelDefines = IfcSpace->get_inverse("IsDefinedBy"); }

                for (const auto &IfcRelDefine: *IfcRelDefines) {
                    if (!IfcRelDefine->declaration().is("IfcRelDefinesByType")) continue;
                    auto IfcRelDefinesByType = IfcRelDefine->template as<typename Schema::IfcRelDefinesByType>();
                    auto IfcTypeObject = IfcRelDefinesByType->RelatingType();
                    if (!IfcTypeObject->declaration().is("IfcSpaceType")) continue;

                    // remove RepresentationMaps from IfcSpaceType
                    auto IfcSpaceType = IfcTypeObject->template as<typename Schema::IfcSpaceType>();
                    typedef IfcTemplatedEntityList<typename Schema::IfcRepresentationMap> q;
                    typename q::ptr m(new q);
                    IfcSpaceType->setRepresentationMaps(m);

                    if (space.IfcRelDefinesByType == nullptr)
                        space.IfcRelDefinesByType = IfcRelDefine;
                    else
                        std::cout << "[Info] Space " << space.id << " (" << space.is_facade << ") already has IsTypedBy attribute. Maybe mapping is not 1 to 1.\n";
                }
            }

    // store ifcproducts contained in space indirectly by saving the IfcRelContainedInSpatialStructure
    boost::shared_ptr<IfcEntityList> entities;
    try { entities = model->instances_by_type("IfcRelContainedInSpatialStructure"); }
    catch (...) {}

    if (entities->size() > 0)
        for (auto it = entities->end() - 1; it >= entities->begin(); --it) {

            auto IfcRelContainedInSpatialStructure = (*it)->template as<typename Schema::IfcRelContainedInSpatialStructure>();
            auto RelatingStructure = IfcRelContainedInSpatialStructure->RelatingStructure();

            if (RelatingStructure->declaration().is("IfcSpace")) {

                std::string guid = RelatingStructure->get("GlobalId")->toString();
                guid = guid.substr(1, guid.size() - 2); // remove quote sign that ifcopenshell returns: 'dasfsfssf'

                for (auto &space: spaces)
                    for (const auto &old: space.old_space_guids)
                        if (old == guid)
                            space.IfcRelContainedInSpatialStructures.insert(IfcRelContainedInSpatialStructure);
            }
        }
}

template<typename Schema>
typename Schema::IfcRelSpaceBoundary *Kernel::cface_to_ifc_rel_space_boundary(cFace *cface, typename Schema::IfcSpaceBoundarySelect *IfcSpaceBoundarySelect, typename Schema::IfcOwnerHistory *IfcOwnerHistory, gp_Trsf trsf, bool complement_faces, bool second_lvl, const std::string &schema) {

    // IfcElement
    typename Schema::IfcElement *IfcElement = nullptr;
    if (!cface->IsOpening()) {
        if (cface->IfcProduct() != nullptr) {
            IfcElement = cface->IfcProduct()->as<typename Schema::IfcElement>();

            if (IfcElement->declaration().is("IfcVirtualElement"))
                ice::add_product_to_storey_aggregation_via_related_elements<Schema>(cface->RelSpace()->storey->as<typename Schema::IfcBuildingStorey>(), IfcElement->template as<typename Schema::IfcVirtualElement>());
        } else
            std::cerr << "[Info] Space boundary is related to non-existent product with guid " << cface->RelProduct()->guid << "." << std::endl;
    } else
        IfcElement = cface->Ancestor()->Opening()->as<typename Schema::IfcElement>();

    // IfcConnectionGeometry
    typename Schema::IfcConnectionGeometry *IfcConnectionGeometry = create_connection_geometry<Schema>(cface, trsf, complement_faces, schema);

    // IfcInternalOrExternalEnumValue
    typename Schema::IfcInternalOrExternalEnum::Value IfcInternalOrExternalEnumValue;
    if (cface->InternalOrExternal() == SB_IE_INTERNAL) IfcInternalOrExternalEnumValue = Schema::IfcInternalOrExternalEnum::IfcInternalOrExternal_INTERNAL;
    else if (cface->InternalOrExternal() == SB_IE_EXTERNAL) IfcInternalOrExternalEnumValue = Schema::IfcInternalOrExternalEnum::IfcInternalOrExternal_EXTERNAL;
    else if (cface->InternalOrExternal() == SB_IE_EXTERNAL_EARTH) IfcInternalOrExternalEnumValue = Schema::IfcInternalOrExternalEnum::IfcInternalOrExternal_EXTERNAL_EARTH;
    else if (cface->InternalOrExternal() == SB_IE_EXTERNAL_WATER) IfcInternalOrExternalEnumValue = Schema::IfcInternalOrExternalEnum::IfcInternalOrExternal_EXTERNAL_WATER;
    else if (cface->InternalOrExternal() == SB_IE_EXTERNAL_FIRE) IfcInternalOrExternalEnumValue = Schema::IfcInternalOrExternalEnum::IfcInternalOrExternal_EXTERNAL_FIRE;
    else IfcInternalOrExternalEnumValue = Schema::IfcInternalOrExternalEnum::IfcInternalOrExternal_NOTDEFINED;

    // IfcPhysicalOrVirtualEnumValue
    typename Schema::IfcPhysicalOrVirtualEnum::Value IfcPhysicalOrVirtualEnumValue;
    if (cface->PhysicalOrVirtual() == SB_PV_PHYSICAL) IfcPhysicalOrVirtualEnumValue = Schema::IfcPhysicalOrVirtualEnum::IfcPhysicalOrVirtual_PHYSICAL;
    else if (cface->PhysicalOrVirtual() == SB_PV_VIRTUAL) IfcPhysicalOrVirtualEnumValue = Schema::IfcPhysicalOrVirtualEnum::IfcPhysicalOrVirtual_VIRTUAL;
    else IfcPhysicalOrVirtualEnumValue = Schema::IfcPhysicalOrVirtualEnum::IfcPhysicalOrVirtual_NOTDEFINED;

    // Description
    std::string Description = second_lvl ? (cface->SBType() == SB_TYPE_2B ? "2b" : "2a") : "1";

    // Name
    std::string Name = "SB" + std::to_string(cface->ID());
    if (second_lvl) return ice::IfcRelSpaceBoundary2ndLevel<Schema>(IfcOwnerHistory, Name, Description, IfcSpaceBoundarySelect, IfcElement, IfcConnectionGeometry, IfcPhysicalOrVirtualEnumValue, IfcInternalOrExternalEnumValue);
    else return ice::IfcRelSpaceBoundary1stLevel<Schema>(IfcOwnerHistory, Name, Description, IfcSpaceBoundarySelect, IfcElement, IfcConnectionGeometry, IfcPhysicalOrVirtualEnumValue, IfcInternalOrExternalEnumValue);
    //  std::string Name = "2ndLevel";
    //  return ice::IfcRelSpaceBoundary<Schema>(IfcOwnerHistory, Name, Description, IfcSpaceBoundarySelect, IfcElement, IfcConnectionGeometry, IfcPhysicalOrVirtualEnumValue, IfcInternalOrExternalEnumValue);
}

template<typename Schema>
typename Schema::IfcRelSpaceBoundary *Kernel::cface_to_ifc_rel_space_boundary_shading(cFace *cface, typename Schema::IfcSpaceBoundarySelect *IfcSpaceBoundarySelect, typename Schema::IfcOwnerHistory *IfcOwnerHistory, gp_Trsf trsf, const std::string &schema, std::string description) {

    typename Schema::IfcElement *IfcElement = nullptr;
    typename Schema::IfcConnectionGeometry *IfcConnectionGeometry = create_connection_geometry<Schema>(cface, trsf, false, schema);
    typename Schema::IfcInternalOrExternalEnum::Value IfcInternalOrExternalEnumValue = Schema::IfcInternalOrExternalEnum::IfcInternalOrExternal_EXTERNAL;
    typename Schema::IfcPhysicalOrVirtualEnum::Value IfcPhysicalOrVirtualEnumValue = Schema::IfcPhysicalOrVirtualEnum::IfcPhysicalOrVirtual_PHYSICAL;
    std::string Description = description;
    std::string Name = "SB" + std::to_string(cface->ID());
    return ice::IfcRelSpaceBoundary<Schema>(IfcOwnerHistory, Name, Description, IfcSpaceBoundarySelect, IfcElement, IfcConnectionGeometry, IfcPhysicalOrVirtualEnumValue, IfcInternalOrExternalEnumValue);
}

bool Kernel::create_connection_geometry_loop(cFace *cface, const TopoDS_ListOfShape &vertices, std::vector<gp_Pnt> &Pnts, double tol) const {

    if (vertices.Size() < 3) {
        std::cerr << "[Warning] Wire of " << cface->Info() << " consists of less than three vertices (" << vertices.Size() << ") !" << std::endl;
        return false;
    }

    for (const auto &vertex: vertices) {
        gp_Pnt temp_P = BRep_Tool::Pnt(TopoDS::Vertex(vertex));
        Pnts.emplace_back(temp_P.X() * conv_fctr, temp_P.Y() * conv_fctr, temp_P.Z() * conv_fctr); // convert e.g. from metre (standard) to desired millimetre in ifc file
    }

    if (cface->Corresponding() != nullptr) Pnts = sort_points_from_polygon_by_min(Pnts, cface); // sort so point removal will be the same for corresponding boundaries because of same starting point
    Pnts = remove_coincident_points_from_polygon(Pnts, tol);

    if (Pnts.size() < 3) {
        std::cerr << "[Warning] Wire of " << cface->Info() << " consists of less than three points after removing coincident points (" << Pnts.size() << ")!" << std::endl;
        return false;
    }

    Pnts = remove_colinear_points_from_polygon(Pnts, tol);

    if (Pnts.size() < 3) {
        std::cerr << "[Warning] Wire of " << cface->Info() << " consists of less than three points after removing collinear points (" << Pnts.size() << ")!" << std::endl;
        return false;
    }
    return true;
}

template<typename Schema>
typename IfcTemplatedEntityList<typename Schema::IfcCartesianPoint>::ptr Kernel::create_connection_geometry_cartesian_points(const std::vector<gp_Pnt> &Pnts, gp_Pnt O, gp_Dir u_plane, gp_Dir v_plane) {

    // IfcCartesianPoints
    typedef IfcTemplatedEntityList<typename Schema::IfcCartesianPoint> points_t;
    typename points_t::ptr Points(new points_t);

    // calculate u,v, doubles of points on plane, relative to origin point. Construct IfcCartesianPoints
    for (const auto &P: Pnts) {
        double u, v;
        calc_point_uv_parameters_on_plane(P, O, u_plane, v_plane, u, v);
        typename Schema::IfcCartesianPoint *IfcCartesianPoint_ = ice::IfcCartesianPoint<Schema>(round_double_for_ifc_write(u), round_double_for_ifc_write(v));
        Points->push(IfcCartesianPoint_);
    }

    // close loop by adding first point as last point
    Points->push(*(Points->begin()));

    return Points;
}

template<typename Schema>
typename Schema::IfcConnectionGeometry *Kernel::create_connection_geometry(cFace *cface, gp_Trsf trsf, bool complement_faces, const std::string &schema) {

    if (cface->face.IsNull()) {
        std::cerr << "[Warning] Failed to create connection geometry of null face!" << std::endl;
        return ice::IfcConnectionSurfaceGeometry<Schema>();
    }

    TopoDS_Face Face = TopoDS::Face(BRepBuilderAPI_Transform(cface->face, trsf).Shape()); // translate/rotate face to have coordinates relative to IfcSpace
    double tol_vertex = 1.0e-5;

    if (face_is_polygon(Face)) {

        //*****************************************************************
        // outer Wire
        TopoDS_Wire outerWire = BRepTools::OuterWire(Face);
        Topo top(outerWire);
        TopoDS_ListOfShape vertices = top.ordered_vertices_of_wire();
        if (complement_faces) vertices.Reverse(); // reverse order so that normals is not pointing into the space but into the wall

        // vertices to points
        std::vector<gp_Pnt> Pnts;
        if (!create_connection_geometry_loop(cface, vertices, Pnts, tol_vertex)) {
            std::cerr << "Skip space boundary geometry creation because of bad outer wire!" << std::endl;
            return ice::IfcConnectionSurfaceGeometry<Schema>();
        }

        // define plane resp. IfcPlane
        gp_Pnt O = Pnts[0]; // Origin of the IfcPlane
        gp_Pnt T = Pnts[1]; // Succeeding point of Origin in OuterWire
        bool was_checked;
        if (!cface->CheckFaceNormal(was_checked)) std::cerr << "[Warning] Wrong orientation? " << cface->Info() << "\t" << std::boolalpha << (cface->SBType() == SB_TYPE_1) << "\t" << (cface->NormalStatus() == FACE_NORMAL_KNOWN) << std::endl;
        gp_Dir n = complement_faces ? face_normal(Face).Reversed() : face_normal(Face); // cface->FaceNormal().Reversed() : cface->FaceNormal(); use normal of transformed face. not necessary for graph approach but for clip(space) apporach where original ifcspaces are used
        gp_Dir u_plane, v_plane; // u and v parameters (unit vectors) of the face. Vectors putting up plane
        calc_plane_uv_vectors(O, T, n, u_plane, v_plane);
        typename Schema::IfcCartesianPoint *IfcCartesianPoint = ice::IfcCartesianPoint<Schema>(round_double_for_ifc_write(Pnts[0].X()), round_double_for_ifc_write(Pnts[0].Y()), round_double_for_ifc_write(Pnts[0].Z()));     // Origin point on IfcPlane
        typename Schema::IfcDirection *IfcDirection_Axis = ice::IfcDirection<Schema>(round_double_for_ifc_write(n.X()), round_double_for_ifc_write(n.Y()), round_double_for_ifc_write(n.Z()));                                          // Face normal
        typename Schema::IfcDirection *IfcDirection_RefDirection = ice::IfcDirection<Schema>(round_double_for_ifc_write(u_plane.X()), round_double_for_ifc_write(u_plane.Y()), round_double_for_ifc_write(u_plane.Z())); // Tangent vector, vector starting from Origin to next vertex in Wire (u_plane)
        typename Schema::IfcAxis2Placement3D *IfcAxis2Placement3D = ice::IfcAxis2Placement3D<Schema>(IfcCartesianPoint, IfcDirection_Axis, IfcDirection_RefDirection);
        typename Schema::IfcPlane *IfcPlane = ice::IfcPlane<Schema>(IfcAxis2Placement3D);

        // IfcCartesianPoints
        auto Points = create_connection_geometry_cartesian_points<Schema>(Pnts, O, u_plane, v_plane);

        // IfcPolyline
        typename Schema::IfcPolyline *IfcPolyline_Outer = ice::IfcPolyline<Schema>(Points);

        // IfcCompositeCurveSegment
        // typename Schema::IfcCompositeCurveSegment *IfcCompositeCurveSegment = ice::IfcCompositeCurveSegment<Schema>(true, IfcPolyline);
        // IfcCompositeCurve
        // typedef IfcTemplatedEntityList<typename Schema::IfcCompositeCurveSegment> curve_segments;
        // typename curve_segments::ptr CurveSegments(new curve_segments);
        // CurveSegments->push(IfcCompositeCurveSegment);
        // typename Schema::IfcCompositeCurve *IfcCompositeCurve = ice::IfcCompositeCurve<Schema>(CurveSegments, false);
        //*****************************************************************

        //*****************************************************************
        // Inner Wire(s)
        typedef IfcTemplatedEntityList<typename Schema::IfcCurve> curves;
        typename curves::ptr InnerCurves(new curves);

        TopoDS_ListOfShape innerWires = Topo(Face).wires();
        innerWires.Remove(outerWire);

        for (const auto &innerWire: innerWires) {

            Topo top2(innerWire);
            TopoDS_ListOfShape vertices2 = top2.ordered_vertices_of_wire();
            if (complement_faces) vertices2.Reverse();

            std::vector<gp_Pnt> Pnts2;
            if (!create_connection_geometry_loop(cface, vertices2, Pnts2, tol_vertex)) {
                ShapeChecker(cface->face).dump_topology();
                std::cerr << "[Warning] Skip space boundary geometry creation because of bad inner wire!" << std::endl;
                return ice::IfcConnectionSurfaceGeometry<Schema>();
            }

            auto Points2 = create_connection_geometry_cartesian_points<Schema>(Pnts2, O, u_plane, v_plane);
            typename Schema::IfcPolyline *IfcPolyline = ice::IfcPolyline<Schema>(Points2);
            InnerCurves->push(IfcPolyline);

            //typename Schema::IfcCompositeCurveSegment *IfcCompositeCurveSegment = ice::IfcCompositeCurveSegment<Schema>(true, IfcPolyline);
            //typedef IfcTemplatedEntityList<typename Schema::IfcCompositeCurveSegment> curve_segments;
            //typename curve_segments::ptr CurveSegments(new curve_segments);
            //CurveSegments->push(IfcCompositeCurveSegment);
            //typename Schema::IfcCompositeCurve *IfcCompositeCurve = ice::IfcCompositeCurve<Schema>(CurveSegments, false);
            //InnerCurves->push(IfcCompositeCurve);
        }
        //*****************************************************************

        typename Schema::IfcCurveBoundedPlane *IfcCurveBoundedPlane = ice::IfcCurveBoundedPlane<Schema>(IfcPlane, IfcPolyline_Outer, InnerCurves); // IfcCurveBoundedPlane
        return ice::IfcConnectionSurfaceGeometry<Schema>(IfcCurveBoundedPlane); // IfcConnectionSurfaceGeometry

    } else {
        std::cerr << "[Info] Face is not planar or edges are curved!" << std::endl;

        // serialize face
        typename Schema::IfcProductDefinitionShape *IfcProductDefinitionShape = nullptr;
        IfcUtil::IfcBaseClass *serial = IfcGeom::serialise(schema, Face, true);

        if (serial == nullptr) {
            std::cerr << "[Warning] Failed to serialize face!" << std::endl;
            return ice::IfcConnectionSurfaceGeometry<Schema>();
        }

        IfcProductDefinitionShape = serial->as<typename Schema::IfcProductDefinitionShape>();

        if ((*IfcProductDefinitionShape->Representations()).size() > 0) {

            auto IfcRepresentation = *(*IfcProductDefinitionShape->Representations()).begin();

            if ((*IfcRepresentation->Items()).size() > 0) {
                auto IfcRepresentationItem = *(*IfcRepresentation->Items()).begin();

                if (IfcRepresentationItem->declaration().is("IfcOpenShell")) {
                    // Library "IfcOpenShell" creates an IfcOpenShell if there is a single face without shell to serialize.
                    // This IfcOpenShell can be connected to the IfcConnectionSurfaceGeometry because
                    // an IfcOpenShell can be linked to an IfcFaceBasedSurfaceModel as FaceSet
                    auto IfcOpenShell = IfcRepresentationItem->template as<typename Schema::IfcOpenShell>();

                    typedef IfcTemplatedEntityList<typename Schema::IfcConnectedFaceSet> q;
                    typename q::ptr FaceSet(new q);
                    FaceSet->push(IfcOpenShell);

                    typename Schema::IfcFaceBasedSurfaceModel *IfcFaceBasedSurfaceModel = ice::IfcFaceBasedSurfaceModel<Schema>(FaceSet);
                    return ice::IfcConnectionSurfaceGeometry<Schema>(IfcFaceBasedSurfaceModel);
                }
            }
        }

        std::cerr << "[Warning] Failed to create connection geometry of curved face!" << std::endl;
        return ice::IfcConnectionSurfaceGeometry<Schema>();
    }

}

template<typename Schema>
void Kernel::add_entities_to_model_worker(std::unique_ptr<IfcParse::IfcFile> &model, std::list<Space> &spaces, const std::set<IfcUtil::IfcBaseEntity *> &virtual_products, const std::unordered_map<std::string, std::list<TopoDS_Face>> &shadings, IfcUtil::IfcBaseClass *IfcGeomReprContext_Model) {

    // set IfcGeometricRepresentationContext *************************
    typename Schema::IfcGeometricRepresentationContext *IfcGeomeReprCont_Copy = IfcGeomReprContext_Model->as<typename Schema::IfcGeometricRepresentationContext>();
    //*****************************************************************

    // calculate rotation matrices for each storey ********************
    // Because spaces are located relative to storey,
    // calculate the position by taken the storey component from the world coordinates
    IfcGeom::Kernel kernel(model.get());
    std::unordered_map<typename Schema::IfcBuildingStorey *, gp_Trsf> storeyTrsf;
    for (const auto &entity: *model->instances_by_type("IfcBuildingStorey")) {
        typename Schema::IfcBuildingStorey *IfcBuildingStorey = entity->as<typename Schema::IfcBuildingStorey>();
        gp_Trsf trsf;
        kernel.convert_placement(IfcBuildingStorey->ObjectPlacement(), trsf);
        trsf.Invert();
        storeyTrsf[IfcBuildingStorey] = trsf;
    }
    //*****************************************************************

    // calculate rotation matrices for each building ******************
    std::unordered_map<typename Schema::IfcBuilding *, gp_Trsf> buildingTrsf;
    for (const auto &entity: *model->instances_by_type("IfcBuilding")) {
        typename Schema::IfcBuilding *IfcBuilding = entity->as<typename Schema::IfcBuilding>();
        gp_Trsf trsf;
        kernel.convert_placement(IfcBuilding->ObjectPlacement(), trsf);
        trsf.Invert();
        buildingTrsf[IfcBuilding] = trsf;
    }
    //*****************************************************************

    // add meta data **************************************************
    auto CreationDate = time(nullptr);
    typename Schema::IfcPerson *IfcPerson = ice::IfcPerson<Schema>("Eric", "Fichter");
    typename Schema::IfcOrganization *IfcOrganization = ice::IfcOrganization<Schema>(IFC2SB_ORGANIZATION);
    typename Schema::IfcPersonAndOrganization *IfcPersonAndOrganization = ice::IfcPersonAndOrganization<Schema>(IfcPerson, IfcOrganization);
    typename Schema::IfcApplication *IfcApplication = ice::IfcApplication<Schema>(IfcOrganization, IFC2SB_VERSION, IFC2SB_FULLNAME, IFC2SB_NAME);
    typename Schema::IfcOwnerHistory *IfcOwnerHistory = ice::IfcOwnerHistory<Schema>(IfcPersonAndOrganization, IfcApplication, IfcPersonAndOrganization, IfcApplication, CreationDate);
    model->addEntity(IfcOwnerHistory);
    //*****************************************************************

    // attributes for virtual elements ********************************
    for (const auto &v: virtual_products) {
        auto IfcVirtualElement = v->template as<typename Schema::IfcVirtualElement>();
        IfcVirtualElement->setOwnerHistory(IfcOwnerHistory);
        IfcVirtualElement->setDescription("created by original IfcSpaces");
    }
    //*****************************************************************

    // reusable data **************************************************
    std::string schema = Schema::get_schema().name();
    typename Schema::IfcCartesianPoint *IfcCartesianPoint = ice::IfcCartesianPoint<Schema>(0, 0, 0);
    typename Schema::IfcDirection *IfcDirection_Axis = ice::IfcDirection<Schema>(0, 0, 1);
    typename Schema::IfcDirection *IfcDirection_RefDirection = ice::IfcDirection<Schema>(1, 0, 0);
    typename Schema::IfcAxis2Placement *IfcAxis2Placement = ice::IfcAxis2Placement3D<Schema>(IfcCartesianPoint, IfcDirection_Axis, IfcDirection_RefDirection);
    //typename Schema::IfcAxis2Placement3D *IfcAxis2Placement3D = ice::IfcAxis2Placement3D<Schema>(IfcCartesianPoint, IfcDirection_Axis, IfcDirection_RefDirection);
    typename Schema::IfcGeometricRepresentationContext *IfcGeometricRepresentationContext = ice::IfcGeometricRepresentationContext<Schema>("Model", IfcGeomeReprCont_Copy->WorldCoordinateSystem(), 0);
    if (IfcGeomeReprCont_Copy->hasTrueNorth()) IfcGeometricRepresentationContext->setTrueNorth(IfcGeomeReprCont_Copy->TrueNorth());
    IfcGeometricRepresentationContext->setPrecision(round_ifc_write_inv);
    typename Schema::IfcGeometricRepresentationSubContext *IfcGeometricRepresentationSubContext = ice::IfcGeometricRepresentationSubContext<Schema>("Body", "Model", IfcGeometricRepresentationContext);
    typename Schema::IfcSIUnit *VolumeUnit = ice::IfcSIUnit<Schema>(Schema::IfcUnitEnum::IfcUnit_VOLUMEUNIT, Schema::IfcSIUnitName::IfcSIUnitName_CUBIC_METRE);
    typename Schema::IfcSIUnit *AreaUnit = ice::IfcSIUnit<Schema>(Schema::IfcUnitEnum::IfcUnit_AREAUNIT, Schema::IfcSIUnitName::IfcSIUnitName_SQUARE_METRE);
    //*****************************************************************

    // add spaces *****************************************************
    for (auto &space: spaces) {

        typename Schema::IfcSpatialElement *IfcSpatialElement;
        bool complement_faces;

        if (space.is_facade) {
            IfcSpatialElement = space_to_ifc_external_spatial_element<Schema>(space, model, IfcOwnerHistory, buildingTrsf, schema, IfcAxis2Placement, IfcGeometricRepresentationSubContext);
            complement_faces = false;
        } else {
            IfcSpatialElement = space_to_ifc_space<Schema>(space, model, IfcOwnerHistory, storeyTrsf, schema, IfcAxis2Placement, IfcGeometricRepresentationSubContext, VolumeUnit, AreaUnit);
            complement_faces = true;
        }

        model->addEntity(IfcSpatialElement);

        gp_Trsf trsf;
        kernel.convert_placement(IfcSpatialElement->ObjectPlacement(), trsf);
        trsf.Invert();

        // 2nd Level
        for (auto &cface: space.SecondLvl) {
            auto *IfcRelSpaceBoundary = cface_to_ifc_rel_space_boundary<Schema>(cface, IfcSpatialElement, IfcOwnerHistory, trsf, complement_faces, true, schema);
            model->addEntity(IfcRelSpaceBoundary);
            cface->SetIfcRelSpaceBoundary(IfcRelSpaceBoundary);
        }

        // 1st Level
        for (auto &cface: space.FirstLvl) {
            auto *IfcRelSpaceBoundary = cface_to_ifc_rel_space_boundary<Schema>(cface, IfcSpatialElement, IfcOwnerHistory, trsf, complement_faces, false, schema);
            model->addEntity(IfcRelSpaceBoundary);
            cface->SetIfcRelSpaceBoundary(IfcRelSpaceBoundary);
        }

        // Shadings
        if (space.is_facade) {
            unsigned scount = 0;
            for (auto &shadings_type: shadings)
                for (auto &shading: shadings_type.second) {
                    cFace cshading(shading, nullptr, scount);
                    auto *IfcRelSpaceBoundary = cface_to_ifc_rel_space_boundary_shading<Schema>(&cshading, IfcSpatialElement, IfcOwnerHistory, trsf, schema, shadings_type.first);
                    model->addEntity(IfcRelSpaceBoundary);
                    scount++;
                }
        }
    }
    //*****************************************************************

    //*****************************************************************
    // Add ParentBoundary attribute for IfcRelSpaceBoundaryFirstLevel
    for (auto &space: spaces) {

        // link space boundaries (windows) to parent boundaries (walls) using the ParentBoundary attribute.
        for (auto &cface: space.FirstLvl) {

            if (cface->Parent() == nullptr || cface->IfcRelSpaceBoundary() == nullptr) continue;

            auto parent = cface->Parent();

            if (parent->IfcRelSpaceBoundary() == nullptr) continue;

            auto A = cface->IfcRelSpaceBoundary()->template as<typename Schema::IfcRelSpaceBoundary1stLevel>();
            auto B = parent->IfcRelSpaceBoundary()->template as<typename Schema::IfcRelSpaceBoundary1stLevel>();
            A->setParentBoundary(B);
        }
    }
    //*****************************************************************

    //*****************************************************************
    // Add ParentBoundary- und CorrespondingBoundary attribute for IfcRelSpaceBoundarySecondLevel
    for (auto &space: spaces) {
        for (auto &cface: space.SecondLvl) {

            // link space boundaries to each other using the CorrespondingBoundary attribute.
            if (cface->Corresponding() == nullptr || cface->IfcRelSpaceBoundary() == nullptr) continue;

            auto corr = cface->Corresponding();

            if (corr->IfcRelSpaceBoundary() == nullptr) continue;

            auto A = cface->IfcRelSpaceBoundary()->template as<typename Schema::IfcRelSpaceBoundary2ndLevel>();
            auto B = corr->IfcRelSpaceBoundary()->template as<typename Schema::IfcRelSpaceBoundary2ndLevel>();

            A->setCorrespondingBoundary(B);
            B->setCorrespondingBoundary(A);
        }

        // link space boundaries (windows) to parent boundaries (walls) using the ParentBoundary attribute.
        for (auto &cface: space.SecondLvl) {

            if (cface->Parent() == nullptr || cface->IfcRelSpaceBoundary() == nullptr) continue;

            auto parent = cface->Parent();

            if (parent->IfcRelSpaceBoundary() == nullptr) continue;

            auto A = cface->IfcRelSpaceBoundary()->template as<typename Schema::IfcRelSpaceBoundary1stLevel>();
            auto B = parent->IfcRelSpaceBoundary()->template as<typename Schema::IfcRelSpaceBoundary1stLevel>();

            A->setParentBoundary(B);
        }
    }
    //*****************************************************************

    //*****************************************************************
    // Additional geometric sb data in an ifc table
    typename Schema::IfcTable *IfcTable = sb_ifctable<Schema>(spaces, AreaUnit);
    if (IfcTable->hasColumns() && IfcTable->hasRows())
        model->addEntity(IfcTable);
    else
        std::cerr << "[Warning] IfcTable is empty: " << IfcTable->data().toString() << "!" << std::endl;
    //*****************************************************************

    //*****************************************************************
    // Additional material sb data in an ifc table
    for (auto &table: sb_material_ifctable<Schema>(spaces, AreaUnit))
        model->addEntity(table);
    //*****************************************************************
}

template<typename Schema>
typename Schema::IfcTable *Kernel::sb_ifctable(std::list<Space> &spaces, typename Schema::IfcSIUnit *AreaUnit) {

    typedef IfcTemplatedEntityList<typename Schema::IfcTableRow> r;
    typename r::ptr Rows(new r);

    typedef IfcTemplatedEntityList<typename Schema::IfcTableColumn> c;
    typename c::ptr Columns(new c);

    // columns
    typename Schema::IfcReference *IfcReference = ice::IfcReference<Schema>();
    Columns->push(ice::IfcTableColumn_noUnit<Schema>("1", "guid", "guid of space boundary", IfcReference));
    Columns->push(ice::IfcTableColumn<Schema>("2", "area", "surface area", AreaUnit, IfcReference));
    Columns->push(ice::IfcTableColumn_noUnit<Schema>("3", "normal_x", "x component of normal vector", IfcReference));
    Columns->push(ice::IfcTableColumn_noUnit<Schema>("4", "normal_y", "y component of normal vector", IfcReference));
    Columns->push(ice::IfcTableColumn_noUnit<Schema>("5", "normal_z", "z component of normal vector", IfcReference));
    Columns->push(ice::IfcTableColumn_noUnit<Schema>("6", "tilt to xy-plane", "angle in degree between 0 and 180", IfcReference));

    // rows
    for (auto &space: spaces) {

        bool complement_faces = !space.is_facade;

        for (auto &cface: space.SecondLvl) {
            IfcEntityList::ptr Values(new IfcEntityList);

            if (cface->IfcRelSpaceBoundary() == nullptr) Values->push(ice::IfcText<Schema>("UNKNOWN"));
            else {
                auto A = cface->IfcRelSpaceBoundary()->template as<typename Schema::IfcRelSpaceBoundary>();
                Values->push(ice::IfcText<Schema>(A->GlobalId()));
            }

            Values->push(ice::IfcReal<Schema>(round_double_two_digits(area(cface->face))));

            gp_Dir n = complement_faces ? face_normal(cface->face).Reversed() : face_normal(cface->face);
            Values->push(ice::IfcReal<Schema>(round_double_two_digits(n.X())));
            Values->push(ice::IfcReal<Schema>(round_double_two_digits(n.Y())));
            Values->push(ice::IfcReal<Schema>(round_double_two_digits(n.Z())));

            gp_Dir n_xy(0, 0, 1);
            double a = round_double_two_digits(n_xy.Angle(n) * 180 / M_PI);
            Values->push(ice::IfcReal<Schema>(a));

            typename Schema::IfcTableRow *Row = ice::IfcTableRow<Schema>(Values, false);
            Rows->push(Row);
        }
    }

    return ice::IfcTable<Schema>("BIM2SIM sb quantities", Rows, Columns);
}

template<typename Schema>
std::list<typename Schema::IfcTable *> Kernel::sb_material_ifctable(std::list<Space> &spaces, typename Schema::IfcSIUnit *AreaUnit) {

    std::list<typename Schema::IfcTable *> tables;

    // columns
    typename Schema::IfcReference *IfcReference = ice::IfcReference<Schema>();
    typename Schema::IfcSIUnit *LengthUnit = ice::IfcSIUnit<Schema>(Schema::IfcUnitEnum::IfcUnit_LENGTHUNIT, Schema::IfcSIUnitName::IfcSIUnitName_METRE);
    typedef IfcTemplatedEntityList<typename Schema::IfcTableColumn> c;
    typename c::ptr Columns(new c);
    Columns->push(ice::IfcTableColumn_noUnit<Schema>("1", "guid", "guid of product", IfcReference));
    Columns->push(ice::IfcTableColumn<Schema>("2", "d", "thickness", LengthUnit, IfcReference));

    for (auto &space: spaces)
        for (auto &cface: space.SecondLvl)
            cface->SetIsTrash(false);

    std::set<std::pair<Product *, Product *>> collisions;

    for (auto &space: spaces)
        for (auto &cface: space.SecondLvl) {

            if (cface->IfcRelSpaceBoundary() == nullptr || cface->Corresponding() == nullptr) // cface->sb_type == SB_TYPE_2B ||
                continue;
            if (cface->Corresponding()->IfcRelSpaceBoundary() == nullptr)
                continue;
            if (cface->IsTrash() || cface->Corresponding()->IsTrash())
                continue;

            auto A = cface->IfcRelSpaceBoundary()->template as<typename Schema::IfcRelSpaceBoundary2ndLevel>();
            auto B = cface->Corresponding()->IfcRelSpaceBoundary()->template as<typename Schema::IfcRelSpaceBoundary2ndLevel>();

            // create map
            std::vector<std::pair<Product *, double>> M = cface->MaterialLayers(collisions);

            if (M.empty()) continue;

            // rows
            typedef IfcTemplatedEntityList<typename Schema::IfcTableRow> r;
            typename r::ptr Rows(new r);

            // ... add d and guid
            for (const auto &e: M) {
                IfcEntityList::ptr Values(new IfcEntityList);
                Values->push(ice::IfcText<Schema>(e.first->guid));
                Values->push(ice::IfcReal<Schema>(round_double_two_digits(e.second)));
                typename Schema::IfcTableRow *Row = ice::IfcTableRow<Schema>(Values, false);
                Rows->push(Row);
            }

            // name
            std::string name = "BIM2SIM_Layers_" + std::to_string(A->data().id()) + "_" + std::to_string(B->data().id());

            // table
            typename Schema::IfcTable *IfcTable = ice::IfcTable<Schema>(name, Rows, Columns);
            tables.push_back(IfcTable);

            cface->SetIsTrash(true);
            cface->Corresponding()->SetIsTrash(true);
        }

    if (!collisions.empty()) {
        std::cout << "[Info] Collision of products in original IFC file between:" << "\n";
        for (const auto &pair: collisions)
            std::cout << "\t" << pair.first->guid << "\t" << pair.second->guid << "\n";
    } else std::cout << "[Info] No collisions of products in original IFC file" << "\n";

    return tables;
}

template<typename Schema>
typename Schema::IfcProductDefinitionShape *Kernel::createRepresentation(const TopoDS_Solid &solid, const std::string &schema, typename Schema::IfcGeometricRepresentationSubContext *IfcGeometricRepresentationSubContext) {

    typename Schema::IfcProductDefinitionShape *IfcProductDefinitionShape = nullptr;
    IfcUtil::IfcBaseClass *serial = IfcGeom::serialise(schema, solid, false);

    if (serial == nullptr) {
        std::cerr << "[Warning] Failed to serialize space solid with normal options!" << std::endl;
        serial = IfcGeom::serialise(schema, solid, true); // try again with advanced faces
    }

    if (serial == nullptr) {
        std::cerr << "[Warning] Failed to serialize space solid with advanced options!" << std::endl;
        serial = IfcGeom::tesselate(schema, solid, 0.5); // try again with tesselation
    }

    if (serial == nullptr)
        std::cerr << "[Error] Failed to serialize space solid with tesselation!" << std::endl;
    else {
        IfcProductDefinitionShape = serial->as<typename Schema::IfcProductDefinitionShape>();

        // link IfcGeometricRepresentationSubContext to IfcRepresentation
        typename Schema::IfcRepresentation *IfcRepresentation = *IfcProductDefinitionShape->Representations()->begin();
        IfcRepresentation->setContextOfItems(IfcGeometricRepresentationSubContext);
    }

    return IfcProductDefinitionShape;
}

template<typename Schema>
typename Schema::IfcSpace *
Kernel::space_to_ifc_space(Space &space, std::unique_ptr<IfcParse::IfcFile> &model, typename Schema::IfcOwnerHistory *IfcOwnerHistory, std::unordered_map<typename Schema::IfcBuildingStorey *, gp_Trsf> &storeyTrsf, const std::string &schema, typename Schema::IfcAxis2Placement *IfcAxis2Placement,
                           typename Schema::IfcGeometricRepresentationSubContext *IfcGeometricRepresentationSubContext, typename Schema::IfcSIUnit *VolumeUnit, typename Schema::IfcSIUnit *AreaUnit) {

    std::string Name = std::to_string(space.id);
    std::string LongName = "Space_" + Name;
    double ElevationWithFlooring = round_double_two_digits(space.floor_elevation) * conv_fctr;

    // create IfcSpace
    typename Schema::IfcSpace *IfcSpace = ice::IfcSpace<Schema>(IfcOwnerHistory, Name, LongName, ElevationWithFlooring);

    // find relating IfcBuildingStorey
    space.storey = space_storey<Schema>(space.FirstLvl, *model->instances_by_type("IfcBuildingStorey")->begin());
    typename Schema::IfcBuildingStorey *IfcBuildingStorey = space.storey->as<typename Schema::IfcBuildingStorey>();

    // add IfcSpace to IfcBuildingStorey's IfcRelAggregates
    ice::add_product_to_storey_aggregation_via_related_objects<Schema>(IfcBuildingStorey, IfcSpace);

    // IfcObjectPlacement (IfcLocalPlacement)
    typename Schema::IfcLocalPlacement *IfcLocalPlacement_Storey = IfcBuildingStorey->ObjectPlacement()->template as<typename Schema::IfcLocalPlacement>();
    typename Schema::IfcObjectPlacement *IfcObjectPlacement = ice::IfcLocalPlacement<Schema>(IfcLocalPlacement_Storey, IfcAxis2Placement);
    IfcSpace->setObjectPlacement(IfcObjectPlacement); // link to IfcSpace. If the PlacementRelTo is not given, then the IfcProduct is placed absolutely within the world coordinate system. When commenting this line, also comment trsf of shell

    //****************************************************************************************************
    // create IfcProductDefinitionShape (and IfcShapeRepresentation)
    TopoDS_Shell shell = TopoDS::Shell(BRepBuilderAPI_Transform(space.shell, storeyTrsf[IfcBuildingStorey]).Shape()); // translate/rotate shape to have coordinates relative to IfcBuildingStorey
    TopoDS_Solid solid = BRepBuilderAPI_MakeSolid(shell).Solid(); // create solid so it's not exported as class IfcOpenShell
    typename Schema::IfcProductDefinitionShape *IfcProductDefinitionShape = createRepresentation<Schema>(solid, schema, IfcGeometricRepresentationSubContext);

    if (IfcProductDefinitionShape != nullptr) {
        IfcSpace->setRepresentation(IfcProductDefinitionShape); // link to IfcSpace
        round_ifc_cartesian_points_traversal<Schema>(model, IfcSpace->Representation()); // round coordinates of IfcCartesianPoints linked to Representation of IfcSpace
    }
    //****************************************************************************************************

    //****************************************************************************************************
    // create QuantitySet (IfcElementQuantity and IfcRelDefinesByProperties)
    space.CalcQuantities();
    typename Schema::IfcQuantityVolume *NetVolume = ice::IfcQuantityVolume<Schema>("NetVolume", VolumeUnit, round_double_two_digits(space.NetVolume));
    typename Schema::IfcQuantityArea *NetFloorArea = ice::IfcQuantityArea<Schema>("NetFloorArea", AreaUnit, round_double_two_digits(space.NetFloorArea));
    typename Schema::IfcQuantityArea *NetWallArea = ice::IfcQuantityArea<Schema>("NetWallArea", AreaUnit, round_double_two_digits(space.NetWallArea));
    typename Schema::IfcQuantityArea *NetCeilingArea = ice::IfcQuantityArea<Schema>("NetCeilingArea", AreaUnit, round_double_two_digits(space.NetCeilingArea));

    typedef IfcTemplatedEntityList<typename Schema::IfcPhysicalQuantity> q;
    typename q::ptr Quantities(new q);
    Quantities->push(NetVolume);
    Quantities->push(NetFloorArea);
    Quantities->push(NetWallArea);
    Quantities->push(NetCeilingArea);

    typename Schema::IfcElementQuantity *IfcElementQuantity = ice::IfcElementQuantity<Schema>(IfcOwnerHistory, "Qto_SpaceBaseQuantities", Quantities);
    typedef IfcTemplatedEntityList<typename Schema::IfcObjectDefinition> r;
    typename r::ptr RelatedObjects(new r);
    RelatedObjects->push(IfcSpace);
    typename Schema::IfcRelDefinesByProperties *IfcRelDefinesByProperties = ice::IfcRelDefinesByProperties<Schema>(IfcOwnerHistory, RelatedObjects, IfcElementQuantity);
    model->addEntity(IfcRelDefinesByProperties);
    //****************************************************************************************************

    //****************************************************************************************************
    // assign information of old space
    set_old_space_info<Schema>(space, IfcSpace);

    // PredefinedType currently skipped because can only be INTERNAL
    // ...

    // link relationships inherited from old spaces
    link_relationships<Schema>(space, IfcSpace);
    //****************************************************************************************************

    return IfcSpace;
}

template<typename Schema>
typename Schema::IfcSpace *
Kernel::space_to_ifc_space_Ifc2x3(Space &space, std::unique_ptr<IfcParse::IfcFile> &model, typename Schema::IfcOwnerHistory *IfcOwnerHistory, std::unordered_map<typename Schema::IfcBuildingStorey *, gp_Trsf> &storeyTrsf, const std::string &schema,
                                  typename Schema::IfcAxis2Placement *IfcAxis2Placement, typename Schema::IfcGeometricRepresentationSubContext *IfcGeometricRepresentationSubContext, typename Schema::IfcSIUnit *VolumeUnit, typename Schema::IfcSIUnit *AreaUnit) {

    // find relating IfcBuildingStorey
    space.storey = space_storey<Schema>(space.FirstLvl, *model->instances_by_type("IfcBuildingStorey")->begin());
    typename Schema::IfcBuildingStorey *IfcBuildingStorey = space.storey->as<typename Schema::IfcBuildingStorey>();

    // create IfcProductDefinitionShape (and IfcShapeRepresentation)
    TopoDS_Shell shell = TopoDS::Shell(BRepBuilderAPI_Transform(space.shell, storeyTrsf[IfcBuildingStorey]).Shape()); // translate/rotate shape to have coordinates relative to IfcBuildingStorey
    TopoDS_Solid solid = BRepBuilderAPI_MakeSolid(shell).Solid(); // create solid so it's not exported as class IfcOpenShell
    typename Schema::IfcProductDefinitionShape *IfcProductDefinitionShape = createRepresentation<Schema>(solid, schema, IfcGeometricRepresentationSubContext);

    // IfcObjectPlacement (IfcLocalPlacement)
    typename Schema::IfcLocalPlacement *IfcLocalPlacement_Storey = IfcBuildingStorey->ObjectPlacement()->template as<typename Schema::IfcLocalPlacement>();
    typename Schema::IfcObjectPlacement *IfcObjectPlacement = ice::IfcLocalPlacement<Schema>(IfcLocalPlacement_Storey, IfcAxis2Placement);

    // create IfcSpace
    std::string Name = std::to_string(space.id);
    std::string LongName = "Space_" + Name;
    double ElevationWithFlooring = round_double_two_digits(space.floor_elevation) * conv_fctr;

    typename Schema::IfcSpace *IfcSpace = ice::IfcSpace_Ifc2x3<Schema>(IfcOwnerHistory, Name, IfcObjectPlacement, IfcProductDefinitionShape, LongName, ElevationWithFlooring);

    // add IfcSpace to IfcBuildingStorey's IfcRelAggregates
    ice::add_product_to_storey_aggregation_via_related_objects<Schema>(IfcBuildingStorey, IfcSpace);

    // round coordinates of IfcCartesianPoints linked to Representation of IfcSpace
    round_ifc_cartesian_points_traversal<Schema>(model, IfcSpace->Representation());

    //****************************************************************************************************
    // create QuantitySet (IfcElementQuantity and IfcRelDefinesByProperties)
    space.CalcQuantities();
    typename Schema::IfcQuantityVolume *NetVolume = ice::IfcQuantityVolume_2x3<Schema>("NetVolume", VolumeUnit, round_double_two_digits(space.NetVolume));
    typename Schema::IfcQuantityArea *NetFloorArea = ice::IfcQuantityArea_2x3<Schema>("NetFloorArea", AreaUnit, round_double_two_digits(space.NetFloorArea));
    typename Schema::IfcQuantityArea *NetWallArea = ice::IfcQuantityArea_2x3<Schema>("NetWallArea", AreaUnit, round_double_two_digits(space.NetWallArea));
    typename Schema::IfcQuantityArea *NetCeilingArea = ice::IfcQuantityArea_2x3<Schema>("NetCeilingArea", AreaUnit, round_double_two_digits(space.NetCeilingArea));

    typedef IfcTemplatedEntityList<typename Schema::IfcPhysicalQuantity> q;
    typename q::ptr Quantities(new q);
    Quantities->push(NetVolume);
    Quantities->push(NetFloorArea);
    Quantities->push(NetWallArea);
    Quantities->push(NetCeilingArea);

    typename Schema::IfcElementQuantity *IfcElementQuantity = ice::IfcElementQuantity<Schema>(IfcOwnerHistory, "Qto_SpaceBaseQuantities", Quantities);
    typedef IfcTemplatedEntityList<typename Schema::IfcObject> r;
    typename r::ptr RelatedObjects(new r);
    RelatedObjects->push(IfcSpace);
    typename Schema::IfcRelDefinesByProperties *IfcRelDefinesByProperties = ice::IfcRelDefinesByProperties_2x3<Schema>(IfcOwnerHistory, RelatedObjects, IfcElementQuantity);
    model->addEntity(IfcRelDefinesByProperties);
    //****************************************************************************************************

    //****************************************************************************************************
    // assign information of old space
    set_old_space_info<Schema>(space, IfcSpace);

    // PredefinedType currently skipped because can only be INTERNAL
    // ...

    // link relationships inherited from old spaces
    link_relationships<Schema>(space, IfcSpace);
    //****************************************************************************************************

    return IfcSpace;
}

template<typename Schema>
void Kernel::round_ifc_cartesian_points_traversal(std::unique_ptr<IfcParse::IfcFile> &model, IfcUtil::IfcBaseClass *E) {

    // round coordinates of IfcCartesianPoints linked to Representation
    auto entities = model->traverse(E);
    for (const auto &entity: *entities) {
        if (entity->declaration().is("IfcCartesianPoint")) {
            typename Schema::IfcCartesianPoint *C = entity->template as<typename Schema::IfcCartesianPoint>();
            std::vector<double> coords(3);
            for (int i = 0; i < 3; i++)
                coords[i] = round_double_for_ifc_write(C->Coordinates()[i] * conv_fctr);
            C->setCoordinates(coords);
        }
    }
}

template<typename Schema>
typename Schema::IfcExternalSpatialElement *
Kernel::space_to_ifc_external_spatial_element(Space &space, std::unique_ptr<IfcParse::IfcFile> &model, typename Schema::IfcOwnerHistory *IfcOwnerHistory, std::unordered_map<typename Schema::IfcBuilding *, gp_Trsf> &buildingTrsf, const std::string &schema,
                                              typename Schema::IfcAxis2Placement *IfcAxis2Placement, typename Schema::IfcGeometricRepresentationSubContext *IfcGeometricRepresentationSubContext) {

    // create IfcExternalSpatialElement
    typename Schema::IfcExternalSpatialElement *IfcExternalSpatialElement = ice::IfcExternalSpatialElement<Schema>(IfcOwnerHistory, "facade", "building envelope");

    // get IfcBuilding
    space.storey = *model->instances_by_type("IfcBuilding")->begin();
    typename Schema::IfcBuilding *IfcBuilding = space.storey->as<typename Schema::IfcBuilding>();

    // add to IfcBuilding's IfcRelAggregates
    ice::add_product_to_building_aggregation_via_related_objects<Schema>(IfcBuilding, IfcExternalSpatialElement);

    // IfcObjectPlacement (IfcLocalPlacement)
    typename Schema::IfcLocalPlacement *IfcLocalPlacement_Building = IfcBuilding->ObjectPlacement()->template as<typename Schema::IfcLocalPlacement>();
    typename Schema::IfcObjectPlacement *IfcObjectPlacement = ice::IfcLocalPlacement<Schema>(IfcLocalPlacement_Building, IfcAxis2Placement);
    IfcExternalSpatialElement->setObjectPlacement(IfcObjectPlacement); // link to IfcExternalSpatialElement

    //****************************************************************************************************
    // create IfcProductDefinitionShape (and IfcShapeRepresentation)
    TopoDS_Shell shell = TopoDS::Shell(BRepBuilderAPI_Transform(space.shell.Complemented(), buildingTrsf[IfcBuilding]).Shape()); // translate/rotate shape to have coordinates relative to IfcBuilding
    TopoDS_Solid solid = BRepBuilderAPI_MakeSolid(shell).Solid(); // create solid so it's not exported as class IfcOpenShell
    typename Schema::IfcProductDefinitionShape *IfcProductDefinitionShape = createRepresentation<Schema>(solid, schema, IfcGeometricRepresentationSubContext);

    if (IfcProductDefinitionShape != nullptr) {
        IfcExternalSpatialElement->setRepresentation(IfcProductDefinitionShape); // link to IfcExternalSpatialElement
        round_ifc_cartesian_points_traversal<Schema>(model, IfcExternalSpatialElement->Representation()); // round coordinates of IfcCartesianPoints linked to Representation
    }
    //****************************************************************************************************

    return IfcExternalSpatialElement;
}

template<typename Schema>
void Kernel::set_old_space_info(Space &space, typename Schema::IfcSpace *IfcSpace) {

    //std::string LongName = "Space_" + Name;
    //    LongName = "oldSpaceGuids_"; // write guids of old IfcSpaces that were linked to new space
    //    for (const auto &guid: space.old_space_guids)
    //        if (guid == *space.old_space_guids.begin())
    //            LongName += guid;
    //        else
    //            LongName += "_" + guid;
    //    IfcSpace->setLongName(LongName);

    // attributes
    std::string Name = space.old_space_info["Name"];
    if (!Name.empty())
        IfcSpace->setName(Name);

    std::string Description = space.old_space_info["Description"];
    if (!Description.empty())
        IfcSpace->setDescription(Description);

    std::string ObjectType = space.old_space_info["ObjectType"];
    if (!ObjectType.empty())
        IfcSpace->setObjectType(ObjectType);

    std::string LongName = space.old_space_info["LongName"];
    if (!LongName.empty())
        IfcSpace->setLongName(LongName);
}

template<typename Schema>
void Kernel::link_relationships(const Space &space, typename Schema::IfcSpace *IfcSpace) {

    //TODO stuff should be saved as IfcRelDefines, IfcRelAssigns, ... instead of IfcRelDefinesByProperties, ...
    // space should be linked by general attribute setting (see below)
    // should always be the fifth argument RelatedObjects

    // link IfcRelAssigns
    auto BaseEntitiesAssigns = space.IfcRelAssigns;
    for (const auto &baseEntity: BaseEntitiesAssigns) {
        auto IfcRelAssigns = baseEntity->template as<typename Schema::IfcRelAssigns>();
        auto L = IfcRelAssigns->RelatedObjects();
        L->push(IfcSpace);
        IfcRelAssigns->setRelatedObjects(L);
    }

    // link IfcRelAssociates
    auto BaseEntitiesAssociates = space.IfcRelAssociates;
    for (const auto &baseEntity: BaseEntitiesAssociates) {
        auto IfcRelAssociates = baseEntity->template as<typename Schema::IfcRelAssociates>();
        auto L = IfcRelAssociates->RelatedObjects();
        L->push(IfcSpace);
        IfcRelAssociates->setRelatedObjects(L);
    }

    // link property sets
    auto BaseEntitiesProp = space.IfcRelDefinesByProperties;
    std::set<std::string> pset_names;
    for (const auto &baseEntity: BaseEntitiesProp) {
        auto IfcRelDefinesByProperties = baseEntity->template as<typename Schema::IfcRelDefinesByProperties>();

        // renaming
        if (IfcRelDefinesByProperties->RelatingPropertyDefinition()->declaration().is("IfcPropertySet")) {
            auto IfcPropertySet = IfcRelDefinesByProperties->RelatingPropertyDefinition()->template as<typename Schema::IfcPropertySet>();
            std::string name = IfcPropertySet->Name();
            unsigned int cname = 1;
            while (pset_names.find(name) != pset_names.end())
                name = IfcPropertySet->Name() + "_" + std::to_string(cname++);
            IfcPropertySet->setName(name);
            pset_names.insert(IfcPropertySet->Name());
        }

        auto L = IfcRelDefinesByProperties->RelatedObjects();
        L->push(IfcSpace);
        IfcRelDefinesByProperties->setRelatedObjects(L);
        // general approach somehow doesnt work. doesnt accept the templated list as argument:
/*        auto* attr = new IfcWrite::IfcWriteArgument();
        attr->set(L);
        IfcRelDefinesByProperties->data().setArgument(4, attr);*/
    }

    // link type sets
    // auto BaseEntitiesType = space.IfcRelDefinesByType;
    // for (const auto &baseEntity: BaseEntitiesType) {
    //   auto IfcRelDefinesByType = baseEntity->template as<typename Schema::IfcRelDefinesByType>();
    if (space.IfcRelDefinesByType != nullptr) {
        auto IfcRelDefinesByType = space.IfcRelDefinesByType->template as<typename Schema::IfcRelDefinesByType>();
        auto L = IfcRelDefinesByType->RelatedObjects();
        L->push(IfcSpace);
        IfcRelDefinesByType->setRelatedObjects(L);
    }
    //}

    // link containment in spatial structure
    auto BaseEntitiesContain = space.IfcRelContainedInSpatialStructures;
    for (const auto &baseEntity: BaseEntitiesContain) {
        auto IfcRelContainedInSpatialStructure = baseEntity->template as<typename Schema::IfcRelContainedInSpatialStructure>();
        IfcRelContainedInSpatialStructure->setRelatingStructure(IfcSpace);
    }
}

template<typename Schema>
IfcUtil::IfcBaseClass *Kernel::space_storey(const std::set<cFace *> &cfaces, IfcUtil::IfcBaseClass *backup_storey) {

    //TODO  multi-storey space is contained (or belongs to) the building storey at which its ground level is, but it is referenced by all the other building storeys, in which it spans.
    // get all IfcWalls in space
    std::set<typename Schema::IfcElement *> walls;
    std::set<std::string> classes;
    for (const auto &cface: cfaces) {
        IfcUtil::IfcBaseEntity *product = cface->IfcProduct();

        if (product == nullptr) {
            std::cerr << "[Info] Space boundary is related to non-existent product with guid " << cface->RelProduct()->guid << "." << std::endl;
            continue;
        }

        if (product->declaration().is("IfcWall") || product->declaration().is("IfcColumn"))
            walls.insert(product->as<typename Schema::IfcElement>());

        classes.insert(cface->IfcClass());
    }

    if (walls.size() == 0) {
        std::cerr << "[Warning] No walls found for space!";
        for (auto &c: classes) std::cerr << " " << c;
        std::cerr << std::endl;
        return backup_storey;
    }

    // get linked IfcBuildingStoreys of IfcWalls
    std::set<typename Schema::IfcBuildingStorey *> space_storeys;
    for (const auto &W: walls) {
        if ((*W->ContainedInStructure()).size() == 0) continue;
        auto IfcRelContainedInSpatialStructure = *(*W->ContainedInStructure()).begin(); // IfcRelContainedInSpatialStructure
        try { auto t = IfcRelContainedInSpatialStructure->RelatingStructure(); }
        catch (...) {
            std::cout << "[Info] No RelatingStructure found for instance! " << W->data().toString() << "\n";
            continue;
        }
        typename Schema::IfcProduct *SE = IfcRelContainedInSpatialStructure->RelatingStructure(); // In Ifc4 it's an IfcSpatialElement, in 2x3 anIfcSpatialStructureElement, use IfcProduct for schema independence
        //typename Schema::IfcSpatialElement *SE = (*W->ContainedInStructure()->begin())->RelatingStructure();
        if (SE->declaration().is("IfcBuildingStorey"))
            space_storeys.insert(SE->template as<typename Schema::IfcBuildingStorey>());
    }

    // if no storey is linked, something is strange in the ifc file
    if (space_storeys.size() == 0) {
        std::cerr << "[Warning] No Storey found for space!" << std::endl;
        return backup_storey;
    }

        // if one storey is linked, link space to this storey
    else if (space_storeys.size() == 1)
        return *space_storeys.begin();

        // if more than one storey is linked, link space with lowest z coordinate to this storey
    else {

        bool all_have_elevation = true;
        for (const auto &space_storey: space_storeys)
            if (!space_storey->hasElevation()) {
                all_have_elevation = false;
                break;
            }

        if (all_have_elevation) {
            IfcUtil::IfcBaseClass *storey = *space_storeys.begin();
            float zmin = 1e6;

            for (const auto &space_storey: space_storeys)
                if (space_storey->Elevation() < zmin) {
                    storey = space_storey;
                    zmin = space_storey->Elevation();
                }

            return storey;
        } else {
            std::cerr << "[Error] Not all walls have elevation. Not implemented yet!" << std::endl;
            return *space_storeys.begin();
        }
    }
}

template<typename Schema>
void Kernel::get_length_conversion_factor_worker(std::unique_ptr<IfcParse::IfcFile> &model) {

    // IFC2SB tool calculates everything in metres. The standard precision of OCC is 1e-6m or 1e-3mm.
    // If Ifc unit is e.g. mm, it would be incorrect to state 1e-6 in the Context and round Points to 6 digits.
    // Instead 1e-3 and 3 digits should be used, hence the if else conditions.

    conv_fctr = 1; // conversion factor from metre to <prefix>-metre

    for (auto &c: *model->instances_by_type("IfcSIUnit")) {

        typename Schema::IfcSIUnit *IfcSIUnit = c->as<typename Schema::IfcSIUnit>();

        if (IfcSIUnit->UnitType() == Schema::IfcUnitEnum::IfcUnit_LENGTHUNIT && IfcSIUnit->Name() == 15) { // 15 is "METRE"

            if (IfcSIUnit->hasPrefix()) {
                if (IfcSIUnit->Prefix() == Schema::IfcSIPrefix::IfcSIPrefix_EXA) {
                    conv_fctr = 1e-18;
                    n_digits_round_double_ifc_write = 12;
                } else if (IfcSIUnit->Prefix() == Schema::IfcSIPrefix::IfcSIPrefix_PETA) {
                    conv_fctr = 1e-15;
                    n_digits_round_double_ifc_write = 12;
                } else if (IfcSIUnit->Prefix() == Schema::IfcSIPrefix::IfcSIPrefix_TERA) {
                    conv_fctr = 1e-12;
                    n_digits_round_double_ifc_write = 12;
                } else if (IfcSIUnit->Prefix() == Schema::IfcSIPrefix::IfcSIPrefix_GIGA) {
                    conv_fctr = 1e-9;
                    n_digits_round_double_ifc_write = 12;
                } else if (IfcSIUnit->Prefix() == Schema::IfcSIPrefix::IfcSIPrefix_MEGA) {
                    conv_fctr = 1e-6;
                    n_digits_round_double_ifc_write = 12;
                } else if (IfcSIUnit->Prefix() == Schema::IfcSIPrefix::IfcSIPrefix_KILO) {
                    conv_fctr = 1e-3;
                    n_digits_round_double_ifc_write = 9;
                } else if (IfcSIUnit->Prefix() == Schema::IfcSIPrefix::IfcSIPrefix_HECTO) {
                    conv_fctr = 1e-2;
                    n_digits_round_double_ifc_write = 8;
                } else if (IfcSIUnit->Prefix() == Schema::IfcSIPrefix::IfcSIPrefix_DECA) {
                    conv_fctr = 1e-1;
                    n_digits_round_double_ifc_write = 7;
                } else if (IfcSIUnit->Prefix() == Schema::IfcSIPrefix::IfcSIPrefix_DECI) {
                    conv_fctr = 1e1;
                    n_digits_round_double_ifc_write = 5;
                } else if (IfcSIUnit->Prefix() == Schema::IfcSIPrefix::IfcSIPrefix_CENTI) {
                    conv_fctr = 1e2;
                    n_digits_round_double_ifc_write = 4;
                } else if (IfcSIUnit->Prefix() == Schema::IfcSIPrefix::IfcSIPrefix_MILLI) {
                    conv_fctr = 1e3;
                    n_digits_round_double_ifc_write = 3;
                } else if (IfcSIUnit->Prefix() == Schema::IfcSIPrefix::IfcSIPrefix_MICRO) {
                    conv_fctr = 1e6;
                    n_digits_round_double_ifc_write = 0;
                } else if (IfcSIUnit->Prefix() == Schema::IfcSIPrefix::IfcSIPrefix_NANO) {
                    conv_fctr = 1e9;
                    n_digits_round_double_ifc_write = 0;
                } else if (IfcSIUnit->Prefix() == Schema::IfcSIPrefix::IfcSIPrefix_PICO) {
                    conv_fctr = 1e12;
                    n_digits_round_double_ifc_write = 0;
                } else if (IfcSIUnit->Prefix() == Schema::IfcSIPrefix::IfcSIPrefix_FEMTO) {
                    conv_fctr = 1e15;
                    n_digits_round_double_ifc_write = 0;
                } else if (IfcSIUnit->Prefix() == Schema::IfcSIPrefix::IfcSIPrefix_ATTO) {
                    conv_fctr = 1e18;
                    n_digits_round_double_ifc_write = 0;
                } else std::cerr << "[Error] Unknown SIPrefix!" << std::endl;

                // recalculate factors
                round_double_ifc_write = pow(10, n_digits_round_double_ifc_write);
                round_ifc_write_inv = 1.0 / round_double_ifc_write;
            }

        }
        break; // stop after first occurance of length unit declaration
    }
}


bool Kernel::add_entities_to_model_space_approach(std::unique_ptr<IfcParse::IfcFile> &model, std::list<Space> &spaces, const std::set<IfcUtil::IfcBaseEntity *> &virtual_products, const std::unordered_map<std::string, std::list<TopoDS_Face>> &shadings,
                                                  std::unordered_map<Space *, ifcspaceInfo *> &space_map) {

    auto start = std::chrono::high_resolution_clock::now();

    // check if there is at least one IfcBuildingStorey in ifc file
    if (model->instances_by_type("IfcBuildingStorey")->size() == 0) {
        std::cerr << "[Error] No storeys found in ifc file!" << std::endl;
        return false;
    }

    // check if there is at least one IfcGeometricRepresentationContext for 3D models in ifc file
    IfcUtil::IfcBaseClass *IfcGeomReprContext_Model = nullptr;
    for (auto &c: *model->instances_by_type_excl_subtypes("IfcGeometricRepresentationContext")) {
        if (c->data().getArgument(1)->toString() == "'Model'" && c->data().getArgument(2)->toString() == "3") {
            IfcGeomReprContext_Model = c;
            break;
        }
    }
    if (IfcGeomReprContext_Model == nullptr) {
        std::cerr << "[Error] No valid representation context found in ifc file!" << std::endl;
        return false;
    }

    // start model enrichment based on schema
    if (ifcSchema == IFC4)
        add_entities_to_model_space_approach_worker<Ifc4>(model, spaces, virtual_products, shadings, IfcGeomReprContext_Model, space_map);
    else if (ifcSchema == IFC4X1)
        add_entities_to_model_space_approach_worker<Ifc4x1>(model, spaces, virtual_products, shadings, IfcGeomReprContext_Model, space_map);
    else if (ifcSchema == IFC4X2)
        add_entities_to_model_space_approach_worker<Ifc4x2>(model, spaces, virtual_products, shadings, IfcGeomReprContext_Model, space_map);
    else if (ifcSchema == IFC4X3_RC1)
        add_entities_to_model_space_approach_worker<Ifc4x3_rc1>(model, spaces, virtual_products, shadings, IfcGeomReprContext_Model, space_map);
    else {
        std::cerr << "[Error] Schema not implemented yet! " << model->schema()->name() << std::endl;
        return false;
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Add entities to model", "");

    return true;
}

template<typename Schema>
void Kernel::add_entities_to_model_space_approach_worker(std::unique_ptr<IfcParse::IfcFile> &model, std::list<Space> &spaces, const std::set<IfcUtil::IfcBaseEntity *> &virtual_products, const std::unordered_map<std::string, std::list<TopoDS_Face>> &shadings,
                                                         IfcUtil::IfcBaseClass *IfcGeomReprContext_Model, std::unordered_map<Space *, ifcspaceInfo *> &space_map) {

    // set IfcGeometricRepresentationContext *************************
    typename Schema::IfcGeometricRepresentationContext *IfcGeomeReprCont_Copy = IfcGeomReprContext_Model->as<typename Schema::IfcGeometricRepresentationContext>();
    //*****************************************************************

    // calculate rotation matrices for each storey ********************
    // Because spaces are located relative to storey,
    // calculate the position by taken the storey component from the world coordinates
    IfcGeom::Kernel kernel(model.get());
    std::unordered_map<typename Schema::IfcBuildingStorey *, gp_Trsf> storeyTrsf;
    for (const auto &entity: *model->instances_by_type("IfcBuildingStorey")) {
        typename Schema::IfcBuildingStorey *IfcBuildingStorey = entity->as<typename Schema::IfcBuildingStorey>();
        gp_Trsf trsf;
        kernel.convert_placement(IfcBuildingStorey->ObjectPlacement(), trsf);
        trsf.Invert();
        storeyTrsf[IfcBuildingStorey] = trsf;
    }
    //*****************************************************************

    // calculate rotation matrices for each building ******************
    std::unordered_map<typename Schema::IfcBuilding *, gp_Trsf> buildingTrsf;
    for (const auto &entity: *model->instances_by_type("IfcBuilding")) {
        typename Schema::IfcBuilding *IfcBuilding = entity->as<typename Schema::IfcBuilding>();
        gp_Trsf trsf;
        kernel.convert_placement(IfcBuilding->ObjectPlacement(), trsf);
        trsf.Invert();
        buildingTrsf[IfcBuilding] = trsf;
    }
    //*****************************************************************

    // add meta data **************************************************
    auto CreationDate = time(nullptr);
    typename Schema::IfcPerson *IfcPerson = ice::IfcPerson<Schema>("Eric", "Fichter");
    typename Schema::IfcOrganization *IfcOrganization = ice::IfcOrganization<Schema>(IFC2SB_ORGANIZATION);
    typename Schema::IfcPersonAndOrganization *IfcPersonAndOrganization = ice::IfcPersonAndOrganization<Schema>(IfcPerson, IfcOrganization);
    typename Schema::IfcApplication *IfcApplication = ice::IfcApplication<Schema>(IfcOrganization, IFC2SB_VERSION, IFC2SB_FULLNAME, IFC2SB_NAME);
    typename Schema::IfcOwnerHistory *IfcOwnerHistory = ice::IfcOwnerHistory<Schema>(IfcPersonAndOrganization, IfcApplication, IfcPersonAndOrganization, IfcApplication, CreationDate);
    model->addEntity(IfcOwnerHistory);
    //*****************************************************************

    // attributes for virtual elements ********************************
    for (const auto &v: virtual_products) {
        auto IfcVirtualElement = v->template as<typename Schema::IfcVirtualElement>();
        IfcVirtualElement->setOwnerHistory(IfcOwnerHistory);
        IfcVirtualElement->setDescription("created by original IfcSpaces");
    }
    //*****************************************************************

    // reusable data **************************************************
    std::string schema = Schema::get_schema().name();
    typename Schema::IfcCartesianPoint *IfcCartesianPoint = ice::IfcCartesianPoint<Schema>(0, 0, 0);
    typename Schema::IfcDirection *IfcDirection_Axis = ice::IfcDirection<Schema>(0, 0, 1);
    typename Schema::IfcDirection *IfcDirection_RefDirection = ice::IfcDirection<Schema>(1, 0, 0);
    typename Schema::IfcAxis2Placement *IfcAxis2Placement = ice::IfcAxis2Placement3D<Schema>(IfcCartesianPoint, IfcDirection_Axis, IfcDirection_RefDirection);
    typename Schema::IfcGeometricRepresentationContext *IfcGeometricRepresentationContext = ice::IfcGeometricRepresentationContext<Schema>("Model", IfcGeomeReprCont_Copy->WorldCoordinateSystem(), 0);
    if (IfcGeomeReprCont_Copy->hasTrueNorth()) IfcGeometricRepresentationContext->setTrueNorth(IfcGeomeReprCont_Copy->TrueNorth());
    IfcGeometricRepresentationContext->setPrecision(round_ifc_write_inv);
    typename Schema::IfcGeometricRepresentationSubContext *IfcGeometricRepresentationSubContext = ice::IfcGeometricRepresentationSubContext<Schema>("Body", "Model", IfcGeometricRepresentationContext);
    typename Schema::IfcSIUnit *VolumeUnit = ice::IfcSIUnit<Schema>(Schema::IfcUnitEnum::IfcUnit_VOLUMEUNIT, Schema::IfcSIUnitName::IfcSIUnitName_CUBIC_METRE);
    typename Schema::IfcSIUnit *AreaUnit = ice::IfcSIUnit<Schema>(Schema::IfcUnitEnum::IfcUnit_AREAUNIT, Schema::IfcSIUnitName::IfcSIUnitName_SQUARE_METRE);
    //*****************************************************************

    // add spaces *****************************************************
    for (auto &space: spaces) {

        if (space.is_facade)
            continue;

        typename Schema::IfcSpatialElement *IfcSpatialElement = space_map[&space]->product->as<typename Schema::IfcSpatialElement>();
        space.storey = ice::storey_of_object<Schema>(IfcSpatialElement);
        //std::cerr << space.storey->data().toString() << std::endl;
        bool complement_faces = true;

        gp_Trsf trsf;
        kernel.convert_placement(IfcSpatialElement->ObjectPlacement(), trsf);
        trsf.Invert();

        // 2nd Level
        for (auto &cface: space.SecondLvl) {
            auto *IfcRelSpaceBoundary = cface_to_ifc_rel_space_boundary<Schema>(cface, IfcSpatialElement, IfcOwnerHistory, trsf, complement_faces, true, schema);
            model->addEntity(IfcRelSpaceBoundary);
            cface->SetIfcRelSpaceBoundary(IfcRelSpaceBoundary);
        }
    }
    //*****************************************************************


    //*****************************************************************
    // Add ParentBoundary- und CorrespondingBoundary attribute for IfcRelSpaceBoundarySecondLevel
    for (auto &space: spaces) {
        for (auto &cface: space.SecondLvl) {

            // link space boundaries to each other using the CorrespondingBoundary attribute.
            if (cface->Corresponding() == nullptr || cface->IfcRelSpaceBoundary() == nullptr) continue;

            auto corr = cface->Corresponding();

            if (corr->IfcRelSpaceBoundary() == nullptr) continue;

            auto A = cface->IfcRelSpaceBoundary()->template as<typename Schema::IfcRelSpaceBoundary2ndLevel>();
            auto B = corr->IfcRelSpaceBoundary()->template as<typename Schema::IfcRelSpaceBoundary2ndLevel>();

            A->setCorrespondingBoundary(B);
            B->setCorrespondingBoundary(A);
        }

        // link space boundaries (windows) to parent boundaries (walls) using the ParentBoundary attribute.
        for (auto &cface: space.SecondLvl) {

            if (cface->Parent() == nullptr || cface->IfcRelSpaceBoundary() == nullptr) continue;

            auto parent = cface->Parent();

            if (parent->IfcRelSpaceBoundary() == nullptr) continue;

            auto A = cface->IfcRelSpaceBoundary()->template as<typename Schema::IfcRelSpaceBoundary1stLevel>();
            auto B = parent->IfcRelSpaceBoundary()->template as<typename Schema::IfcRelSpaceBoundary1stLevel>();

            A->setParentBoundary(B);
        }
    }
    //*****************************************************************

    //*****************************************************************
    // Additional geometric sb data in an ifc table
    typename Schema::IfcTable *IfcTable = sb_ifctable<Schema>(spaces, AreaUnit);
    if (IfcTable->hasColumns() && IfcTable->hasRows())
        model->addEntity(IfcTable);
    else
        std::cerr << "[Warning] IfcTable is empty: " << IfcTable->data().toString() << "!" << std::endl;
    //*****************************************************************

    //*****************************************************************
    // Additional material sb data in an ifc table
    for (auto &table: sb_material_ifctable<Schema>(spaces, AreaUnit))
        model->addEntity(table);
    //*****************************************************************
}
