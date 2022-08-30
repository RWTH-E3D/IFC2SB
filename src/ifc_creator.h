// Copyright 2022 Eric Fichter
#ifndef IFC_CREATOR_H
#define IFC_CREATOR_H

#include "headers.h"

namespace ice {

    typedef IfcParse::IfcGlobalId guid;
    typedef std::string S;
    boost::none_t const null = boost::none;

    template<typename Schema>
    typename Schema::IfcPerson *IfcPerson(std::string GivenName, std::string FamilyName) {

        return new typename Schema::IfcPerson(
                null,       // Identification
                FamilyName, // FamilyName
                GivenName,  // GivenName
                null,       // MiddleNames
                null,       // PrefixTitles
                null,       // SuffixTitles
                null,       // Roles
                null        // Addresses
        );

    }

    template<typename Schema>
    typename Schema::IfcOrganization *IfcOrganization(const std::string &Name) {

        return new typename Schema::IfcOrganization(
                null,           // Identification
                Name,           // Name
                null,           // Description
                null,           // Roles
                null            // Addresses
        );

    }

    template<typename Schema>
    typename Schema::IfcPersonAndOrganization *IfcPersonAndOrganization(typename Schema::IfcPerson *ThePerson, typename Schema::IfcOrganization *TheOrganization) {

        return new typename Schema::IfcPersonAndOrganization(
                ThePerson,          // ThePerson
                TheOrganization,    // TheOrganization
                null                // Roles
        );

    }

    template<typename Schema>
    typename Schema::IfcApplication *IfcApplication(typename Schema::IfcOrganization *ApplicationDeveloper, const std::string &Version, const std::string &ApplicationFullName, const std::string &ApplicationIdentifier) {

        return new typename Schema::IfcApplication(
                ApplicationDeveloper,       // ApplicationDeveloper
                Version,                    // Version
                ApplicationFullName,        // ApplicationFullName
                ApplicationIdentifier       // ApplicationIdentifier
        );

    }

    template<typename Schema>
    typename Schema::IfcOwnerHistory *IfcOwnerHistory(typename Schema::IfcPersonAndOrganization *OwningUser, typename Schema::IfcApplication *OwningApplication, typename Schema::IfcPersonAndOrganization *LastModifyingUser, typename Schema::IfcApplication *LastModifyingApplication,
                                                      int CreationDate) {

        return new typename Schema::IfcOwnerHistory(
                OwningUser,                 // OwningUser
                OwningApplication,          // OwningApplication
                null,                       // State
                null,                       // ChangeAction
                CreationDate,               // LastModifiedDate
                LastModifyingUser,          // LastModifyingUser
                LastModifyingApplication,   // LastModifyingApplication
                CreationDate                // CreationDate
        );

    }

    template<typename Schema>
    typename Schema::IfcOwnerHistory *IfcOwnerHistory_Ifc2x3(typename Schema::IfcPersonAndOrganization *OwningUser, typename Schema::IfcApplication *OwningApplication, typename Schema::IfcPersonAndOrganization *LastModifyingUser, typename Schema::IfcApplication *LastModifyingApplication,
                                                             int CreationDate) {

        return new typename Schema::IfcOwnerHistory(
                OwningUser,                                         // OwningUser
                OwningApplication,                                  // OwningApplication
                null,                                               // State
                Schema::IfcChangeActionEnum::IfcChangeAction_ADDED, // ChangeAction
                CreationDate,                                       // LastModifiedDate
                LastModifyingUser,                                  // LastModifyingUser
                LastModifyingApplication,                           // LastModifyingApplication
                CreationDate                                        // CreationDate
        );

    }

    template<typename Schema>
    typename Schema::IfcSpace *IfcSpace(typename Schema::IfcOwnerHistory *OwnerHistory, const std::string &Name, const std::string &LongName, const double &ElevationWithFlooring) {

        return new typename Schema::IfcSpace(
                guid(),                                                             // GlobalId
                OwnerHistory,                                                       // OwnerHistory
                Name,                                                               // Name
                null,                                                               // Description
                null,                                                               // ObjectType
                0,                                                                  // ObjectPlacement
                0,                                                                  // Representation
                LongName,                                                           // LongName
                Schema::IfcElementCompositionEnum::IfcElementComposition_ELEMENT,   // CompositionType
                Schema::IfcSpaceTypeEnum::IfcSpaceType_INTERNAL,                    // PredefinedType
                ElevationWithFlooring                                               // ElevationWithFlooring
        );

    }

    template<typename Schema>
    typename Schema::IfcSpace *IfcSpace_Ifc2x3(typename Schema::IfcOwnerHistory *OwnerHistory, const std::string &Name, typename Schema::IfcObjectPlacement *ObjectPlacement, typename Schema::IfcProductRepresentation *Representation, const std::string &LongName, const double &ElevationWithFlooring) {

        return new typename Schema::IfcSpace(
                guid(),                                                             // GlobalId
                OwnerHistory,                                                       // OwnerHistory
                Name,                                                               // Name
                null,                                                               // Description
                null,                                                               // ObjectType
                ObjectPlacement,                                                    // ObjectPlacement
                Representation,                                                     // Representation
                LongName,                                                           // LongName
                Schema::IfcElementCompositionEnum::IfcElementComposition_ELEMENT,   // CompositionType
                Schema::IfcInternalOrExternalEnum::IfcInternalOrExternal_INTERNAL,  // InteriorOrExteriorSpace
                ElevationWithFlooring                                               // ElevationWithFlooring
        );

    }

    template<typename Schema>
    typename Schema::IfcGeometricRepresentationContext *IfcGeometricRepresentationContext(const std::string &ContextType, typename Schema::IfcAxis2Placement *WorldCoordinateSystem, typename Schema::IfcDirection *TrueNorth) {

        return new typename Schema::IfcGeometricRepresentationContext(
                null,                                                   // ContextIdentifier
                ContextType,                                            // ContextType
                3,                                                      // CoordinateSpaceDimension
                null,                                                   // Precision
                WorldCoordinateSystem,                                  // WorldCoordinateSystem
                TrueNorth                                               // TrueNorth
        );

    }

    template<typename Schema>
    typename Schema::IfcGeometricRepresentationSubContext *IfcGeometricRepresentationSubContext(const std::string &ContextIdentifier, const std::string &ContextType, typename Schema::IfcGeometricRepresentationContext *ParentContext) {

        return new typename Schema::IfcGeometricRepresentationSubContext(
                ContextIdentifier,                                                      // ContextIdentifier
                ContextType,                                                            // ContextType
                ParentContext,                                                          // ParentContext
                null,                                                                   // TargetScale
                Schema::IfcGeometricProjectionEnum::IfcGeometricProjection_MODEL_VIEW,  // TargetView
                null                                                                    // UserDefinedTargetView
        );

    }

    template<typename Schema>
    typename Schema::IfcLocalPlacement *IfcLocalPlacement(typename Schema::IfcObjectPlacement *PlacementRelTo, typename Schema::IfcAxis2Placement *RelativePlacement) {

        return new typename Schema::IfcLocalPlacement(
                PlacementRelTo,     // PlacementRelTo
                RelativePlacement   // RelativePlacement
        );

    }

    template<typename Schema>
    typename Schema::IfcCartesianPoint *IfcCartesianPoint(double x, double y, double z) {

        std::vector<double> v = {x, y, z};
        return new typename Schema::IfcCartesianPoint(v);

    }

    template<typename Schema>
    typename Schema::IfcCartesianPoint *IfcCartesianPoint(double x, double y) {

        std::vector<double> v = {x, y};
        return new typename Schema::IfcCartesianPoint(v);

    }

    template<typename Schema>
    typename Schema::IfcDirection *IfcDirection(double x, double y, double z) {

        std::vector<double> v = {x, y, z};
        return new typename Schema::IfcDirection(v);

    }

    template<typename Schema>
    typename Schema::IfcAxis2Placement3D *IfcAxis2Placement3D(typename Schema::IfcCartesianPoint *Location, typename Schema::IfcDirection *Axis, typename Schema::IfcDirection *RefDirection) {

        return new typename Schema::IfcAxis2Placement3D(
                Location,       // Location
                Axis,           // Axis
                RefDirection    // RefDirection
        );

    }

    template<typename Schema>
    typename Schema::IfcRelSpaceBoundary *
    IfcRelSpaceBoundary(typename Schema::IfcOwnerHistory *OwnerHistory, const std::string &Name, const std::string &Description, typename Schema::IfcSpaceBoundarySelect *RelatingSpace, typename Schema::IfcElement *RelatedBuildingElement, typename Schema::IfcConnectionGeometry *ConnectionGeometry,
                        typename Schema::IfcPhysicalOrVirtualEnum::Value PhysicalOrVirtualBoundary, typename Schema::IfcInternalOrExternalEnum::Value InternalOrExternalBoundary) {

        return new typename Schema::IfcRelSpaceBoundary(
                guid(),                         // GlobalId
                OwnerHistory,                   // OwnerHistory
                Name,                           // Name
                Description,                    // Description
                RelatingSpace,                  // RelatingSpace
                RelatedBuildingElement,         // RelatedBuildingElement
                ConnectionGeometry,             // ConnectionGeometry
                PhysicalOrVirtualBoundary,      // PhysicalOrVirtualBoundary
                InternalOrExternalBoundary      // InternalOrExternalBoundary
        );

    }

    template<typename Schema>
    typename Schema::IfcRelSpaceBoundary *
    IfcRelSpaceBoundary_Ifc2x3(typename Schema::IfcOwnerHistory *OwnerHistory, const std::string &Name, const std::string &Description, typename Schema::IfcSpace *RelatingSpace, typename Schema::IfcElement *RelatedBuildingElement, typename Schema::IfcConnectionGeometry *ConnectionGeometry,
                               typename Schema::IfcPhysicalOrVirtualEnum::Value PhysicalOrVirtualBoundary, typename Schema::IfcInternalOrExternalEnum::Value InternalOrExternalBoundary) {

        return new typename Schema::IfcRelSpaceBoundary(
                guid(),                         // GlobalId
                OwnerHistory,                   // OwnerHistory
                Name,                           // Name
                Description,                    // Description
                RelatingSpace,                  // RelatingSpace
                RelatedBuildingElement,         // RelatedBuildingElement
                ConnectionGeometry,             // ConnectionGeometry
                PhysicalOrVirtualBoundary,      // PhysicalOrVirtualBoundary
                InternalOrExternalBoundary      // InternalOrExternalBoundary
        );

    }

    template<typename Schema>
    typename Schema::IfcRelSpaceBoundary2ndLevel *
    IfcRelSpaceBoundary2ndLevel(typename Schema::IfcOwnerHistory *OwnerHistory, const std::string &Name, const std::string &Description, typename Schema::IfcSpaceBoundarySelect *RelatingSpace, typename Schema::IfcElement *RelatedBuildingElement,
                                typename Schema::IfcConnectionGeometry *ConnectionGeometry,
                                typename Schema::IfcPhysicalOrVirtualEnum::Value PhysicalOrVirtualBoundary, typename Schema::IfcInternalOrExternalEnum::Value InternalOrExternalBoundary) {

        return new typename Schema::IfcRelSpaceBoundary2ndLevel(
                guid(),                         // GlobalId
                OwnerHistory,                   // OwnerHistory
                Name,                           // Name
                Description,                    // Description
                RelatingSpace,                  // RelatingSpace
                RelatedBuildingElement,         // RelatedBuildingElement
                ConnectionGeometry,             // ConnectionGeometry
                PhysicalOrVirtualBoundary,      // PhysicalOrVirtualBoundary
                InternalOrExternalBoundary,     // InternalOrExternalBoundary
                0,                              // ParentBoundary
                0                               // CorrespondingBoundary
        );

    }

    template<typename Schema>
    typename Schema::IfcRelSpaceBoundary1stLevel *
    IfcRelSpaceBoundary1stLevel(typename Schema::IfcOwnerHistory *OwnerHistory, const std::string &Name, const std::string &Description, typename Schema::IfcSpaceBoundarySelect *RelatingSpace, typename Schema::IfcElement *RelatedBuildingElement,
                                typename Schema::IfcConnectionGeometry *ConnectionGeometry,
                                typename Schema::IfcPhysicalOrVirtualEnum::Value PhysicalOrVirtualBoundary, typename Schema::IfcInternalOrExternalEnum::Value InternalOrExternalBoundary) {

        return new typename Schema::IfcRelSpaceBoundary1stLevel(
                guid(),                         // GlobalId
                OwnerHistory,                   // OwnerHistory
                Name,                           // Name
                Description,                    // Description
                RelatingSpace,                  // RelatingSpace
                RelatedBuildingElement,         // RelatedBuildingElement
                ConnectionGeometry,             // ConnectionGeometry
                PhysicalOrVirtualBoundary,      // PhysicalOrVirtualBoundary
                InternalOrExternalBoundary,     // InternalOrExternalBoundary
                0                              // ParentBoundary
        );

    }

    template<typename Schema>
    typename Schema::IfcFaceBasedSurfaceModel *IfcFaceBasedSurfaceModel(typename IfcTemplatedEntityList<typename Schema::IfcConnectedFaceSet>::ptr s) {

        return new typename Schema::IfcFaceBasedSurfaceModel(
                s   // v1_FbsmFaces
        );

    }

    template<typename Schema>
    typename Schema::IfcConnectionSurfaceGeometry *IfcConnectionSurfaceGeometry() {

        return new typename Schema::IfcConnectionSurfaceGeometry(0, 0);

    }


    template<typename Schema>
    typename Schema::IfcConnectionSurfaceGeometry *IfcConnectionSurfaceGeometry(typename Schema::IfcSurfaceOrFaceSurface *SurfaceOnRelatingElement) {

        return new typename Schema::IfcConnectionSurfaceGeometry(
                SurfaceOnRelatingElement,   // SurfaceOnRelatingElement
                0                           // SurfaceOnRelatedElement
        );

    }

    template<typename Schema>
    typename Schema::IfcPolyline *IfcPolyline(typename IfcTemplatedEntityList<typename Schema::IfcCartesianPoint>::ptr Points) {

        return new typename Schema::IfcPolyline(
                Points       // Points
        );

    }

    template<typename Schema>
    typename Schema::IfcCompositeCurveSegment *IfcCompositeCurveSegment(bool SameSense, typename Schema::IfcCurve *ParentCurve) {

        return new typename Schema::IfcCompositeCurveSegment(
                Schema::IfcTransitionCode::IfcTransitionCode_CONTINUOUS,    // Transition
                SameSense,                                                  // SameSense
                ParentCurve                                                 // ParentCurve
        );

    }

    template<typename Schema>
    typename Schema::IfcCompositeCurve *IfcCompositeCurve(typename IfcTemplatedEntityList<typename Schema::IfcCompositeCurveSegment>::ptr Segments, bool SelfIntersect) {

        return new typename Schema::IfcCompositeCurve(
                Segments,       // Segments
                SelfIntersect   // SelfIntersect
        );

    }

    template<typename Schema>
    typename Schema::IfcCurveBoundedPlane *IfcCurveBoundedPlane(typename Schema::IfcPlane *BasisSurface, typename Schema::IfcCurve *OuterBoundary, typename IfcTemplatedEntityList<typename Schema::IfcCurve>::ptr InnerBoundaries) {

        return new typename Schema::IfcCurveBoundedPlane(
                BasisSurface,       // BasisSurface
                OuterBoundary,   // OuterBoundary
                InnerBoundaries   // InnerBoundaries
        );

    }

    template<typename Schema>
    typename Schema::IfcPlane *IfcPlane(typename Schema::IfcAxis2Placement3D *Position) {

        return new typename Schema::IfcPlane(
                Position      // Position
        );

    }

    template<typename Schema>
    typename Schema::IfcExternalSpatialElement *IfcExternalSpatialElement(typename Schema::IfcOwnerHistory *OwnerHistory, const std::string &Name, const std::string &Description) {

        return new typename Schema::IfcExternalSpatialElement(
                guid(),                                                                             // GlobalId
                OwnerHistory,                                                                       // OwnerHistory
                Name,                                                                               // Name
                Description,                                                                        // Description
                null,                                                                               // ObjectType
                0,                                                                                  // ObjectPlacement
                0,                                                                                  // Representation
                null,                                                                               // LongName
                Schema::IfcExternalSpatialElementTypeEnum::IfcExternalSpatialElementType_EXTERNAL   // PredefinedType
        );

    }

    template<typename Schema>
    typename Schema::IfcVirtualElement *IfcVirtualElement() {

        return new typename Schema::IfcVirtualElement(
                guid(),                                             // GlobalId
                0,                                                  // OwnerHistory
                null,                                               // Name
                null,                                               // Description
                null,                                               // ObjectType
                0,                                                  // ObjectPlacement
                0,                                                  // Representation
                null                                                // Tag
        );

    }

    template<typename Schema>
    typename Schema::IfcSIUnit *IfcSIUnit(typename Schema::IfcUnitEnum::Value UnitType, typename Schema::IfcSIUnitName::Value Name) {

        return new typename Schema::IfcSIUnit(
                UnitType,   // UnitType
                null,       // Prefix
                Name        // Name
        );

    }

    template<typename Schema>
    typename Schema::IfcQuantityVolume *IfcQuantityVolume(const std::string &Name, typename Schema::IfcSIUnit *Unit, const double &VolumeValue) {

        return new typename Schema::IfcQuantityVolume(
                Name,           // Name
                null,           // Description
                Unit,           // Unit
                VolumeValue,    // VolumeValue
                null            // Formula
        );

    }

    template<typename Schema>
    typename Schema::IfcQuantityVolume *IfcQuantityVolume_2x3(const std::string &Name, typename Schema::IfcSIUnit *Unit, const double &VolumeValue) {

        return new typename Schema::IfcQuantityVolume(
                Name,           // Name
                null,           // Description
                Unit,           // Unit
                VolumeValue    // VolumeValue
        );

    }

    template<typename Schema>
    typename Schema::IfcQuantityArea *IfcQuantityArea(const std::string &Name, typename Schema::IfcSIUnit *Unit, const double &VolumeValue) {

        return new typename Schema::IfcQuantityArea(
                Name,           // Name
                null,           // Description
                Unit,           // Unit
                VolumeValue,    // VolumeValue
                null            // Formula
        );

    }

    template<typename Schema>
    typename Schema::IfcQuantityArea *IfcQuantityArea_2x3(const std::string &Name, typename Schema::IfcSIUnit *Unit, const double &VolumeValue) {

        return new typename Schema::IfcQuantityArea(
                Name,           // Name
                null,           // Description
                Unit,           // Unit
                VolumeValue     // VolumeValue
        );

    }

    template<typename Schema>
    typename Schema::IfcRelDefinesByProperties *IfcRelDefinesByProperties(typename Schema::IfcOwnerHistory *OwnerHistory, typename IfcTemplatedEntityList<typename Schema::IfcObjectDefinition>::ptr RelatedObjects, typename Schema::IfcPropertySetDefinitionSelect *RelatingPropertyDefinition) {

        return new typename Schema::IfcRelDefinesByProperties(
                guid(),                                             // GlobalId
                OwnerHistory,                                       // OwnerHistory
                null,                                               // Name
                null,                                               // Description
                RelatedObjects,                                     // RelatedObjects
                RelatingPropertyDefinition                          // RelatingPropertyDefinition
        );

    }

    template<typename Schema>
    typename Schema::IfcRelDefinesByProperties *IfcRelDefinesByProperties_2x3(typename Schema::IfcOwnerHistory *OwnerHistory, typename IfcTemplatedEntityList<typename Schema::IfcObject>::ptr RelatedObjects, typename Schema::IfcPropertySetDefinition *RelatingPropertyDefinition) {

        return new typename Schema::IfcRelDefinesByProperties(
                guid(),                                             // GlobalId
                OwnerHistory,                                       // OwnerHistory
                null,                                               // Name
                null,                                               // Description
                RelatedObjects,                                     // RelatedObjects
                RelatingPropertyDefinition                          // RelatingPropertyDefinition
        );

    }

    template<typename Schema>
    typename Schema::IfcElementQuantity *IfcElementQuantity(typename Schema::IfcOwnerHistory *OwnerHistory, const std::string &Name, typename IfcTemplatedEntityList<typename Schema::IfcPhysicalQuantity>::ptr Quantities) {

        return new typename Schema::IfcElementQuantity(
                guid(),                             // GlobalId
                OwnerHistory,                       // OwnerHistory
                Name,                               // Name
                null,                               // Description
                null,                               // MethodOfMeasurement
                Quantities                          // Quantities
        );

    }

    template<typename Schema>
    typename Schema::IfcTable *IfcTable(const std::string &Name, typename IfcTemplatedEntityList<typename Schema::IfcTableRow>::ptr Rows, typename IfcTemplatedEntityList<typename Schema::IfcTableColumn>::ptr Columns) {

        return new typename Schema::IfcTable(
                Name,   // Name
                Rows,   // Rows
                Columns // Columns
        );
    }

    template<typename Schema>
    typename Schema::IfcTableRow *IfcTableRow(IfcEntityList::ptr RowCells, bool IsHeading) {

        return new typename Schema::IfcTableRow(
                RowCells,   // RowCells
                IsHeading   // IsHeading
        );
    }

    template<typename Schema>
    typename Schema::IfcTableColumn *IfcTableColumn(const std::string &Identifier, const std::string &Name, const std::string &Description, typename Schema::IfcUnit *Unit, typename Schema::IfcReference *ReferencePath) {

        return new typename Schema::IfcTableColumn(
                Identifier,     // Identifier
                Name,           // Name
                Description,    // Description
                Unit,           // Unit
                ReferencePath   // ReferencePath
        );
    }

    template<typename Schema>
    typename Schema::IfcTableColumn *IfcTableColumn_noUnit(const std::string &Identifier, const std::string &Name, const std::string &Description, typename Schema::IfcReference *ReferencePath) {

        return new typename Schema::IfcTableColumn(
                Identifier,     // Identifier
                Name,           // Name
                Description,    // Description
                0,           // Unit
                ReferencePath   // ReferencePath
        );
    }

    template<typename Schema>
    typename Schema::IfcReference *IfcReference() {

        return new typename Schema::IfcReference(
                null,        // TypeIdentifier
                null,        // AttributeIdentifier
                null,        // InstanceName
                null,        // ListPositions
                0            // InnerReference
        );
    }

    template<typename Schema>
    typename Schema::IfcText *IfcText(const std::string &v) {
        return new typename Schema::IfcText(v);
    }

    template<typename Schema>
    typename Schema::IfcIdentifier *IfcIdentifier(const std::string &v) {
        return new typename Schema::IfcIdentifier(v);
    }

    template<typename Schema>
    typename Schema::IfcLabel *IfcLabel(const std::string &v) {
        return new typename Schema::IfcLabel(v);
    }

    template<typename Schema>
    typename Schema::IfcReal *IfcReal(double v) {
        return new typename Schema::IfcReal(v);
    }

    template<typename Schema>
    typename Schema::IfcBoolean *IfcBoolean(bool v) {
        return new typename Schema::IfcBoolean(v);
    }

//    template<typename Schema>
//    void add_product_to_storey_aggregation_via_related_objects(typename Schema::IfcBuildingStorey *IfcBuildingStorey, typename Schema::IfcProduct *IfcProduct) {
//        //typename Schema::IfcRelAggregates *IfcRelAggregates = *(IfcBuildingStorey->IsDecomposedBy()->begin());
//        auto IfcRelAggregates = *(IfcBuildingStorey->IsDecomposedBy()->begin()); // in 2x3 it's IfcRelDecomposes
//        auto RelatedObjects = IfcRelAggregates->RelatedObjects();
//        RelatedObjects->push(IfcProduct);
//        IfcRelAggregates->setRelatedObjects(RelatedObjects);
//    }

    template<typename Schema>
    void add_product_to_storey_aggregation_via_related_objects(typename Schema::IfcBuildingStorey *IfcBuildingStorey, typename Schema::IfcProduct *IfcProduct) {
        //typename Schema::IfcRelAggregates *IfcRelAggregates = *(IfcBuildingStorey->IsDecomposedBy()->begin());
        auto L = *(IfcBuildingStorey->IsDecomposedBy()); // in 2x3 it's IfcRelDecomposes

        if (L.size() == 0) {
            std::cerr << "[Error] Building storey " << IfcBuildingStorey->data().toString() << " has no IfcRelAggregate!" << std::endl;
            return;
        }

        auto IfcRelAggregates = *L.begin();
        auto RelatedObjects = IfcRelAggregates->RelatedObjects();
        RelatedObjects->push(IfcProduct);
        IfcRelAggregates->setRelatedObjects(RelatedObjects);
    }

    template<typename Schema>
    void add_product_to_building_aggregation_via_related_objects(typename Schema::IfcBuilding *IfcBuilding, typename Schema::IfcProduct *IfcProduct) {

        auto L = *(IfcBuilding->IsDecomposedBy()); // in 2x3 it's IfcRelDecomposes

        if (L.size() == 0) {
            std::cerr << "[Error] Building storey " << IfcBuilding->data().toString() << " has no IfcRelAggregate!" << std::endl;
            return;
        }

        auto IfcRelAggregates = *L.begin();
        auto RelatedObjects = IfcRelAggregates->RelatedObjects();
        RelatedObjects->push(IfcProduct);
        IfcRelAggregates->setRelatedObjects(RelatedObjects);
    }

    template<typename Schema>
    void add_product_to_storey_aggregation_via_related_elements(typename Schema::IfcBuildingStorey *IfcBuildingStorey, typename Schema::IfcProduct *IfcProduct) {
        auto IfcRelContainedInSpatialStructure = *(IfcBuildingStorey->ContainsElements()->begin());
        auto RelatedElements = IfcRelContainedInSpatialStructure->RelatedElements();

        if (RelatedElements->contains(IfcProduct)) return;

        if (IfcProduct->declaration().is("IfcElement"))
            if (IfcProduct->template as<typename Schema::IfcElement>()->ContainedInStructure()->size() != 0) return;

        RelatedElements->push(IfcProduct);
        IfcRelContainedInSpatialStructure->setRelatedElements(RelatedElements);
    }

    template<typename Schema>
    typename Schema::IfcObjectDefinition *storey_of_object(typename Schema::IfcObject *IfcObject) {

        auto Decomposes = IfcObject->Decomposes();

        if (Decomposes->size() == 0)
            std::cerr << "[Error] No Decomposes!" << std::endl;

        auto IfcRelAggregates = *(Decomposes->begin());
        auto RelatingObject = IfcRelAggregates->RelatingObject();

        if (!RelatingObject->declaration().is("IfcBuildingStorey"))
            std::cerr << "[Error] No storey!" << std::endl;

        return RelatingObject;
    }

    inline void remove_entities(std::unique_ptr<IfcParse::IfcFile> &model, const boost::shared_ptr<IfcEntityList> &entities) {

        if (entities == nullptr) return;

        for (auto it = entities->end() - 1; it >= entities->begin(); --it) {

            unsigned int id = (*it)->data().id();
            try { auto e = model->instance_by_id(id); }
            catch (...) {
                std::cerr << "[Warning] Entity with id " << id << " not found." << std::endl;
                continue;
            }
            IfcUtil::IfcBaseClass *const inst = *it;
            //std::cout << "[Info] Remove " << inst->data().toString() << std::endl;
            model->removeEntity(inst);
        }
    }

    inline void remove_entity(std::unique_ptr<IfcParse::IfcFile> &model, IfcUtil::IfcBaseClass *entity) {

        if (entity == nullptr) return;

        unsigned int id = entity->data().id();
        try { auto e = model->instance_by_id(id); }
        catch (...) {
            std::cerr << "[Warning] Entity with id " << id << " not found." << std::endl;
            return;
        }
        //std::cout << "[Info] Remove " << entity->data().toString() << std::endl;
        model->removeEntity(entity);
    }

    template<typename Schema>
    void remove_duplicate_guids(std::unique_ptr<IfcParse::IfcFile> &model) {

        boost::shared_ptr<IfcEntityList> IfcEntityList = nullptr;
        try { IfcEntityList = model->instances_by_type("IfcProduct"); }
        catch (...) {}

        if (IfcEntityList == nullptr) return;

        std::set<std::string> unique_guids;

        for (auto E: *IfcEntityList) {

            auto P = E->as<typename Schema::IfcProduct>();

            if (unique_guids.find(P->GlobalId()) != unique_guids.end()) {
                std::string new_guid = guid().formatted();

                while (unique_guids.find(new_guid) != unique_guids.end())
                    new_guid = guid().formatted();

                std::cerr << "[Info] Change guid of product with line number " << P->data().id() << " from " << P->GlobalId() << " to " << new_guid << std::endl;

                P->setGlobalId(new_guid);
            }
            unique_guids.insert(P->GlobalId());
        }
    }

}
#endif //IFC_CREATOR_H
