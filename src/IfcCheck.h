// Copyright 2022 Eric Fichter
#ifndef IFCCHECK_H
#define IFCCHECK_H

#include "headers.h"

class IfcCheck {

public:
    IfcCheck(IfcParse::IfcFile *_model, unsigned int _num_threads);

private:
    IfcParse::IfcFile *model;
    IfcGeom::IteratorSettings settings;
    const unsigned int num_threads;

    enum space_errors {
        NEGATIVE_VOLUME,
        NON_CONSISTENT_ORIENTATION,
        NOT_CLOSED,
        SMALL_AREA,
        WRONG_NUMBER_OF_SHELLS
    };

    std::unordered_map<space_errors, std::string> error_to_string{
            {NEGATIVE_VOLUME,            "has negative volume"},
            {NON_CONSISTENT_ORIENTATION, "is not consistently oriented"},
            {NOT_CLOSED,                 "is not closed"},
            {SMALL_AREA,                 "has small area"},
            {WRONG_NUMBER_OF_SHELLS,     "has wrong number of shells"}
    };

    struct IfcSpace {
        IfcUtil::IfcBaseEntity *product;
        TopoDS_Shape shape;
        std::string guid;
        std::list<space_errors> errors;

        IfcSpace(IfcUtil::IfcBaseEntity *_product, TopoDS_Shape _shape, std::string _guid) : product(_product), shape(std::move(_shape)), guid(std::move(_guid)) {}
    };

    struct IfcExternalSpatialElement {
        IfcUtil::IfcBaseEntity *product;
        TopoDS_Shape shape;
        std::string guid;
        std::list<space_errors> errors;

        IfcExternalSpatialElement(IfcUtil::IfcBaseEntity *_product, TopoDS_Shape _shape, std::string _guid) : product(_product), shape(std::move(_shape)), guid(std::move(_guid)) {}
    };

    struct IfcSpaceBoundary {
        IfcUtil::IfcBaseEntity *product;
        IfcUtil::IfcBaseClass *space;
        TopoDS_Shape shape;
        std::string guid;
        std::list<space_errors> errors;

        IfcSpaceBoundary(IfcUtil::IfcBaseEntity *_product, IfcUtil::IfcBaseClass *_space, TopoDS_Shape _shape, std::string _guid) : product(_product), space(_space), shape(std::move(_shape)), guid(std::move(_guid)) {}
    };

    void check_model();

    void check_ifc_spaces();

    template<typename T>
    void check_ifc_spaces_shell_number(std::list<T> &spaces);

    template<typename T>
    void check_ifc_spaces_shell(std::list<T> &spaces);

    void check_ifc_external_spatial_elements();

    template<typename T>
    void evaluation(const std::list<T> &L, const std::string &name);

    void check_ifc_space_boundaries(const std::string &ifcclass);

    static void check_ifc_spaceboundaries_shell(std::list<IfcSpaceBoundary> &boundaries);
};

#endif //IFCCHECK_H