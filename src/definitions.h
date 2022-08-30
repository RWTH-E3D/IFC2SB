// Copyright 2022 Eric Fichter
#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#define VISUALIZATION
#define PARALLEL_PROCESSING
#define IFCCONVERT_DOUBLE_PRECISION
#define IFC2SB_VERSION "0.9"
#define IFC2SB_FULLNAME "IFC-based Space Boundary Generator"
#define IFC2SB_NAME "IFC2SB"
#define IFC2SB_AUTHOR "Eric Fichter"
#define IFC2SB_ORGANIZATION "Institute of Energy Efficiency and Sustainable Building at RWTH Aachen University"

#ifdef IFCCONVERT_DOUBLE_PRECISION
typedef double real_t;
#else
typedef float real_t;
#endif

//! Template switch for IfcOpenShell.
enum IFC_SCHEMA {
    IFC2X3, IFC4, IFC4X1, IFC4X2, IFC4X3_RC1
};

//! According to IFC4 space boundary definition.
enum SB_Type {
    SB_TYPE_1, SB_TYPE_2A, SB_TYPE_2B
};

//! According to IFC4 space boundary definition.
enum SB_internal_or_external_boundary {
    SB_IE_INTERNAL,           //	The space boundary faces a physical or virtual element where there is an internal space on the other side.
    SB_IE_EXTERNAL,           //	The space boundary faces a physical or virtual element where there is an external space on the other side.
    SB_IE_EXTERNAL_EARTH,     //	The space boundary faces a physical or virtual element where there is earth (or terrain) on the other side.
    SB_IE_EXTERNAL_WATER,     //	The space boundary faces a physical or virtual element where there is water (water component of terrain) on the other side.
    SB_IE_EXTERNAL_FIRE,      //	The space boundary faces a physical or virtual element where there is another building on the other side.
    SB_IE_NOTDEFINED
};

//! According to IFC4 space boundary definition.
enum SB_physical_or_virtual_boundary {
    SB_PV_PHYSICAL,             //	The space boundary is provided physically (by a physical element).
    SB_PV_VIRTUAL,              //	The space boundary is provided virtually (by a logical divider that has no physical manifestation).
    SB_PV_NOTDEFINED            //	No information available.
};

//! Defines if the topological element is used according to orientation of its geometry.
enum used_orientation {
    USEDFORWARD,
    USEDREVERSED,
};

//! Defines if the correct orientation of the face is known.
enum face_normal_status {
    FACE_NORMAL_KNOWN,
    FACE_NORMAL_UNKNOWN,
};

//! Defines if entity belongs to outer wire or hole or is within a face.
enum typeEdge2d {
    TYPE2D_OUTER_ORIG,
    TYPE2D_OUTER,
    TYPE2D_HOLE,
    TYPE2D_INNER,
};

//! Collects several IfcSpace data
struct ifcspaceInfo {
    IfcUtil::IfcBaseEntity *product;
    std::string guid;
    TopoDS_Solid shape;
    Bnd_Box bbox;
    double volume;

    ifcspaceInfo(IfcUtil::IfcBaseEntity *product, std::string guid, TopoDS_Solid shape, const Bnd_Box &bbox, double volume) : product(product), guid(std::move(std::move(guid))), shape(std::move(shape)), bbox(bbox), volume(volume) {}
};

typedef std::list<ifcspaceInfo> ifcspaceInfoList;

//! Structure holding extended cface data for space search.
class cFace;
class iFace;

struct oriFace {
    iFace *iface;
    used_orientation orient;

    oriFace(iFace *iface, used_orientation o) : iface(iface), orient(o) {
    }

    oriFace(){
        iface = nullptr;
        orient = USEDFORWARD;
    }

};

inline bool operator <(const oriFace& lhs, const oriFace& rhs) {
    return std::tie(lhs.iface, lhs.orient) < std::tie(rhs.iface, rhs.orient);
}

inline bool operator ==(const oriFace& lhs, const oriFace& rhs) {
    return std::tie(lhs.iface, lhs.orient) == std::tie(rhs.iface, rhs.orient);
}

struct iFace {
    cFace *cface;
    bool isUsedForward;
    bool isUsedReversed;
    bool bad;
    bool used;
    bool duplicate;
    std::unordered_map<used_orientation, std::unordered_map<unsigned int, oriFace>> nb;

    iFace(cFace *cface) : cface(cface) {
        isUsedForward = false;
        isUsedReversed = false;
        bad = false;
        used = false;
        duplicate = false;
    }

};

typedef std::map<cFace *, std::set<cFace *>> cfaceSetMap;

#endif //DEFINITIONS_H