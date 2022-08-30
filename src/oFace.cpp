// Copyright 2022 Eric Fichter
#include "oFace.h"

oFace::oFace(TopoDS_Face _face, Product *_product, unsigned short int _shell_id) : face(_face), product(_product), shell_id(_shell_id) {

    n = Kernel::face_normal(face);
    status_normal = FACE_NORMAL_UNKNOWN;
    isOffset = false;
    opening = nullptr;
    parent_face_id = 0;
}

IfcUtil::IfcBaseClass *oFace::Opening() const { return opening; }

unsigned int oFace::ShellID() const { return shell_id; }

unsigned int oFace::ParentID() const { return parent_face_id; }

bool oFace::IsOffset() const { return isOffset; }

bool oFace::IsVirtual() const {
    if (RelProduct() == nullptr) return false;
    return RelProduct()->ifcproduct->declaration().is("IfcVirtualElement");
}

bool oFace::IsOpening() const { return opening != nullptr; }

gp_Dir oFace::Normal() const {

    if (status_normal == FACE_NORMAL_UNKNOWN)
        std::cerr << "[Error] IFC Face (" << RelProduct()->ifcproduct->declaration().name() << ") has no fixed face normal. Returned value could be wrong." << std::endl;
    return n;
}

Product *oFace::RelProduct() const { return product; }

void oFace::SetParentID(unsigned int v) { parent_face_id = v; }

void oFace::SetOpening(IfcUtil::IfcBaseClass *o) { opening = o; }

void oFace::SetIsOffset(bool b) { isOffset = b; }

void oFace::SetShellID(unsigned int v) { shell_id = v; }

void oFace::SetNormalStatus(face_normal_status s) { status_normal = s; }

face_normal_status oFace::NormalStatus() const { return status_normal; }

std::string oFace::IfcClass() const {

    if (RelProduct() == nullptr) {
        std::cerr << "[Warning] No product." << std::endl;
        return "No product";
    } else return RelProduct()->IfcClass();
}

std::string oFace::IfcGuid() const {

    if (RelProduct() == nullptr) {
        std::cerr << "[Warning] No product." << std::endl;
        return "No product";
    } else return RelProduct()->guid;
}

bool oFace::IsIfcClass(std::string s) const {

    if (RelProduct() == nullptr) {
        std::cerr << "[Warning] No product." << std::endl;
        return false;
    } else return RelProduct()->IsIfcClass(std::move(s));
}

std::string oFace::Info() const { return IfcClass() + ", " + IfcGuid() + ", " + std::to_string(Kernel::hash(face)); }
