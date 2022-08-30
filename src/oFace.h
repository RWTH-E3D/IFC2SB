// Copyright 2022 Eric Fichter
#ifndef ORIGFACE_H
#define ORIGFACE_H

#include "headers.h"

class Product;

//! Describes a face, taken from the shape of an ifc product
class oFace {

public:
    oFace(TopoDS_Face _face, Product *_product, unsigned short int _shell_id);

    bool operator==(const oFace &O) const { return (product == O.product && shell_id == O.shell_id && face == O.face); }

    //! TopoDS_Face of this. TopoDS_Face references a TShape (geometry). Has its own Location and attributes (e.g. orientation).
    //! TopoDS_Face is not a pointer!
    TopoDS_Face face;

    //! Returns the related ifc product, if this is part of an opening (e.g. window). Else returns nullptr.
    IfcUtil::IfcBaseClass *Opening() const;

    //! Returns if face is an offset face and therefore introduced by this tool by extension of ifc faces.
    bool IsOffset() const;

    //! Returns if face's RelProduct is an IfcVirtualElement.
    bool IsVirtual() const;

    //! Returns true, if face is linked to an opening (e.g. window).
    bool IsOpening() const;

    //! Returns surface normal.
    gp_Dir Normal() const;

    //! Returns attribute shell_id. Faces from an ifc product belonging to same shell share an id.
    unsigned int ShellID() const;

    //! Returns attribute parent_face_id. Only faces from walls that were cut by an opening have a parent_face_id differing from 0
    unsigned int ParentID() const;

    //! Returns related product.
    Product *RelProduct() const;

    //! Returns attribute status_normal.
    face_normal_status NormalStatus() const;

    //! Sets opening product.
    void SetOpening(IfcUtil::IfcBaseClass *o);

    //! Sets attribute parent_face_id.
    void SetParentID(unsigned int v);

    //! Sets attribute isOffset.
    void SetIsOffset(bool b);

    //! Sets attribute shell_id.
    void SetShellID(unsigned int v);

    //! Sets attribute status_normal.
    void SetNormalStatus(face_normal_status s);

    //! Returns name of the IfcClass.
    std::string IfcClass() const;

    //! Returns guid of the relating product.
    std::string IfcGuid() const;

    //! Returns true if product is of specified IfcClass.
    bool IsIfcClass(std::string s) const;

    //! String with basic info.
    std::string Info() const;

private:

    bool isOffset;
    gp_Dir n; // surface normal
    face_normal_status status_normal;
    unsigned short int shell_id; // id of the shell of the product, the face belonged to. it's an id generated using a counter, it's not the hash of the TopoDS_Shell
    unsigned int parent_face_id;
    Product *product;
    IfcUtil::IfcBaseClass *opening; // faces of openings e.g. windows are linked to the wall by product variable, but using OpeningProduct to the ifcwindow

};

#endif //ORIGFACE_H