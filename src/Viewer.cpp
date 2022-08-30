// Copyright 2022 Eric Fichter
#ifdef VISUALIZATION
#include "Viewer.h"

void Viewer::visualize_products(const std::list<Product> &products) {

    std::cout << "Visualize products\n\n\n";
    std::list<viewerHelper::DisplayShapes> ds;
    for (const auto &product: products) {
        ds.emplace_back();
        ds.back().shape = product.shape;
        ds.back().clr_random = true;
    }
    ViewerMain::start_viewer(ds);
}

void Viewer::visualize_products_openings(const std::list<Product> &products) {

    std::cout << "Visualize opening products\n\n\n";
    std::list<viewerHelper::DisplayShapes> ds;
    for (const auto &product: products) {
        ds.emplace_back();
        ds.back().shape = product.shape;

        if (product.IsIfcClass("IfcDoor") || product.IsIfcClass("IfcWindow")) {
            ds.back().clr_by_string = true;
            ds.back().transparency = 0.0;
            ds.back().clr_string = "LIGHTGREEN";
        }
        if (!product.VoidFillingElements.empty()) {
            ds.back().clr_by_string = true;
            ds.back().transparency = 0.5;
            ds.back().clr_string = "BLUE";
        } else ds.back().transparency = 1;

    }
    ViewerMain::start_viewer(ds);
}

void Viewer::visualize_orig_faces(const std::list<oFace> &orig_faces) {

    std::cout << "Visualize original faces\n\n\n";
    std::list<viewerHelper::DisplayShapes> ds;
    for (const auto &f: orig_faces) {
        if (f.IsOffset()) continue;
        ds.emplace_back();
        ds.back().shape = f.face;
        ds.back().clr_random = true;
    }
    ViewerMain::start_viewer(ds);
}

void Viewer::visualize_orig_faces_fixed_normals(const std::list<oFace> &orig_faces) {

    std::cout << "Visualize original faces with fixed normals\n\n\n";
    std::list<viewerHelper::DisplayShapes> ds;
    for (const auto &f: orig_faces) {
        // if (f.IsOffset()) continue;
        ds.emplace_back();
        ds.back().shape = f.face;
        ds.back().clr_by_string = true;
        ds.back().transparency = 0.5;
        f.NormalStatus() == FACE_NORMAL_KNOWN ? ds.back().clr_string = "LIGHTGREEN" : ds.back().clr_string = "RED";
    }
    ViewerMain::start_viewer(ds);
}

void Viewer::visualize_orig_faces_openings(const std::list<oFace> &orig_faces) {

    std::cout << "Visualize original opening faces\n\n\n";
    std::list<viewerHelper::DisplayShapes> ds;
    for (const auto &f: orig_faces) {
        ds.emplace_back();
        ds.back().shape = f.face;
        ds.back().clr_by_string = true;
        ds.back().transparency = 0.5;
        if (f.IsOpening()) {
            ds.back().clr_string = "LIGHTGREEN";
            ds.back().transparency = 0.0;
        } else if (f.IsIfcClass("IfcDoor") || f.IsIfcClass("IfcWindow")) {
            ds.back().clr_string = "BLUE";
            ds.back().transparency = 0.0;
        } else {
            //ds.back().clr_string = "RED";
            ds.back().transparency = 0.8;
        }
    }
    ViewerMain::start_viewer(ds);
}

void Viewer::visualize_orig_faces_virtual_product(const std::list<oFace> &orig_faces) {

    std::cout << "Visualize original faces with virtual product\n\n\n";
    std::list<viewerHelper::DisplayShapes> ds;
    for (const auto &f: orig_faces) {
        ds.emplace_back();
        ds.back().shape = f.face;
        ds.back().clr_by_string = true;
        ds.back().transparency = 0.5;

        if (f.IsVirtual()) {
            ds.back().clr_string = "RED";
            ds.back().transparency = 0.5;
        } else {
            ds.back().clr_string = "LIGHTGREEN";
            ds.back().transparency = 1.0;
        }
    }
    ViewerMain::start_viewer(ds);
}

void Viewer::visualize_cFaces(const std::list<cFace> &cFaces) {

    std::cout << "Visualize cFaces\n\n\n";
    std::list<viewerHelper::DisplayShapes> ds;
    for (const auto &cface: cFaces) {

        if (cface.RelSpace() != nullptr)
            if (cface.RelSpace()->is_facade) continue;

        ds.emplace_back();
        ds.back().shape = cface.face;
//        ds.back().show_adaptor_face_normal = true;
//        ds.back().adaptor_face_normal = Kernel::face_normal(TopoDS::Face(ds.back().shape));
//        ds.back().center_of_face = Kernel::face_center(TopoDS::Face(ds.back().shape));
        ds.back().clr_by_string = true;
        ds.back().transparency = 0.6;

        if (cface.RelSpace() != nullptr)
            if (cface.RelSpace()->is_facade)
                ds.back().line_clr_string = "DARKGREY";

        if (cface.SBType() == SB_TYPE_2A && cface.InternalOrExternal() == SB_IE_INTERNAL)
            ds.back().clr_string = "LIGHTGREEN";
        else if (cface.SBType() == SB_TYPE_2A && cface.InternalOrExternal() == SB_IE_EXTERNAL)
            ds.back().clr_string = "BLUE";
        else if (cface.SBType() == SB_TYPE_2A && cface.InternalOrExternal() == SB_IE_EXTERNAL_EARTH)
            ds.back().clr_string = "BROWN";
        else if (cface.SBType() == SB_TYPE_2A && cface.InternalOrExternal() == SB_IE_EXTERNAL_FIRE)
            ds.back().clr_string = "ORANGE";
        else if (cface.SBType() == SB_TYPE_2A && cface.InternalOrExternal() == SB_IE_EXTERNAL_WATER)
            ds.back().clr_string = "CYAN";
        else if (cface.SBType() == SB_TYPE_2B && cface.InternalOrExternal() == SB_IE_INTERNAL)
            ds.back().clr_string = "LIGHTRED";
        else if (cface.SBType() == SB_TYPE_2B && cface.InternalOrExternal() == SB_IE_EXTERNAL)
            ds.back().clr_string = "RED";
        else if (cface.SBType() == SB_TYPE_2B && cface.InternalOrExternal() == SB_IE_EXTERNAL_EARTH)
            ds.back().clr_string = "DARKRED";
        //else
        //    ds.back().clr_string = "BLACK";

    }
    ViewerMain::start_viewer(ds);
}

void Viewer::visualize_component_search(const std::list<cFace> &cFaces, const std::set<std::pair<iFace *, used_orientation>> &L, const std::pair<iFace *, used_orientation> &oface) {

    std::list<viewerHelper::DisplayShapes> ds;

    for (const auto &cface: cFaces) {
        ds.emplace_back();
        ds.back().shape = cface.face;
        ds.back().transparency = 1;
    }

    for (auto l: L) {
        ds.emplace_back();
        ds.back().shape = l.first->cface->face;
        ds.back().clr_by_string = true;
        ds.back().transparency = 0.0;
        l.first->cface->IsIfcClass("IfcVirtualElement") ? ds.back().clr_string = "BLUE" : ds.back().clr_string = "LIGHTGREEN";
        ds.back().show_adaptor_face_normal = true;
        ds.back().adaptor_face_normal = l.second == USEDFORWARD ? Kernel::face_normal(TopoDS::Face(ds.back().shape)) : Kernel::face_normal(TopoDS::Face(ds.back().shape)).Reversed();
        ds.back().center_of_face = Kernel::face_center(TopoDS::Face(ds.back().shape));
    }

    ds.emplace_back();
    ds.back().shape = oface.first->cface->face;
    ds.back().clr_by_string = true;
    ds.back().transparency = 0.0;
    oface.first->cface->IsIfcClass("IfcVirtualElement") ? ds.back().clr_string = "ORANGE" : ds.back().clr_string = "RED";
    ds.back().show_adaptor_face_normal = true;
    ds.back().adaptor_face_normal = oface.second == USEDFORWARD ? Kernel::face_normal(TopoDS::Face(ds.back().shape)) : Kernel::face_normal(TopoDS::Face(ds.back().shape)).Reversed();
    ds.back().center_of_face = Kernel::face_center(TopoDS::Face(ds.back().shape));

    ViewerMain::start_viewer(ds);
}

void Viewer::visualize_cFaces_fixed_normals(const std::list<cFace> &cFaces) {

    std::cout << "Visualize cFaces with fixed normal\n\n\n";
    std::list<viewerHelper::DisplayShapes> ds;
    for (const auto &cface: cFaces) {
        ds.emplace_back();
        ds.back().shape = cface.face;
        ds.back().clr_by_string = true;
        ds.back().transparency = 0.5;

        if (cface.NormalStatus() == FACE_NORMAL_KNOWN) {
            ds.back().clr_string = "LIGHTGREEN";
            ds.back().transparency = 0.5;
        } else {
            ds.back().clr_string = "RED";
            ds.back().transparency = 1.0;
        }
    }
    ViewerMain::start_viewer(ds);
}

void Viewer::visualize_cFaces_coplanar_and_inner(const std::list<cFace> &cFaces) {

    std::cout << "Visualize cFaces flagged as coplanar or inner\n\n\n";
    std::list<viewerHelper::DisplayShapes> ds;
    for (const auto &cface: cFaces) {
        ds.emplace_back();
        ds.back().shape = cface.face;
        ds.back().clr_by_string = true;

        if (cface.IsCoplanar()) {
            ds.back().clr_string = "LIGHTGREEN";
            ds.back().transparency = 0.5;
        } else if (cface.IsInner()) {
            ds.back().clr_string = "BLUE";
            ds.back().transparency = 0.5;
        } else {
            ds.back().clr_string = "RED";
            ds.back().transparency = 1.0;
        }

    }
    ViewerMain::start_viewer(ds);
}

void Viewer::visualize_cFaces_hanging(const std::list<cFace> &cFaces) {

    std::cout << "Visualize cFaces flagged as hanging\n\n\n";
    std::list<viewerHelper::DisplayShapes> ds;
    unsigned int n = 0;
    for (const auto &cface: cFaces) {
        ds.emplace_back();
        ds.back().shape = cface.face;
        ds.back().clr_by_string = true;

        if (cface.IsHanging()) {
            ds.back().clr_string = "RED";
            ds.back().transparency = 0.5;
            n++;
        } else {
            ds.back().clr_string = "LIGHTGREEN";
            ds.back().transparency = 1.0;
        }
    }
    std::cout << n << " of " << cFaces.size() << " are hanging\n";
    ViewerMain::start_viewer(ds);
}

void Viewer::visualize_cFaces_bad_corresponding(std::list<cFace> &cFaces) {

    std::cout << "Visualize cFaces flagged as hanging\n\n\n";
    std::list<viewerHelper::DisplayShapes> ds;

    for (auto &cface: cFaces) {
        ds.emplace_back();
        ds.back().shape = cface.face;
        ds.back().clr_by_string = true;

        if (cface.Corresponding() != nullptr) {
            if (!cface.CheckCorresponding(0.001, 3.05)) {
                ds.back().clr_string = "RED";
                ds.back().transparency = 0;
            } else {
                ds.back().clr_string = "LIGHTGREEN";
                ds.back().transparency = 0.5;
            }
        } else
            ds.back().transparency = 1.0;
    }

    ViewerMain::start_viewer(ds);
}

void Viewer::visualize_cFaces_materials(const std::list<cFace> &cFaces) {

    std::cout << "Visualize cFaces\n\n\n";
    std::list<viewerHelper::DisplayShapes> ds;
    for (const auto &cface: cFaces) {
        ds.emplace_back();
        ds.back().shape = cface.face;
        ds.back().clr_by_string = true;
        ds.back().transparency = 0.0;

        size_t layers = 0;
        size_t products = cface.Materials().size();
        for (auto &m: cface.Materials()) layers += m.second.size();

        if (layers == 0) ds.back().clr_string = "LIGHTGREEN";
        else if (layers == 1) {
            if (products == 1) ds.back().clr_string = "DARKBLUE";
            else ds.back().clr_string = "CYAN";
        } else if (layers == 2) {
            if (products == 1) ds.back().clr_string = "MIDGREY";
            else ds.back().clr_string = "DARKERGREY";
        }
        else if (layers == 3) ds.back().clr_string = "YELLOW";
        else if (layers == 4) ds.back().clr_string = "DARKRED";
        else ds.back().clr_string = "PURPLE";

        //std::cerr << "M " << products << "\t" << layers << std::endl;
    }
    ViewerMain::start_viewer(ds);
}

void Viewer::visualize_cFaces_openings(const std::list<cFace> &cFaces) {

    std::cout << "Visualize opening cFaces\n\n\n";
    std::list<viewerHelper::DisplayShapes> ds;
    for (const auto &cface: cFaces) {
        ds.emplace_back();
        ds.back().shape = cface.face;
        ds.back().clr_by_string = true;

        if (cface.IsOpening()) {
            ds.back().clr_string = "LIGHTGREEN";
            ds.back().transparency = 0.0;
        } else if (cface.IsIfcClass("IfcDoor") || cface.IsIfcClass("IfcWindow")) {
            ds.back().clr_string = "BLUE";
            ds.back().transparency = 0.0;
        } else {
            ds.back().transparency = 0.9;
        }
    }
    ViewerMain::start_viewer(ds);
}

void Viewer::visualize_cFaces_non_openings(const std::list<cFace> &cFaces) {

    std::cout << "Visualize non-opening cFaces\n\n\n";
    std::list<viewerHelper::DisplayShapes> ds;
    for (const auto &cface: cFaces) {
        if (cface.IsOpening() || cface.IsIfcClass("IfcDoor") || cface.IsIfcClass("IfcWindow")) continue;
        ds.emplace_back();
        ds.back().shape = cface.face;
        ds.back().transparency = 0.0;
    }
    ViewerMain::start_viewer(ds);
}

void Viewer::visualize_cFaces_holes(const std::list<cFace> &cFaces) {

    std::cout << "Visualize cFaces with holes\n\n\n";
    std::list<viewerHelper::DisplayShapes> ds;
    for (const auto &cface: cFaces) {
        ds.emplace_back();
        ds.back().shape = cface.face;
        if (cface.HasHoles()) {
            ds.back().transparency = 0.0;
            ds.back().clr_by_string = true;
            ds.back().clr_string = "RED";
        } else ds.back().transparency = 1;
    }
    ViewerMain::start_viewer(ds);
}

void Viewer::visualize_cFaces_physical_virtual(const std::list<cFace> &cFaces) {

    std::cout << "Visualize cFaces phys/virt\n\n\n";
    std::list<viewerHelper::DisplayShapes> ds;
    for (const auto &cface: cFaces) {
        ds.emplace_back();
        ds.back().shape = cface.face;
        ds.back().clr_by_string = true;
        ds.back().transparency = 0.5;
        cface.PhysicalOrVirtual() == SB_PV_VIRTUAL ? ds.back().clr_string = "RED" : cface.PhysicalOrVirtual() == SB_PV_PHYSICAL ? ds.back().clr_string = "LIGHTGREEN" : ds.back().clr_string = "BLACK";
    }
    ViewerMain::start_viewer(ds);
}

void Viewer::visualize_cFaces_internal_external(const std::list<cFace> &cFaces) {

    std::cout << "Visualize cFaces int/ext\n\n\n";
    std::list<viewerHelper::DisplayShapes> ds;
    for (const auto &cface: cFaces) {
        ds.emplace_back();
        ds.back().shape = cface.face;
        ds.back().clr_by_string = true;
        ds.back().transparency = 0.5;
        cface.InternalOrExternal() == SB_IE_INTERNAL ? ds.back().clr_string = "LIGHTGREEN" : cface.InternalOrExternal() == SB_IE_NOTDEFINED ? ds.back().clr_string = "BLACK" : ds.back().clr_string = "RED";
    }
    ViewerMain::start_viewer(ds);
}

void Viewer::visualize_cFaces_type(const std::list<cFace> &cFaces) {

    std::cout << "Visualize cFaces type\n\n\n";
    std::list<viewerHelper::DisplayShapes> ds;
    for (const auto &cface: cFaces) {
        ds.emplace_back();
        ds.back().shape = cface.face;
        ds.back().clr_by_string = true;
        ds.back().transparency = 0.5;
        cface.SBType() == SB_TYPE_2A ? ds.back().clr_string = "LIGHTGREEN" : cface.SBType() == SB_TYPE_2B ? ds.back().clr_string = "RED" : ds.back().clr_string = "BLACK";
    }
    ViewerMain::start_viewer(ds);
}

void Viewer::visualize_cFaces_and_inners(const std::list<cFace> &cFaces, const std::list<cFace> &inners) {

    std::cout << "Visualize cFaces and inner faces\n\n\n";
    std::list<viewerHelper::DisplayShapes> ds;
    for (const auto &cface: cFaces) {
        ds.emplace_back();
        ds.back().shape = cface.face;
        ds.back().transparency = 1;
    }
    for (const auto &cface: inners) {
        ds.emplace_back();
        ds.back().shape = cface.face;
        ds.back().clr_by_string = true;
        ds.back().transparency = 0.0;
        ds.back().clr_string = "LIGHTGREEN";
    }
    ViewerMain::start_viewer(ds);
}

void Viewer::visualize_cFace_and_shape(const cFace &cface, const TopoDS_Shape &S) {

    std::cout << "Visualize cFace and shape\n\n\n";
    std::list<viewerHelper::DisplayShapes> ds;
    ds.emplace_back();
    ds.back().shape = cface.face;
    ds.back().transparency = 1;

    for (const auto &f: Topo(S).faces()) {
        ds.emplace_back();
        ds.back().shape = f;
        ds.back().transparency = 0.5;
        ds.back().clr_by_string = true;
        ds.back().clr_string = "LIGHTGREEN";
    }
    ViewerMain::start_viewer(ds);
}

void Viewer::visualize_cFace(const cFace &cface) {

    std::cout << "Visualize cFace\n\n\n";
    std::list<viewerHelper::DisplayShapes> ds;
    ds.emplace_back();
    ds.back().shape = cface.face;
    ds.back().show_adaptor_face_normal = true;
    ds.back().adaptor_face_normal = Kernel::face_normal(TopoDS::Face(ds.back().shape));
    ds.back().center_of_face = Kernel::face_center(TopoDS::Face(ds.back().shape));
    ViewerMain::start_viewer(ds);
}

void Viewer::visualize_cFace(const cFace &cface, const std::list<cFace> &cFaces) {

    std::cout << "Visualize cFace\n\n\n";
    std::list<viewerHelper::DisplayShapes> ds;

    for (const auto &f: cFaces) {
        if (f == cface) continue;
        ds.emplace_back();
        ds.back().shape = f.face;
        ds.back().transparency = 0.95;
    }

    ds.emplace_back();
    ds.back().shape = cface.face;
    ds.back().show_adaptor_face_normal = true;
    ds.back().adaptor_face_normal = Kernel::face_normal(TopoDS::Face(ds.back().shape));
    ds.back().center_of_face = Kernel::face_center(TopoDS::Face(ds.back().shape));
    ds.back().clr_by_string = true;
    ds.back().clr_string = "RED";
    ds.back().transparency = 0;
    ViewerMain::start_viewer(ds);
}

void Viewer::visualize_face_graph(std::list<cFace> &cFaces) {

    std::cout << "Visualize face graph\n\n\n";
    std::list<viewerHelper::DisplayShapes> ds;

    for (auto &cface: cFaces)
        cface.SetPointOnFace(Kernel::point_on_face(cface.face, cface.FixedFaceNormal(), cface.Info()));

    for (const auto &cface: cFaces) {

        if (Kernel::aabb(cface.face, 0).CornerMin().Z() > 2.4) continue;

        ds.emplace_back();
        ds.back().shape = cface.face;
        ds.back().transparency = 0.3;

        auto v = BRepBuilderAPI_MakeVertex(cface.PointOnFace()).Vertex();
        ds.emplace_back();
        ds.back().shape = v;
        ds.back().clr_by_string = true;
        ds.back().line_clr_string = true;
        ds.back().transparency = 0.0;
        ds.back().clr_string = "DARKRED";
        ds.back().line_clr_string = "DARKRED";

        for (const auto &a: cface.AdjacentFacesList()) {
            if (Kernel::aabb(a->face, 0).CornerMin().Z() > 2.4) continue;
            auto v2 = BRepBuilderAPI_MakeVertex(a->PointOnFace()).Vertex();
            auto e = BRepBuilderAPI_MakeEdge(v, v2).Edge();
            ds.emplace_back();
            ds.back().shape = e;
            ds.back().edge_thickness = 0.015;
            ds.back().clr_by_string = true;
            ds.back().transparency = 0.0;
            ds.back().clr_string = "BLACK";
        }
    }
    ViewerMain::start_viewer(ds);
}

void Viewer::visualize_cFaces_shadings(const std::list<cFace> &cFaces) {

    std::cout << "Visualize cFaces shadings\n\n\n";
    std::list<viewerHelper::DisplayShapes> ds;

    std::set<Product *> prod;
    for (auto &cface: cFaces)
        if (!cface.RelSpace()->is_facade)
            prod.insert(cface.RelProduct());

    for (const auto &cface: cFaces) {
        if (!cface.RelSpace()->is_facade) continue;
        ds.emplace_back();
        ds.back().shape = cface.RelSpace()->Shell();
        ds.back().transparency = 1;
        break;
    }

    for (const auto &cface: cFaces) {

        if (!cface.IsPotentiallyShading(prod)) continue;

        ds.emplace_back();
        ds.back().shape = cface.face;
        ds.back().clr_by_string = true;
        ds.back().transparency = 0.0;
        ds.back().clr_string = "PURPLE";
    }
    ViewerMain::start_viewer(ds);
}

void Viewer::visualize_shadings(const std::list<cFace> &cFaces, const std::unordered_map<std::string, std::list<TopoDS_Face>> &shadings) {

    std::cout << "Visualize shadings\n\n\n";
    std::list<viewerHelper::DisplayShapes> ds;

    for (const auto &cface: cFaces) {
        if (!cface.RelSpace()->is_facade) continue;
        ds.emplace_back();
        ds.back().shape = cface.face;
        ds.back().transparency = 0.7;
    }

    for (const auto &shading: shadings)
        for (const auto &s: shading.second) {
            ds.emplace_back();
            ds.back().shape = s;
            ds.back().clr_by_string = true;
            ds.back().transparency = 0.0;
            ds.back().clr_string = "PURPLE";
        }

    ViewerMain::start_viewer(ds);
}

void Viewer::visualize_external_walls(const std::list<cFace> &cFaces) {

    std::cout << "Visualize external walls\n\n\n";
    std::list<viewerHelper::DisplayShapes> ds;

    std::set<Product *> prod;
    for (auto &cface: cFaces)
        if (cface.RelSpace()->is_facade && cface.IsIfcClass("IfcWall"))
            prod.insert(cface.RelProduct());

    for (const auto &p: prod) {
        ds.emplace_back();
        ds.back().shape = p->shape; // product.shape = TopoDS_Shape(); must be removed
        ds.back().clr_string = "GREY";
    }

    ViewerMain::start_viewer(ds);
}

void Viewer::visualize_collisions(const std::set<std::pair<Product *, Product *>> &collisions) {

    std::cout << "Visualize collisions\n\n\n";
    std::list<viewerHelper::DisplayShapes> ds;

    std::set<Product *> products;
    for (const auto &c: collisions)
        for (const auto &pair: collisions) {
            products.insert(pair.first);
            products.insert(pair.second);
        }

    for (const auto &p: products) {
        ds.emplace_back();
        ds.back().shape = p->shape; // product.shape = TopoDS_Shape(); must be removed
        ds.back().transparency = 0.5;
    }

    ViewerMain::start_viewer(ds);
}

void Viewer::visualize_spaces(const std::list<Space> &spaces, bool show_shells, bool firstLvl) {

    std::cout << "Visualize spaces\n\n\n";

    std::list<viewerHelper::DisplayShapes> ds;

    if (show_shells)
        for (const auto &space: spaces) {
            if (space.is_facade) continue;
            ds.emplace_back();
            ds.back().shape = space.Shell();
            ds.back().clr_random = true;
            ds.back().transparency = 0.8;
        }

    else
        for (const auto &space: spaces) {
            if (space.is_facade) continue;

            BRep_Builder B;
            TopoDS_Compound C;
            B.MakeCompound(C);

            if (firstLvl)
                for (const auto &cface: space.FirstLevel())
                    B.Add(C, cface->face);
            else
                for (const auto &cface: space.SecondLevel())
                    B.Add(C, cface->face);

            ds.emplace_back();
            ds.back().shape = C;
            ds.back().clr_random = true;
            ds.back().transparency = 0.0;
        }
    ViewerMain::start_viewer(ds);
}

void Viewer::visualize_spaces_facade(const std::list<Space> &spaces, bool show_shells, bool firstLvl) {

    std::cout << "Visualize facade space\n\n\n";

    std::list<viewerHelper::DisplayShapes> ds;

    if (show_shells)
        for (const auto &space: spaces) {
            if (!space.is_facade) continue;
            ds.emplace_back();
            ds.back().shape = space.Shell();
            ds.back().clr_random = true;
            ds.back().transparency = 0.8;
        }

    else
        for (const auto &space: spaces) {
            if (!space.is_facade) continue;

            BRep_Builder B;
            TopoDS_Compound C;
            B.MakeCompound(C);

            if (firstLvl)
                for (const auto &cface: space.FirstLevel())
                    B.Add(C, cface->face);
            else
                for (const auto &cface: space.SecondLevel())
                    B.Add(C, cface->face);

            ds.emplace_back();
            ds.back().shape = C;
            ds.back().clr_random = true;
            ds.back().transparency = 0.0;
        }
    ViewerMain::start_viewer(ds);
}

void Viewer::visualize_spaces_to_delete(std::list<Space> &spaces, std::list<cFace> &cFaces, bool show_shells, const std::set<Space *> &del) {

    std::cout << "Visualize spaces to delete\n\n\n";

    std::list<viewerHelper::DisplayShapes> ds;

    if (show_shells)
        for (auto &space: spaces) {

            if (space.is_facade) {
                ds.emplace_back();
                ds.back().shape = space.Shell();
                ds.back().line_clr_string = "LIGHTGREY";
                ds.back().transparency = 1.0;
            }

            if (del.find(&space) != del.end()) {
                ds.emplace_back();
                ds.back().shape = space.Shell();
                ds.back().clr_by_string = true;
                ds.back().clr_string = "RED";
                ds.back().transparency = 0;
            }
        }
    else
        for (auto &space: spaces) {

            if (del.find(&space) != del.end()) {
                BRep_Builder B;
                TopoDS_Compound C;
                B.MakeCompound(C);
                for (const auto &cface: space.SecondLevel())
                    B.Add(C, cface->face);
                ds.emplace_back();
                ds.back().shape = C;
                ds.back().clr_by_string = true;
                ds.back().clr_string = "RED";
                ds.back().transparency = 0;
            }

            if (space.is_facade) {
                BRep_Builder B;
                TopoDS_Compound C;
                B.MakeCompound(C);
                for (const auto &cface: space.SecondLevel())
                    B.Add(C, cface->face);
                ds.emplace_back();
                ds.back().shape = C;
                ds.back().line_clr_string = "LIGHTGREY";
                ds.back().transparency = 1.0;
            }
        }
    ViewerMain::start_viewer(ds);
}

void Viewer::visualize_cFaces_as_space_boundaries(const std::list<cFace> &cFaces) {

    std::cout << "Visualize cFaces as space boundaries\n\n\n";
    std::list<viewerHelper::DisplayShapes> ds;
    for (const auto &cface: cFaces) {

        if (cface.RelSpace()->is_facade) continue;

        ds.emplace_back();
        ds.back().transparency = 0.7;
        ds.back().clr_by_string = true;
        ds.back().shape = cface.face;

        ds.back().show_adaptor_face_normal = true;
        ds.back().adaptor_face_normal = Kernel::face_normal(TopoDS::Face(ds.back().shape));
        ds.back().center_of_face = Kernel::face_center(TopoDS::Face(ds.back().shape));

        if (cface.RelSpace() != nullptr)
            if (cface.RelSpace()->is_facade)
                ds.back().line_clr_string = "DARKGREY";

        if (cface.SBType() == SB_TYPE_2A && cface.InternalOrExternal() == SB_IE_INTERNAL)
            ds.back().clr_string = "LIGHTGREEN";
        else if (cface.SBType() == SB_TYPE_2A && cface.InternalOrExternal() == SB_IE_EXTERNAL)
            ds.back().clr_string = "BLUE";
        else if (cface.SBType() == SB_TYPE_2A && cface.InternalOrExternal() == SB_IE_EXTERNAL_EARTH)
            ds.back().clr_string = "BROWN";
        else if (cface.SBType() == SB_TYPE_2B && cface.InternalOrExternal() == SB_IE_INTERNAL)
            ds.back().clr_string = "LIGHTRED";
        else if (cface.SBType() == SB_TYPE_2B && cface.InternalOrExternal() == SB_IE_EXTERNAL)
            ds.back().clr_string = "RED";
        else if (cface.SBType() == SB_TYPE_2B && cface.InternalOrExternal() == SB_IE_EXTERNAL_EARTH)
            ds.back().clr_string = "DARKRED";
        // else
        //    ds.back().clr_string = "BLACK";
    }
    ViewerMain::start_viewer(ds);
}

void Viewer::visualize_zones(const std::vector<std::vector<std::pair<gp_Pnt, std::set<oFace *>>>> &zones) {

    for (const auto &zone: zones) {

        std::list<viewerHelper::DisplayShapes> ds;

        for (const auto &space: zone) {

            // fluid point
            ds.emplace_back();
            ds.back().shape = BRepPrimAPI_MakeSphere(gp_Pnt(space.first.X(), space.first.Y(), space.first.Z()), 0.05);
            ds.back().clr_by_string = true;
            ds.back().clr_string = "GREEN";

            TopoDS_ListOfShape L;
            for (const auto &face: space.second)
                L.Append(face->face);

            ds.emplace_back();
            ds.back().shape = Kernel::compound_from_shape_list(L);
            ds.back().clr_by_string = true;
            ds.back().clr_string = "RED";
            ds.back().transparency = 0.9;

        }
        ViewerMain::start_viewer(ds);
    }
}

void Viewer::visualize_octree_products(std::list<oFace> &orig_faces, const std::vector<std::pair<std::array<double, 3>, std::set<std::string>>> &spaces_guids, const std::vector<std::pair<std::array<double, 3>, std::set<unsigned int>>> &spaces_triangles,
                                       const std::vector<std::set<unsigned int>> &zones, const std::vector<oFace *> &tri2orig) {

    for (unsigned int i = 0; i < zones.size(); i++) {
        std::cout << "Zone  " << i << ":\t";
        for (const auto &space: zones[i])
            std::cout << "\t" << space;
        std::cout << std::endl;
    }

    for (unsigned int i = 0; i < spaces_guids.size(); i++) {
        auto mid = spaces_guids[i].first;
        std::cout << "Space " << i << ":\t" << mid[0] << ",\t" << mid[1] << ",\t" << mid[2] << "  |\t" << spaces_guids[i].second.size() << "\t" << spaces_triangles[i].second.size() << std::endl;
    }

    std::map<std::string, std::list<oFace *>> M;

    for (auto &orig_face: orig_faces)
        M[orig_face.RelProduct()->guid].push_back(&orig_face);

    std::cout << "Visualize octree products\n\n\n";

    for (const auto &zone: zones) {

        std::list<viewerHelper::DisplayShapes> ds;

        for (const auto &space: zone) {

            auto spaces_tri = spaces_triangles[space];

            ds.emplace_back();
            ds.back().shape = BRepPrimAPI_MakeSphere(gp_Pnt(spaces_tri.first[0], spaces_tri.first[1], spaces_tri.first[2]), 0.05);
            ds.back().clr_by_string = true;
            ds.back().clr_string = "GREEN";

            std::set<oFace *> O;
            std::set<std::string> guids;

            // add every face that lies in space and store in guid list
            for (const auto &index: spaces_tri.second) {
                auto orig_face = tri2orig[index];
                O.insert(orig_face);
                guids.insert(orig_face->RelProduct()->guid);
            }

            // now add orig faces that belong to guid and are no offset faces
            for (const auto &guid: guids)
                for (const auto &orig_face: M[guid])
                    if (!orig_face->IsOffset())
                        O.insert(orig_face);

            for (const auto &f: O) {
                ds.emplace_back();
                ds.back().shape = f->face;
                ds.back().clr_by_string = true;
                ds.back().clr_string = "RED";
                ds.back().transparency = 0.9;
            }
        }
        ViewerMain::start_viewer(ds);
    }

    std::cout << "Visualize octree products\n\n\n";

    for (const auto &spaces_tri: spaces_triangles) {

        std::list<viewerHelper::DisplayShapes> ds;

        ds.emplace_back();
        ds.back().shape = BRepPrimAPI_MakeSphere(gp_Pnt(spaces_tri.first[0], spaces_tri.first[1], spaces_tri.first[2]), 0.05);
        ds.back().clr_by_string = true;
        ds.back().clr_string = "GREEN";

        std::set<oFace *> O;
        std::set<std::string> guids;

        // add every face that lies in space and store in guid list
        for (const auto &index: spaces_tri.second) {
            auto orig_face = tri2orig[index];
            O.insert(orig_face);
            guids.insert(orig_face->RelProduct()->guid);
        }

        // now add orig faces that belong to guid and are no offset faces
        for (const auto &guid: guids)
            for (const auto &orig_face: M[guid])
                if (!orig_face->IsOffset())
                    O.insert(orig_face);

        for (const auto &f: O) {
            ds.emplace_back();
            ds.back().shape = f->face;
            ds.back().clr_by_string = true;
            ds.back().clr_string = "RED";
            ds.back().transparency = 0.9;
        }

        ViewerMain::start_viewer(ds);
    }

    std::cout << "Visualize octree products\n\n\n";

    for (const auto &space_guids: spaces_guids) {

        std::list<viewerHelper::DisplayShapes> ds;

        ds.emplace_back();
        ds.back().shape = BRepPrimAPI_MakeSphere(gp_Pnt(space_guids.first[0], space_guids.first[1], space_guids.first[2]), 0.05);
        ds.back().clr_by_string = true;
        ds.back().clr_string = "GREEN";

        for (const auto &guid: space_guids.second)

            for (const auto &orig_face: M[guid]) {

                ds.emplace_back();
                ds.back().shape = orig_face->face;
                ds.back().clr_by_string = true;
                ds.back().clr_string = "RED";
                ds.back().transparency = 0.9;

            }

        ViewerMain::start_viewer(ds);
    }

    std::cout << "Visualize octree products\n\n\n";

    for (const auto &spaces_tri: spaces_triangles) {

        std::list<viewerHelper::DisplayShapes> ds;

        ds.emplace_back();
        ds.back().shape = BRepPrimAPI_MakeSphere(gp_Pnt(spaces_tri.first[0], spaces_tri.first[1], spaces_tri.first[2]), 0.05);
        ds.back().clr_by_string = true;
        ds.back().clr_string = "GREEN";

        std::set<oFace *> O;

        for (const auto &index: spaces_tri.second)
            O.insert(tri2orig[index]);

        for (const auto &f: O) {
            ds.emplace_back();
            ds.back().shape = f->face;
            ds.back().clr_by_string = true;
            ds.back().clr_string = "RED";
            ds.back().transparency = 0.9;
        }

        ViewerMain::start_viewer(ds);
    }
}

void Viewer::visualize_shape(const TopoDS_Shape &S) {

    std::cout << "Visualize shape\n\n\n";
    std::list<viewerHelper::DisplayShapes> ds;
    ds.emplace_back();
    ds.back().shape = S;
    ds.back().clr_by_string = true;
    ds.back().transparency = 0.0;
    ViewerMain::start_viewer(ds);
}

void Viewer::visualize_shape_triangulation(const TopoDS_Shape &S) {

    std::cout << "Visualize shape triangulation\n\n\n";

    std::list<viewerHelper::DisplayShapes> ds;

    for (auto &elem: Topo(S).faces()) {

        TopoDS_Face F = TopoDS::Face(elem);

        TopLoc_Location L;
        Handle_Poly_Triangulation Poly_Triangulation = BRep_Tool::Triangulation(F, L);

        if (Poly_Triangulation.IsNull()) {
            std::cerr << "[Warning] Null triangulation." << std::endl;
            continue;
        }

        const Poly_Array1OfTriangle &Poly_Array1OfTriangle = Poly_Triangulation->Triangles();

        for (int i = 1; i <= Poly_Triangulation->NbTriangles(); ++i) {

            auto Poly_Triangle = Poly_Array1OfTriangle.Value(i);
            int i1, i2, i3;
            Poly_Triangle.Get(i1, i2, i3);

            gp_Pnt p1 = Poly_Triangulation->Node(i1).Transformed(L.Transformation());
            gp_Pnt p2 = Poly_Triangulation->Node(i2).Transformed(L.Transformation());
            gp_Pnt p3 = Poly_Triangulation->Node(i3).Transformed(L.Transformation());

            std::cerr << p1.X() << "\t" << p1.Y() << "\t" << p1.Z() << std::endl;
            std::cerr << p2.X() << "\t" << p2.Y() << "\t" << p2.Z() << std::endl;
            std::cerr << p3.X() << "\t" << p3.Y() << "\t" << p3.Z() << std::endl;

            TopoDS_Wire wire = BRepBuilderAPI_MakePolygon(p1, p2, p3, true).Wire();
            TopoDS_Face poly = BRepBuilderAPI_MakeFace(wire).Face();

            ds.emplace_back();
            ds.back().shape = poly;
            ds.back().clr_random = true;
            ds.back().transparency = 0.3;
        }
    }

    ViewerMain::start_viewer(ds);
}

void Viewer::visualize_polygon(const std::vector<std::vector<gp_Pnt2d>> &ClipWires, const std::vector<std::vector<gp_Pnt2d>> &Holes2D) {

    std::cout << "Visualize polygon\n\n\n";
    std::list<viewerHelper::DisplayShapes> ds;
    for (auto &W: ClipWires) {
        BRepBuilderAPI_MakePolygon m;
        for (auto &P: W) m.Add(gp_Pnt(P.X(), P.Y(), 0));
        m.Close();
        TopoDS_Face F = BRepBuilderAPI_MakeFace(m.Wire()).Face();
        ds.emplace_back();
        ds.back().shape = m.Wire();
        ds.back().clr_by_string = true;
        ds.back().transparency = 0.5;
        ds.back().clr_by_string = true;
        ds.back().clr_string = "GREEN";
    }
    for (auto &W: Holes2D) {
        BRepBuilderAPI_MakePolygon m;
        for (auto &P: W) m.Add(gp_Pnt(P.X(), P.Y(), 0));
        m.Close();
        TopoDS_Face F = BRepBuilderAPI_MakeFace(m.Wire()).Face();
        ds.emplace_back();
        ds.back().shape = m.Wire();
        ds.back().clr_by_string = true;
        ds.back().transparency = 0.5;
        ds.back().clr_by_string = true;
        ds.back().clr_string = "RED";
    }
    ViewerMain::start_viewer(ds);
}

void Viewer::visualize_shapelist(const TopoDS_ListOfShape &L) {

    std::cout << "Visualize shape\n\n\n";
    std::list<viewerHelper::DisplayShapes> ds;
    for (auto &S: L) {
        ds.emplace_back();
        ds.back().shape = S;
        ds.back().clr_by_string = true;
        ds.back().transparency = 0.5;
    }
    ViewerMain::start_viewer(ds);
}

void Viewer::visualize_clip(const ClipperLib::Paths &poly1, const ClipperLib::Paths &poly2, const ClipperLib::Paths &sol) {

    std::cout << "Visualize clip\n\n\n";
    std::list<viewerHelper::DisplayShapes> ds;

    ds.emplace_back();
    ds.back().shape = visualize_clip(poly1, 0);
    ds.back().clr_by_string = true;
    ds.back().transparency = 0.5;
    ds.back().clr_by_string = true;
    ds.back().clr_string = "GREEN";

    ds.emplace_back();
    ds.back().shape = visualize_clip(poly2, 10000);
    ds.back().clr_by_string = true;
    ds.back().transparency = 0.5;
    ds.back().clr_by_string = true;
    ds.back().clr_string = "RED";

    ds.emplace_back();
    ds.back().shape = visualize_clip_sol(sol, 100000);
    ds.back().clr_by_string = true;
    ds.back().transparency = 0.5;
    ds.back().clr_by_string = true;
    ds.back().clr_string = "RED";

    ViewerMain::start_viewer(ds);
}

TopoDS_Face Viewer::visualize_clip(const ClipperLib::Paths &poly, double z) {

    auto outerWire1 = poly[0];
    BRepBuilderAPI_MakePolygon m;
    for (auto &P: outerWire1) m.Add(gp_Pnt(P.X, P.Y, z));
    m.Close();
    TopoDS_Wire W1 = m.Wire();
    TopoDS_Face F1 = BRepBuilderAPI_MakeFace(W1);
    BRepBuilderAPI_MakeFace MF(F1);

    for (unsigned int i = 1; i < poly.size(); i++) {
        BRepBuilderAPI_MakePolygon n;
        for (auto &P: poly[i]) n.Add(gp_Pnt(P.X, P.Y, z));
        n.Close();
        MF.Add(n.Wire());
    }

    return MF.Face();
}

TopoDS_Face Viewer::visualize_clip_sol(const ClipperLib::Paths &poly, double z) {

    TopoDS_Wire W1;
    TopoDS_ListOfShape iW;

    for (const auto &s: poly) {

        if (fabs(ClipperLib::Area(s)) < 1e6) continue; // skip tiny polygons

        if (ClipperLib::Orientation(s)) {
            BRepBuilderAPI_MakePolygon m;
            for (auto &P: s) m.Add(gp_Pnt(P.X, P.Y, z));
            m.Close();
            W1 = m.Wire();
        } else {
            BRepBuilderAPI_MakePolygon n;
            for (auto &P: s) n.Add(gp_Pnt(P.X, P.Y, z));
            n.Close();
            iW.Append(n.Wire());
        }
    }

    TopoDS_Face F1 = BRepBuilderAPI_MakeFace(W1);
    BRepBuilderAPI_MakeFace MF(F1);

    for (auto &w: iW)
        MF.Add(TopoDS::Wire(w));

    return MF.Face();
}

void Viewer::visualize_sFaces(const std::list<sFace> &sFaces) {

    std::cout << "Visualize sFaces\n\n\n";
    std::list<viewerHelper::DisplayShapes> ds;
    for (const auto &sface: sFaces) {
        ds.emplace_back();
        ds.back().shape = sface.face;
        ds.back().show_adaptor_face_normal = true;
        ds.back().adaptor_face_normal = Kernel::face_normal(TopoDS::Face(ds.back().shape));
        ds.back().center_of_face = Kernel::face_center(TopoDS::Face(ds.back().shape));
        ds.back().transparency = 0.0;
    }
    ViewerMain::start_viewer(ds);
}
#endif