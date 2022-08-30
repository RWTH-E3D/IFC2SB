// Copyright 2022 Eric Fichter
#include "Space.h"

#include <memory>

Space::Space(unsigned int _id, std::set<cFace *> &_faces, bool update_faces_after_sewing) : FirstLvl(_faces), id(_id) {
    wasFlipped = false;
    shell = build_shell_from_cFaces(update_faces_after_sewing);
    Bnd_Box aabb = Kernel::aabb(shell, 0);
    aabb.Get(aabb_xmin, aabb_ymin, aabb_zmin, aabb_xmax, aabb_ymax, aabb_zmax);
    volume = Kernel::volume(shell);
    calc_floor_elevation();
    is_facade = false;
    classifier = nullptr;
    //intersector = nullptr;
    storey = nullptr;
    nbs.clear();
    SecondLvl.clear();
    old_space_guids.clear();
    IfcRelContainedInSpatialStructures.clear();
    IfcRelDefinesByProperties.clear();
    IfcRelDefinesByType = nullptr;
    IfcRelAssigns.clear();
    IfcRelAssociates.clear();
    NetVolume = 0;
    NetFloorArea = 0;
    NetWallArea = 0;
    NetCeilingArea = 0;

    old_space_info["Name"];
    old_space_info["Description"];
    old_space_info["ObjectType"];
    old_space_info["LongName"];
    old_space_info["PredefinedType"];
}

TopoDS_Shell Space::build_shell_from_cFaces(bool update_faces_after_sewing) {

    BRepBuilderAPI_Sewing sew;

    for (auto &cface: FirstLvl)
        sew.Add(cface->face);

    sew.SetNonManifoldMode(true); // Non-manifold mode enables sewing of non-manifold topology.
    sew.Perform();
    TopoDS_Shape shape = sew.SewedShape();
    TopExp_Explorer Ex;
    TopoDS_Shell S;

    if (update_faces_after_sewing)
        for (auto &cface: FirstLvl) {
            auto M = sew.Modified(cface->face);
            if (!cface->face.IsSame(M))
                cface->face = TopoDS::Face(M);
        }

    for (Ex.Init(shape, TopAbs_SHELL); Ex.More(); Ex.Next()) {
        S = TopoDS::Shell(Ex.Current());
        break;
    }

    if (S.IsNull()) std::cerr << "\t[Warning] No shell for space " << id << std::endl;

    BRepCheck_Shell c = BRepCheck_Shell(S);
    BRepCheck_Status Closed = c.Closed();
    BRepCheck_Status Orientation = c.Orientation();

    if (Closed != BRepCheck_NoError)
        std::cerr << "\t[Warning] Shell is not closed for space " << id << std::endl;

    if (Orientation != BRepCheck_NoError)
        std::cerr << "\t[Warning] Shell is not consistently oriented for space " << id << std::endl;

    if (Kernel::volume(S) < 0) {
        S.Complement();
        wasFlipped = true;
        //std::cerr << "\tWarning: Shell normals pointing inwards" << std::endl;
    }

    if (Kernel::volume(S) < 0)
        std::cerr << "\t[Warning] Shell normals pointing inwards again for space " << id << std::endl;

    ShapeAnalysis_Shell a;
    a.LoadShells(TopoDS::Shell(S));
    a.CheckOrientedShells(S, true, true);

    //! Edges which do not fulfill these conditions are bad:
    //! edge is, either present once (free edge) or twice (connected
    //! edge) but with different orientations (FORWARD/REVERSED)
    if (a.HasBadEdges()) {

        std::set<unsigned int> bad_hashes;

        std::cerr << "\t[Warning] Shell of space " << id << " has \"bad\" edges (indicates non-manifold space (is ok) or bad face orientations (not ok)):";
        for (auto &b: Topo(a.BadEdges()).edges()) {
            bad_hashes.insert(Kernel::hash(b));
            //std::cerr << " " << Kernel::hash(b);
        }
        std::cerr << std::endl;

        // std::cerr << "\t          Bad edges:\n";

        std::unordered_map<unsigned int, std::list<unsigned int>> Me;
        for (auto &e: Topo(S).edges())
            if(bad_hashes.find(Kernel::hash(e)) != bad_hashes.end())
                Me[Kernel::hash(e)].push_back(e.Orientation());

        for (auto &e: Me) {
            std::cerr << "\t             " << e.first << ":";
            for (auto &h: e.second)
                std::cerr << " " << h;
            std::cerr << std::endl;
        }

        //std::cerr << "\tProducts in space:";
        //for (auto &p: ProductsInSpace()) std::cerr << " " << p->data().getArgument(0)->toString();
        //std::cerr << std::endl;
        //ShapeChecker(S).dump_topology();
    }
    if (a.HasFreeEdges()) std::cerr << "\t[Warning] Shell free edges " << id << std::endl;

    return S;
}

bool Space::is_point_in_aabb(gp_Pnt P) const {
    if (P.X() > aabb_xmax) return false;
    if (P.Y() > aabb_ymax) return false;
    if (P.Z() > aabb_zmax) return false;
    if (P.X() < aabb_xmin) return false;
    if (P.Y() < aabb_ymin) return false;
    if (P.Z() < aabb_zmin) return false;
    return true;
}

bool Space::is_point_in_aabb(gp_Pnt P, double tol) const {
    if (P.X() > aabb_xmax + tol) return false;
    if (P.Y() > aabb_ymax + tol) return false;
    if (P.Z() > aabb_zmax + tol) return false;
    if (P.X() < aabb_xmin - tol) return false;
    if (P.Y() < aabb_ymin - tol) return false;
    if (P.Z() < aabb_zmin - tol) return false;
    return true;
}

bool Space::is_point_in_shell(gp_Pnt P) {
    classifier->Perform(P, 0.001);
    return classifier->State() == TopAbs_IN;
}

gp_Pnt Space::random_point_on_cface_shell() {
    auto f = (*FirstLvl.begin())->face;
    auto v = Topo(f).vertices().First();
    return BRep_Tool::Pnt(TopoDS::Vertex(v));
}

void Space::create_solid_classifier_from_shell() { classifier = std::make_unique<BRepClass3d_SolidClassifier>(shell); }

/*void Space::create_intersector_from_cfaces(double tol) {

    TopoDS_ListOfShape L;
    for (auto &cface: cfaces)
        L.Append(cface->face);
    auto Comp = Kernel::compound_from_shape_list(L);

    intersector = std::make_unique<IntCurvesFace_ShapeIntersector>();
    intersector->Load(Comp, tol);
}*/

void Space::get_second_level() {

    for (auto &cface: FirstLvl)
        for (auto &subface: cface->Subfaces())
            SecondLvl.insert(subface);
}

void Space::calc_floor_elevation() { floor_elevation = aabb_zmin; }

void Space::CalcQuantities() {

    NetVolume = volume;
    NetFloorArea = 0;
    NetWallArea = 0;
    NetCeilingArea = 0;

    double tol = 0.3; // angular tolerance
    double inv = 0.5 * M_PI - tol;

    gp_Dir up(0, 0, 1);
    gp_Dir down(0, 0, -1);

    for (auto &cface: SecondLvl) {

        gp_Dir n = Kernel::face_normal(cface->face);
        double a = Kernel::area(cface->face);

        if (n.Z() > 0.1 && n.Angle(up) < tol && !cface->IsIfcClass("IfcVirtualElement")) {
            NetFloorArea += a;
        } else if (n.Z() < 0.1 && n.Angle(down) < tol && !cface->IsIfcClass("IfcVirtualElement")) {
            NetCeilingArea += a;
        } else if (n.Angle(up) > inv && n.Angle(down) > inv && !cface->IsOpening() && !cface->IsIfcClass("IfcVirtualElement")) {
            NetWallArea += a;
        }
    }

}

double Space::Volume() const { return volume; }

void Space::AddSecondLvlFace(cFace *cface) { SecondLvl.insert(cface); }

std::set<cFace *> Space::SecondLevel() const { return SecondLvl; }

std::set<cFace *> Space::FirstLevel() const { return FirstLvl; }

TopoDS_Shell Space::Shell() const { return shell; }

void Space::Log() const {
    std::cout << "\t" << id;
    std::cout << "\t" << std::boolalpha << is_facade;
    std::cout << "\t" << std::boolalpha << wasFlipped;
    std::cout << "\t" << std::setfill('0') << std::setw(5) << FirstLvl.size();
    std::cout << "\t" << std::setfill('0') << std::setw(3) << old_space_guids.size();
    std::cout << "\t| " << std::setfill('0') << std::setw(3) << old_space_info.size();
    std::cout << "\t" << std::setfill('0') << std::setw(3) << IfcRelDefinesByProperties.size();
    std::cout << "\t" << std::setfill('0') << std::setw(3) << (IfcRelDefinesByType == nullptr);
    std::cout << "\t" << std::setfill('0') << std::setw(3) << IfcRelContainedInSpatialStructures.size();
    std::cout << "\t" << std::setfill('0') << std::setw(10) << std::fixed << std::setprecision(3) << Volume();
    std::cout << "\t| " << std::setfill('0') << std::setw(5) << std::fixed << std::setprecision(2) << aabb_zmin;
    std::cout << "\t " << std::setfill('0') << std::setw(5) << std::fixed << std::setprecision(2) << aabb_zmax;
    std::cout << "\t" << std::setfill('0') << std::setw(4) << std::fixed << std::setprecision(2) << nbs.size();
    std::cout << "\t" << std::setfill('0') << std::setw(5) << std::fixed << std::setprecision(2) << floor_elevation;
    std::cout << "\t" << std::setfill('0') << std::setw(5) << std::fixed << std::setprecision(2) << NetVolume;
    std::cout << "\t| " << std::setfill('0') << std::setw(5) << std::fixed << std::setprecision(2) << NetFloorArea;
    std::cout << "\t" << std::setfill('0') << std::setw(5) << std::fixed << std::setprecision(2) << NetWallArea;
    std::cout << "\t" << std::setfill('0') << std::setw(5) << std::fixed << std::setprecision(2) << NetCeilingArea;
    std::cout << "\n";
}

void Space::ClearFirstLevel() { FirstLvl.clear(); }

std::set<IfcUtil::IfcBaseClass *> Space::ProductsInSpace() {

    std::set<IfcUtil::IfcBaseClass *> s;

    for (auto &cface: FirstLvl) {
        if (cface->RelProduct() != nullptr)
            s.insert(cface->RelProduct()->ifcproduct);
        if (cface->IsOpening() && cface->Ancestor()->Opening() != nullptr)
            s.insert(cface->Ancestor()->Opening());
    }

    return s;
}
