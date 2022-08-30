// Copyright 2022 Eric Fichter
#include "Kernel.h"

// Concerns cFace objects only. Small functions with independent processing of faces.

void Kernel::correct_face_normals(std::list<cFace> &cFaces) {

    auto start = std::chrono::high_resolution_clock::now();

#ifdef PARALLEL_PROCESSING
    tbb::parallel_for_each(cFaces.begin(), cFaces.end(), [&](cFace &cface) { cface.CheckSetFaceNormal(); });
#else
    for (auto &cface: cFaces)
        cface.CheckSetFaceNormal();
#endif

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Correct face normals", std::to_string(cFaces.size()));
}

bool Kernel::check_faces(std::list<cFace> &cFaces, bool mark_as_trash, double l_warn, double l_crit) {

    auto start = std::chrono::high_resolution_clock::now();

    bool b = true;

    for (auto &cface: cFaces) {
        std::string guid = cface.Ancestor() == nullptr ? "No Ancestor" : cface.IfcGuid();
        if (!check_face(cface.face, cface.ID(), cface.Info(), l_warn, l_crit)) {
            b = false;
            if (mark_as_trash) cface.SetIsTrash(true);
        }
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Check faces", std::to_string(cFaces.size()));

    return b;
}

void Kernel::calculate_points_on_faces(std::list<cFace> &cFaces_2ndLvl) const {

    auto start = std::chrono::high_resolution_clock::now();

#ifdef PARALLEL_PROCESSING
    std::vector<cFace *> V;
    V.reserve(cFaces_2ndLvl.size());
    for (auto &cface: cFaces_2ndLvl)
        V.push_back(&cface);

#pragma omp parallel for default(none) shared(V) schedule (static) num_threads(num_threads)
    for (auto cface: V)
        cface->SetPointOnFace(point_on_face(cface->face, cface->FixedFaceNormal(), cface->Info()));
#else
    for (auto &cface: cFaces_2ndLvl)
        cface.SetPointOnFace(point_on_face(cface.face, cface.FixedFaceNormal()));
#endif

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Calculate points on faces", std::to_string(cFaces_2ndLvl.size()));
}

void Kernel::create_missing_parent_faces(std::list<cFace> &cFaces, std::list<oFace> &orig_faces, unsigned int &fid) {

    // If opening is positioned in embrasure of wall and couldn't be unified by a wall face, the opening face is missing a parent face on the one hand and the space is not closed by non-opening faces.
    // A wall face must be added to have parent relationship.

    auto start = std::chrono::high_resolution_clock::now();

    for (auto &cface: cFaces) {
        if (!cface.IsOpening() || cface.Parent() != nullptr) continue;

        std::cout << "[Info] Opening face has no parent face " << cface.Info() << "\t | " << cface.Ancestor()->Opening()->data().getArgument(0)->toString() << "\t" << cface.Ancestor()->Opening()->declaration().name();

        if (cface.RelProduct()->ifcproduct == cface.Ancestor()->Opening()) { // happens when ifcwindow and ifcdoor are not integrated into wall faces but stay their own object
            std::cout << "\n";
            continue;
        }

        orig_faces.push_back(*cface.Ancestor());
        orig_faces.back().SetOpening(nullptr);

        cFaces.push_back(cface);
        cFaces.back().face = TopoDS::Face(shape_copy(cface.face)); // copy, otherwise the intersector of the space will only detect original face
        cFaces.back().SetAncestor(&orig_faces.back());
        cFaces.back().SetID(fid);
        fid++;

        cFaces.back().RelSpace()->FirstLvl.insert(&cFaces.back());

        cface.SetParent(&cFaces.back());

        std::cout << "\t | Add face " << cFaces.back().Info() << ".\n";
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Create missing parent faces", std::to_string(orig_faces.size()) + ", " + std::to_string(cFaces.size()));
}

void Kernel::trash_2b_opening_faces(std::list<cFace> &cFaces) {

    auto start = std::chrono::high_resolution_clock::now();

    for (auto &cface: cFaces)
        if (cface.IsOpening() && cface.SBType() != SB_TYPE_2A) {
            std::cout << "[Info] Unpaired opening face: " << cface.Info() << "\n";
            cface.SetIsTrash(true);

            if (cface.Superface()->InternalOrExternal() == SB_IE_NOTDEFINED)
                cface.Superface()->SetInternalOrExternal(cface.RelSpace()->is_facade ? SB_IE_EXTERNAL : SB_IE_INTERNAL);
        }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Trash 2b opening faces", std::to_string(cFaces.size()));
}

void Kernel::reduce_edges_cFaces(std::list<cFace> &cFaces) {

    // Despite unify of adjacent cfaces, there are colinear adjecent edges that should be combined to reduce edge number for output
    auto start = std::chrono::high_resolution_clock::now();

#ifdef PARALLEL_PROCESSING
    tbb::parallel_for_each(cFaces.begin(), cFaces.end(), [&](cFace &cface) { cface.ReduceEdges(); });
#else
    for (auto &cface: cFaces)
        cface.ReduceEdges();
#endif

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Reduce edges in faces ", std::to_string(cFaces.size()));
}

void Kernel::identify_decoupled_offset_faces(std::list<cFace> &cFaces) {

    // identify offset faces, that are not connected to the shell of the face, they were offsetted from
    auto start = std::chrono::high_resolution_clock::now();

    for (auto &cface: cFaces)
        if (cface.IsOffset() && !cface.IsConnectedToShell()) cface.SetIsTrash(true);

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Identify decoupled offset faces", std::to_string(cFaces.size()));
}

void Kernel::clear_cface_maps(std::list<cFace> &cFaces) {

    auto start = std::chrono::high_resolution_clock::now();

    for (auto &cface: cFaces)
        cface.ClearMaps();

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Clear face maps", std::to_string(cFaces.size()));
}

void Kernel::update_adjacence_map_after_unification(std::list<cFace> &cFaces) {

    // update adjacence info for all non-trash (no unifyed) faces. because of inconsistent merging between faces, half edge data cannot be used anymore.
    cfaceSetMap M;

    for (auto &cface: cFaces) {

        cFace *key = cface.UnifyingFace() == nullptr ? &cface : cface.UnifyingFace();

        for (auto &nb: cface.AdjacentFacesList()) {
            cFace *u = nb->UnifyingFace() == nullptr ? nb : nb->UnifyingFace();
            M[key].insert(u);
        }
    }

    for (auto &m: M) {
        auto cface = m.first;
        auto &nbs = m.second;
        nbs.erase(cface);
        cface->ClearMaps();
        for (auto &nb: nbs)
            cface->adjacentfaces[0].push_back(nb);
    }
}

void Kernel::check_corresponding_face_pairs(std::list<cFace> &cFaces) {

    auto start = std::chrono::high_resolution_clock::now();

    unsigned int i = 0;

    for (auto &cface: cFaces)
        if (cface.Corresponding() != nullptr) {
            if (cface.Corresponding()->IsTrash()) std::cerr << "[Warning] Corresponding face " << cface.Corresponding()->Info() << " of face " << cface.Info() << " is trash." << std::endl;

            if (!cface.CheckCorresponding(0.001, 3.05)) {
                i++;
                // nullify corresponding's attributes if the corresponding face points to cface
                if (cface.Corresponding() == &cface) {
                    cface.Corresponding()->SetCorresponding(nullptr);
                    cface.Corresponding()->SetSBType(SB_TYPE_2B);
                }
                // nullify own attributes
                cface.SetCorresponding(nullptr);
                cface.SetSBType(SB_TYPE_2B);
            }
        } else if (cface.SBType() == SB_TYPE_2A) {
            std::cerr << "[Warning] No corresponding face for 2a face " << cface.Info() << "." << std::endl;
            cface.SetSBType(SB_TYPE_2B);
            i++;
        }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Check corresponding faces ", std::to_string(i) + "/" + std::to_string(cFaces.size()));
}

void Kernel::check_fixed_normal(std::list<cFace> &cFaces) {

    auto start = std::chrono::high_resolution_clock::now();

    unsigned int n = 0;

    for (auto &cface: cFaces)
        if (cface.NormalStatus() == FACE_NORMAL_UNKNOWN) {
            std::cerr << "[Warning] Face normal is unknown " << cface.Info() << "." << std::endl;
            n++;
        }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Check fixed normal", std::to_string(n) + "/" + std::to_string(cFaces.size()));
}

void Kernel::check_redundant_faces(std::list<cFace> &cFaces, std::list<cFace> &innerFaces) {

    auto start = std::chrono::high_resolution_clock::now();

    unsigned int n = 0;
    std::unordered_map<unsigned int, cFace *> m;

    for (auto &innerFace: innerFaces)
        m[innerFace.FaceID()] = &innerFace;

    for (auto &cface: cFaces) {

        if (m.find(cface.FaceID()) == m.end()) continue;

        const auto &innerFace = m[cface.FaceID()];
        // std::cout << "[Info] Duplicate face " << cface.ID() << "\t" << cface.IfcClass() << "\t" << innerFace->ID() << "\t" << innerFace->IfcClass() << "." << std::endl;

        if (innerFace->IfcProduct() == cface.IfcProduct()) {
            if (cface.IsOpening() || innerFace->IsOpening())
                std::cout << "[Info] Duplicate face from same product " << cface.Info() << "\t" << innerFace->Info() << ".\n";
            else {
                std::cerr << "[Warning] Duplicate face from same product " << cface.Info() << "\t" << innerFace->Info() << "." << std::endl;
                n++;
            }
        }
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Check redundant faces", std::to_string(n) + "/" + std::to_string(cFaces.size()));
}

void Kernel::check_closed_space_edge_id(std::list<Space> &spaces) {

    auto start = std::chrono::high_resolution_clock::now();

    for (auto &space: spaces) {

        // closed definition in this fct: every edge must occur in even numbers
        std::unordered_map<unsigned int, std::list<TopAbs_Orientation>> ids;

        for (auto &cface: space.FirstLvl)
            for (auto &edge: Topo(cface->face).edges())
                ids[hash(edge)].push_back(edge.Orientation());

        for (auto &id: ids)
            if (id.second.size() % 2 != 0) {
                if (id.second.size() == 1 && id.second.front() == TopAbs_INTERNAL)
                    std::cout << "[Info] Space " << space.id << "\t (" << std::round(space.Volume()) << ") contains seam edge (will be removed later). " << id.second.size() << " " << id.first << "\n";
                else {
                    std::cerr << "[Warning] Space " << space.id << "\t (" << std::round(space.Volume()) << ") is not closed. " << id.second.size() << " " << id.first << std::endl;
                    break;
                }
            }
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Check space closed", std::to_string(spaces.size()));
}

void Kernel::check_adjacency_self_reference(std::list<cFace> &cFaces) {

    auto start = std::chrono::high_resolution_clock::now();

    unsigned int n = 0;

    for (auto &cface: cFaces) {
        auto L = cface.AdjacentFacesList();
        if (L.find(&cface) != L.end()) {
            std::cerr << "[Warning] Face has self-reference in adjacent faces " << cface.Info() << "." << std::endl;
            n++;
        }
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Check adjacency self reference", std::to_string(n) + "/" + std::to_string(cFaces.size()));
}

void Kernel::check_holes(std::list<cFace> &cFaces) {

    auto start = std::chrono::high_resolution_clock::now();

    unsigned int n = std::count_if(cFaces.begin(), cFaces.end(), [&](cFace &cface) { return cface.HasHoles(); });

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Check for holes in faces ", std::to_string(n) + "/" + std::to_string(cFaces.size()));
}

void Kernel::check_loops(std::list<cFace> &cFaces, double tol) {

    auto start = std::chrono::high_resolution_clock::now();

    unsigned int n = 0;

    for (auto &cface: cFaces)
        for (auto &wire: Topo(cface.face).wires()) {
            std::vector<gp_Pnt> Pnts;
            if (!create_connection_geometry_loop(&cface, Topo(wire).ordered_vertices_of_wire(), Pnts, tol)) {
                std::cerr << "[Warning] Face " << cface.Info() << " has bad wires!" << std::endl;
                n++;
            }
        }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Check loops of faces ", std::to_string(n) + "/" + std::to_string(cFaces.size()));
}

void Kernel::check_duplicate_vertices(std::list<cFace> &cFaces) {

    auto start = std::chrono::high_resolution_clock::now();

    for (auto &cface: cFaces) {

        TopoDS_ListOfShape verts = Topo(cface.face).vertices();
        std::set<unsigned int> unique;
        std::list<unsigned int> all, seam;

        for (auto &v: verts) {
            unique.insert(hash(v));
            all.push_back(hash(v));
        }

        // vertices of seam edges
        for (auto &e: Topo(cface.face).seam_edges())
            for (auto &v: Topo(e).vertices())
                seam.push_back(hash(v));

        if (verts.Size() > 2 * unique.size())
            std::cerr << "[Warning] Face " << cface.Info() << " contains duplicate vertices! " << verts.Size() << " " << unique.size() << std::endl;

        for (auto &u: unique) {
            auto n = std::count(all.begin(), all.end(), u);
            if (n > 2) {
                auto s = std::count(seam.begin(), seam.end(), u);
                if (n - s > 2)
                    std::cerr << "[Warning] Face " << cface.Info() << " seems to have pointwise self-intersection (maybe inner wire intersects with outer wire in one vertex) " << u << std::endl;
                // else it's just the seam edge's vertices that are contributing to this disbalance
            }
        }
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Check for duplicate vertices", std::to_string(cFaces.size()));
}

void Kernel::check_superface(std::list<cFace> &cFaces) {

    auto start = std::chrono::high_resolution_clock::now();

    unsigned int n = 0;

    for (auto &cface: cFaces)
        if (cface.Superface() != nullptr)
            if (cface.Superface()->IsTrash()) {
                std::cerr << "[Warning] Superface of face " << cface.Info() << " is null or trash!" << std::endl;
                n++;
            }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Check superface", std::to_string(n) + "/" + std::to_string(cFaces.size()));
}

void Kernel::check_super_and_subfaces(std::list<cFace> &cFaces) {

    auto start = std::chrono::high_resolution_clock::now();

    for (auto &cface: cFaces) {

        if (cface.Superface() != nullptr)
            if (cface.Superface()->IsTrash())
                std::cerr << "[Warning] Superface of face " << cface.Info() << " is null or trash!" << std::endl;

        if (cface.Subfaces().empty())
            std::cout << "[Info] No subfaces for face " << cface.Info() << " (Probably because its subface was a 2b opening face which was removed)!\n";
        else
            for (auto &sub: cface.Subfaces())
                if (sub == nullptr)
                    std::cerr << "[Warning] Subface of face " << cface.Info() << " is null!" << std::endl;
                else if (sub->IsTrash())
                    std::cerr << "[Warning] Subface of face " << cface.Info() << " is trash!" << std::endl;
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Check super and subfaces", std::to_string(cFaces.size()));
}

void Kernel::check_convex(std::list<cFace> &cFaces) {

    auto start = std::chrono::high_resolution_clock::now();

    unsigned int n = std::count_if(cFaces.begin(), cFaces.end(), [&](cFace &cface) { return cface.IsConvex(); });

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Check for convexity of faces", std::to_string(n) + "/" + std::to_string(cFaces.size()));
}

void Kernel::change_cface_sbtype(std::list<cFace> &cFaces_1stLvl, std::list<cFace> &cFaces_2ndLvl) {

    auto start = std::chrono::high_resolution_clock::now();

    // inherit internal/external attribute to first level
    for (auto &cface: cFaces_2ndLvl) {
        if (cface.Superface() == nullptr) continue;
        if (cface.SBType() != SB_TYPE_2B) // 2b don't know a thing about internal/external
            cface.Superface()->SetInternalOrExternal(cface.InternalOrExternal());
    }

    for (auto &cface: cFaces_1stLvl) {
        cface.SetSBType(SB_TYPE_1);

        if (!cface.Subfaces().empty()) {
            cface.SetPhysicalOrVirtual(cface.Subfaces().front()->PhysicalOrVirtual());

            if (cface.InternalOrExternal() == SB_IE_NOTDEFINED)
                cface.SetInternalOrExternal(cface.Subfaces().front()->InternalOrExternal());
        }
        //else
        // just a problem of my order of functions. if 2nd lvl face is removed e.g. because its an unpaired opening its information about enums is lost for this function
        // std::cout << "[Info] Face has no subfaces. " << cface.Info() << "!" << std::endl;

        if (cface.IsIfcClass("IfcVirtualElement") && cface.PhysicalOrVirtual() != SB_PV_VIRTUAL)
            std::cerr << "[Warning] Non-virtual face related to IfcVirtualElement " << cface.Info() << "!" << std::endl;
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Set enumerations", std::to_string(cFaces_1stLvl.size()));
}

void Kernel::check_parent_faces(std::list<cFace> &cFaces) {

    auto start = std::chrono::high_resolution_clock::now();

    unsigned int i = 0;

    for (auto &cface: cFaces)
        if (cface.Parent() != nullptr) {
            if (!cface.CheckParent(1.0e-4, 0.001)) {
                cface.SetParent(nullptr);
                i++;
            }
        } else {
            if (cface.IsOpening()) std::cerr << "[Warning] Opening face has no parent " << cface.Info() << "!" << std::endl;
            else if (cface.IsIfcClass("IfcDoor") || cface.IsIfcClass("IfcWindow")) std::cout << "[Info] Door/opening face has no parent " << cface.Info() << "!\n";
        }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Check parent faces", std::to_string(i) + "/" + std::to_string(cFaces.size()));
}

void Kernel::check_adjacency_null(std::list<cFace> &cFaces) {

    auto start = std::chrono::high_resolution_clock::now();

    for (auto &cface: cFaces)
        for (auto &nb: cface.AdjacentFacesList()) {
            unsigned int id = nb->ID();
            bool trash = nb->IsTrash();
            if (trash) std::cerr << "[Warning] Neighbor " << id << " of face " << cface.Info() << " is trash!" << std::endl;
        }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Check correct pointers in adjacency", std::to_string(cFaces.size()));
}

void Kernel::check_space_behind(std::list<cFace> &cFaces) {

    auto start = std::chrono::high_resolution_clock::now();

    unsigned int n = std::count_if(cFaces.begin(), cFaces.end(), [&](cFace &cface) { return !cface.CheckSpaceBehind(); });

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Check space behind", std::to_string(n) + "/" + std::to_string(cFaces.size()));
}

void Kernel::check_pointers(std::list<cFace> &cFaces, bool secondlvl) {

    auto start = std::chrono::high_resolution_clock::now();

    for (auto &cface: cFaces)
        cface.CheckPointers(secondlvl);

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Check pointers", std::to_string(cFaces.size()));
}

void Kernel::check_enum_intext(std::list<cFace> &cFaces) {

    auto start = std::chrono::high_resolution_clock::now();

    unsigned int i = 0;

    for (auto &cface: cFaces)
        if (cface.InternalOrExternal() == SB_IE_NOTDEFINED)
            std::cerr << "[Warning] Undefined internal/external enum. " << cface.Info() << std::endl;

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Check internal/external enum", std::to_string(i) + "/" + std::to_string(cFaces.size()));
}

void Kernel::check_redundant_cFaces(std::list<cFace> &cFaces) {

    auto start = std::chrono::high_resolution_clock::now();

    std::map<std::tuple<Product *, unsigned int, unsigned int, bool>, std::list<cFace *>> M; // product, shell id, face id
    for (auto &cface: cFaces)
        M[std::make_tuple(cface.RelProduct(), cface.Ancestor()->ShellID(), cface.FaceID(), cface.IsOffset())].push_back(&cface);

    for (auto &m: M)
        if (m.second.size() > 1) // can occur, when two faces of the shell are close to each other, resulting in a coplanar face paur while fusing
            std::cout << "[Warning] Duplicate cface within product and shell " << std::get<0>(m.first)->guid << "\t" << std::get<1>(m.first) << "\t" << std::get<2>(m.first) << "\t" << std::get<3>(m.first) << "\tSIZE " << m.second.size() << "\n";

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Check for redundant cFaces", std::to_string(cFaces.size()));
}

void Kernel::update_half_edges(std::list<cFace> &cFaces) {

    auto start = std::chrono::high_resolution_clock::now();

    unsigned int count = 0;
    for (auto &cface: cFaces) {
        cface.UpdateHalfEdges();
        count++;
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Save half-edge information (skip seam edges)", std::to_string(count));
}

void Kernel::remove_seam_edges_from_faces(std::list<cFace> &cFaces) {

    auto start = std::chrono::high_resolution_clock::now();

#ifdef PARALLEL_PROCESSING
    tbb::parallel_for_each(cFaces.begin(), cFaces.end(), [&](cFace &cface) { cface.RemoveSeamEdgesFromFace(); });
#else
    for (auto &cface: cFaces)
        cface.RemoveSeamEdgesFromFace();
#endif

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Remove seam edges from faces", std::to_string(cFaces.size()));
}

void Kernel::log_cFaces(const std::list<cFace> &cFaces) {
    unsigned int i = 0;
    for (const auto &cface: cFaces) {
        cface.Log(i);
        i++;
    }
}

unsigned int Kernel::number_of_faces_of_unknown_orientation(const std::list<cFace> &cFaces) {
    unsigned int n = std::count_if(cFaces.begin(), cFaces.end(), [&](const cFace &cface) { return cface.NormalStatus() != FACE_NORMAL_KNOWN; });
    std::cout << "[Info] " << n << " of " << cFaces.size() << " faces have an unknown orientation." << std::endl;
    return n;
}

void Kernel::remove_inner_window_and_door_faces(std::list<cFace> &cFaces) {

    auto start = std::chrono::high_resolution_clock::now();

    auto n = cFaces.size();

    cFaces.remove_if([](cFace &cface) { return cface.IsOpening() || cface.IsIfcClass("IfcWindow") || cface.IsIfcClass("IfcDoor"); });

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Remove inner and coplanar window and door faces", std::to_string(n) + "/" + std::to_string(cFaces.size()));
}

TopoDS_Shape Kernel::mold(const std::list<Space> &spaces) {

    auto start = std::chrono::high_resolution_clock::now();

    TopoDS_ListOfShape shells;

    for (const auto &space: spaces) {
        std::set < cFace * > L;
        for (const auto &cface: space.SecondLevel())
            if (!cface->IsVirtual()) {
                L.insert(cface);

//                if (cface->SBType() == SB_TYPE_2A) {
//                    TopoDS_Shape P = BRepPrimAPI_MakePrism(cface->face, gp_Vec(cface->FaceNormal()).Scaled(cface->DistanceToSpaceBoundaryBehind() + 0.01).Reversed()).Shape();
//                    TopoDS_Shell p = TopoDS::Shell(Topo(P).shells().First());
//                    if (volume(p) > 0) p.Complement();
//                    prisms.Append(p);
//                }
            }

        BRepBuilderAPI_Sewing sew;
        for (const auto &cface: L)
            sew.Add(cface->face);
        sew.SetNonManifoldMode(true);
        sew.SetTolerance(1.0e-3);
        sew.Perform();

        for (const auto &s: Topo(sew.SewedShape()).shells()) {
            TopoDS_Shape shell = TopoDS::Shell(s);
            if (space.is_facade && volume(shell) < 0) shell.Complement();
            else if (!space.is_facade && volume(shell) > 0) shell.Complement();
            shells.Append(shell);
        }
    }

    BRepBuilderAPI_MakeSolid m;
    for (const auto &s: shells)
        m.Add(TopoDS::Shell(s));
    TopoDS_Solid solid = m.Solid();

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Create mold", std::to_string(Topo(solid).shells().Size()));

    return solid;
}
