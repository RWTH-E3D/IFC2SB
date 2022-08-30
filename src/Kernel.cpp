// Copyright 2022 Eric Fichter
#include "Kernel.h"
#include "remove_holes_SB.h"

// Concerns mostly cFace and Space objects only. Longer functions with deeper geometric processing.

Kernel::Kernel(const unsigned int _num_threads) : num_threads(_num_threads) {

    n_digits_round_double_ifc_write = 8; // if number is too low, rounding, processing of sb connection geometry and recreation from ifc file can create wrong numbers and therefore holes and overlaps
    round_double_ifc_write = pow(10, n_digits_round_double_ifc_write);
    round_ifc_write_inv = 1.0 / round_double_ifc_write;
    conv_fctr = 1.0;

    ifcSchema = IFC4;

    shapeEnum_to_string = {{TopAbs_COMPOUND,  "TopAbs_COMPOUND"},
                           {TopAbs_COMPSOLID, "TopAbs_COMPSOLID"},
                           {TopAbs_SOLID,     "TopAbs_SOLID"},
                           {TopAbs_SHELL,     "TopAbs_SHELL"},
                           {TopAbs_FACE,      "TopAbs_FACE"},
                           {TopAbs_WIRE,      "TopAbs_WIRE"},
                           {TopAbs_EDGE,      "TopAbs_EDGE"},
                           {TopAbs_VERTEX,    "TopAbs_VERTEX"},
                           {TopAbs_SHAPE,     "TopAbs_SHAPE"}};
}

bool Kernel::fuse_original_faces(TopoDS_Shape &fuse, std::list<oFace> &orig_faces, std::list<cFace> &cFaces, double fuzzy_tol, unsigned int &fid) {

    auto start = std::chrono::high_resolution_clock::now();

    // Attention: Creating copies of shapes, seems to raise complexity of shape (see dump) and therefore lowers speed of fusing!
    // Same goes for shape that have a triangulation

    // populate builder
    BOPAlgo_Builder builder;
    for (const auto &f: orig_faces)
        builder.AddArgument(f.face);

    builder.SetNonDestructive(false); // Safe input shapes option allows preventing modification of the input shapes
    builder.SetCheckInverted(false);  // Enables/Disables the check of the input solids for inverted status.
    // print(precision_Confusion()); // Returns the recommended precision value when checking coincidence of two points in real space.
    builder.SetFuzzyValue(fuzzy_tol); // Fuzzy option allows setting the additional tolerance for the operation
    builder.SetRunParallel(true);  // Parallel processing option allows running the algorithm in parallel mode
    //builder.SetParallelMode(true);
    builder.SetUseOBB(true);  // Usage of Oriented Bounding Boxes in the operation
    // from OCC.Core.BOPAlgo import BOPAlgo_GlueShift
    // builder.SetGlue(BOPAlgo_GlueShift);  // Gluing option allows speeding-up the intersection of the arguments
    builder.SetToFillHistory(true);  // Allows disabling the history collection.

    builder.Perform();

    if (builder.HasWarnings()) {
        std::cerr << "Warnings:" << std::endl;
        builder.DumpWarnings(std::cerr);
    }
    if (builder.HasErrors()) {
        std::cerr << "Errors:" << std::endl;
        builder.DumpErrors(std::cerr);
    }

    if (!builder.HasModified()) {
        std::cerr << "Nothing was modified." << std::endl;
        return false;
    }

    for (auto &orig_face: orig_faces) {

        bool isDeleted = builder.IsDeleted(orig_face.face); // in case face was deleted, don't add to cfaces (maybe add the generated?). Without ignoring, the cface.face will not be present is fuse shape

        if (isDeleted) {
            auto Generated = builder.Generated(orig_face.face);
            std::cerr << "[Error] Face was deleted! " << hash(orig_face.face) << "\t" << orig_face.IfcGuid() << "\t" << orig_face.IfcClass() << "\t" << isDeleted << "\t" << Generated.Size() << std::endl;
            continue;
        }

        // because faces are the arguments in fuse, only entities TopAbs_VERTEX and TopAbs_EDGE can be generated
        // auto Generated = builder.Generated(orig_face.face);
        // if (!Generated.IsEmpty())
        //     for (const TopoDS_Shape &shape : Generated)
        //         std::cerr << hash(orig_face.face) << "\t" << hash(shape) << "\t" << shapeEnum_to_string[shape.ShapeType()] << "\t" << shape.Orientation() << std::endl;

        auto L = builder.Modified(orig_face.face);

        if (L.IsEmpty()) {  // face was not modified
            cFaces.emplace_back(orig_face.face, &orig_face, fid); // original TopoDS_Face, pointer on oFace
            fid++;
        } else // face was modified into one or more new faces
            for (const TopoDS_Shape &shape: L) {
                cFaces.emplace_back(TopoDS::Face(shape), &orig_face, fid); // new face, pointer on oFace
                fid++;
            }
    }

    fuse = builder.Shape();

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Fuse original faces", std::to_string(Topo(fuse).faces().Size()) + ", " + std::to_string(cFaces.size()));

    return true;

}

void Kernel::identify_enclosed_faces(const TopoDS_Shape &fuse, std::list<cFace> &cFaces) {

    // enclosed faces will be flagged as trash. Attention: enclosed faces are faces that are inside another shell and adjacent to at least one of the shell faces
    // all other faces have to be found using hanging methods

    auto start = std::chrono::high_resolution_clock::now();

    std::unordered_map<unsigned int, std::list<cFace *>> id2cFaces;
    for (auto &cFace: cFaces) {
        if (id2cFaces.find(cFace.FaceID()) != id2cFaces.end())
            id2cFaces[cFace.FaceID()].emplace_back(&cFace);
        else
            id2cFaces[cFace.FaceID()] = {&cFace};
    }

    TopTools_IndexedDataMapOfShapeListOfShape edgeFaceMap;
    TopExp::MapShapesAndAncestors(fuse, TopAbs_EDGE, TopAbs_FACE, edgeFaceMap);

    TopExp_Explorer Ex;
    for (auto &cface: cFaces)
        if (cface.IsEnclosed(edgeFaceMap, id2cFaces, cFaces)) {
            cface.SetIsTrash(true);
            cface.SetIsInner(true);
        }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Identify enclosed faces", std::to_string(cFaces.size()));

}

void Kernel::remove_trash(std::list<cFace> &cFaces, std::list<cFace> &cFaces_trash) {

    auto start = std::chrono::high_resolution_clock::now();

    cFaces.remove_if([](cFace &cface) { return cface.IsTrash(); });

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Remove trash faces", std::to_string(cFaces.size()));
}

void Kernel::remove_corresponding_trash(std::list<cFace> &cFaces) {

    // if no angle criterion is used for corresponding face pairs, after some algorithms (e.g. when split_corresponding_equivalent() is used)
    // there can be a mismatch for situations where there is no exact 1 to 1 but 1 to 2 correspondence.

    auto start = std::chrono::high_resolution_clock::now();

    unsigned int i = 0;

    for (auto &cface: cFaces) {
        if (cface.IsTrash()) continue;

        if (cface.Corresponding() != nullptr)
            if (cface.Corresponding()->IsTrash()) {
                i++;
                std::cerr << "[Info] Corresponding face " << cface.Corresponding()->Info() << " of face " << cface.Info() << " is trash. Remove link." << std::endl;
                cface.SetCorresponding(nullptr);
                cface.SetSBType(SB_TYPE_2B);
            }
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Remove trash faces from corresponding attribute", std::to_string(i) + "/" + std::to_string(cFaces.size()));
}

void Kernel::identify_duplicate_faces(std::list<cFace> &cFaces, std::unordered_map<std::string, unsigned int> ranks) {

    // identify mat-mat faces by checking face normals of faces with same id (which were created during fuse, when to faces were coplanar).

    auto start = std::chrono::high_resolution_clock::now();

    // find duplicate faces
    std::unordered_map<unsigned int, std::list<cFace *>> id2pairs;

    for (auto &cface: cFaces) {
        unsigned int id = cface.FaceID();
        if (id2pairs.find(id) == id2pairs.end()) {
            std::list<cFace *> l = {&cface};
            id2pairs[id] = l;
        } else
            id2pairs[id].push_back(&cface);
    }

    // check for mat-mat faces using face normals
    for (auto &id2pair: id2pairs) {

        // skip faces that have no duplicate ids
        if (id2pair.second.size() == 1) continue;

        // get all mat faces of the duplicate faces
        std::list<cFace *> mat_faces, off_faces;

        for (auto &it_face: id2pair.second)
            it_face->IsOffset() ? off_faces.push_back(it_face) : mat_faces.push_back(it_face);

        if (mat_faces.size() == 1) { // if only one mat face paired with offset faces, then delete all offsets and keep mat
            for (auto &off_face: off_faces)
                off_face->SetIsTrash(true);

        } else if (mat_faces.size() > 1) { // if more than one mat face paired with offset faces, then ...

            // delete all offsets and ...
            for (auto &off_face: off_faces)
                off_face->SetIsTrash(true);

            // check for same direction of mat faces
            bool all_same = compare_normals(mat_faces);

            if (all_same) { // if all same normal, keep one

                // with ranks
                std::multimap<int, cFace *> rank_to_mat_face = sort_mat_mat_faces(mat_faces, ranks);
                auto it_face = rank_to_mat_face.begin();
                for (std::advance(it_face, 1); it_face != rank_to_mat_face.end(); ++it_face) { // skip first face
                    it_face->second->SetIsTrash(true);
                    it_face->second->SetIsCoplanar(true);
                }

                // without ranks
                // std::list<cFace *>::iterator it_face = mat_faces.begin();
                //for (std::advance(it_face, 1); it_face != mat_faces.end(); ++it_face) // skip first face
                //    (*it_face)->isTrash = true;

            } else // else, delete all of them
                for (auto &mat_face: mat_faces) {
                    mat_face->SetIsTrash(true);
                    mat_face->SetIsCoplanar(true);
                }

        } else { // if no mat faces and only offset faces

            // check for same direction of off faces
            bool all_same = compare_normals(off_faces);

            if (all_same) { // if all same normal, keep one
                auto it_face = off_faces.begin();
                for (std::advance(it_face, 1); it_face != off_faces.end(); ++it_face) // skip first face
                    (*it_face)->SetIsTrash(true);
            } else // else, delete all of them
                for (auto &off_face: off_faces)
                    off_face->SetIsTrash(true);

        }
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Identify duplicate faces", std::to_string(cFaces.size()));
}

bool Kernel::compare_normals(const std::list<cFace *> &cFaces) {

    gp_Dir n1, n2;

    for (const auto &face_1: cFaces)
        for (const auto &face_2: cFaces) {
            if (face_1 == face_2) continue;
            n1 = face_1->FixedFaceNormal();
            n2 = face_2->FixedFaceNormal();
            if (n1.Angle(n2) > 0.01) {
                return false;
            }
        }

    return true; // returns true, if normals are all same
}

std::multimap<int, cFace *> Kernel::sort_mat_mat_faces(const std::list<cFace *> &mat_faces, std::unordered_map<std::string, unsigned int> &ranks) {

    std::multimap<int, cFace *> rank_to_mat_face;

    for (const auto &mat_face: mat_faces) {

        std::string ifcClass = mat_face->IfcClass();
        unsigned int rank;

        if (ranks.find(ifcClass) != ranks.end())
            rank = ranks[ifcClass];
        else
            rank = 999;

        rank_to_mat_face.insert(std::pair<int, cFace *>(rank, mat_face));
    }

    return rank_to_mat_face;
}

void Kernel::check_duplicate_faces(std::list<cFace> &cFaces) {

    auto start = std::chrono::high_resolution_clock::now();

    std::set<unsigned int> unique_ids;
    for (auto &cface: cFaces) {
        if (unique_ids.find(cface.FaceID()) != unique_ids.end())
            std::cerr << "[Error] TopoDS_Face ID twice in cFaces! Indicates OCC fuse error!" << cface.Info() << std::endl;
        unique_ids.insert(cface.FaceID());
    }

    unique_ids.clear();
    for (auto &cface: cFaces) {
        if (unique_ids.find(cface.ID()) != unique_ids.end())
            std::cerr << "[Error] cFace ID twice in cFaces!" << cface.Info() << std::endl;
        unique_ids.insert(cface.ID());
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Check duplicate face ids", std::to_string(cFaces.size()));
}

void Kernel::update_face_adjacencies(std::list<cFace> &cFaces, const TopoDS_Shape &fuse) {

    // Attention: Seam edges were not saved by update_half_edges.
    // The seam edge from one face, could be a normal edge for another face.
    // So a face could refer to another face but not the other way, because edge isn't existent in second face
    // save adjency will consider this and will remove adjacent from edge

    auto start = std::chrono::high_resolution_clock::now();

    std::unordered_map<unsigned int, cFace *> id2cFace = setup_face_cFace_map(cFaces);

    TopTools_IndexedDataMapOfShapeListOfShape edgeFaceMap;
    TopExp::MapShapesAndAncestors(fuse, TopAbs_EDGE, TopAbs_FACE, edgeFaceMap);

    for (auto &cface: cFaces)
        cface.UpdateFaceAdjacencies(edgeFaceMap, id2cFace);

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Save face adjacency information", std::to_string(cFaces.size()));

}

std::unordered_map<unsigned int, cFace *> Kernel::setup_face_cFace_map(std::list<cFace> &cFaces) {

    std::unordered_map<unsigned int, cFace *> id2cFace;
    for (auto &cface: cFaces) {
        if (id2cFace.find(cface.FaceID()) != id2cFace.end())
            std::cerr << "Error: Two faces with same ID in cFaces list! " << cface.Info() << std::endl;
        id2cFace[cface.FaceID()] = &(cface);
    }
    return id2cFace;
}

void Kernel::remove_adjacency_by_orientation(std::list<cFace> &cFaces) {

    auto start = std::chrono::high_resolution_clock::now();

    for (auto &cface: cFaces)
        cface.RemoveWrongOrientedAdjacentFaces();

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Remove adjacence information between faces by orientation", std::to_string(cFaces.size()));
}

void Kernel::identify_hanging_faces(std::list<cFace> &cFaces) {

    auto start = std::chrono::high_resolution_clock::now();

    std::string s;

    while (identify_hanging_faces_while(cFaces))
        s += "#";

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Identify hanging faces " + s, std::to_string(cFaces.size()));

}

bool Kernel::identify_hanging_faces_while(std::list<cFace> &cFaces) {

    bool found_hanging_faces = false;
    for (auto &cface: cFaces)
        if (cface.SetHanging())
            found_hanging_faces = true;

    return found_hanging_faces;
}

void Kernel::remove_trash_and_face_adjacency(std::list<cFace> &cFaces) {

    auto start = std::chrono::high_resolution_clock::now();

    std::list<cFace *> del;

    for (auto &cface: cFaces) { // remove adjacency information pointing to soon deleted cfaces, before removing trash cfaces from cFaces preventing invalid memory read
        if (cface.IsTrash()) continue;
        cface.RemoveTrashedAdjacentFaces();
    }

    cFaces.remove_if([](cFace &cface) { return cface.IsTrash(); });

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Remove trash faces and adjacency information of trashed faces", std::to_string(cFaces.size()));
}

void Kernel::find_spaces(std::list<cFace> &cFaces, std::list<Space> &spaces, unsigned int &space_id) {

    auto start = std::chrono::high_resolution_clock::now();

    while (true) {

        // get start face
        cFace *start_face = nullptr;
        for (auto &cface: cFaces)
            if (!cface.WasVisited()) {
                start_face = &cface;
                break;
            }

        if (start_face == nullptr)
            break;

        std::set<cFace *> spaceFaces;
        spaceFaces.insert(start_face);

        std::stack<cFace *> face_stack;
        face_stack.push(start_face);

        while (!face_stack.empty()) {

            cFace *pF = face_stack.top();
            face_stack.pop();

            if (!pF->WasVisited()) // if not visited, add adjacent faces to stack
                for (auto &it: pF->adjacentfaces)
                    for (auto &adjacent_cface: it.second)
                        face_stack.push(adjacent_cface);

            pF->SetWasVisited(true);
            spaceFaces.insert(pF);
        }

        spaces.emplace_back(space_id, spaceFaces, false);
        for (auto &spaceFace: spaceFaces) spaceFace->SetSpace(&(spaces.back()));
        space_id++;

    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Search components", std::to_string(spaces.size()));

}

void Kernel::identify_facade_space(std::list<Space> &spaces) {

    auto start = std::chrono::high_resolution_clock::now();

    double v = 0;
    Space *s = nullptr;

    for (auto &space: spaces)
        if (space.Volume() > v) {
            v = space.Volume();
            s = &space;
        }

    if (s == nullptr)
        std::cerr << "[Warning] No facade space found!" << std::endl;
    else {
        s->is_facade = true;
        s->shell.Complement();
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Identify facade space", std::to_string(spaces.size()));
}

void Kernel::remove_trash_from_spaces(std::list<cFace> &cFaces, std::list<Space> &spaces) {

    auto start = std::chrono::high_resolution_clock::now();

    for (auto &space: spaces) {

        std::list<cFace *> del;

        for (auto &cface: space.FirstLvl)
            if (cface->IsTrash())
                del.push_back(cface);

        for (auto &d: del)
            space.FirstLvl.erase(d);

        del.clear();

        for (auto &cface: space.SecondLvl)
            if (cface->IsTrash())
                del.push_back(cface);

        for (auto &d: del)
            space.SecondLvl.erase(d);
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Remove trash faces from spaces", std::to_string(spaces.size()));

}

bool Kernel::write_spaces_to_stl(std::list<cFace> &cFaces, std::list<Space> &spaces, const std::string &input_ifc_filepath, const std::string &output_ifc_file, const std::string &path, const std::string &mode, bool ascii, bool ignore_facade_space,
                                 const std::map<std::string, TopoDS_Shape> &additional_shapes, double fuzzy_tol) {

    auto start = std::chrono::high_resolution_clock::now();

    // start model enrichment based on schema
    if (ifcSchema == IFC2X3)
        write_spaces_to_stl_worker<Ifc2x3>(cFaces, spaces, input_ifc_filepath, output_ifc_file, path, mode, ascii, ignore_facade_space, additional_shapes, fuzzy_tol);
    else if (ifcSchema == IFC4)
        write_spaces_to_stl_worker<Ifc4>(cFaces, spaces, input_ifc_filepath, output_ifc_file, path, mode, ascii, ignore_facade_space, additional_shapes, fuzzy_tol);
    else if (ifcSchema == IFC4X1)
        write_spaces_to_stl_worker<Ifc4x1>(cFaces, spaces, input_ifc_filepath, output_ifc_file, path, mode, ascii, ignore_facade_space, additional_shapes, fuzzy_tol);
    else if (ifcSchema == IFC4X2)
        write_spaces_to_stl_worker<Ifc4x2>(cFaces, spaces, input_ifc_filepath, output_ifc_file, path, mode, ascii, ignore_facade_space, additional_shapes, fuzzy_tol);
    else if (ifcSchema == IFC4X3_RC1)
        write_spaces_to_stl_worker<Ifc4x3_rc1>(cFaces, spaces, input_ifc_filepath, output_ifc_file, path, mode, ascii, ignore_facade_space, additional_shapes, fuzzy_tol);
    else {
        std::cerr << "[Error] Schema not implemented yet! " << std::endl;
        return false;
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Write spaces to stl", std::to_string(cFaces.size()));

    return true;
}

template<typename Schema>
void Kernel::write_spaces_to_stl_worker(std::list<cFace> &cFaces, std::list<Space> &spaces, const std::string &input_ifc_filepath, const std::string &output_ifc_file, std::string path, const std::string &mode, bool ascii, bool ignore_facade_space,
                                        const std::map<std::string, TopoDS_Shape> &additional_shapes, double fuzzy_tol) {

    auto start = std::chrono::high_resolution_clock::now();

    SB_StlAPI_Writer E;

    // ASCII or binary output
    Standard_Boolean &is_ascii = E.ASCIIMode();
    if (ascii)
        is_ascii = true;
    else {
        is_ascii = false;
        std::cerr << "binary not implemented" << std::endl;
        return;
    }

    auto in_time_t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::stringstream ss;
    ss << std::put_time(std::localtime(&in_time_t), "%Y-%m-%d-%X");

    if (mode == "productWiseInOne") {

        remove(output_ifc_file.c_str()); // all products are appended to one single file, so check if file exists beforehand

        std::ifstream file(output_ifc_file.c_str());
        if (file.good()) {
            std::cerr << "[Error] Could not overwrite file." << std::endl;
            return;
        }

        bool fuse_instead_of_sew = true;

        for (auto &space: spaces) {

            if (space.is_facade && ignore_facade_space) continue;

            // sewing should be default approach, but sometimes occ has a bug: there are hanging faces despite shell check says no free edges.
            // after fuse they are gone

            std::unordered_map<std::string, TopoDS_ListOfShape> M; // collect solids by summarizing faces from same product

            if (fuse_instead_of_sew) {

                // fuse
                BOPAlgo_Builder fuse;
                for (auto &cface: space.FirstLevel())
                    fuse.AddArgument(cface->face);
                fuse.SetFuzzyValue(fuzzy_tol);
                fuse.Perform();
                TopoDS_Shape shell = fuse.Shape();

                // triangulate
                BRepMesh_IncrementalMesh(shell, 5.0, false, 5.0, true);

                for (auto &cface: space.FirstLvl) {

                    std::string n = cface->IsOpening() ?
                                    cface->Ancestor()->Opening()->declaration().name() + "_" + cface->Ancestor()->Opening()->as<typename Schema::IfcProduct>()->GlobalId() :
                                    cface->IfcClass() + "_" + cface->IfcGuid();

                    auto L = fuse.Modified(cface->face);

                    if (L.IsEmpty())
                        M[n].Append(cface->face);
                    else
                        for (const TopoDS_Shape &shape: L)
                            M[n].Append(TopoDS::Face(shape));
                }
            } else {

                // sew
                BRepBuilderAPI_Sewing sew;
                for (auto &cface: space.FirstLvl)
                    sew.Add(cface->face);
                sew.Perform();
                TopoDS_Shape shell = sew.SewedShape();

                // triangulate
                BRepMesh_IncrementalMesh(shell, 5.0, false, 5.0, true);

                for (auto &cface: space.FirstLvl) {

                    std::string n = cface->IsOpening() ?
                                    cface->Ancestor()->Opening()->declaration().name() + "_" + cface->Ancestor()->Opening()->as<typename Schema::IfcProduct>()->GlobalId() :
                                    cface->IfcClass() + "_" + cface->IfcGuid();

                    M[n].Append(sew.Modified(cface->face));
                }
            }

            for (const auto &L: M) {
                TopoDS_Compound c = compound_from_shape_list(L.second);
                std::string solid_name = std::to_string(space.id) + "_" + L.first;
                if (!E.Write(c, output_ifc_file.c_str(), solid_name.c_str(), false)) std::cerr << "[Error] Compound could not be exported!" << std::endl;
            }
        }

        for (auto &additional_shape: additional_shapes)
            if (!E.Write(additional_shape.second, output_ifc_file.c_str(), additional_shape.first.c_str(), false))
                std::cerr << "[Error] Additional shape could not be exported!" << std::endl;

    } else if (mode == "productWiseInOne2") {

        //std::string filename = path + boost::filesystem::path(input_ifc_filepath).stem().string() + "_" + ss.str() + ".stl";
        //std::cout << " (target path: " << filename << ") ";
        const std::string &filename = output_ifc_file;

        // all products are appended to one single file, so check if file exists beforehand
        remove(filename.c_str());
        //if (remove(filename.c_str()) != 0)
        //perror("[Info] File not existing or could not be deleted.");

        std::ifstream file(filename.c_str());
        if (file.good()) {
            std::cerr << "[Error] Could not overwrite file." << std::endl;
            return;
        }

        for (auto &space: spaces) {
            if (space.is_facade && ignore_facade_space) continue;

            for (auto &cface: space.FirstLvl)
                if (!E.Write(cface->face, filename.c_str(), cface->STLSolidName().c_str(), false))
                    std::cerr << "[Error] Face could not be exported!" << std::endl;
        }

        for (auto &additional_shape: additional_shapes)
            if (!E.Write(additional_shape.second, filename.c_str(), additional_shape.first.c_str(), false))
                std::cerr << "[Error] Additional shape could not be exported!" << std::endl;

    } else if (mode == "productWiseInProduct") {

        // create directory if not existing
        boost::filesystem::path dir(path);
        if (boost::filesystem::create_directory(dir))
            std::cout << "... directory created: " << path << " ... ";

        for (auto &space: spaces) {
            if (ignore_facade_space && space.is_facade) continue;

            for (auto &cface: space.FirstLvl) {

                std::string solidname = cface->STLSolidName();
                std::string filename = path + boost::filesystem::path(input_ifc_filepath).stem().string() + "_" + ss.str() + "_" + solidname + ".stl";

                if (!E.Write(cface->face, filename.c_str(), solidname.c_str(), true))
                    std::cerr << "[Error] Face could not be exported!" << std::endl;
            }
        }

        for (auto &additional_shape: additional_shapes) {

            std::string solidname = "add_" + additional_shape.first;
            std::string filename = path + boost::filesystem::path(input_ifc_filepath).stem().string() + "_" + ss.str() + "_" + solidname + ".stl";

            if (!E.Write(additional_shape.second, filename.c_str(), solidname.c_str(), true))
                std::cerr << "[Error] Additional shape could not be exported!" << std::endl;
        }

    } else
        std::cerr << "not implemented" << std::endl;
}

void Kernel::process_space_shells_for_export(std::list<Space> &spaces, bool skip_facade) {

    auto start = std::chrono::high_resolution_clock::now();

#ifdef PARALLEL_PROCESSING
    if (skip_facade)
        tbb::parallel_for_each(spaces.begin(), spaces.end(), [&](Space &space) { if (!space.is_facade) space.shell = TopoDS::Shell(process_shape_for_export(space.shell)); });
    else
        tbb::parallel_for_each(spaces.begin(), spaces.end(), [&](Space &space) { space.shell = TopoDS::Shell(process_shape_for_export(space.shell)); });
#else
    if (skip_facade)
        for (auto &space: spaces) {
            if (space.is_facade) continue;
            space.shell = TopoDS::Shell(process_shape_for_export(space.shell));
        }
    else
        for (auto &space: spaces)
            space.shell = TopoDS::Shell(process_shape_for_export(space.shell));
#endif

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Process space shells for export", std::to_string(spaces.size()));
}

std::set<std::string> Kernel::find_openings_in_products(std::list<Product> &products) {

    auto start = std::chrono::high_resolution_clock::now();

    std::set<std::string> guids;

    if (ifcSchema == IFC2X3)
        find_openings_in_products_worker<Ifc2x3>(products, guids);
    else if (ifcSchema == IFC4)
        find_openings_in_products_worker<Ifc4>(products, guids);
    else if (ifcSchema == IFC4X1)
        find_openings_in_products_worker<Ifc4x1>(products, guids);
    else if (ifcSchema == IFC4X2)
        find_openings_in_products_worker<Ifc4x2>(products, guids);
    else if (ifcSchema == IFC4X3_RC1)
        find_openings_in_products_worker<Ifc4x3_rc1>(products, guids);
    else
        std::cerr << "[Error] Schema not implemented yet! " << std::endl;

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Find not filled openings in products", std::to_string(products.size()));

    return guids;
}

template<typename Schema>
void Kernel::find_openings_in_products_worker(std::list<Product> &products, std::set<std::string> &guids) {

    /*
     IfcElement (IfcWall)
     IfcElement.HasOpenings -> IfcRelVoidsElements
     IfcRelVoidsElement.RelatedOpeningElement -> IfcOpeningElement
     IfcOpeningElement.HasFillings -> IfcRelFillsElements
     IfcRelFillsElement.RelatedBuildingElement -> IfcElement (IfcWindow)
     */

    for (auto &product: products) {

        if (!product.IsIfcClass("IfcElement")) continue;

        auto IfcElement = product.ifcproduct->as<typename Schema::IfcElement>();

        if (IfcElement->HasOpenings()->size() > 0) {

            auto IfcRelVoidsElements = IfcElement->HasOpenings();

            for (const auto &IfcRelVoidsElement: *IfcRelVoidsElements) {

                // IfcFeatureElementSubtraction -> IfcOpeningElement
                auto IfcFeatureElementSubtraction = IfcRelVoidsElement->RelatedOpeningElement();
                if (!IfcFeatureElementSubtraction->declaration().is("IfcOpeningElement")) continue;
                auto *IfcOpeningElement = IfcFeatureElementSubtraction->template as<typename Schema::IfcOpeningElement>();

                guids.insert(IfcOpeningElement->GlobalId()); // add IfcOpening to the guid list. Products in the shape will get a shape later

                if (IfcOpeningElement->HasFillings()->size() == 0)
                    product.IfcOpeningElements.push_back(IfcOpeningElement); // hole is not filled by another product
                else { // get the entity filling the hole, eg. window
                    auto IfcRelFillsElements = IfcOpeningElement->HasFillings();
                    for (const auto &IfcRelFillsElement: *IfcRelFillsElements)
                        product.VoidFillingElements.emplace_back(IfcOpeningElement, IfcRelFillsElement->RelatedBuildingElement());
                }

            }
        }
    }
}

void Kernel::subtract_openings_from_products(std::list<Product> &products, std::unordered_map<IfcUtil::IfcBaseClass *, TopoDS_Shape> &opening_shapes) {

    auto start = std::chrono::high_resolution_clock::now();

#ifdef PARALLEL_PROCESSING
    tbb::parallel_for_each(products.begin(), products.end(), [&](Product &product) {

        if (product.IfcOpeningElements.empty()) return;

        TopoDS_ListOfShape L;

        for (auto &OpeningElement: product.IfcOpeningElements)
            if (opening_shapes.find(OpeningElement) != opening_shapes.end())
                L.Append(opening_shapes[OpeningElement]);

        if (L.IsEmpty()) return;

        TopoDS_Shape Tool;

        if (L.Size() == 1) Tool = L.First();
        else
            Tool = compound_from_shape_list(L);

        product.shape = BRepAlgoAPI_Cut(product.shape, Tool).Shape();
    });
#else
    for (auto &product: products) {
        if (product.IfcOpeningElements.empty()) continue;

        TopoDS_ListOfShape L;

        for (auto &OpeningElement: product.IfcOpeningElements)
            if (opening_shapes.find(OpeningElement) != opening_shapes.end())
                L.Append(opening_shapes[OpeningElement]);

        if (L.IsEmpty()) continue;

        TopoDS_Shape Tool;

        if (L.Size() == 1) Tool = L.First();
        else {
            BRep_Builder builder;
            TopoDS_Compound C;
            builder.MakeCompound(C);
            for (const auto &S: L)
                builder.Add(C, S);
            Tool = C;
        }

        product.shape = BRepAlgoAPI_Cut(product.shape, Tool).Shape();
//        BRepAlgoAPI_Cut Cutter;
//        TopTools_ListOfShape L1, L2;
//        L1.Append(product.shape);
//        L2.Append(Tool);
//        Cutter.SetArguments(L1);
//        Cutter.SetTools(L2);
//        Cutter.SetFuzzyValue(1.0e-5);
//        Cutter.Build();
//        product.shape = Cutter.Shape();
    }
#endif

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Subtract openings from products", std::to_string(opening_shapes.size()));
}

void Kernel::identify_blocked_faces(std::list<cFace> &cFaces) {

    // After deleting enclosed offset faces (enclosed means in a product) and deleting hanging faces (faces that do not have at least one adjacent consistent oriented SHELL face on all half-edges),
    // there are now remaining offset faces that are closing gaps (e.g. beam in the roof of FZK-Haus Ifc4). It can be correct to keep them and delete shell faces that
    // are behind the newly created boundary (blocked faces), for this:
    //TODO again remove inner faces but dont take into account offset/shell classification or linked product
    // Current quickfix to handle situation is instead of identifying and removing the blocked faces, the blocking offset faces are removed

    auto start = std::chrono::high_resolution_clock::now();

    for (auto &cface: cFaces)
        if (cface.IsOffset()) // && !cface.ifcproduct()->declaration().is("IfcWall")
            cface.SetIsTrash(true);

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Identify blocked faces (by joined offset faces)", std::to_string(cFaces.size()));
}

void Kernel::update_space_neighbourships(std::list<Space> &spaces, double s) {

    rtree_lib::RTree<Space *, double, 3, double> tree;

    // populate rtree and remove old nb relationships
    for (auto &space: spaces) {
        space.nbs.clear();
        double min[3] = {space.aabb_xmin - s, space.aabb_ymin - s, space.aabb_zmin - s};
        double max[3] = {space.aabb_xmax + s, space.aabb_ymax + s, space.aabb_zmax + s};
        tree.Insert(min, max, &space);
    }

    // find intersections (spaces will find themselves, too)
    for (auto &space: spaces) {
        double min[3] = {space.aabb_xmin - s, space.aabb_ymin - s, space.aabb_zmin - s};
        double max[3] = {space.aabb_xmax + s, space.aabb_ymax + s, space.aabb_zmax + s};
        tree.Search(min, max, [&space](Space *found_space) {
            space.nbs.insert(found_space);
            return true;
        });
    }

    // delete self-reference
    for (auto &space: spaces) {
        space.nbs.erase(&space);
        if (space.nbs.empty()) std::cerr << "[Warning] Space has no space neighbour!" << std::endl;
    }
}

void Kernel::remove_invalid_spaces(std::list<Space> &spaces, std::list<cFace> &cFaces) {

    auto start = std::chrono::high_resolution_clock::now();

    std::set<Space *> del;

    // find neighbouring spaces of space using rtree
    update_space_neighbourships(spaces, 0.05);

    for (auto &space: spaces) {

        if (space.is_facade)
            continue;

        // Criteria 1: space in space **************
        gp_Pnt c = space.random_point_on_cface_shell();

        for (auto &nb_space: space.nbs) {

            if (nb_space->is_facade || !nb_space->is_point_in_aabb(c))
                continue;

            if (nb_space->is_point_in_shell(c)) {
                std::cout << "[Warning] Space " << space.id << " is located within space " << nb_space->id << ". Space will be removed\n";
                del.insert(&space);
                break;
            }
        }
        // *****************************************

        // Criteria 2: small volume ****************
        if (space.Volume() < 0.1) {
            std::cout << "[Warning] Space " << space.id << " has small volume " << space.Volume() << ". Space will be removed\n";
            del.insert(&space);
        }
        // *****************************************

        // Criteria 3: surrounds a building element and not an air volume
        if (!space.wasFlipped && !space.is_facade) {
            std::cout << "[Warning] Space " << space.id << " seems to surround material instead of air. Space will be removed\n";
            del.insert(&space);
        }
        // *****************************************

        // Criteria 4: consisting of one product only
        std::set<Product *> l;
        for (auto &cface: space.FirstLvl)
            l.insert(cface->RelProduct());
        if (l.size() == 1) {
            std::cout << "[Warning] Space " << space.id << " has only faces from one product " << (*l.begin())->guid << ", " << (*l.begin())->IfcClass() << ". Space will be removed\n";
            del.insert(&space);
        }
        // *****************************************
    }

    // Setting cfaces belonging to space to trash
    for (auto &cface: cFaces)
        if (del.find(cface.RelSpace()) != del.end()) {
            cface.SetIsTrash(true);
            cface.SetIsInner(true); // also add to inner faces for ray tracing
        }

    // Remove space
    for (auto &space: del)
        spaces.remove(*space);

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Remove invalid spaces", std::to_string(spaces.size()));
}

void Kernel::initialize_space_solid_classifiers(std::list<Space> &spaces) {

    auto start = std::chrono::high_resolution_clock::now();

    for (auto &space: spaces)
        space.create_solid_classifier_from_shell();

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Calculate solid classifiers of spaces", std::to_string(spaces.size()));
}

std::map<std::string, TopoDS_Shape> Kernel::add_additional_shapes(std::unique_ptr<IfcParse::IfcFile> &model, IfcGeom::IteratorSettings settings, const std::set<std::string> &include_entities, const std::set<std::string> &cfd_entities, bool heal_shapes, bool triangulate_shapes) const {

    auto start = std::chrono::high_resolution_clock::now();

    std::map<std::string, TopoDS_Shape> M;
    std::set<std::string> guids;

    std::set<std::string> ignore = {
            "IfcAnnotation",
            "IfcBuilding",
            "IfcBuildingStorey",
            "IfcCivilElement",
            "IfcOpeningElement",
            "IfcSite",
            "IfcSpace",
            "IfcVirtualElement",
            "IfcBuildingElementProxy"
    };

    for (const auto &c: include_entities)
        ignore.insert(c);


    //***********************************************************
    // Get guids of products to be used
    std::set<std::string> classes;
    for (const auto &c: cfd_entities) {

        boost::shared_ptr<IfcEntityList> IfcEntityList = nullptr;
        try { IfcEntityList = model->instances_by_type(c); }
        catch (...) {}

        if (IfcEntityList != nullptr)
            for (auto E: *IfcEntityList) { // for the entities, check if they are child of IfcElement class, only those can have or be a void filling product

                bool to_ignore = false;

                for (const auto &i: ignore)
                    if (E->declaration().is(i)) {
                        to_ignore = true;
                        break;
                    }

                if (to_ignore) continue;

                std::string s = E->data().getArgument(0)->toString();
                s.pop_back();
                s.erase(s.begin());
                guids.insert(s);
                classes.insert(E->declaration().name());
            }
    }

    if (guids.empty()) {
        std::cout << "[Info] No additional cfd entities found.\n";
        return M;
    } else {
        std::cout << "[Info] Add " << guids.size() << " cfd entities.\n";
        for (const auto &c: classes) std::cout << "\t" << c << "\n";
    }
    //***********************************************************


    //***********************************************************
    // Create shapes
    IfcGeom::attribute_filter attribute_filter;
    attribute_filter.include = true;
    attribute_filter.traverse = false;
    attribute_filter.attribute_name = "GlobalId";
    attribute_filter.populate(guids);
    std::vector<IfcGeom::filter_t> filter_funcs;
    filter_funcs.emplace_back(boost::ref(attribute_filter));

    IfcGeom::Iterator<real_t> geom_iterator(settings, model.get(), filter_funcs, num_threads); // construct iterator

    if (!geom_iterator.initialize()) {
        std::cout << "[Info] No geometrical entities found. Initialization failed.\n";
        return M;
    }

    do {
        IfcGeom::Element<real_t> *geom_object = geom_iterator.get();
        if (!check_ifc_context(geom_object)) continue;
        std::string s = "X_" + geom_object->product()->declaration().name() + "_" + geom_object->guid();
        TopoDS_Shape S = geom_object_to_shape(geom_object);
        if (!S.IsNull())
            M[s] = S;
    } while (geom_iterator.next());
    //***********************************************************


    //***********************************************************
    // Healing
    if (heal_shapes) {
#ifdef PARALLEL_PROCESSING
        std::vector<std::pair<std::string, TopoDS_Shape>> V;
        for (auto &m: M) V.emplace_back(m.first, m.second);

#pragma omp parallel for default(none) shared(V, M, std::cerr) schedule(dynamic, 1) num_threads(num_threads)
        for (unsigned int i = 0; i < V.size(); i++)
            V[i].second = ShapeHealing(shape_copy(V[i].second)).heal_shape();

        M.clear();
        for (auto &item: V)
            if (item.second.IsNull()) std::cerr << "[Warning] Shape healing failed." << item.first << std::endl;
            else M[item.first] = item.second;
#else
        std::set <std::string> del;
            for (auto &item: M) {
                item.second = ShapeHealing(item.second).heal_shape();
                if (item.second.IsNull()) {
                    std::cerr << "[Warning] Shape healing failed." << "\n";
                    del.insert(item.first);
                }
            }
            // remove null shapes
            for (auto &d: del)
                M.erase(d);
#endif
    }
    //***********************************************************


    //***********************************************************
    // Triangulation
    if (triangulate_shapes)
        for (auto &item: M)
            BRepMesh_IncrementalMesh(item.second, 5.0, false, 5.0, true);
    //***********************************************************

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Add additional shapes from IFC", std::to_string(M.size()));

    return M;
}

void Kernel::add_attribut_strings(Argument *arg, const std::string &name, Space &space) {

    if (space.old_space_info.find(name) == space.old_space_info.end()) return;

    if (!arg->isNull()) {
        std::string str = arg->toString();
        str.erase(remove(str.begin(), str.end(), '\''), str.end());

        if (space.old_space_info.find(name)->second.empty())
            space.old_space_info[name] = str;
        else
            space.old_space_info[name].append("_" + str);
    }
}

void Kernel::split_spaces_by_IfcSpaces(std::list<Space> &spaces, std::list<cFace> &cFaces, std::list<oFace> &orig_faces, std::list<Product> &products, ifcspaceInfoList &old_spaces, unsigned int &fid) {

    auto start = std::chrono::high_resolution_clock::now();

    if (old_spaces.empty()) return;

    // create rtree and save intersections
    auto comp = IfcSpaces_rtree(old_spaces, spaces);
    if (comp.empty()) return;

    // evaluate and filter according to overlap volumes
    link_spaces_and_IfcSpaces(spaces, comp, false);
    if (comp.empty()) return;

    // split spaces by old spaces
    if (ifcSchema == IFC2X3)
        old_spaces_split_and_create_new_spaces<Ifc2x3>(cFaces, orig_faces, products, spaces, comp, fid);
    else if (ifcSchema == IFC4)
        old_spaces_split_and_create_new_spaces<Ifc4>(cFaces, orig_faces, products, spaces, comp, fid);
    else if (ifcSchema == IFC4X1)
        old_spaces_split_and_create_new_spaces<Ifc4x1>(cFaces, orig_faces, products, spaces, comp, fid);
    else if (ifcSchema == IFC4X2)
        old_spaces_split_and_create_new_spaces<Ifc4x2>(cFaces, orig_faces, products, spaces, comp, fid);
    else if (ifcSchema == IFC4X3_RC1)
        old_spaces_split_and_create_new_spaces<Ifc4x3_rc1>(cFaces, orig_faces, products, spaces, comp, fid);
    else {
        std::cerr << "[Error] Schema not implemented yet! " << std::endl;
        return;
    }

/*    // print matching info
    link_spaces_and_IfcSpaces(spaces, comp, true);
    for (auto &space: spaces) {
        if (space.is_facade) continue;
        const auto &matches = comp[&space];
        std::cout << "[Info] Space " << space.id << " with volume of " << space.Volume() << " has " << matches.size() << " matches with old spaces:\n";

        for (auto &found_space: matches) {
            std::cout << "\tGUID: " << found_space->guid << " with volume of " << found_space->volume;
            if (fabs(found_space->volume - space.Volume()) > 1.0e-4) std::cout << "\tvolume difference\n";
            else std::cout << ".\n";
        }
    }*/

/*    // evaluate and filter according to overlap volumes
    link_spaces_and_IfcSpaces(spaces, comp, true);

    // save guids of old spaces in space
    for (auto &space: spaces)
        if (comp.find(&space) != comp.end())
            for (const auto &subspace: comp[&space])
                space.old_space_guids.insert(std::get<1>(*subspace));

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
    }*/

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Split spaces by IfcSpaces from original file", std::to_string(old_spaces.size()));
}

std::unordered_map<Space *, std::list<ifcspaceInfoList::iterator>> Kernel::IfcSpaces_rtree(ifcspaceInfoList &old_spaces, std::list<Space> &spaces) {

    rtree_lib::RTree<ifcspaceInfoList::iterator, double, 3, double> tree;

    // fill rtree
    for (auto it = old_spaces.begin(); it != old_spaces.end(); ++it) {
        Bnd_Box bnd = it->bbox;
        double min[3] = {bnd.CornerMin().X(), bnd.CornerMin().Y(), bnd.CornerMin().Z()};
        double max[3] = {bnd.CornerMax().X(), bnd.CornerMax().Y(), bnd.CornerMax().Z()};
        tree.Insert(min, max, it);
    }

    // find intersections
    std::unordered_map<Space *, std::list<ifcspaceInfoList::iterator>> comp;

    for (auto &space: spaces) {
        if (space.is_facade) continue;
        double min[3] = {space.aabb_xmin, space.aabb_ymin, space.aabb_zmin};
        double max[3] = {space.aabb_xmax, space.aabb_ymax, space.aabb_zmax};
        tree.Search(min, max, [&space, &comp](ifcspaceInfoList::iterator found_space) {
            comp[&space].push_back(found_space);
            return true;
        });
    }

    return comp;
}

std::list<ifcspaceInfoList::iterator>
Kernel::link_spaces_and_IfcSpaces_worker(const Space &space, const std::list<ifcspaceInfoList::iterator> &comp, const TopoDS_Solid &new_solid) {

    std::list<ifcspaceInfoList::iterator> del;

    double V_space = volume(new_solid);
    double Vnew = (space.aabb_xmax - space.aabb_xmin) * (space.aabb_ymax - space.aabb_ymin) * (space.aabb_zmax - space.aabb_zmin);

    for (auto &found_space: comp) {

        auto Pmin = found_space->bbox.CornerMin();
        auto Pmax = found_space->bbox.CornerMax();

        // aabb overlaps
        double dV = overlapping_volume_aabbs(space.aabb_xmin, space.aabb_xmax, space.aabb_ymin, space.aabb_ymax, space.aabb_zmin, space.aabb_zmax, Pmin.X(), Pmax.X(), Pmin.Y(), Pmax.Y(), Pmin.Z(), Pmax.Z());
        double Vold = (Pmax.X() - Pmin.X()) * (Pmax.Y() - Pmin.Y()) * (Pmax.Z() - Pmin.Z());

        // filter
        if (dV / Vnew < 0.4 && dV / Vold < 0.4) {
            del.push_back(found_space);
            continue;
        }

        // shape overlaps
        auto common = BRepAlgoAPI_Common(new_solid, found_space->shape).Shape();
        double dV_bool = volume(common);
        double v1 = dV_bool / V_space;
        double v2 = dV_bool / found_space->volume;

        // filter
        if (v1 < 0.93 && v2 < 0.93) {
            del.push_back(found_space);
            continue;
        }
    }

    return del;

}

void Kernel::link_spaces_and_IfcSpaces(std::list<Space> &spaces, std::unordered_map<Space *, std::list<ifcspaceInfoList::iterator>> &comp, bool reconsider) {

    std::map<Space *, TopoDS_Solid> solids;
    for (auto &space: spaces) {
        if (space.is_facade) continue;
        TopoDS_Solid new_solid = BRepBuilderAPI_MakeSolid(space.shell).Solid();
        if (volume(new_solid) < 0)
            new_solid.Complement();
        solids[&space] = new_solid;
    }
//
//#ifdef PARALLEL_PROCESSING
//    std::vector<std::pair<Space *, std::list<std::list<std::tuple<IfcUtil::IfcBaseEntity *, std::string, TopoDS_Solid, Bnd_Box>>::iterator>>> V;
//    V.reserve(spaces.size());
//    for (auto &space: spaces) {
//        std::list<std::list<std::tuple<IfcUtil::IfcBaseEntity *, std::string, TopoDS_Solid, Bnd_Box>>::iterator> t;
//        V.emplace_back(&space, t);
//    }
//
//#pragma omp parallel for default(none) shared(V, reconsider, comp, solids) schedule(static,10) num_threads(num_threads)
//    for (unsigned int i = 0; i < V.size(); i++) {
//        auto space = V[i].first;
//        if (reconsider && comp[space].size() < 2) continue;
//        // find wrong spaces and delete them
//        V[i].second = link_spaces_and_IfcSpaces_worker(*space, comp[space], solids[space]);
//    }
//
//    for (auto &i: V)
//        for (const auto &d: i.second)
//            comp[i.first].remove(d);
//#else
    for (auto &space: spaces) {
        if (space.is_facade) continue;
        if (reconsider && comp[&space].size() < 2) continue;

        // find wrong spaces and delete them
        auto del = link_spaces_and_IfcSpaces_worker(space, comp[&space], solids[&space]);
        for (const auto &d: del)
            comp[&space].remove(d);
    }
//#endif

}

template<typename Schema>
void Kernel::old_spaces_split_and_create_new_spaces(std::list<cFace> &cFaces, std::list<oFace> &orig_faces, std::list<Product> &products, std::list<Space> &spaces, std::unordered_map<Space *, std::list<ifcspaceInfoList::iterator>> &comp, unsigned int &fid) {

    // spaces to be processed
    std::list<Space *> stack;
    for (auto &space: spaces) if (comp[&space].size() > 1) stack.push_back(&space);

    // id for space creation
    unsigned int space_id = 0;
    for (const auto &space: spaces) if (space.id > space_id) space_id = space.id;
    space_id++;

    std::set<cFace *> del_cface; // cfaces to delete. adding of faces and spaces is done within the loop
    std::set<Space *> del_space; // spaces to delete

    for (auto &space: stack)
        old_spaces_split_space<Schema>(space, cFaces, orig_faces, spaces, products, space_id, fid, comp, del_cface, del_space);

    // delete replaced cfaces in cFaces
    for (auto &d: del_cface)
        cFaces.remove(*d);

    // delete old spaces
    for (auto &d: del_space)
        spaces.remove(*d);

    std::cout << "[Info] Number of removed spaces: " << del_cface.size() << "\n";
}

template<typename Schema>
void Kernel::old_spaces_split_space(Space *space, std::list<cFace> &cFaces, std::list<oFace> &orig_faces, std::list<Space> &spaces, std::list<Product> &products, unsigned int &space_id, unsigned int &fid, std::unordered_map<Space *, std::list<ifcspaceInfoList::iterator>> &comp,
                                    std::set<cFace *> &del_cface, std::set<Space *> &del_space) {

    // setup cut faces
    TopoDS_ListOfShape cutFaces;
    for (auto &old_space: comp[space])
        for (auto &f: Topo(old_space->shape).faces())
            cutFaces.Append(f);

    // List of new subspaces created from new space
    TopoDS_ListOfShape subspaces;
    TopoDS_Solid new_solid = BRepBuilderAPI_MakeSolid(space->shell).Solid();
    subspaces.Append(new_solid);

    // link TopoDS_Faces of new space solid to cfaces
    std::unordered_map<unsigned int, cFace *> id2cFaces;
    for (auto &cface: space->FirstLvl)
        id2cFaces[cface->FaceID()] = cface;

    // split space
    auto H = split_spaces_by_faces(cutFaces, subspaces, 4, 0.4);

    // skip if no valid solids can be created
    if (subspaces.Size() == 1) return;

    // mark to delete
    del_space.insert(space);

    // find virtual face pairs
    auto virtual_pairs = virtual_face_pairs(subspaces, H); // within virtual pair function spaces' faces can be split. Subspaces are updated. Only faces without history are changed.

    // create virtual elements and link them to virtual pair faces
    std::unordered_map<unsigned int, Product *> id2VirtualElement;
    for (const auto &pair: virtual_pairs) {
        auto IfcVirtualElement = ice::IfcVirtualElement<Schema>();
        products.emplace_back(IfcVirtualElement, IfcVirtualElement->GlobalId(), *pair.begin());

        for (const auto &p: pair) // save in map
            id2VirtualElement[hash(p)] = &products.back();
    }

    for (const auto &subspace: subspaces) { // each of the subspaces becomes a new Space

        // list of new cfaces forming new subspace
        std::set<cFace *> spaceFaces;

        // Iterate through new space faces
        for (auto &f: Topo(subspace).faces()) {

            TopoDS_Face F = TopoDS::Face(f); // Use this face as key for the maps
            TopoDS_Face C = TopoDS::Face(F.Complemented()); // Complement face normal because cFaces should point into space. Solid space faces are pointing outside. Use this face for new cfaces and origfaces.

            if (H.IsBound(F)) { // face is from original space solid ...

                const TopoDS_Face &origin_face = H.Find(F);
                cFace *origin_cface = id2cFaces[hash(origin_face)];

                if (hash(origin_face) == hash(F))  // face was not altered at all
                    spaceFaces.insert(origin_cface);

                else { // face was altered, create new cface and add old cface to delete

                    // new cface
                    cFaces.emplace_back(C, origin_cface->Ancestor(), fid);
                    fid++;
                    del_cface.insert(origin_cface);
                    spaceFaces.insert(&cFaces.back());
                }
            } else { // ... or face is subface of the cutting faces. new cface is created linked to virtualelement.

                // new orig face
                auto P = id2VirtualElement[hash(F)];
                orig_faces.emplace_back(C, P, 0);
                orig_faces.back().SetIsOffset(false);

                // new cface
                P->orig_faces.push_back(&orig_faces.back());
                cFaces.emplace_back(C, &orig_faces.back(), fid);
                cFaces.back().SetNormalStatus(FACE_NORMAL_KNOWN); // because splitting is operation on solids, orientation of face is correct
                fid++;
                spaceFaces.insert(&cFaces.back());
            }
        }

        // add new space and link cfaces to space
        spaces.emplace_back(space_id, spaceFaces, true); // ATTENTION: Somehow during sewing, the edges of faces that were created during space spliting and not by the general fuse are altered. So the sew.modified faces replace them
        std::cout << "[Info] Add space " << spaces.back().id << " derived from space " << space->id << "\n";
        for (auto &spaceFace: spaceFaces)
            spaceFace->SetSpace(&(spaces.back()));
        space_id++;

        // save old_spaces
        comp[&spaces.back()] = comp[space];

        // update adjacency information
        for (auto &spaceFace: spaceFaces) {
            spaceFace->halfedges.clear();
            spaceFace->adjacentfaces.clear();
            spaceFace->UpdateHalfEdges();
        }
        std::unordered_map<unsigned int, cFace *> id2cFace;
        for (auto &spaceFace: spaceFaces) id2cFace[spaceFace->FaceID()] = spaceFace;
        TopTools_IndexedDataMapOfShapeListOfShape edgeFaceMap;
        TopExp::MapShapesAndAncestors(spaces.back().shell, TopAbs_EDGE, TopAbs_FACE, edgeFaceMap);
        for (auto &spaceFace: spaceFaces) spaceFace->UpdateFaceAdjacencies(edgeFaceMap, id2cFace);
    }
}

bool Kernel::accept_solid_split(const TopoDS_ListOfShape &solids, double V_crit, double h_crit) {

    for (const auto &solid: solids) // accept solution only for solids big enough
        if (volume(solid) < V_crit || aabb(solid, 0).IsZThin(h_crit))
            return false;

    for (const auto &solid: solids) { // shells must be closed and ok

        TopoDS_ListOfShape shells = Topo(solid).shells();

        for (const auto &shell: shells) {

            BRepCheck_Shell c = BRepCheck_Shell(TopoDS::Shell(shell));
            BRepCheck_Status Closed = c.Closed();
            BRepCheck_Status Orientation = c.Orientation();

            if (Topo(shell).faces().Size() < 6)
                return false;

            if (Closed != BRepCheck_NoError)
                return false;

            if (Orientation != BRepCheck_NoError)
                return false;

            if (Kernel::volume(shell) < 0)
                return false;
        }
    }
    return true;
}

NCollection_DataMap<TopoDS_Face, TopoDS_Face> Kernel::split_spaces_by_faces(const TopoDS_ListOfShape &cutFaces, TopoDS_ListOfShape &subspaces, double V_crit, double h_crit) {

    // V_crit: Critical volume. No new subspaces are created with volumes lower than this (prevent small spaces)
    // h_crit: Critical height. No new subspaces are created with an delta z of the aabb lower than this (prevent spaces within slab opening)

    // inheritence of original faces to sub and subsub faces created by split
    NCollection_DataMap<TopoDS_Face, TopoDS_Face> H;

    // all original faces of the new space are saved
    for (auto &f: Topo(subspaces.First()).faces())
        H.Bind(TopoDS::Face(f), TopoDS::Face(f));

    for (const auto &cutFace: cutFaces) {

        TopoDS_ListOfShape add; // all new created solids
        TopoDS_ListOfShape del; // all solids that are replaced

        for (const auto &subspace: subspaces) {

            // Split solids
            BOPAlgo_MakerVolume mv;
            TopTools_ListOfShape args;
            args.Append(subspace);
            args.Append(cutFace);
            mv.SetArguments(args);
            mv.SetRunParallel(true);
            mv.SetParallelMode(true);
            mv.Perform();
            if (mv.HasWarnings()) {
                std::cerr << "Warnings:" << std::endl;
                mv.DumpWarnings(std::cerr);
            }
            if (mv.HasErrors()) {
                std::cerr << "Errors:" << std::endl;
                mv.DumpErrors(std::cerr);
            }

            // Extract new solids
            TopoDS_ListOfShape solids = Topo(mv.Shape()).solids();

            if (solids.Size() > 1 && !mv.HasWarnings() && !mv.HasErrors()) {

                if (accept_solid_split(solids, V_crit, h_crit)) {

                    // mark for add and delete
                    del.Append(subspace);
                    for (const auto &solid: solids)
                        add.Append(solid);

                    // identify changes done to subshape
                    for (auto &F: Topo(subspace).faces()) {

                        if (!H.IsBound(TopoDS::Face(F))) continue; // skip faces that are not part of the subspace but that were introduced by the cutting face

                        // save link between modified and new faces
                        for (const TopoDS_Shape &itr: mv.Modified(F))
                            H.Bind(TopoDS::Face(itr), H.Find(TopoDS::Face(F)));
                    }
                }
            }
        }

        for (const auto &d: del)
            subspaces.Remove(d);

        for (const auto &a: add)
            subspaces.Append(a);
    }

    return H;
}

std::list<std::list<TopoDS_Face>> Kernel::virtual_face_pairs(TopoDS_ListOfShape &subspaces, const NCollection_DataMap<TopoDS_Face, TopoDS_Face> &H) {

    // collect virtual faces in map
    std::unordered_map<unsigned int, std::list<TopoDS_Face>> raw_pairs;
    NCollection_DataMap<TopoDS_Face, TopoDS_Shape> shapes_faces;
    for (auto &subspace: subspaces)
        for (auto &f: Topo(subspace).faces()) {
            shapes_faces.Bind(TopoDS::Face(f), subspace);
            if (!H.IsBound(TopoDS::Face(f)))
                raw_pairs[hash(f)].push_back(TopoDS::Face(f));
        }


    // find face pairs by HashCode
    std::list<std::list<TopoDS_Face>> virtual_pairs;
    std::vector<TopoDS_Face> S;
    for (auto &p: raw_pairs)
        if (p.second.size() > 1)
            virtual_pairs.push_back(p.second);
        else
            S.push_back(*p.second.begin());

    if (S.empty()) return virtual_pairs;
    //**********************************

    // find face pairs by face center and normal
    std::map<std::vector<int>, std::list<TopoDS_Face>> temp;
    int s = 100;

    for (const auto &p: S) {

        gp_Pnt c = face_center(p);
        gp_Dir n = face_normal(p);

        std::vector<int> t = {
                double_to_scaled_int(c.X(), s),
                double_to_scaled_int(c.Y(), s),
                double_to_scaled_int(c.Z(), s),
                abs(double_to_scaled_int(n.X(), s)),
                abs(double_to_scaled_int(n.Y(), s)),
                abs(double_to_scaled_int(n.Z(), s)),
        };
        temp[t].push_back(p);
    }
    S.clear();
    for (auto &p: temp)
        if (p.second.size() > 1)
            virtual_pairs.push_back(p.second);
        else
            S.push_back(*p.second.begin());

    if (S.empty()) return virtual_pairs;
    //**********************************

    // generate face pairs by boolean
    NCollection_DataMap<TopoDS_Face, TopoDS_ListOfShape> replacing_faces;
    for (int i = 0; i < S.size(); i++) {
        for (int j = 0; j < S.size(); j++) {

            // skip redundant comparisons, faces from same space and faces with same normal
            gp_Dir n_i = face_normal(S[i]);
            gp_Dir n_j = face_normal(S[j]);
            if (j >= i || shapes_faces.Find(S[i]).IsSame(shapes_faces.Find(S[j])) || n_i.Angle(n_j) < 0.1) continue;

            auto common = BRepAlgoAPI_Common(S[i], S[j]).Shape();

            // skip if common does not yield to a face
            auto L = Topo(common).faces();
            if (L.Size() > 1)
                std::cerr << "[Warning] More than one face via common." << std::endl;
            if (common.IsNull() || L.Size() != 1) continue;

            TopoDS_Face C = TopoDS::Face(L.First());
            double A_i = area(S[i]);
            double A_j = area(S[j]);
            double A_c = area(common);
            double tol = 0.001;
            gp_Dir n_c = face_normal(C);

            if (fabs(A_i - A_c) > tol) {
                TopoDS_ListOfShape l;
                if (replacing_faces.IsBound(S[i])) l = replacing_faces.Find(S[i]);
                if (n_i.Angle(n_c) < 0.1) l.Append(C);
                else l.Append(C.Complemented());
                replacing_faces.Bind(S[i], l);
            }

            if (fabs(A_j - A_c) > tol) {
                TopoDS_ListOfShape l;
                if (replacing_faces.IsBound(S[j])) l = replacing_faces.Find(S[j]);
                if (n_j.Angle(n_c) < 0.1) l.Append(C);
                else l.Append(C.Complemented());
                replacing_faces.Bind(S[j], l);
            }

        }
    }

    // reshape spaces
    TopoDS_ListOfShape subspaces_new;

    for (auto &subspace: subspaces) {

        TopoDS_ListOfShape L;

        for (auto &f: Topo(subspace).faces())

            if (!replacing_faces.IsBound(TopoDS::Face(f)))
                L.Append(f);
            else {
                S.erase(std::remove(S.begin(), S.end(), TopoDS::Face(f)), S.end());
                for (auto &new_face: replacing_faces.Find(TopoDS::Face(f))) {
                    L.Append(new_face);
                    S.push_back(TopoDS::Face(new_face));
                }
            }

        subspaces_new.Append(compound_from_shape_list(L));

    }

    // find face pairs by face center and normal
    temp.clear();
    for (const auto &p: S) {

        gp_Pnt c = face_center(p);
        gp_Dir n = face_normal(p);

        std::vector<int> t = {
                double_to_scaled_int(c.X(), s),
                double_to_scaled_int(c.Y(), s),
                double_to_scaled_int(c.Z(), s),
                abs(double_to_scaled_int(n.X(), s)),
                abs(double_to_scaled_int(n.Y(), s)),
                abs(double_to_scaled_int(n.Z(), s)),
        };
        temp[t].push_back(p);
    }
    for (auto &p: temp)
        virtual_pairs.push_back(p.second);

    subspaces = subspaces_new;

    return virtual_pairs;

}

std::string Kernel::print_time(double t, const std::string &s1, const std::string &s2) {

    int c = 68; // tabstop
    int places = 4; // number of decimal places

    std::string r;

    if (s2.empty())
        r = s1;
    else
        r = s1 + " (" + s2 + ")";

    std::string a;
    int n = c - r.size();
    if (n > 0)
        a = std::string(n, ' ');

    double t_rounded = round_double_to_n_decimal_places(t, places);
    std::string time = std::to_string(t_rounded);
    std::vector<std::string> token = split_string_by_character(time, '.');

    if (token.size() == 1) // add zeros
        time = token[0] + "." + std::string(places, '0');
    else if (token.size() == 2)
        if (token[1].size() < places) // add zeros
            time = token[0] + "." + token[1] + std::string(places - token[1].size(), '0');
        else  // cut
            time = token[0] + "." + token[1].substr(0, places);

    r += a + " ... \tElapsed time: " + time + " s\n";
    return r;
}

void Kernel::check_containment_in_fuse(const TopoDS_Shape &fuse, std::list<cFace> &cFaces, bool check_edges) {

    auto start = std::chrono::high_resolution_clock::now();

    TopoDS_ListOfShape S = Topo(fuse).faces();

    // TopoDS_ListOfShape.Contains is orientation sensitive! if face has 0 and in fuse it's 1 it will return false!

#ifdef PARALLEL_PROCESSING
    if (check_edges) {

        TopoDS_ListOfShape E = Topo(fuse).edges();

        tbb::parallel_for_each(cFaces.begin(), cFaces.end(), [&](cFace &cface) {

            // face in fuse
            if (!S.Contains(cface.face))
                if (!S.Contains(cface.face.Complemented()))
                    std::cerr << "[Error] cFace not contained in fuse! " << hash(cface.face) << "\t" << cface.Info() << std::endl;

            // edges in fuse
            for (auto &e: Topo(cface.face).edges())
                if (!E.Contains(e))
                    if (!E.Contains(e.Complemented()))
                        std::cerr << "[Error] Edge not contained in fuse! " << hash(cface.face) << "\t" << cface.Info() << "\t" << hash(e) << std::endl;
        });

    } else {
        tbb::parallel_for_each(cFaces.begin(), cFaces.end(), [&](cFace &cface) {
            if (!S.Contains(cface.face))
                if (!S.Contains(cface.face.Complemented()))
                    std::cerr << "[Error] cFace not contained in fuse! " << hash(cface.face) << "\t" << cface.Info() << std::endl;
        });
    }
#else
    if (check_edges) {

        TopoDS_ListOfShape E = Topo(fuse).edges();

        for (auto &cface: cFaces) {

            // face in fuse
            if (!S.Contains(cface.face))
                if (!S.Contains(cface.face.Complemented()))
                    std::cerr << "[Error] cFace not contained in fuse! " << hash(cface.face) << "\t" << cface.Info() << std::endl;

            // edges in fuse
            for (auto &e: Topo(cface.face).edges())
                if (!E.Contains(e))
                    if (!E.Contains(e.Complemented()))
                        std::cerr << "[Error] Edge not contained in fuse! " << hash(cface.face) << "\t" << cface.Info() << "\t" << hash(e) << std::endl;
        }

    } else {
        for (auto &cface: cFaces) {
            if (!S.Contains(cface.face))
                if (!S.Contains(cface.face.Complemented()))
                    std::cerr << "[Error] cFace not contained in fuse! " << hash(cface.face) << "\t" << cface.Info() << std::endl;
        }
    }
#endif

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Check containment in fuse", std::to_string(S.Size()) + ", " + std::to_string(cFaces.size()));
}

void Kernel::get_zones_from_octree(std::vector<std::pair<std::array<double, 3>, std::set<std::string
>>> &spaces_guids, std::vector<std::pair<std::array<double, 3>, std::set<unsigned int>>> &spaces_triangles,
                                   std::vector<std::set<unsigned int>> &zones, std::list<oFace>
                                   &orig_faces,
                                   double maximum_voxel_size,
                                   int maximum_octree_levels,
                                   double write_vtk, std::vector<oFace *>
                                   &tri2orig,
                                   bool apply_refinement_octrees
) {
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<std::tuple<double, double, double>> vertices;
    std::vector<std::tuple<int, int, int>> faces;
    std::vector<std::string> attrs;
    std::vector<double> X_min; // calculating bbox
    std::vector<double> Y_min;
    std::vector<double> Z_min;
    std::vector<double> X_max;
    std::vector<double> Y_max;
    std::vector<double> Z_max;
    double tol = 1.0e-6;
    for (auto &origface: orig_faces) {
        auto CopyFace = shape_copy(origface.face);
        BRepMesh_IncrementalMesh(CopyFace, 12.0, false, 12.0, true);
        TopExp_Explorer Ex;
        for (Ex.Init(CopyFace, TopAbs_FACE); Ex.More(); Ex.Next()) {
            TopoDS_Face Face = TopoDS::Face(Ex.Current());
            TopLoc_Location L;
            Handle_Poly_Triangulation Poly_Triangulation = BRep_Tool::Triangulation(Face, L);
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
                if (p1.IsEqual(p2, tol) || p2.IsEqual(p3, tol) || p3.IsEqual(p1, tol)) {
                    std::cerr << "[Warning] Overlapping points." << std::endl;
                    continue;
                }
                if (are_points_colinear(p1, p2, p3, 1.0e-5)) {
                    std::cerr << "[Warning] Colinear points." << std::endl;
                    continue;
                }
                vertices.emplace_back(p1.X(), p1.Y(), p1.Z());
                vertices.emplace_back(p2.X(), p2.Y(), p2.Z());
                vertices.emplace_back(p3.X(), p3.Y(), p3.Z());
                faces.emplace_back(vertices.size() - 3, vertices.size() - 2, vertices.size() - 1);
                attrs.push_back(origface.IfcGuid());
                tri2orig.push_back(&origface);
                X_min.push_back(std::min({p1.X(), p2.X(), p3.X()}));
                X_max.push_back(std::max({p1.X(), p2.X(), p3.X()}));
                Y_min.push_back(std::min({p1.Y(), p2.Y(), p3.Y()}));
                Y_max.push_back(std::max({p1.Y(), p2.Y(), p3.Y()}));
                Z_min.push_back(std::min({p1.Z(), p2.Z(), p3.Z()}));
                Z_max.push_back(std::max({p1.Z(), p2.Z(), p3.Z()}));
            }
        }
    }
    double xmin = *std::min_element(X_min.begin(), X_min.end());
    double xmax = *std::max_element(X_max.begin(), X_max.end());
    double ymin = *std::min_element(Y_min.begin(), Y_min.end());
    double ymax = *std::max_element(Y_max.begin(), Y_max.end());
    double zmin = *std::min_element(Z_min.begin(), Z_min.end());
    double zmax = *std::max_element(Z_max.begin(), Z_max.end());
    double dx = xmax - xmin;
    double dy = ymax - ymin;
    double dz = zmax - zmin;
    int octree_depth = std::ceil(log2(std::max({dx, dy, dz}) / maximum_voxel_size));
    if (octree_depth > maximum_octree_levels)
        octree_depth = maximum_octree_levels;  // std::cerr << "\t" << xmin << "\t" << xmax << "\t" << ymin << "\t" << ymax << "\t" << zmin << "\t" << zmax << "\t" << dx << "\t" << dy << "\t" << dz << "\t" << std::max({dx, dy, dz}) / pow(2, 5) << "\t" << octree_depth << std::endl;

    std::set<std::array<double, 3>> fluid_points;
    octree_interface::process_mesh(spaces_guids, spaces_triangles, zones, vertices, faces, attrs, fluid_points, octree_depth, write_vtk
    );  // For each space, run refined octree to minimize number of triangles
    if (apply_refinement_octrees)
        for (const auto &zone: zones) {
            for (const auto &space_index: zone) {
                auto space_tri = &spaces_triangles[space_index];
                auto space_gui = &spaces_guids[space_index];
                refined_octree(space_tri->first, space_tri->second, space_gui->second, vertices, faces, attrs, maximum_voxel_size, maximum_octree_levels);
            }
        }
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Octree", std::to_string(orig_faces.size()) + "/" + std::to_string(spaces_guids.size()));
}

void Kernel::refined_octree(const std::array<double, 3> &fluid_point, std::set<unsigned int> &old_space_triangles, std::set<std::string> &old_spaces_guids, const std::vector<std::tuple<double, double, double>> &global_vertices, const std::vector<std::tuple<int, int, int>> &global_faces,
                            const std::vector<std::string> &global_attrs, double maximum_voxel_size, int maximum_octree_levels) {

    std::vector<std::tuple<double, double, double>> vertices;
    std::vector<std::tuple<int, int, int>> faces;
    std::vector<std::string> attrs;
    std::map<unsigned int, unsigned int> M; //  lookup from new triangle index to old triangle index

    // calculating bbox
    std::vector<double> X_min;
    std::vector<double> Y_min;
    std::vector<double> Z_min;
    std::vector<double> X_max;
    std::vector<double> Y_max;
    std::vector<double> Z_max;

    for (const auto &triangle_index: old_space_triangles) {

        M[faces.size()] = triangle_index;

        auto face = global_faces[triangle_index];
        auto v1 = global_vertices[std::get<0>(face)];
        auto v2 = global_vertices[std::get<1>(face)];
        auto v3 = global_vertices[std::get<2>(face)];
        const auto &attr = global_attrs[triangle_index];

        // fill vertices, faces and attributes
        vertices.emplace_back(v1);
        vertices.emplace_back(v2);
        vertices.emplace_back(v3);
        faces.emplace_back(vertices.size() - 3, vertices.size() - 2, vertices.size() - 1);
        attrs.push_back(attr);

        X_min.push_back(std::min({std::get<0>(v1), std::get<0>(v2), std::get<0>(v3)}));
        X_max.push_back(std::max({std::get<0>(v1), std::get<0>(v2), std::get<0>(v3)}));
        Y_min.push_back(std::min({std::get<1>(v1), std::get<1>(v2), std::get<1>(v3)}));
        Y_max.push_back(std::max({std::get<1>(v1), std::get<1>(v2), std::get<1>(v3)}));
        Z_min.push_back(std::min({std::get<2>(v1), std::get<2>(v2), std::get<2>(v3)}));
        Z_max.push_back(std::max({std::get<2>(v1), std::get<2>(v2), std::get<2>(v3)}));
    }

    double xmin = *std::min_element(X_min.begin(), X_min.end());
    double xmax = *std::max_element(X_max.begin(), X_max.end());
    double ymin = *std::min_element(Y_min.begin(), Y_min.end());
    double ymax = *std::max_element(Y_max.begin(), Y_max.end());
    double zmin = *std::min_element(Z_min.begin(), Z_min.end());
    double zmax = *std::max_element(Z_max.begin(), Z_max.end());
    double dx = xmax - xmin;
    double dy = ymax - ymin;
    double dz = zmax - zmin;

    double finer_maximum_voxel_size = 0.20 * maximum_voxel_size;
    int octree_depth = std::ceil(log2(std::max({dx, dy, dz}) / finer_maximum_voxel_size));
    if (octree_depth > maximum_octree_levels) octree_depth = maximum_octree_levels;

    std::vector<std::pair<std::array<double, 3>, std::set<std::string>>> spaces_guids;
    std::vector<std::pair<std::array<double, 3>, std::set<unsigned int>>> spaces_triangles;
    std::vector<std::set<unsigned int>> zones;

    std::set<std::array<double, 3>> fluid_points = {fluid_point};
    octree_interface::process_mesh(spaces_guids, spaces_triangles, zones, vertices, faces, attrs, fluid_points, octree_depth, false);

    if (spaces_guids.size() != 1 || spaces_triangles.size() != 1) {
        std::cerr << "[Warning] Refined octree has incorrect number of spaces (" << spaces_guids.size() << ", " << spaces_triangles.size() << ")." << std::endl;
        return;
    }

    if (old_space_triangles == spaces_triangles[0].second) { // nothing changed
        std::cerr << "[Info] Refined octree has no new results." << std::endl;
        return;
    }

    std::set<unsigned int> t;
    for (const auto &tri: spaces_triangles[0].second)
        t.insert(M[tri]);

    // update lists
    old_space_triangles = t;
    old_spaces_guids = spaces_guids[0].second;
}

std::vector<std::tuple<std::array<double, 3>, TopoDS_Shape, std::list<cFace>>>
Kernel::perform_fuse_on_zone_builders(std::vector<BOPAlgo_Builder>
                                      &B,
                                      std::vector<std::list<oFace *>> &O, std::vector<NCollection_DataMap<TopoDS_Face, TopoDS_Face>>
                                      &N) const {

    auto start = std::chrono::high_resolution_clock::now();

    std::vector<std::tuple<std::array<double, 3>, TopoDS_Shape, std::list<cFace>>>
            zones(B
                          .

                                  size()

    );

#pragma omp parallel for default(none) shared(B, O, N, zones, std::cerr) schedule(dynamic, 1) num_threads(num_threads)
    for (
            unsigned int i = 0;
            i < B.

                    size();

            i++) {

        auto &builder = B[i];
        std::list<cFace> cFaces;
        unsigned int id_counter = 10000000 * i; // so that there is global unique counting

        builder.

                Perform();

        if (builder.

                HasWarnings()

                ) {
            std::cerr << "Warnings:" <<
                      std::endl;
            builder.
                    DumpWarnings(std::cerr);
        }
        if (builder.

                HasErrors()

                ) {
            std::cerr << "Errors:" <<
                      std::endl;
            builder.
                    DumpErrors(std::cerr);
        }

        if (!builder.

                HasModified()

                )
            std::cerr << "Nothing was modified." <<
                      std::endl;

        for (
            auto &orig_face
                : O[i]) {

            TopoDS_Face F = N[i].Find(orig_face->face);

            bool isDeleted = builder.IsDeleted(F);

            if (isDeleted) {
                auto Generated = builder.Generated(F);
                std::cerr << "[Error] Face was deleted! " <<
                          hash(F)
                          << "\t" << orig_face->IfcGuid() << "\t" << orig_face->IfcClass() << "\t" << isDeleted << "\t" << Generated.

                        Size()

                          <<
                          std::endl;
            }

            auto L = builder.Modified(F);

            if (L.

                    IsEmpty()

                    ) {  // face was not modified
                cFaces.
                        emplace_back(F, orig_face, id_counter
                ); // original TopoDS_Face, pointer on oFace
                id_counter++;
            } else // face was modified into one or more new faces
                for (
                    const TopoDS_Shape &shape
                        : L) {
                    cFaces.
                            emplace_back(TopoDS::Face(shape), orig_face, id_counter
                    ); // new face, pointer on oFace
                    id_counter++;
                }
        }

        std::array<double, 3> mid{0, 0, 0};
        zones[i] =
                std::make_tuple(mid, builder
                        .

                                Shape(), cFaces

                );

    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout <<
              print_time(elapsed
                                 .

                                         count(),

                         "Perform builder fusing",
                         std::to_string(zones
                                                .

                                                        size()

                         ));

    return
            zones;
}

std::vector<std::pair<std::vector<gp_Pnt>, std::set<oFace *>>>
Kernel::get_origfaces_from_octree(std::list<oFace>
                                  &orig_faces,
                                  double maximum_voxel_size,
                                  int maximum_octree_levels,
                                  double write_vtk,
                                  bool apply_refinement_octrees
) {

    auto start = std::chrono::high_resolution_clock::now();

// get space and zone information from octree ****************************************************************
    std::vector<std::pair<std::array<double, 3>, std::set<std::string>>>
            octree_spaces_guids;       // guids of each space, populated by octree
    std::vector<std::pair<std::array<double, 3>, std::set<unsigned int>>>
            octree_spaces_triangles;  // triangles of each space, populated by octree
    std::vector<std::set<unsigned int>> zone_spaces;                                                // points to spaces in a zone
    std::vector<oFace *> tri2orig;                                                               // link triangles to origface

    get_zones_from_octree(octree_spaces_guids, octree_spaces_triangles, zone_spaces, orig_faces, maximum_voxel_size, maximum_octree_levels, write_vtk, tri2orig, apply_refinement_octrees
    ); // list of guids belonging to the zones
//************************************************************************************************************

// get gp_Pnt (random air point) and relevant oFace* for each space ***************************************
// add all faces found via triangle to zone and also all faces belonging to touched products and not being offset
    std::vector<std::pair<gp_Pnt, std::set<oFace * >>>
            spaces_orig;
    std::map<std::string, std::list<oFace *>> M; // lookup table, finding orig_face by guid

    for (
        auto &orig_face
            : orig_faces)
        M[orig_face.IfcGuid()].
                push_back(&orig_face);

    for (
        const auto &space
            : octree_spaces_triangles) {

        spaces_orig.

                emplace_back();

        spaces_orig.

                        back()

                .
                        first = gp_Pnt(space.first[0], space.first[1], space.first[2]);

        std::set<std::string> guids;

// add every face that lies in space and store in guid list
        for (
            const auto &index
                : space.second) {
            auto orig_face = tri2orig[index];
            spaces_orig.

                            back()

                    .second.
                            insert(orig_face);
            guids.
                    insert(orig_face
                                   ->IfcGuid());
        }

// now add orig faces that belong to guid and are no offset faces
        for (
            const auto &guid
                : guids)
            for (
                const auto &orig_face
                    : M[guid])
                if (!orig_face->IsOffset())
                    spaces_orig.

                                    back()

                            .second.
                                    insert(orig_face);
    }
//************************************************************************************************************

// store all information in zone vector **********************************************************************
    std::vector<std::vector<std::pair<gp_Pnt, std::set<oFace * >>>>
            zones;
    zones.
            reserve(zone_spaces
                            .

                                    size()

    );
// zone 0: space 0, space 1, ...

    for (
        const auto &spaces
            : zone_spaces) {
        zones.

                emplace_back();

        for (
            const auto &space
                : spaces)
            zones.

                            back()

                    .
                            push_back(spaces_orig[space]);
    }
//************************************************************************************************************

    for (
            unsigned int i = 0;
            i < zones.

                    size();

            i++) {
        std::cout << "Zone  " << i << ": (" << zones[i].

                size()

                  << ")\n";
        for (
            const auto &space
                : zones[i])
            std::cout << "\t\t\tSpace " << "(" << space.second.

                    size()

                      << "\t | " <<
                      round_double_to_n_decimal_places(space
                                                               .first.

                                                               X(),

                                                       2) << "\t" <<
                      round_double_to_n_decimal_places(space
                                                               .first.

                                                               Y(),

                                                       2) << "\t" <<
                      round_double_to_n_decimal_places(space
                                                               .first.

                                                               Z(),

                                                       2) << ")\n";
    }

#ifdef VISUALIZATION
// visualize_zones(zones);
#endif

// for further processing combine original faces of a zone's spaces ******************************************
    std::vector<std::pair<std::vector<gp_Pnt>, std::set<oFace *>>>
            Z; // zone has complete list of faces and all fluid points
    Z.
            reserve(zones
                            .

                                    size()

    );

    for (
        const auto &zone
            : zones) {
        Z.

                emplace_back();

        for (
            const auto &space
                : zone) {
            Z.

                            back()

                    .first.
                            push_back(space
                                              .first);
            for (
                const auto &face
                    : space.second)
                Z.

                                back()

                        .second.
                                insert(face);
        }
    }
//************************************************************************************************************

    for (
        const auto &z
            : Z)
        std::cout << "Zone: " << z.first.

                size()

                  << "\t" << z.second.

                size()

                  << "\n";

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout <<
              print_time(elapsed
                                 .

                                         count(),

                         "Get zones from octree",
                         std::to_string(zones
                                                .

                                                        size()

                         ) + ", " +
                         std::to_string(spaces_orig
                                                .

                                                        size()

                         ));

    return
            Z;
}

void Kernel::find_spaces_filter_start_face_by_fluid_points(std::list<cFace> &GlobalcFaces, std::list<cFace> &cFaces, std::list<Space> &spaces, unsigned int &space_id, const std::vector<gp_Pnt> &Points_air) {

    auto start = std::chrono::high_resolution_clock::now();

    unsigned int found_spaces = 0;
    unsigned int n_fluid_points = Points_air.size();
    const double tol = 0.00001;
    const double point_movement = 0.00010;

    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(-1, 1);

    // create intersector
    TopoDS_ListOfShape L;
    NCollection_DataMap<TopoDS_Shape, cFace *> N;
    for (auto &cface: cFaces) {
        L.Append(cface.face);
        N.Bind(L.Last(), &cface);
    }
    TopoDS_Shape Comp = compound_from_shape_list(L);
    IntCurvesFace_ShapeIntersector intersector;
    intersector.Load(Comp, tol);

    for (const auto P_air: Points_air) {

        cFace *start_face = nullptr;

        //**********************************************************
        // Find start face. first guess is random face center (helps for facade space)
        gp_Pnt C = face_center(cFaces.begin()->face);
        gp_Lin Line1(P_air, gp_Vec(P_air, C));
        intersector.PerformNearest(Line1, 0, RealLast());

        if (intersector.NbPnt() != 0) {
            const auto &hit_face = intersector.Face(1);
            auto P_isct = intersector.Pnt(1);
            start_face = N.Find(hit_face);
        }
        //**********************************************************

        //**********************************************************
        // Find start face
        if (start_face == nullptr)
            while (true) {
                double a = distribution(generator);
                double b = distribution(generator);
                double c = distribution(generator);
                gp_Lin Line(P_air, gp_Vec(a, b, c));
                intersector.PerformNearest(Line, 0, RealLast());

                if (intersector.NbPnt() == 0) continue;

                const auto &hit_face = intersector.Face(1);
                auto P_isct = intersector.Pnt(1);

                // not necessary in "fine approach" if (!is_point_in_face(hit_face, P_isct, tol)) continue; // only check normal condition, if no edge was hit
                start_face = N.Find(hit_face);
                break;
            }
        //**********************************************************

        //**********************************************************
        // Search space component
        std::set<cFace *> spaceFaces;
        spaceFaces.insert(start_face);

        std::stack<cFace *> face_stack;
        face_stack.push(start_face);

        while (!face_stack.empty()) {

            cFace *pF = face_stack.top();
            face_stack.pop();

            if (!pF->WasVisited()) // if not visited, add adjacent faces to stack
                for (auto &it: pF->adjacentfaces)
                    for (auto &adjacent_cface: it.second)
                        face_stack.push(adjacent_cface);

            pF->SetWasVisited(true);
            spaceFaces.insert(pF);
        }
        //**********************************************************

        //**********************************************************
        // Copy faces from components to graph's global cFaces list. Build space using the global cFaces, instead of the local ones, that are deleted when leaving the scope
        // cFace objects point to adjacent cFaces, so this info has to be updated, too
        std::set<cFace *> G; // spaceFaces stored in global cFace list
        std::map<cFace *, cFace *> M; // store map to later update adjacence info of cfaces
        for (const auto &spaceface: spaceFaces) {
            //GlobalcFaces.push_back(std::move(*spaceface));
            GlobalcFaces.push_back(*spaceface);
            G.insert(&GlobalcFaces.back());
            M[spaceface] = &GlobalcFaces.back();
        }

        // Update adjacence info
        for (const auto &global_spaceface: G) {
            for (auto &adjacent_cfaces: global_spaceface->adjacentfaces) {
                std::list<cFace *> updated;
                for (auto &adjacent_cface: adjacent_cfaces.second)
                    updated.push_back(M[adjacent_cface]);
                adjacent_cfaces.second = updated;
            }
        }

        // create space objects and link cfaces
        spaces.emplace_back(space_id, G, false);
        for (auto &F: G) F->SetSpace(&(spaces.back()));
        space_id++;
        found_spaces++;
        //**********************************************************

    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Search components", std::to_string(n_fluid_points) + "/" + std::to_string(found_spaces));

}

std::vector<std::tuple<std::vector<gp_Pnt>, TopoDS_Shape, std::list<cFace>>>
Kernel::fuse_zones_in_parallel(std::vector<std::pair<std::vector<gp_Pnt>, std::set<oFace *>>
> zones_orig,
                               double fuzzy_tol
) const {

    auto start = std::chrono::high_resolution_clock::now();

    std::vector<std::tuple<std::vector<gp_Pnt>, TopoDS_Shape, std::list<cFace>>>
            zones(zones_orig
                          .

                                  size()

    );

#pragma omp parallel for default(none) shared(fuzzy_tol, zones_orig, zones, std::cerr) schedule(dynamic, 1) num_threads(std::min(int(zones_orig.size()), int(num_threads)))
    for (
            int i = 0;
            i < zones_orig.

                    size();

            i++) {

        TopoDS_Shape fuse;
        std::list<cFace> cFaces;

// populate builder
        BOPAlgo_Builder builder;
        for (
            const auto &f
                : zones_orig[i].second)
            builder.
                    AddArgument(f
                                        ->face);

        builder.SetNonDestructive(true);
        builder.SetCheckInverted(false);
        builder.SetFuzzyValue(fuzzy_tol);
        builder.SetRunParallel(true);
        builder.SetUseOBB(true);
        builder.SetToFillHistory(true);

        builder.

                Perform();

        if (builder.

                HasWarnings()

                ) {
            std::cerr << "Warnings:" <<
                      std::endl;
            builder.
                    DumpWarnings(std::cerr);
        }
        if (builder.

                HasErrors()

                ) {
            std::cerr << "Errors:" <<
                      std::endl;
            builder.
                    DumpErrors(std::cerr);
        }

        if (!builder.

                HasModified()

                )
            std::cerr << "Nothing was modified." <<
                      std::endl;

        unsigned int id_counter = 0;

        for (
            auto &orig_face
                : zones_orig[i].second) {

            bool isDeleted = builder.IsDeleted(orig_face->face); // TODO in case face was deleted, don't add to cfaces (maybe add the generated?). Without ignoring, the cface.face will not be present is fuse shape

            if (isDeleted) {
                auto Generated = builder.Generated(orig_face->face);
                std::cerr << "[Error] Face was deleted! " <<
                          hash(orig_face
                                       ->face) << "\t" << orig_face->IfcGuid() << "\t" << orig_face->IfcClass() << "\t" << isDeleted << "\t" << Generated.

                        Size()

                          <<
                          std::endl;
            }

            auto L = builder.Modified(orig_face->face);

            if (L.

                    IsEmpty()

                    ) {  // face was not modified
                cFaces.
                        emplace_back(orig_face
                                             ->face, orig_face, id_counter); // original TopoDS_Face, pointer on oFace
                id_counter++;
            } else // face was modified into one or more new faces
                for (
                    const TopoDS_Shape &shape
                        : L) {
                    cFaces.
                            emplace_back(TopoDS::Face(shape), orig_face, id_counter
                    ); // new face, pointer on oFace
                    id_counter++;
                }
        }

        std::get<0>(zones[i]) = zones_orig[i].first;
        std::get<1>(zones[i]) = builder.Shape();
        std::get<2>(zones[i]) = cFaces;

    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout <<
              print_time(elapsed
                                 .

                                         count(),

                         "Fuse original faces",
                         std::to_string(zones
                                                .

                                                        size()

                         ));

    return
            zones;
}

std::vector<std::tuple<std::vector<gp_Pnt>, TopoDS_Shape, std::list<cFace>>>
Kernel::get_cfaces_from_octree(std::list<oFace>
                               &orig_faces,
                               double maximum_voxel_size,
                               int maximum_octree_levels,
                               double write_vtk,
                               double fuzzy_tol,
                               bool apply_refinement_octrees
) {

    auto zones_orig = get_origfaces_from_octree(orig_faces, maximum_voxel_size, maximum_octree_levels, write_vtk, apply_refinement_octrees);
    return
            fuse_zones_in_parallel(zones_orig, fuzzy_tol
            );
}

void Kernel::check_manifoldness(const std::list<cFace> &cFaces) {

    // after fuse, non-manifold edges can occur if an offset faces cuts t2o adjacent faces in the common edge -> edge exists three times
    // also if two objects are touching on an edge, there are for faces connected -> edge exists four times
    // 4 offset faces touching, too
    // At the last step, manifold edges are still present, but each faces has only one adjacent face connected on this edge

    /*
     // examples (id, guid, offset, orientation, adjacent ids)
     - after fuse and first adjacency update
    1009077553	3LJODRGPbDdfXHXShFL8tL	false	0	1545389761	2030732324
	2030732324	23Np8uMAvEN9H6Ds49debS	true	0	1545389761	1009077553
	1545389761	3LJODRGPbDdfXHXShFL8tL	false	1	1009077553	2030732324

    808124657	3LJODRGPbDdfXHXShFLBXo	true	0	807065761	1662567540	1328434820
	807065761	3LJODRGPbDdfXHXShFLBXo	true	1	808124657	1662567540	1328434820
	1662567540	3LJODRGPbDdfXHXShFL8tL	true	0	808124657	807065761	1328434820
	1328434820	3LJODRGPbDdfXHXShFL8tL	true	1	808124657	807065761	1662567540

	741006385	23Np8uMAvEN9H6Ds49debS	false	1	740736801	1478880625	1143262737
	740736801	23Np8uMAvEN9H6Ds49debS	false	0	741006385	1478880625	1143262737
	1478880625	0CWDnGbzX4Sf1BesL9Tj2L	false	0	741006385	740736801	1143262737
	1143262737	0CWDnGbzX4Sf1BesL9Tj2L	false	1	741006385	740736801	1478880625

     - after deletion of adjacency by edge orientation
    1009077553	3LJODRGPbDdfXHXShFL8tL	false	0	1545389761
	2030732324	23Np8uMAvEN9H6Ds49debS	true	0	1545389761
	1545389761	3LJODRGPbDdfXHXShFL8tL	false	1	1009077553	2030732324

    808124657	3LJODRGPbDdfXHXShFLBXo	true	0	807065761	1328434820
	807065761	3LJODRGPbDdfXHXShFLBXo	true	1	808124657	1662567540
	1662567540	3LJODRGPbDdfXHXShFL8tL	true	0	807065761	1328434820
	1328434820	3LJODRGPbDdfXHXShFL8tL	true	1	808124657	1662567540

	741006385	23Np8uMAvEN9H6Ds49debS	false	1	740736801	1478880625
	740736801	23Np8uMAvEN9H6Ds49debS	false	0	741006385	1143262737
	1478880625	0CWDnGbzX4Sf1BesL9Tj2L	false	0	741006385	1143262737
	1143262737	0CWDnGbzX4Sf1BesL9Tj2L	false	1	740736801	1478880625

     - after deletion of offset faces
     solved
     solved

	741006385	23Np8uMAvEN9H6Ds49debS	false	1	740736801	1478880625
	740736801	23Np8uMAvEN9H6Ds49debS	false	0	741006385	1143262737
	1478880625	0CWDnGbzX4Sf1BesL9Tj2L	false	0	741006385	1143262737
	1143262737	0CWDnGbzX4Sf1BesL9Tj2L	false	1	740736801	1478880625

    - after deletion of adjacency by face angle
    solved
    solved

	741006385	23Np8uMAvEN9H6Ds49debS	false	1	1478880625
	740736801	23Np8uMAvEN9H6Ds49debS	false	0	1143262737
	1478880625	0CWDnGbzX4Sf1BesL9Tj2L	false	0	741006385
	1143262737	0CWDnGbzX4Sf1BesL9Tj2L	false	1	740736801
     */

    auto start = std::chrono::high_resolution_clock::now();

    std::unordered_map<unsigned int, std::set<cFace *>> all_edges, non_manifold_edges;

    bool non_manifoldness = false;

    for (const auto &cface: cFaces)
        for (const auto &adj: cface.adjacentfaces) {
            if (adj.second.size() != 1) {
                non_manifoldness = true;
                //std::cout << "Face: " << cface.Info() << "\t" << adj.first << "\t" << adj.second.size() << "\n";
            }
            for (const auto &adj_cface: adj.second)
                all_edges[adj.first].insert(adj_cface);
        }

    if (non_manifoldness) {

        std::cout << "[Info] Non-manifoldness appeared\n";

        for (const auto &e: all_edges)
            if (e.second.size() > 2)
                non_manifold_edges.insert(e);

        // Within fuse shape there can be edges, that have the same hash code but are not the same (they are different locations and have different vertices). The manifold check will also display those "pairings":
        // because the edge face map from OCC, only give the correct adjacent faces on the edge, those manifold pairings can be ignored. See Location hashcode (should be different)
        for (const auto &n: non_manifold_edges) {

            unsigned int id = n.first;
            std::cout << "\tEdge: " << id << "\t" << n.second.size() << "\n";

            for (const auto &cface: n.second) {
                TopoDS_Edge Edge = TopoDS::Edge(cface->halfedges[id]);
                auto n1 = cface->FaceNormal();
                std::cout << std::boolalpha << "\t\t" << cface->Info() << "\t" << cface->IsOffset() << "\tO: " << cface->halfedges[id].Orientation() << "\tL: " << cface->halfedges[id].Location().HashCode(INT_MAX);
                std::cout << "\tV: " << TopExp::FirstVertex(Edge).HashCode(INT_MAX) << "\t" << TopExp::LastVertex(Edge).HashCode(INT_MAX) << "\tn:" << Kernel::round_double_two_digits(n1.X()) << "\t" << Kernel::round_double_two_digits(n1.Y())
                          << "\t" << Kernel::round_double_two_digits(n1.Z()) << "\tAdj:";
                for (const auto &adj: cface->adjacentfaces[id])
                    std::cout << "\t" << adj->FaceID();
                std::cout << "\n";
            }
        }
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Non-manifold edges", std::to_string(non_manifold_edges.size()));
}

void Kernel::nullify_orig_faces(std::list<oFace> &orig_faces) {

    for (auto &orig_face: orig_faces)
        orig_face.face = TopoDS_Face();

}

void Kernel::unify_cFaces(std::list<cFace> &cFaces, bool close_filled_wall_holes) {

    // Unify faces with adjacent faces from the original face, (if they are from the same sb type)
    // faces that were already merged with another face are set to trash and then skipped
    // Attention: when wanting to parallelize. Adjacent faces are set to trash, so this info has to be shared. Otherwise unified faces appear more than one

    auto start = std::chrono::high_resolution_clock::now();

    // just unify all faces. holes in wall faces and window faces will be kept. perfect for cfd
    if (!close_filled_wall_holes) {

        for (auto &cface: cFaces) {
            if (cface.IsTrash()) continue;
            std::set<cFace *> faces_to_unify = cface.UnifyableFacesMaintainHoles();
            cface.UnifyFaces(faces_to_unify);
        }

    } else { // unify wall faces with all faces. window faces are "overtaken" by wall faces. reboot window faces afterwards and unify them with other wall faces

        // skip opening faces
        for (auto &cface: cFaces) {
            if (cface.IsTrash() || cface.IsOpening()) continue;
            std::set<cFace *> faces_to_unify = cface.UnifyableFacesOvertakeOpenings();
            cface.UnifyFaces(faces_to_unify);
        }

        // reboot the opening faces that were marked as trash after unifying with wall faces
        for (auto &cface: cFaces)
            if (cface.IsTrash() && cface.IsOpening()) {
                cface.SetIsTrash(false);
                cface.SetParent(cface.UnifyingFace());
                cface.SetUnifyingFace(nullptr);
            }

        // unify opening faces with other opening faces
        for (auto &cface: cFaces) {
            if (cface.IsTrash() || !cface.IsOpening()) continue; // skip wall faces
            std::set<cFace *> faces_to_unify = cface.UnifyableFacesOpeningsOnly();
            cface.UnifyFaces(faces_to_unify);
        }

        // update corresponding faces because trash faces are removed later
        for (auto &cface: cFaces) {
            if (cface.Corresponding() != nullptr) {
                if (cface.Corresponding()->IsTrash())
                    cface.SetCorresponding(cface.Corresponding()->UnifyingFace());
            }
            if (cface.Parent() != nullptr) {
                if (cface.Parent()->IsTrash())
                    cface.SetParent(cface.Parent()->UnifyingFace());
            }
        }
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Unify faces", std::to_string(cFaces.size()));
}

void Kernel::remove_non_manifold_adjacency_by_angle(std::list<cFace> &cFaces) {

    auto start = std::chrono::high_resolution_clock::now();

    for (auto &cface: cFaces)
        cface.RemoveNonManifoldAdjacencyByAngle();

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Remove adjacency of faces on non-manifold edges by angle", std::to_string(cFaces.size()));
}

TopoDS_ListOfShape Kernel::create_extension_faces(const TopoDS_Edge &edge, gp_Dir n1, gp_Dir n2, double length) {

    TopoDS_ListOfShape L;

    double a = n1.Angle(n2);

    if (a < 0.1) {
        // if faces are coplanar only one extended face is needed
        gp_Dir n_m(0.5 * (n1.X() + n2.X()), 0.5 * (n1.Y() + n2.Y()), 0.5 * (n1.Z() + n2.Z()));
        gp_Vec v(n_m);
        v.Scale(length);
        TopoDS_Shape S = BRepPrimAPI_MakePrism(edge, v).Shape();
        if (!S.IsNull()) L.Append(S);
    } else {
        gp_Vec v1(n1);
        v1.Scale(length);
        TopoDS_Shape S1 = BRepPrimAPI_MakePrism(edge, v1).Shape();
        if (!S1.IsNull()) L.Append(S1);

        gp_Vec v2(n2);
        v2.Scale(length);
        TopoDS_Shape S2 = BRepPrimAPI_MakePrism(edge, v2).Shape();
        if (!S2.IsNull()) L.Append(S2);
    }

    return L;
}

bool Kernel::fuse_cFaces(TopoDS_Shape &fuse, const std::list<cFace> &old_cFaces, std::list<cFace> &new_cFaces, double fuzzy_tol, unsigned int &fid) {

    auto start = std::chrono::high_resolution_clock::now();

    // populate builder
    BOPAlgo_Builder builder;
    for (const auto &cface: old_cFaces)
        builder.AddArgument(cface.face);

    builder.SetNonDestructive(false);
    builder.SetCheckInverted(false);
    builder.SetFuzzyValue(fuzzy_tol);
    builder.SetRunParallel(true);
    builder.SetUseOBB(true);
    builder.SetToFillHistory(true);

    builder.Perform();

    if (builder.HasWarnings()) {
        std::cerr << "Warnings:" << std::endl;
        builder.DumpWarnings(std::cerr);
    }
    if (builder.HasErrors()) {
        std::cerr << "Errors:" << std::endl;
        builder.DumpErrors(std::cerr);
    }

    if (!builder.HasModified()) {
        std::cerr << "Nothing was modified." << std::endl;
        return false;
    }

    for (auto &cface: old_cFaces) {

        bool isDeleted = builder.IsDeleted(cface.face);

        if (isDeleted) {
            auto Generated = builder.Generated(cface.face);
            std::cerr << "[Error] Face was deleted! " << cface.Info() << "\t" << isDeleted << "\t" << Generated.Size() << std::endl;
            continue;
        }

        //  if (!builder.Generated(orig_face.face).IsEmpty())
        //      std::cerr << "[Error] Face generated new faces! " << hash(orig_face.face) << "\t" << isDeleted << "\t" << builder.Generated(orig_face.face).Size() << "\t" << builder.Modified(orig_face.face).Size() << std::endl;

        auto L = builder.Modified(cface.face);

        if (L.IsEmpty()) {  // face was not modified
            new_cFaces.emplace_back(cface.face, cface.Ancestor(), fid); // original TopoDS_Face, pointer on oFace
            fid++;
        } else // face was modified into one or more new faces
            for (const TopoDS_Shape &shape: L) {
                new_cFaces.emplace_back(TopoDS::Face(shape), cface.Ancestor(), fid); // new face, pointer on oFace
                fid++;
            }
    }

    fuse = builder.Shape();

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Fuse cFaces", std::to_string(Topo(fuse).faces().Size()) + ", " + std::to_string(new_cFaces.size()));

    return true;

}

void Kernel::check_fuse(const TopoDS_Shape &fuse) {

    auto start = std::chrono::high_resolution_clock::now();

    TopTools_IndexedDataMapOfShapeListOfShape M;
    TopExp::MapShapesAndAncestors(fuse, TopAbs_EDGE, TopAbs_FACE, M);

//    for (auto &m : M)
//        if (!m.First().IsSame(m.Last())) // sometimes in non-closed shape same face references to itself
//            std::cout << "F1 " << hash(m.First()) << "\tF2 " << hash(m.Last()) << "\n";

    TopExp_Explorer Ex;
    for (Ex.Init(fuse, TopAbs_EDGE); Ex.More(); Ex.Next()) {
        const TopTools_ListOfShape &faces = M.FindFromKey(Ex.Current());

        if (faces.Size() == 1) {
            if (Ex.Current().Orientation() == TopAbs_INTERNAL) std::cerr << "[Warning] Only one adjacent face (id) on seam edge " << hash(Ex.Current()) << "\t" << hash(faces.First()) << "\n";
            else std::cerr << "[Warning] Only one adjacent face (id) on non-seam edge " << hash(Ex.Current()) << "\t" << hash(faces.First()) << "\n";
        }

//        std::cout << "HE " << hash(Ex.Current()) << " \t F ";
//        for (auto &face : faces)
//            std::cout << " " << hash(face);
//        std::cout << "\n";
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Check fuse", std::to_string(Topo(fuse).faces().Size()));

}

void Kernel::check_adjacency(const TopoDS_Shape &fuse, std::list<cFace> &cFaces, bool b1) {

    auto start = std::chrono::high_resolution_clock::now();

    bool skip1 = b1; // skip checks, which would throw an error for valid fuses where inner faces were deleted

    std::unordered_map<unsigned int, std::list<cFace *>> id2cFaces;
    for (auto &cface: cFaces)
        if (id2cFaces.find(cface.FaceID()) != id2cFaces.end())
            id2cFaces[cface.FaceID()].emplace_back(&cface);
        else
            id2cFaces[cface.FaceID()] = {&cface};

    TopTools_IndexedDataMapOfShapeListOfShape edgeFaceMap;
    TopExp::MapShapesAndAncestors(fuse, TopAbs_EDGE, TopAbs_FACE, edgeFaceMap);

    // TODO check if faces are adjacent at all. same id, different location?

    for (auto &cface: cFaces)
        for (auto &it: cface.halfedges) {

            auto TopoFaces = edgeFaceMap.FindFromKey(it.second);

            unsigned int edge_id = it.first;
            std::map<std::pair<Product *, unsigned int>, std::list<cFace *>> M; // product, shell id
            std::list<cFace *> coplanar_and_self;
            std::list<cFace *> seam_connection;
            std::list<cFace *> offsets;
            std::list<cFace *> flat;
            unsigned int c1(0); // adjacent faces still existing as cFace
            std::list<std::string> T;
            bool same_shell_product = false;
            bool error = false;

            //********************************************************
            for (auto &face: TopoFaces) {

                unsigned int face_id = hash(face);
                if (id2cFaces.find(face_id) == id2cFaces.end()) continue;

                for (auto &c: id2cFaces[face_id]) { // all cfaces containing topods_face of specific id
                    c1++;
                    if (c->ID() == cface.ID()) coplanar_and_self.push_back(c); // self-reference
                    else if (!c->IsIdInHalfEdges(edge_id)) seam_connection.push_back(c); // connection using seam half-edge
                    else if (c->IsOffset()) offsets.push_back(c);// offset
                    else M[std::make_pair(c->RelProduct(), c->Ancestor()->ShellID())].push_back(c);

                    if (c->FaceID() == cface.FaceID() && (c->IfcGuid() != cface.IfcGuid() || c->Ancestor()->ShellID() != cface.Ancestor()->ShellID())) coplanar_and_self.push_back(c); // coplanar face, can be redundant to M
                    if (c->ID() != cface.ID() && c->FaceID() == cface.FaceID() && c->IfcGuid() == cface.IfcGuid() && c->Ancestor()->ShellID() == cface.Ancestor()->ShellID())
                        flat.push_back(c); // if 3D object was flattened by fusing, there can be 3 adjacent faces from same product and shell. One of the adjacent faces has same faceid (flat face pair)
                }
            }
            //********************************************************

            //********************************************************
            // evaluation

            // Face info
            T.push_back(
                    "\nFace: " + cface.Info() + "\t" + std::to_string(cface.Ancestor()->ShellID()) + "\t" + std::to_string(edge_id) + "\t" +
                    std::to_string(TopoFaces.Size()) + "\t" + std::to_string(Topo(cface.halfedges[edge_id]).vertices().Size()) + "\n");
            T.emplace_back("\t+++ Faces on edge");
            for (const auto &g: TopoFaces)
                if (id2cFaces.find(hash(g)) == id2cFaces.end()) T.push_back("\t" + std::to_string(hash(g)) + " (" + "[Case] No cFace!" + ")");
                else for (auto &c: id2cFaces[hash(g)]) T.push_back("\t" + std::to_string(hash(g)) + " (" + std::to_string(c->ID()) + ")");
            T.emplace_back("\n");

            // Hanging
            if (TopoFaces.Size() < 2) {
                T.emplace_back("\t[Case] hanging1\n");
                error = true;
            }//TODO filter by orientation
            if (c1 < 2) {
                T.emplace_back("\t[Case] hanging2\n");
                error = true;
            }

            // Flat
            if (!flat.empty()) {
                T.push_back("\tFlat (" + std::to_string(flat.size()) + "): ");
                for (const auto &c: flat)
                    T.push_back("\t(" + std::to_string(c->ID()) + ", " + std::to_string(c->FaceID()) + ", " + std::to_string(c->Ancestor()->ShellID()) + ", " + c->IfcGuid() + ")");
                T.emplace_back("\n");
            }

            // Same id
            T.push_back("\tCoplanar and self (" + std::to_string(coplanar_and_self.size()) + "): ");
            for (const auto &c: coplanar_and_self)
                T.push_back("\t(" + std::to_string(c->ID()) + ", " + std::to_string(c->FaceID()) + ", " + std::to_string(c->Ancestor()->ShellID()) + ", " + c->IfcGuid() + ")");
            T.emplace_back("\n");

            // No connection to each other
            if (!seam_connection.empty()) {
                T.emplace_back("\t[Case] OneSided\n");
                error = true;
            }
            T.push_back("\tOne-sided connection (" + std::to_string(seam_connection.size()) + "): "); // other face doesnt have adjacwnce info about this face
            for (const auto &c: seam_connection)
                T.push_back("\t(" + std::to_string(c->ID()) + ", " + std::to_string(c->FaceID()) + ", " + std::to_string(c->Ancestor()->ShellID()) + ", " + c->IfcGuid() + ")");
            T.emplace_back("\n");

            // Offset
            if (!offsets.empty()) T.emplace_back("\t[Case] offset\n");
            T.push_back("\tOffsets (" + std::to_string(offsets.size()) + "): ");
            for (const auto &c: offsets)
                T.push_back("\t(" + std::to_string(c->ID()) + ", " + std::to_string(c->FaceID()) + ", " + std::to_string(c->Ancestor()->ShellID()) + ", " + c->IfcGuid() + ")");
            T.emplace_back("\n");

            // All other connections
            T.push_back("\tShell connections (" + std::to_string(M.size()) + "): " + "\n");
            for (const auto &m: M) {
                T.push_back("\t\tM (" + std::to_string(m.second.size()) + ") -> ");

                // Same product, same shell
                if (m.first.first->guid == cface.IfcGuid() && m.first.second == cface.Ancestor()->ShellID()) {
                    T.emplace_back("Same product, same shell:  -> ");
                    same_shell_product = true;
                    if (m.second.size() != 1) { // duplicate faces
                        std::set<unsigned int> temp;
                        for (const auto &t: m.second) temp.insert(t->FaceID());
                        if (temp.size() != m.second.size()) {
                            T.emplace_back("\t[Case] SPSS D -> ");
                            error = true;
                        } else {
                            T.emplace_back("\t[Case] SPSS U -> ");
                            error = true;
                        }
                    } else T.emplace_back("           -> ");

                    // Same product, different shell
                } else if (m.first.first->guid == cface.IfcGuid() && m.first.second != cface.Ancestor()->ShellID()) {
                    T.emplace_back("Same product, diff. shell: -> ");
                    if (m.second.size() != 2) {
                        T.emplace_back("\t[Case] SPDS    -> ");
                        error = true;
                    } else T.emplace_back("           -> ");

                    // Different product, specific shell
                } else if (m.first.first->guid != cface.IfcGuid()) {
                    T.emplace_back("Different product:         -> ");
                    if (m.second.size() == 1) {

                        bool b = false;
                        for (const auto &l: coplanar_and_self)
                            if (m.second.front()->Ancestor()->ShellID() == l->Ancestor()->ShellID() && m.second.front()->IfcGuid() == l->IfcGuid()) {
                                b = true;
                                break;
                            }

                        if (b) T.emplace_back("\tH C       -> "); // At least one of the adjacent faces of the other product of specific shell is coplanar to the face
                        else {
                            m.second.front()->IsIfcClass("IfcVirtualElement") ? T.emplace_back("\t[Case] 1DPV    -> ") : T.emplace_back("\t[Case] 1DP     -> ");
                            if (!skip1) error = true;
                        }

                    } else if (m.second.size() > 2) {
                        T.emplace_back("\t[Case] MDP     -> ");
                        error = true;
                    } // Too many adjacent faces
                    else T.emplace_back("           -> ");
                }

                T.push_back(std::to_string(m.first.second) + ", " + m.first.first->guid + ": ");

                for (const auto &c: m.second)
                    T.push_back("\t(" + std::to_string(c->ID()) + ", " + std::to_string(c->FaceID()) + ", " + std::to_string(c->Ancestor()->ShellID()) + ", " + c->IfcGuid() + ", " + c->IfcClass() + ")");
                T.emplace_back("\n");
            }
            if (!skip1 && !same_shell_product) {
                T.emplace_back("\t[Case] Missing adjacent product/shell face\n");
                error = true;
            }

            if (error) for (const auto &t: T) std::cout << t;
        }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Check adjacency", std::to_string(cFaces.size()));
}

void Kernel::simplify_products(std::list<Product> &products) const {

    auto start = std::chrono::high_resolution_clock::now();

#ifdef PARALLEL_PROCESSING

    std::vector<Product *> V;
    V.reserve(products.size());
    for (auto &product: products)
        V.push_back(&product);

#pragma omp parallel for default(none) shared(V, std::cout) schedule(static, 5) num_threads(num_threads)
    for (unsigned int i = 0; i < V.size(); i++) {
        Product *p = V[i];
        if (p->IsIfcClass("IfcDoor") || p->IsIfcClass("IfcWindow") || p->IsIfcClass("IfcCurtainWall")) {
            TopoDS_ListOfShape L;
            L.Append(p->shape);
            TopoDS_Shape S = best_fitting_bbox(L);
            if (!S.IsNull()) {
                p->shape = S;
                std::cout << "\t[Info] Simplify " << p->guid << "\t" << p->IfcClass() << "\n";
            }
        }
    }
#else
    for (auto &p: products) {
        if (p.IsIfcClass("IfcDoor") || p.IsIfcClass("IfcWindow") || p.IsIfcClass("IfcCurtainWall")) {
            TopoDS_ListOfShape L;
            L.Append(p.shape);
            TopoDS_Shape S = best_fitting_bbox(L);
            if (!S.IsNull()) p.shape = S;
            std::cout << "\t[Info] Simplify " << p.guid << "\t" << p.IfcClass() << "\n";
        }
    }
#endif

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Simplify products", std::to_string(products.size()));
}

bool Kernel::check_shape_tolerances(const TopoDS_Shape &S, double critical) {

    bool b = true;

    ShapeAnalysis_ShapeTolerance ST;
    double s_min = ST.Tolerance(S, -1);
    double s_max = ST.Tolerance(S, 1);

    double v_min = ST.Tolerance(S, -1, TopAbs_VERTEX);
    double v_max = ST.Tolerance(S, 1, TopAbs_VERTEX);

    double e_min = ST.Tolerance(S, -1, TopAbs_EDGE);
    double e_max = ST.Tolerance(S, 1, TopAbs_EDGE);

    double f_min = ST.Tolerance(S, -1, TopAbs_FACE);
    double f_max = ST.Tolerance(S, 1, TopAbs_FACE);

    if (v_max >= critical) {
        std::cout << "Vertex tolerance bigger than critical value: " << v_max << "\t" << critical << "\n";
        b = false;
    }

    if (e_max >= critical) {
        std::cout << "Edge tolerance bigger than critical value:   " << e_max << "\t" << critical << "\n";
        b = false;
    }

    if (f_max >= critical) {
        std::cout << "Face tolerance bigger than critical value:   " << f_max << "\t" << critical << "\n";
        b = false;
    }

    if (v_max + 1.0e-8 < s_max) {
        std::cout << "Vertex tolerance smaller than shape maximum value: " << v_max << "\t" << s_max << "\n";
        b = false;
    }

    if (e_max + 1.0e-8 < s_max) {
        std::cout << "Edge tolerance smaller than shape maximum value:   " << e_max << "\t" << s_max << "\n";
        b = false;
    }

    if (f_max + 1.0e-8 < s_max) {
        std::cout << "Face tolerance smaller than shape maximum value:   " << f_max << "\t" << s_max << "\n";
        b = false;
    }

    if (!b) {
        std::cout << "\tShape:  " << s_min << "\t" << s_max << "\n";
        std::cout << "\tVertex: " << v_min << "\t" << v_max << "\n";
        std::cout << "\tEdge:   " << e_min << "\t" << e_max << "\n";
        std::cout << "\tFace:   " << f_min << "\t" << f_max << "\n";
    }

    return b;
}

void Kernel::check_product_tolerances(std::list<Product> &products, double critical) {

    auto start = std::chrono::high_resolution_clock::now();

    for (auto &p: products)
        bool good = check_shape_tolerances(p.shape, critical);

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Check product tolerances", std::to_string(products.size()));
}

void Kernel::fuse_products(std::list<Product> &products, const double tol) const {

    auto start = std::chrono::high_resolution_clock::now();

#ifdef PARALLEL_PROCESSING

    std::vector<Product *> V;
    V.reserve(products.size());
    for (auto &product: products)
        V.push_back(&product);

#pragma omp parallel for default(none) shared(V, tol) num_threads(num_threads)
    for (unsigned int i = 0; i < V.size(); i++) {
        Product *p = V[i];
        if (p->valid) continue;
        TopoDS_Shape S = fuse_shape(p->shape, tol, p->guid);
        if (!S.IsNull()) p->shape = S;
    }
#else
    for (auto &p: products) {
        if (p.valid) continue;
        TopoDS_Shape S = fuse_shape(p.shape, tol, p.guid);
        if (!S.IsNull()) p.shape = S;
    }
#endif

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Fuse products", std::to_string(products.size()));
}

void Kernel::cut_products(std::list<Product> &products, const std::set<std::string> &C, const double fuzzy_tol) const {

    auto start = std::chrono::high_resolution_clock::now();

    rtree_lib::RTree<Product *, double, 3, double> tree;
    std::map<Product *, std::pair<gp_Pnt, gp_Pnt >> M;

    // populate rtree
    for (auto &p: products) {
        Bnd_Box bnd = aabb(p.shape, fuzzy_tol);
        gp_Pnt Pmin = bnd.CornerMin();
        gp_Pnt Pmax = bnd.CornerMax();
        M[&p] = std::make_pair(Pmin, Pmax);
        double min[3] = {Pmin.X(), Pmin.Y(), Pmin.Z()};
        double max[3] = {Pmax.X(), Pmax.Y(), Pmax.Z()};
        tree.Insert(min, max, &p);
    }

    // get unique pairs
    std::set<std::set<Product *>> Pairs;

    for (auto &p: products) {

        // get all nbs
        gp_Pnt Pmin = M[&p].first;
        gp_Pnt Pmax = M[&p].second;
        double min[3] = {Pmin.X(), Pmin.Y(), Pmin.Z()};
        double max[3] = {Pmax.X(), Pmax.Y(), Pmax.Z()};

        std::set<Product *> comp;

        tree.Search(min, max, [&comp](Product *found_product) {
            comp.insert(found_product);
            return true;
        });

        comp.erase(&p);

        for (auto &p2: comp) {
            std::set<Product *> t = {&p, p2};
            Pairs.insert(t);
        }
    }

    // delete unwanted pairs and store tasks
    std::map<Product *, std::set<Product * >> tasks;
    for (auto &pair: Pairs) {

        Product *Prod1;
        Product *Prod2;

        unsigned int i = 0;
        for (auto &t: pair) {
            if (i == 0) Prod1 = t;
            else if (i == 1) Prod2 = t;
            else break;
            i++;
        }

        // skip
        if (C.find(Prod1->guid) == C.end() && C.find(Prod2->guid) == C.end()) continue;

        double V1 = volume(Prod1->shape);
        double V2 = volume(Prod2->shape);

        if (V1 > V2) tasks[Prod1].insert(Prod2);
        else tasks[Prod2].insert(Prod1);
    }

    // prepare parallel vector
    std::vector<Product *> V;
    V.reserve(products.size());
    for (auto &product: products)
        V.push_back(&product);

#pragma omp parallel for default(none) shared(V, tasks, fuzzy_tol, std::cerr) schedule(dynamic) num_threads(num_threads)
    for (auto p: V) {
        if (tasks.find(p) == tasks.end()) continue;
        if (p->shape.IsNull()) continue;

        std::set<Product *> comp = tasks[p];

        for (auto &p2: comp) {
            if (p->shape.IsNull() || p2->shape.IsNull()) continue;

            BRepAlgoAPI_Cut builder;

            TopTools_ListOfShape aLS;
            aLS.Append(p->shape);
            TopTools_ListOfShape aLT;;
            aLT.Append(p2->shape);

            builder.SetArguments(aLS);
            builder.SetTools(aLT);
            builder.SetFuzzyValue(fuzzy_tol);
            builder.SetCheckInverted(true);
            builder.SetNonDestructive(true);

            builder.Build();

            if (builder.HasWarnings()) {
                Standard_Integer iDMin1, iDMax1, iDMin2, iDMax2;
                BOPTools_AlgoTools::Dimensions(p->shape, iDMin1, iDMax1);
                BOPTools_AlgoTools::Dimensions(p2->shape, iDMin2, iDMax2);
                std::cerr << "Warnings " << "\tArg: " << p->guid << "(" << iDMin1 << ", " << iDMax1 << " )\tTool: " << p2->guid << "(" << iDMin2 << ", " << iDMax2 << ")" << std::endl;
                builder.DumpWarnings(std::cerr);
                std::cerr << "/Warnings " << std::endl;
            }
            if (builder.HasErrors()) { // For Boolean operation CUT the minimal dimension of Tools should not be less than the maximal dimension of Objects. (BOPAlgo/BOPAlgo_BOP.cxx)
                // So it's not allowed that a solid is cut by a shell
                Standard_Integer iDMin1, iDMax1, iDMin2, iDMax2;
                BOPTools_AlgoTools::Dimensions(p->shape, iDMin1, iDMax1);
                BOPTools_AlgoTools::Dimensions(p2->shape, iDMin2, iDMax2);
                std::cerr << "Errors " << "\tArg: " << p->guid << "(" << iDMin1 << ", " << iDMax1 << " )\tTool: " << p2->guid << "(" << iDMin2 << ", " << iDMax2 << ")" << std::endl;
                builder.DumpErrors(std::cerr);
                std::cerr << "/Errors " << std::endl;
            }

            if (builder.Shape().IsNull())
                std::cerr << "Result is null " << "\t" << p->guid << "\t" << p2->guid << std::endl; // is not always bad, e.g. when element is in another element

            if (!builder.HasErrors())
                p->shape = builder.Shape();
        }
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Cut products", std::to_string(products.size()));

    /* Brute force, redundant cuts. If two processes are calculating same pair (visa versa), a hole is created
         rtree_lib::RTree<Product *, double, 3, double> tree;
    std::map<Product *, std::pair<gp_Pnt, gp_Pnt>> M;

    // populate rtree
    for (auto &p : products) {
        Bnd_Box bnd = aabb(p.shape, fuzzy_tol);
        gp_Pnt Pmin = bnd.CornerMin();
        gp_Pnt Pmax = bnd.CornerMax();
        M[&p] = std::make_pair(Pmin, Pmax);
        double min[3] = {Pmin.X(), Pmin.Y(), Pmin.Z()};
        double max[3] = {Pmax.X(), Pmax.Y(), Pmax.Z()};
        tree.Insert(min, max, &p);
    }

    std::vector<Product *> V;
     V.reserve(products.size());
    for (auto &product : products)
        V.push_back(&product);

#pragma omp parallel for default(none) shared(V, C, M, tree, fuzzy_tol, std::cerr) schedule(static, 15) num_threads(num_threads)
    for (unsigned int i = 0; i < V.size(); i++) {
        Product *p = V[i];

        if (p->shape.IsNull()) continue;
        if (C.find(p->guid) == C.end()) continue;

        gp_Pnt Pmin = M[p].first;
        gp_Pnt Pmax = M[p].second;
        double min[3] = {Pmin.X(), Pmin.Y(), Pmin.Z()};
        double max[3] = {Pmax.X(), Pmax.Y(), Pmax.Z()};

        std::set<Product *> comp;

        tree.Search(min, max, [&comp](Product *found_product) {
            comp.insert(found_product);
            return true;
        });

        comp.erase(p);

        for (auto &p2 : comp) {
            if (p->shape.IsNull() || p2->shape.IsNull()) continue;

            BRepAlgoAPI_Cut builder;

            TopTools_ListOfShape aLS;
            aLS.Append(p->shape);
            TopTools_ListOfShape aLT;;
            aLT.Append(p2->shape);
            builder.SetArguments(aLS);

            builder.SetTools(aLT);
            builder.SetFuzzyValue(fuzzy_tol);
            builder.SetCheckInverted(true);
            builder.SetNonDestructive(true);

            builder.Build();

            if (builder.HasWarnings()) {
                std::cerr << "Warnings " << "\t" << p->guid << "\t" << p2->guid << std::endl;
                builder.DumpWarnings(std::cerr);
            }
            if (builder.HasErrors()) {
                std::cerr << "Errors" << "\t" << p->guid << "\t" << p2->guid << std::endl;
                builder.DumpErrors(std::cerr);
            }

            p->shape = builder.Shape();
        }

    }
     */
}

void Kernel::add_prism_faces_by_guid(std::list<oFace> &orig_faces, double offset, const std::set<std::string> &offset_guids) {

    auto start = std::chrono::high_resolution_clock::now();

    std::list<oFace> prism_faces;

    for (const auto &O: orig_faces) {
        //if(area(O.face)>2) continue;
        if (offset_guids.find(O.RelProduct()->guid) == offset_guids.end()) continue;
        if (Topo(O.face).wires().Size() > 1) continue; // don't extrude faces with holes
        if (O.Opening() != nullptr) continue; // skip wall faces that are product of cut with window/door

        TopoDS_Shape P = prism_from_face(O.face, O.Normal(), offset, true, true);

        for (const auto &F: Topo(P).faces()) {
            //if (fabs(Kernel::face_normal(TopoDS::Face(F)).Z()) > 0.9 || face_center(TopoDS::Face(F)).Z() < 0.3) continue;
            prism_faces.emplace_back(TopoDS::Face(F), O.RelProduct(), O.ShellID());
            prism_faces.back().SetIsOffset(true);
        }
    }

    for (const auto &f: prism_faces)
        orig_faces.push_back(f);

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Add prism faces", std::to_string(prism_faces.size()) + "/" + std::to_string(orig_faces.size()));
}

void Kernel::add_prism_faces(std::list<oFace> &orig_faces, double offset, const std::set<std::string> &classes_for_face_extension) {

    auto start = std::chrono::high_resolution_clock::now();

    std::list<oFace> prism_faces;

    for (const auto &O: orig_faces) {

        bool to_extend = false;
        for (const auto &cl: classes_for_face_extension)
            if (O.IsIfcClass(cl)) {
                to_extend = true;
                break;
            }

        if (!to_extend)
            continue;

        if (Topo(O.face).wires().Size() > 1) continue; // don't extrude faces with holes
        if (O.Opening() != nullptr) continue; // skip wall faces that are product of cut with window/door

        TopoDS_Shape P = prism_from_face(O.face, O.Normal(), offset, true, true);

        for (const auto &F: Topo(P).faces()) {
            prism_faces.emplace_back(TopoDS::Face(F), O.RelProduct(), O.ShellID());
            prism_faces.back().SetIsOffset(true);
        }
    }

    for (const auto &f: prism_faces)
        orig_faces.push_back(f);

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Add prism faces", std::to_string(prism_faces.size()) + "/" + std::to_string(orig_faces.size()));
}

void Kernel::add_offset_face_reversed_duplicates(std::list<cFace> &cFaces, unsigned int &fid) {

    auto start = std::chrono::high_resolution_clock::now();

    std::list<cFace> offset_duplicates;

    for (auto &cface: cFaces) {
        if (!cface.IsOffset()) continue;
        offset_duplicates.push_back(cface);
        offset_duplicates.back().SetID(fid);
        offset_duplicates.back().face = TopoDS::Face(cface.face.Complemented());
        offset_duplicates.back().SetNormalStatus(FACE_NORMAL_KNOWN);
        fid++;
    }
    for (auto &f: offset_duplicates) {
        cFaces.push_back(f);

        // add adjacency to other faces
        cFaces.back().UpdateHalfEdges();
        for (auto &it: cFaces.back().adjacentfaces) {
            for (auto &a: it.second)
                a->adjacentfaces[it.first].push_back(&cFaces.back());
        }
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Add reversed duplicates of offset faces", std::to_string(cFaces.size()));
}

void Kernel::nonoffset_inner_or_coplanar_faces(std::list<cFace> &cFaces, std::list<cFace> &L, bool consider_duplicates) {

    if (!consider_duplicates) {

        for (auto &cface: cFaces)
            if (!cface.IsOffset() && (cface.IsInner() || cface.IsCoplanar()))
                L.push_back(cface);

    } else { // do not add, if there is the same TopoDS_Face in cFaces

        std::set<cFace *> del;
        std::set<unsigned int> ids;

        for (auto &cface: cFaces)
            if (!cface.IsOffset() && (cface.IsInner() || cface.IsCoplanar()))
                del.insert(&cface);
            else
                ids.insert(cface.FaceID());

        for (const auto &d: del)
            if (ids.find(d->FaceID()) == ids.end())
                L.push_back(*d);
            else
                std::cout << "[Info] Does not add face " << d->Info() << " to inner faces because of coplanar face.\n";
    }
}

void Kernel::nonoffset_hanging_faces(std::list<cFace> &cFaces, std::list<cFace> &L) {
    for (auto &cface: cFaces)
        if (!cface.IsOffset() && cface.IsHanging())
            L.push_back(cface);
}

void Kernel::add_clipping_box_to_tree(gp_Vec v, double l, cFace *f, cface_tree3D &tree_inner) {
    v.Scale(l);
    Bnd_Box bnd = aabb(BRepPrimAPI_MakePrism(f->face, v).Shape(), 0.001);
    double min[3] = {bnd.CornerMin().X(), bnd.CornerMin().Y(), bnd.CornerMin().Z()};
    double max[3] = {bnd.CornerMax().X(), bnd.CornerMax().Y(), bnd.CornerMax().Z()};
    tree_inner.Insert(min, max, f);
}

void Kernel::find_clipping_pairs(cfaceSetMap &M_planar, cfaceSetMap &M_curved, std::list<cFace> &cFaces, std::list<cFace> &innerFaces, double depth_2b) {

    auto start = std::chrono::high_resolution_clock::now();

    cface_tree3D tree_inner;  // fill rtree
    for (auto &f: innerFaces) {
        if (f.IsVirtual()) continue;

        if (f.NormalStatus() == FACE_NORMAL_KNOWN)
            add_clipping_box_to_tree(f.FixedFaceNormal(), depth_2b, &f, tree_inner);
        else {
            auto n = f.FaceNormal();
            add_clipping_box_to_tree(n, depth_2b, &f, tree_inner);
            add_clipping_box_to_tree(n.Reversed(), depth_2b, &f, tree_inner);
        }
    }

    for (auto &cface: cFaces) {
        if (cface.IsVirtual()) continue;

        Bnd_Box bnd = aabb(cface.face, 0.001);
        double min[3] = {bnd.CornerMin().X(), bnd.CornerMin().Y(), bnd.CornerMin().Z()};
        double max[3] = {bnd.CornerMax().X(), bnd.CornerMax().Y(), bnd.CornerMax().Z()};
        std::list<cFace *> found_inner_faces;
        tree_inner.Search(min, max, [&found_inner_faces](cFace *found_inner_face) {
            found_inner_faces.push_back(found_inner_face);
            return true;
        });

        for (auto &f: found_inner_faces) {

            if (f->FaceID() == cface.FaceID()) {
                std::cout << "[Info] Skip coplanar inner face. " << cface.Info() << "\t" << f->Info() << "\n";
                continue;
            }
            //if (cface.FixedFaceNormal().Angle(f->FixedFaceNormal()) > max_angle) continue; // no need to check here, because normals could be unknown
            if (cface.IsPolygon() && f->IsPolygon()) M_planar[&cface].insert(f);
            else {
                std::cerr << "non_planar" << std::endl;
                M_curved[&cface].insert(f);
            }
        }
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Find clip face pairs using non-sb faces", std::to_string(M_planar.size()) + ", " + std::to_string(M_curved.size()));
}

TopoDS_ListOfShape Kernel::clip_faces(std::list<cFace> &cFaces, cfaceSetMap &M_planar, cfaceSetMap &M_curved, double depth_2b, double fuzzy_tol) {

    auto start = std::chrono::high_resolution_clock::now();

    TopoDS_ListOfShape L;

    for (auto &cface: cFaces) {

        TopoDS_ListOfShape Clips;
        Clipper clip(cface.face, cface.FixedFaceNormal(), true);

        for (auto &f: M_planar[&cface]) {
            clip.clip(f->face, f->FixedFaceNormal());
            if (clip.success) Clips.Append(clip.Result);
        }

        for (
            auto &f
                : M_curved[&cface]) {
            TopoDS_ListOfShape R = clip_faces(cface.face, f->face, depth_2b);
            Clips.
                    Append(R);
        }

        if (Clips.

                Size()

            < 2)
            L.
                    Append(Clips);
        else {
            BOPAlgo_Builder builder;
            for (
                auto &f
                    : Clips)
                builder.
                        AddArgument(f);
            builder.SetCheckInverted(false);
            builder.SetFuzzyValue(0.1 * fuzzy_tol);
            builder.SetRunParallel(true);
            builder.SetUseOBB(false);
            builder.SetToFillHistory(false);
            builder.

                    Perform();

            TopoDS_ListOfShape l = Topo(builder.Shape()).faces();

            TopoDS_ListOfShape del;
            for (
                const auto &t
                    : l)
                if (
                        area(t)
                        < 0.000001)
                    del.
                            Append(t);
            if (!del.

                    IsEmpty()

                    )
                for (
                    const auto &d
                        : del)
                    l.
                            Remove(d);

            L.
                    Append(l);
        }
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout <<
              print_time(elapsed
                                 .

                                         count(),

                         "Clip faces",
                         std::to_string(L
                                                .

                                                        Size()

                         ));

    return
            L;
}

TopoDS_ListOfShape Kernel::clip_faces(const TopoDS_Face &F1, const TopoDS_Face &F2, double depth_2b) {

    TopoDS_ListOfShape R;
    gp_Vec v; // projection direction

    // If F1 is not planar, estimate the projection vector by connection face centers
    if (face_is_polygon(F1))
        v = gp_Vec(face_normal(F1)).Scaled(depth_2b);
    else
        v = gp_Vec(face_center(F2), face_center(F1)).Scaled(depth_2b);

    // if F2 is not planar, a triangulation is needed, because curved wires/faces can't be projected or extruded correctly
    if (face_is_polygon(F2)) {
        TopoDS_Shape prism = BRepPrimAPI_MakePrism(F2, v).Shape();
        TopoDS_ListOfShape L = Topo(BRepAlgoAPI_Common(F1, prism).Shape()).faces();

        if (L.IsEmpty()) return R;

        for (const auto &r: L)
            R.Append(r);
    } else {
        TopoDS_Shape S = polygonize_shape(F2, 2, false, 2, false);
        if (S.IsNull()) return R;
        for (const auto &f: Topo(S).faces()) {
            TopoDS_Shape prism = BRepPrimAPI_MakePrism(f, v).Shape();
            TopoDS_ListOfShape L = Topo(BRepAlgoAPI_Common(F1, prism).Shape()).faces();

            if (L.IsEmpty()) return R;

            for (const auto &r: L)
                R.Append(r);
        }
    }

    return R;
}

void Kernel::check_containment_of_fuse_in_cFaces(const TopoDS_Shape &fuse, std::list<cFace> &cFaces) {

    auto start = std::chrono::high_resolution_clock::now();

    std::set<unsigned int> L;
    for (auto &cface: cFaces)
        L.insert(cface.FaceID());

    TopoDS_ListOfShape S = Topo(fuse).faces();

    for (const auto &f: S)
        if (L.find(hash(f)) == L.end())
            std::cerr << "[Error] Fuse face not contained in cFaces! " << hash(f) << std::endl;

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Check containment in cFaces", std::to_string(S.Size()) + ", " + std::to_string(cFaces.size()));
}

bool Kernel::fuse_clipped_faces(TopoDS_Shape &fuse, std::list<cFace> &cFaces, const TopoDS_ListOfShape &clips, std::list<cFace> &new_cFaces, double fuzzy_tol, unsigned int &fid) {

    auto start = std::chrono::high_resolution_clock::now();

    // populate builder
    BOPAlgo_Builder builder;
    for (const auto &f: cFaces)
        builder.AddArgument(f.face);
    for (const auto &f: clips)
        builder.AddArgument(f);

    // builder.SetNonDestructive; // Safe input shapes option allows preventing modification of the input shapes
    builder.SetCheckInverted(false);
    builder.SetFuzzyValue(fuzzy_tol);
    builder.SetRunParallel(true);
    builder.SetUseOBB(true);
    builder.SetToFillHistory(true);

    builder.Perform();

    if (!builder.HasModified()) {
        std::cerr << "Nothing was modified." << std::endl;
        return false;
    }

    for (auto &cface: cFaces) {

        bool isDeleted = builder.IsDeleted(cface.face); // TODO in case face was deleted, don't add to cfaces (maybe add the generated?). Without ignoring, the cface.face will not be present is fuse shape

        if (isDeleted) {
            auto Generated = builder.Generated(cface.face);
            std::cerr << "[Error] Face was deleted and may have been generated new faces! " << cface.Info() << "\t" << isDeleted << "\t" << Generated.Size() << std::endl;
        }

        auto L = builder.Modified(cface.face);

        if (L.IsEmpty()) {  // face was not modified
            new_cFaces.emplace_back(cface.face, cface.Ancestor(), fid); // TODO wichtige infos mit rübernehmen
            fid++;
        } else // face was modified into one or more new faces
            for (const TopoDS_Shape &shape: L) {
                new_cFaces.emplace_back(TopoDS::Face(shape), cface.Ancestor(), fid); // new face, pointer on oFace
                fid++;
            }
    }

    fuse = builder.Shape();

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Fuse original faces", std::to_string(Topo(fuse).faces().Size()) + ", " + std::to_string(cFaces.size()));

    return true;
}

void Kernel::find_clipping_pairs_2(cfaceSetMap &M_planar, cfaceSetMap &M_curved, std::list<cFace> &cFaces, double depth_2b, double min_angle) {

    auto start = std::chrono::high_resolution_clock::now();

    cface_tree3D tree_1st;  // fill rtree
    for (auto &f: cFaces) {
        if (f.IsVirtual()) continue;
        gp_Vec v2(f.FixedFaceNormal().Reversed());
        v2.Scale(depth_2b);
        Bnd_Box bnd2 = aabb(BRepPrimAPI_MakePrism(f.face, v2).Shape(), 0.001);
        double min2[3] = {bnd2.CornerMin().X(), bnd2.CornerMin().Y(), bnd2.CornerMin().Z()};
        double max2[3] = {bnd2.CornerMax().X(), bnd2.CornerMax().Y(), bnd2.CornerMax().Z()};
        tree_1st.Insert(min2, max2, &f);
    }

    for (auto &cface: cFaces) {
        if (cface.IsVirtual()) continue;

        Bnd_Box bnd = aabb(cface.face, 0.001);
        double min[3] = {bnd.CornerMin().X(), bnd.CornerMin().Y(), bnd.CornerMin().Z()};
        double max[3] = {bnd.CornerMax().X(), bnd.CornerMax().Y(), bnd.CornerMax().Z()};
        std::list<cFace *> found_faces;
        tree_1st.Search(min, max, [&found_faces](cFace *found_inner_face) {
            found_faces.push_back(found_inner_face);
            return true;
        });

        for (auto &f: found_faces) {

            if (cface.ID() == f->ID() || cface.FixedFaceNormal().Angle(f->FixedFaceNormal()) < min_angle) continue;

            if (f->FaceID() == cface.FaceID()) {
                std::cout << "[Info] Skip coplanar space boundary. " << cface.Info() << "\t" << f->Info() << "\n";
                continue;
            }

            if (cface.IsPolygon() && f->IsPolygon()) M_planar[&cface].insert(f);
            else {
                std::cerr << "non_planar" << std::endl;
                M_curved[&cface].insert(f);
            }
        }
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Find clip face pairs between sbs", std::to_string(M_planar.size()) + ", " + std::to_string(M_curved.size()));
}

void Kernel::spaces_introduce_second_level_cfaces(std::list<Space> &spaces) {
    for (auto &space: spaces)
        space.get_second_level();
}

void Kernel::clip_face_add_single_face(cFace &cface, std::list<cFace> &sec, unsigned int &fid) {

    sec.emplace_back(cface.face, cface.Ancestor(), fid, &cface);
    fid++;
    sec.back().SetSpace(cface.RelSpace());
    cFace::AddSuperfaceSubfaceRelationship(&cface, &sec.back());
}

void Kernel::clip_face_add_multiple_faces(cFace &cface, const TopoDS_Shape &C, std::list<cFace> &sec, unsigned int &fid) {

    std::list<cFace> T;

    for (const auto &f: Topo(C).faces()) {
        T.emplace_back(TopoDS::Face(f), cface.Ancestor(), fid, &cface);
        fid++;
        T.back().SetSpace(cface.RelSpace());
    }

    // check_fuse(C);
    // check_redundant_cFaces(T);
    // check_containment_in_fuse(C, T, false);
    // check_containment_of_fuse_in_cFaces(C, T);

    std::list<cFace *> cluster;
    for (auto &t: T) {
        t.UpdateHalfEdges();
        sec.push_back(t);
        cluster.push_back(&sec.back());
        cFace::AddSuperfaceSubfaceRelationship(&cface, &sec.back());
    }

    std::unordered_map<unsigned int, cFace *> id2cFace = setup_face_cFace_map(cluster);
    TopTools_IndexedDataMapOfShapeListOfShape edgeFaceMap;
    TopExp::MapShapesAndAncestors(C, TopAbs_EDGE, TopAbs_FACE, edgeFaceMap);
    for (auto &t: cluster) t->UpdateFaceAdjacencies(edgeFaceMap, id2cFace);
}

TopoDS_Shape Kernel::clip_face_occ(cFace &cface, const std::set<cFace *> &planars, const std::set<cFace *> &curved, double length, double fuzzy_tol) {

    TopoDS_ListOfShape Clips;
    Clipper clip(cface.face, cface.FixedFaceNormal(), true);

    for (auto &f: planars) {
        clip.clip(f->face, f->FixedFaceNormal());
        if (clip.success) Clips.Append(clip.Result);
    }
    for (auto &f: curved) {
        TopoDS_ListOfShape R = clip_faces(cface.face, f->face, length);
        Clips.Append(R);
    }

    // check
    for (auto &f: Clips)
        check_face(TopoDS::Face(f), hash(f), cface.Info(), 1.0e-5, 1.0e-6);

    // store
    if (Clips.IsEmpty()) { // face was not clipped, just copy

        return {};

    } else { // face was clipped. get smallest faces possible by fusing

        BOPAlgo_Builder builder;
        for (auto &f: Clips) builder.AddArgument(f);
        builder.AddArgument(cface.face);
        builder.SetCheckInverted(false);
        builder.SetNonDestructive(true);
        builder.SetFuzzyValue(fuzzy_tol);
        builder.SetRunParallel(true);
        builder.SetUseOBB(false);
        builder.SetToFillHistory(false);
        builder.Perform();

        if (builder.HasWarnings()) {
            std::cerr << "Warnings for " << cface.Info() << ":" << std::endl;
            builder.DumpWarnings(std::cerr);
        }
        if (builder.HasErrors()) {
            std::cerr << "Errors for " << cface.Info() << ":" << std::endl;
            builder.DumpErrors(std::cerr);
        }

        return BRepAlgoAPI_Common(cface.face, builder.Shape()).Shape();
    }
}

TopoDS_Shape Kernel::clip_face(cFace &cface, const std::set<cFace *> &planars, const double length, double fuzzy_tol) {

    typedef std::vector<gp_Pnt2d> chain2D;
    std::vector<chain2D> ClipWires, Holes2D;
    chain2D Wire2D;

    Clipper clip(cface.face, cface.FixedFaceNormal(), false);

    for (auto &f: planars) {
        clip.clip(f->face, f->FaceNormal());
        if (clip.success) {
            ClipWires.insert(ClipWires.end(), clip.Result2DWires.begin(), clip.Result2DWires.end());
            Wire2D = clip.F1_wire2D;
            Holes2D = clip.F1_holes2D;
        }
    }

    // store
    if (ClipWires.empty()) { // face was not clipped, just copy

        return {};

    } else { // face was clipped. get smallest faces possible by fusing

        gp_Pnt P = BRep_Tool::Pnt(TopoDS::Vertex(Topo(cface.face).vertices().First()));
        const gp_Dir n = cface.FixedFaceNormal();
        double a;
        unsigned int proj_axis;
        Clipper::calculate_algebraic_values(P, n, a, proj_axis);

        Topo2D G(fuzzy_tol, proj_axis, n, a);

        G.Add(Wire2D, 0, TYPE2D_OUTER_ORIG);
        int i = 1;
        for (const auto &w: Holes2D) {
            G.Add(w, i, TYPE2D_HOLE);
            i++;
        }
        for (const auto &w: ClipWires) {
            G.Add(w, i, TYPE2D_INNER);
            i++;
        }

        TopoDS_Shape C;
        if (!G.Fuse(C)) {
            std::cerr << "Not fused!" << std::endl;
            return {};
        } else {
            return C;
        }
    }
}

bool Kernel::clip_faces_surface_check(const TopoDS_Shape &C, const cFace &cface, bool level) {

    if (C.IsNull()) return true;

    double A1 = cface.SurfaceArea();
    double A2 = area(C);
    double v = fabs(A1 - A2) / A1;

    if (v < 0.02) return true;

    if (level) std::cerr << "[Warning] Surface area changed! " << cface.Info() << "\t" << A1 << "\t" << A2 << "\t" << v << std::endl;
    else std::cout << "[Info] Surface area changed! " << cface.Info() << "\t" << A1 << "\t" << A2 << "\t" << v << "\n";

    return false;
}

void Kernel::clip_faces(std::list<cFace> &cFaces, std::list<cFace> &sec, cfaceSetMap &M_planar, cfaceSetMap &M_curved, double depth_2b, double fuzzy_tol, unsigned int &fid) const {

    auto start = std::chrono::high_resolution_clock::now();

#ifdef PARALLEL_PROCESSING
    std::vector<cFace *> V;
    V.reserve(cFaces.size());

    std::vector<std::set<cFace *>> planars, curved;
    planars.reserve(cFaces.size());
    curved.reserve(cFaces.size());

    for (auto &f: cFaces) {
        V.push_back(&f);
        planars.push_back(M_planar[&f]);
        curved.push_back(M_curved[&f]);
    }

    std::vector<TopoDS_Shape> Compounds(V.size());

    unsigned int chunk = std::ceil(0.1 * cFaces.size() / num_threads);

#pragma omp parallel for default(none) shared(V, Compounds, planars, curved, depth_2b, fuzzy_tol, chunk) schedule(dynamic, chunk) num_threads(num_threads)
    for (unsigned int i = 0; i < V.size(); i++) {
        cFace *cface = V[i];
        if (curved[i].empty()) {
            Compounds[i] = clip_face(*cface, planars[i], depth_2b, fuzzy_tol);
            if (!clip_faces_surface_check(Compounds[i], *cface))
                Compounds[i] = clip_face_occ(*cface, planars[i], curved[i], depth_2b, fuzzy_tol);
        } else
            Compounds[i] = clip_face_occ(*cface, planars[i], curved[i], depth_2b, fuzzy_tol);
    }

    for (unsigned int i = 0; i < Compounds.size(); i++) {

        const TopoDS_Shape &C = Compounds[i];
        cFace *cface = V[i];

        if (C.IsNull())
            clip_face_add_single_face(*cface, sec, fid);
        else {
            clip_faces_surface_check(C, *cface, true);
            clip_face_add_multiple_faces(*cface, C, sec, fid);
        }
    }
#else
    for (auto &cface: cFaces) {

        TopoDS_Shape C;

        if (M_curved[&cface].empty()) {
            C = clip_face(cface, M_planar[&cface], depth_2b, fuzzy_tol);
            if (!clip_faces_surface_check(C, cface))
                C = clip_face_occ(cface, M_planar[&cface], M_curved[&cface], depth_2b, fuzzy_tol);
        }
        else
            C = clip_face_occ(cface, M_planar[&cface], M_curved[&cface], depth_2b, fuzzy_tol);

        if (C.IsNull())
            clip_face_add_single_face(cface, sec, fid);
        else {
            clip_faces_surface_check(C, cface, true);
            clip_face_add_multiple_faces(cface, C, sec, fid);
        }
    }
#endif

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Clip faces", std::to_string(sec.size()));
}

void Kernel::update_face_adjacencies(std::list<cFace *> &cFaces, const TopoDS_Shape &fuse) {

    auto start = std::chrono::high_resolution_clock::now();

    std::unordered_map<unsigned int, cFace *> id2cFace = setup_face_cFace_map(cFaces);

    TopTools_IndexedDataMapOfShapeListOfShape edgeFaceMap;
    TopExp::MapShapesAndAncestors(fuse, TopAbs_EDGE, TopAbs_FACE, edgeFaceMap);

    for (auto &cface: cFaces)
        cface->UpdateFaceAdjacencies(edgeFaceMap, id2cFace);

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Save face adjacency information", std::to_string(cFaces.size()));
}

std::unordered_map<unsigned int, cFace *> Kernel::setup_face_cFace_map(std::list<cFace *> &cFaces) {

    std::unordered_map<unsigned int, cFace *> id2cFace;
    for (auto &cface: cFaces) {
        if (id2cFace.find(cface->FaceID()) != id2cFace.end())
            std::cerr << "Error: Two faces with same ID in cFaces list! " << cface->Info() << std::endl;
        id2cFace[cface->FaceID()] = cface;
    }
    return id2cFace;
}

void Kernel::site(std::unique_ptr<IfcParse::IfcFile> &model, IfcGeom::IteratorSettings settings, std::list<cFace> &cFaces, std::list<Space> &spaces, unsigned int &fid) const {

    auto start = std::chrono::high_resolution_clock::now();

    TopoDS_ListOfShape L = site_ifcsite(model, settings); // ifcsite shape
    if (L.IsEmpty())
        return;

    // find facade space
    Space *facade_space = nullptr;
    for (auto &space: spaces)
        if (space.is_facade) {
            facade_space = &space;
            break;
        }
    if (facade_space == nullptr) {
        std::cerr << "[Info] No facade space found." << std::endl;
        return;
    }
    Bnd_Box aabb_facade = Kernel::aabb(facade_space->shell, 0.01);

    // create site face list
    std::list<site_face> site_faces;
    double site_min = 1.0e9;
    double site_max = -1.0e9;

    for (auto &F: L) {
        Bnd_Box box = aabb(F, 1.0e-5);
        if (box.IsOut(aabb_facade)) continue;

        site_faces.emplace_back(TopoDS::Face(F), box, face_normal(TopoDS::Face(F)));
        if (site_faces.back().aabb.CornerMin().Z() < site_min) site_min = site_faces.back().aabb.CornerMin().Z();
        if (site_faces.back().aabb.CornerMax().Z() > site_max) site_max = site_faces.back().aabb.CornerMax().Z();
    }

    // create rtree
    rtree_lib::RTree<site_face *, double, 3, double> tree;

    for (auto &f: site_faces) {
        double min2[3] = {f.aabb.CornerMin().X(), f.aabb.CornerMin().Y(), f.aabb.CornerMin().Z()};
        double max2[3] = {f.aabb.CornerMax().X(), f.aabb.CornerMax().Y(), f.aabb.CornerMax().Z()};
        tree.Insert(min2, max2, &f);
    }

    cfaceSetMap oldToNew = site_clip(cFaces, facade_space, tree, site_min, site_max, fid); // clip faces
    site_update_face_relationships(cFaces, facade_space, oldToNew);

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Split external faces by site", std::to_string(cFaces.size()));
}

TopoDS_ListOfShape Kernel::site_ifcsite(std::unique_ptr<IfcParse::IfcFile> &model, IfcGeom::IteratorSettings settings) const {

    // find IfcSites
    std::set<std::string> guids;
    TopoDS_ListOfShape L;

    boost::shared_ptr<IfcEntityList> IfcEntityList = nullptr;
    try { IfcEntityList = model->instances_by_type("IfcSite"); }
    catch (...) {
        std::cout << "[Info] No IfcSites found in ifc file.\n";
        return L;
    }
    for (auto &E: *IfcEntityList) {
        std::string g = E->data().getArgument(0)->toString();
        if (g.size() < 5) continue;
        guids.insert(g.substr(1, g.size() - 2));
    }

    //*****************************************************
    // create shapes
    settings.set(IfcGeom::IteratorSettings::DISABLE_OPENING_SUBTRACTIONS, false);
    IfcGeom::attribute_filter attribute_filter;
    attribute_filter.include = true;
    attribute_filter.traverse = false;
    attribute_filter.attribute_name = "GlobalId";
    attribute_filter.populate(guids);
    std::vector<IfcGeom::filter_t> filter_funcs;
    filter_funcs.emplace_back(boost::ref(attribute_filter));
    IfcGeom::Iterator<real_t> geom_iterator(settings, model.get(), filter_funcs, num_threads);

    if (!geom_iterator.initialize()) {
        std::cout << "[Info] No IfcSite geometry found.\n";
        return L;
    }

    do {
        IfcGeom::Element<real_t> *geom_object = geom_iterator.get();
        if (!check_ifc_context(geom_object)) continue;

        TopoDS_Shape geom = geom_object_to_shape(geom_object);
        ShapeUpgrade_UnifySameDomain U(geom, false, true, false);
        U.SetAngularTolerance(Precision::Angular() * 2);
        U.SetLinearTolerance(Precision::Confusion() * 2);
        U.Build();
        TopoDS_Shape shp = U.Shape();

        for (auto &F: Topo(shp).faces()) {
            if (face_is_polygon(TopoDS::Face(F))) L.Append(F);
            else for (auto &f2: Topo(polygonize_shape(F, 12.0, false, 12.0, true)).faces()) L.Append(f2);
        }
    } while (geom_iterator.next());

    return L;
}

cfaceSetMap Kernel::site_clip(std::list<cFace> &cFaces, Space *facade_space, const rtree_lib::RTree<site_face *, double, 3, double> &tree, double site_min, double site_max, unsigned int &fid) const {

    cfaceSetMap oldToNew; // needed to update parent faces in later step

    for (auto &cface: facade_space->FirstLvl) {

        Bnd_Box box = aabb(cface->face, 1.0e-5);
        if (box.CornerMin().Z() > site_max || box.CornerMax().Z() < site_min) continue;

        std::set<site_face *> nbs;
        double min[3] = {box.CornerMin().X(), box.CornerMin().Y(), box.CornerMin().Z()};
        double max[3] = {box.CornerMax().X(), box.CornerMax().Y(), box.CornerMax().Z()};
        tree.Search(min, max, [&nbs, &cface](site_face *found) {
            if (found->n.Angle(cface->FixedFaceNormal()) < 0.017 || found->n.Angle(cface->FixedFaceNormal()) > 3.124)
                nbs.insert(found);
            return true;
        });

        if (nbs.empty()) continue;

        TopoDS_ListOfShape Clips;

        Clipper clip(cface->face, cface->FixedFaceNormal(), true);

        for (const auto &nb: nbs) {
            clip.clip(nb->face, face_normal(nb->face));
            if (clip.same) cface->SetInternalOrExternal(SB_IE_EXTERNAL_EARTH); // if an earth face is same size, cface is completely below site
            if (clip.success) {
                for (auto &t: clip.Result)
                    for (auto &t2: Topo(BRepAlgoAPI_Common(cface->face, t).Shape()).faces())
                        Clips.Append(t2);
            }
        }

        if (Clips.IsEmpty()) continue;

        BOPAlgo_Builder builder;
        for (auto &f: Clips)
            builder.AddArgument(f);
        builder.AddArgument(cface->face);
        builder.SetCheckInverted(false);
        //builder.SetNonDestructive(true);
        builder.SetFuzzyValue(1.0e-5);
        builder.SetRunParallel(true);
        builder.SetUseOBB(false);
        builder.SetToFillHistory(true);
        builder.Perform();

        if (builder.HasWarnings()) {
            std::cerr << "Warnings for " << cface->Info() << ":" << std::endl;
            builder.DumpWarnings(std::cerr);
        }
        if (builder.HasErrors()) {
            std::cerr << "Errors for " << cface->Info() << ":" << std::endl;
            builder.DumpErrors(std::cerr);
        }

        std::set<unsigned int> earth;

        for (auto &f: builder.Arguments()) {
            if (f.IsSame(cface->face)) continue;
            auto M = builder.Modified(f);
            for (const TopoDS_Shape &shape: M)
                earth.insert(hash(shape));
        }

        for (const auto &f: Topo(builder.Shape()).faces()) {
            cFaces.push_back(*cface);
            cFaces.back().face = TopoDS::Face(f);
            if (!cface->AngleToFaceSmaller90(cFaces.back()))
                cFaces.back().ComplementFace();
            cFaces.back().RelSpace()->FirstLvl.insert(&cFaces.back());
            cFaces.back().SetID(fid);
            fid++;
            oldToNew[cface].insert(&cFaces.back());
            if (earth.find(hash(f)) != earth.end()) cFaces.back().SetInternalOrExternal(SB_IE_EXTERNAL_EARTH);
        }
        cface->SetIsTrash(true);

/*
            std::list<cFace *> ext;
            if (ext.empty() || cface->Corresponding() == nullptr) continue;
            // corresponding
            cFace *corr = ext.front()->Corresponding();
            for (const auto &ex: ext) {
                gp_Vec v(corr->FixedFaceNormal());
                v.Scale(10);
                TopoDS_Shape prism = BRepPrimAPI_MakePrism(ex->face, v).Shape();
                TopoDS_Shape common = BRepAlgoAPI_Common(corr->face, prism).Shape();
                if (common.IsNull()) continue;
                TopoDS_ListOfShape CL = Topo(common).faces();
                if (CL.IsEmpty()) continue;
                TopoDS_Face F = TopoDS::Face(CL.First());
                cFaces.push_back(*corr);
                cFaces.back().face = TopoDS::Face(F);
                cFaces.back().RelSpace()->SecondLvl.insert(&cFaces.back());
                cFaces.back().SetCorresponding(ex);
                cFaces.back().SetID(fid);
                ex->SetCorresponding(&cFaces.back());
                fid++;
            }
            corr->SetIsTrash(true);
            */
    }

    return oldToNew;
}

void Kernel::site_update_face_relationships(std::list<cFace> &cFaces, Space *facade_space, cfaceSetMap &oldToNew) {

    // update parent faces
    // assumes no non-facade face has facade face as parent
    for (auto &cface: facade_space->FirstLvl) {
        if (cface->Parent() == nullptr) continue;
        if (cface->Parent()->IsTrash() && oldToNew.find(cface->Parent()) != oldToNew.end()) {
            std::set<cFace *> &possible_parents = oldToNew[cface->Parent()];
            if (possible_parents.empty()) continue;
            cFace *new_parent = *possible_parents.begin();
            double common = -1;
            for (const auto &c: oldToNew[cface->Parent()]) {
                double A = area(BRepAlgoAPI_Common(cface->face, c->face).Shape());
                if (A > common) {
                    common = A;
                    new_parent = c;
                }
            }
            cface->SetParent(new_parent);
        }
    }

    // update adjacence info, which is only used for shading calculation after this. because the generated faces are probably no shading faces, delete them.
    for (auto &cface: facade_space->FirstLvl) {
        if (cface->IsTrash()) continue;
        for (auto &nbs: cface->adjacentfaces)
            nbs.second.remove_if([](cFace *nb) { return nb->IsTrash(); });
    }

}

void Kernel::triangulate_cfaces_for_ray_tracer(std::list<cFace> &cFaces) const {

    auto start = std::chrono::high_resolution_clock::now();

#ifdef PARALLEL_PROCESSING
    std::vector<cFace *> V;
    V.reserve(cFaces.size());
    for (auto &f: cFaces)
        V.push_back(&f);

#pragma omp parallel for default(none) shared(V) num_threads(num_threads)
//pragma omp parallel for default(none) shared(V, std::cerr) num_threads(num_threads)
    for (auto cface: V) {
        //std::ostringstream address;
        //address << (void const *)&cface->face.Location();
        //std::string pnt = address.str();
        //if (cface->face.Location().HashCode(INT_MAX) != 1)
        //std::cerr << cface->Info() << "\tLocationHash " << cface->face.Location().HashCode(INT_MAX)<< "\tLocationP " << pnt << std::endl;
        TopLoc_Location L;
        Handle(Poly_Triangulation) T = BRep_Tool::Triangulation(cface->face, L);
        if (T.IsNull() || T->NbTriangles() == 0) // only triangulate, if no triangulation exists
            BRepMesh_IncrementalMesh(cface->face, 5.0, false, 5.0, false);
    }
#else
    for (auto &cface: cFaces) {
        TopLoc_Location L;
        Handle(Poly_Triangulation) T = BRep_Tool::Triangulation(cface.face, L);
        if (T.IsNull() || T->NbTriangles() == 0) // only triangulate, if no triangulation exists
            BRepMesh_IncrementalMesh(cface.face, 5.0, false, 5.0, false);
    }
#endif

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Triangulate faces for ray tracer", std::to_string(cFaces.size()));
}

void Kernel::create_mesh_for_ray_tracer(std::list<cFace> &cFaces, std::list<std::array<double, 3>> &vertices, std::list<std::array<unsigned int, 3>> &triangles, std::list<unsigned int> &attributes, std::unordered_map<unsigned int, cFace *> &M) {

    auto start = std::chrono::high_resolution_clock::now();

    double tol = 1.0e-7;

    unsigned int ca = 0;
    unsigned int ct = 0;

    for (auto &cface: cFaces) {

        if (cface.IsOffset()) std::cout << "[Info] Offset face. " << cface.Info() << "\n";

        ca += 1;
        M[ca] = &cface;

        TopLoc_Location L;
        Handle_Poly_Triangulation Poly_Triangulation = BRep_Tool::Triangulation(cface.face, L);

        if (Poly_Triangulation.IsNull()) {
            std::cerr << "[Warning] Null triangulation. " << cface.Info() << std::endl;
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

            if (p1.IsEqual(p2, tol) || p2.IsEqual(p3, tol) || p3.IsEqual(p1, tol)) {
                std::cerr << "[Warning] Overlapping points." << std::endl;
                continue;
            }
            if (are_points_colinear(p1, p2, p3, tol)) {
                std::cerr << "[Warning] Collinear points." << std::endl;
                continue;
            }
            vertices.push_back({p1.X(), p1.Y(), p1.Z()});
            vertices.push_back({p2.X(), p2.Y(), p2.Z()});
            vertices.push_back({p3.X(), p3.Y(), p3.Z()});

            triangles.push_back({ct, ct + 1, ct + 2});
            ct += 3;

            attributes.push_back(ca); // triangles from same cface share id and are therefore not recognized independently by ray tracer
        }
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Create mesh for ray tracer", std::to_string(vertices.size()) + ", " + std::to_string(triangles.size()) + ", " + std::to_string(attributes.size()));
}

void Kernel::create_rays_for_ray_tracer(std::list<cFace> &cFaces_2ndLvl, std::vector<cFace *> &send, std::vector<IntersectionRay> &rays, double tol) {

    auto start = std::chrono::high_resolution_clock::now();

    unsigned int n = cFaces_2ndLvl.size(); // std::count_if(cFaces_2ndLvl.begin(), cFaces_2ndLvl.end(), [&](const cFace &cface) { return cface.sb_type != SB_TYPE_2B; });
    rays.reserve(n);
    send.reserve(n);

    // create ray
    for (auto &cface: cFaces_2ndLvl) {
        if (cface.SBType() == SB_TYPE_2B) continue;
        IntersectionRay ray = cface.RayBehind(tol);
        rays.push_back(ray);
        send.push_back(&cface);
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Create rays", std::to_string(rays.size()));

    /*
    std::vector<cFace *> V;
    V.reserve(cFaces_2ndLvl.size());
    for (auto &cface: cFaces_2ndLvl) {
        if (cface.sb_type == SB_TYPE_2B) continue;
        V.push_back(&cface);
    }

    std::vector<cFace *> send;
    std::vector<IntersectionRay> rays;
    rays.reserve(V.size());
    send.reserve(V.size());

    const double tol = 1.0e-5;
    const double tol2 = 2 * tol;

#pragma omp parallel for default(none) shared(V, tol, tol2, rays, send) schedule (static) num_threads(num_threads)
    for (int i = 0; i < V.size(); i++) {
        cFace *cface = V[i];
        IntersectionRay ray = cface->create_ray_behind(tol, tol2);
        rays[i] = ray;
        send[i] = cface;
    }
*/
}

void Kernel::create_rays_for_ray_tracer(std::list<cFace> &cFaces_2ndLvl, std::vector<cFace *> &send, std::vector<IntersectionRay> &rays, double tol, double l) {

    auto start = std::chrono::high_resolution_clock::now();

    unsigned int n = cFaces_2ndLvl.size();
    rays.reserve(n);
    send.reserve(n);

    // create ray
    for (auto &cface: cFaces_2ndLvl) {
        IntersectionRay ray = cface.RayBehind(tol, l);
        rays.push_back(ray);
        send.push_back(&cface);
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Create rays", std::to_string(n));
}

void Kernel::perform_ray_tracing(IntersectorInterface &intersector, const std::vector<cFace *> &send, std::vector<IntersectionRay> &rays, double tol, std::unordered_map<unsigned int, cFace *> &M) {

    auto start = std::chrono::high_resolution_clock::now();

    for (unsigned int i = 0; i < send.size(); i++) {

        IntersectionRay &ray = rays[i];
        cFace *cface = send[i];

        std::multimap<double, unsigned int> hits = intersector.Perform(ray);

        if (hits.empty()) continue;

        for (const auto &hit: hits) {
            cFace *hit_cface = M[hit.second];
            double dist = hit.first - tol; // correct distance by length origin was moved
            dist = round_double_to_n_decimal_places(dist, 5);
            cface->AppendToMaterial(hit_cface->RelProduct(), dist);
        }
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Perform ray tracing", std::to_string(send.size()));

    /*
#pragma omp parallel for default(none) shared(rays, send, intersector, M, tol) schedule (static) num_threads(num_threads)
    for (int i = 0; i < send.size(); i++) {
        IntersectionRay &ray = rays[i];
        cFace *cface = send[i];

        std::multimap<double, unsigned int> hits = intersector.Perform(ray);

        if (hits.empty()) continue;

        for (const auto &hit: hits) {
            cFace *hit_cface = M[hit.second];
            double dist = hit.first - tol; // correct distance by length origin was moved
            dist = round_double_to_n_decimal_places(dist, 5);
            cface->materials[hit_cface->RelProduct()].insert(dist);
        }
    }
     */
}

void Kernel::complete_material_list_with_inner_faces(std::list<cFace> &inner_cFaces, std::list<cFace> &cFaces_2ndLvl) {

    // empty mesh data containers
    std::list<std::array<double, 3>> vertices;
    std::list<std::array<unsigned int, 3>> triangles;
    std::list<unsigned int> attributes;
    std::unordered_map<unsigned int, cFace *> M;

    // pre-processing
    triangulate_cfaces_for_ray_tracer(inner_cFaces);
    create_mesh_for_ray_tracer(inner_cFaces, vertices, triangles, attributes, M);

    // build intersector
    auto start = std::chrono::high_resolution_clock::now();

    IntersectorInterface intersector(vertices, triangles, attributes, 7, 10, false);
    // IntersectorInterface intersector("/home/fluid/Schreibtisch/a.stl", 7, 10, false);

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Build intersector", std::to_string(triangles.size()));

    // calculate rays
    std::vector<cFace *> send;
    std::vector<IntersectionRay> rays;
    const double tol = 1.0e-5;
    create_rays_for_ray_tracer(cFaces_2ndLvl, send, rays, tol);

    // perform intersections
    perform_ray_tracing(intersector, send, rays, tol, M);

    /* OpenCascade Intersector
     * auto start = std::chrono::high_resolution_clock::now();

    // create compound from inner faces
    TopoDS_ListOfShape L;
    for (const auto &f: inner_cFaces)
        L.Append(f.face);

    TopoDS_Compound Comp = compound_from_shape_list(L);

    // lookup table: face to cface
    // Attention: NCollection_DataMap is a map not a multimap. If there are two coplanar cfaces (sharing same topods_face), intersector will return topods_face twice, but map can only return one cface
    // Store a list as value to get all cfaces
    NCollection_DataMap<TopoDS_Face, std::list<cFace *>> id2cFace;
    for (auto &cface: inner_cFaces) {
        if (id2cFace.IsBound(cface.face)) {
            std::list<cFace *> l = id2cFace.Find(cface.face);
            l.push_back(&cface);
            id2cFace.Bind(cface.face, l);
        } else {
            std::list<cFace *> l = {&cface};
            id2cFace.Bind(cface.face, l);
        }
    }

    // create intersector
    const double tol = 0.00001;
    IntCurvesFace_ShapeIntersector intersector;
    intersector.Load(Comp, tol);

    // perform intersection
    const double ml = -1.0e-6; // ensure self-hit

    for (auto &cface: cFaces_2ndLvl) {
        if (cface.sb_type == SB_TYPE_2B) continue;

        // create line
        gp_Pnt C = cface.P; //point_on_face(cface.face, cface.ancestor->n);
        gp_Pnt P = move_point_along_scaled_unit_vector(C, cface.ancestor->n.Reversed(), 0);
        gp_Lin Line(P, cface.ancestor->n.Reversed());

        // find length of isct segment by evaluating the maximum distance in material list
        double l_line = 0;
        if (!cface.materials.empty())
            for (const auto &e: cface.materials)
                if (*e.second.rbegin() > l_line) l_line = *e.second.rbegin(); // set is ordered, so take last element
        l_line += -ml; // extend length

        // perform intersection
        intersector.Perform(Line, ml, l_line);

        if (intersector.NbPnt() > 0) {

            TopoDS_ListOfShape all_hit_faces;

            for (int i = 1; i <= intersector.NbPnt(); i++) {

                const auto &f = intersector.Face(i);

                // because of the mentioned map issue, topods_faces are checked here, if they were already hit. This approach has one disadvantage: for curved faces, where 2 intersections would occur at different distances would now only know one of the intersections
                // this behaviour is ok for 2nd lvl sb, because it is restricted to planar faces
                if (all_hit_faces.Contains(f))
                    continue;
                else
                    all_hit_faces.Append(f);

                for (const auto &hit_cface: id2cFace.Find(f)) {
                    auto dist = P.Distance(intersector.Pnt(i));
                    dist = round_double_to_n_decimal_places(dist, 5);
                    cface.materials[hit_cface->RelProduct()].insert(dist);
                }
            }
        }
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Calculate intersections with inner and coplanar faces", std::to_string(cFaces_2ndLvl.size()));*/
}

void Kernel::identify_sb_types_ray_tracing(std::list<cFace> &cFaces, gp_XYZ bounds_min, double transmission_length, double min_angle) {

    // empty mesh data containers
    std::list<std::array<double, 3>> vertices;
    std::list<std::array<unsigned int, 3>> triangles;
    std::list<unsigned int> attributes;
    std::unordered_map<unsigned int, cFace *> M;

    // pre-processing
    triangulate_cfaces_for_ray_tracer(cFaces);
    create_mesh_for_ray_tracer(cFaces, vertices, triangles, attributes, M);

    // build intersector
    auto start = std::chrono::high_resolution_clock::now();

    IntersectorInterface intersector(vertices, triangles, attributes, 7, 10, false);

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Build intersector", std::to_string(triangles.size()));

    // calculate rays
    const double tol = 1.0e-5;
    std::vector<cFace *> send;
    std::vector<IntersectionRay> rays;
    create_rays_for_ray_tracer(cFaces, send, rays, tol, transmission_length);

    // perform intersections
    perform_ray_tracing(intersector, send, rays, tol, M, min_angle, bounds_min);
}

void Kernel::identify_sb_types_ray_tracing_first_level(std::list<cFace> &cFaces, gp_XYZ bounds_min, double transmission_length) {

    // empty mesh data containers
    std::list<std::array<double, 3>> vertices;
    std::list<std::array<unsigned int, 3>> triangles;
    std::list<unsigned int> attributes;
    std::unordered_map<unsigned int, cFace *> M;

    // pre-processing
    triangulate_cfaces_for_ray_tracer(cFaces);
    create_mesh_for_ray_tracer(cFaces, vertices, triangles, attributes, M);

    // build intersector
    auto start = std::chrono::high_resolution_clock::now();

    IntersectorInterface intersector(vertices, triangles, attributes, 7, 10, false);

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Build intersector", std::to_string(triangles.size()));

    // calculate rays
    const double tol = 1.0e-5;
    std::vector<cFace *> send;
    std::vector<IntersectionRay> rays;
    create_rays_for_ray_tracer(cFaces, send, rays, tol, transmission_length);

    // perform intersections
    perform_ray_tracing_first_level(intersector, send, rays, tol, M, bounds_min);
}

void Kernel::perform_ray_tracing(IntersectorInterface &intersector, const std::vector<cFace *> &send, std::vector<IntersectionRay> &rays, double tol, std::unordered_map<unsigned int, cFace *> &M, double min_angle, gp_XYZ bounds_min) {

    auto start = std::chrono::high_resolution_clock::now();

    std::unordered_map<Product *, std::vector<cFace * >> virtuals;

    for (unsigned int i = 0; i < send.size(); i++) {

        cFace *cface = send[i];

        if (cface->IsVirtual()) {
            cface->SetPropertiesVirtualElement();
            virtuals[cface->RelProduct()].push_back(cface);
        } else {

            IntersectionRay &ray = rays[i];
            cFace *corr = nullptr;
            bool is_corresponding = false;
            double dist = -1;

            std::multimap<double, unsigned int> hits = intersector.Perform(ray);

            if (hits.empty()) { // there should be at least a self-hit. If not, it is often small faces ("parallel" criteria in intersects_triangle applies because of low values compared to tol)
                std::array<double, 3> PS = ray.get_start_point();
                std::array<double, 3> PE = ray.get_end_point();
                std::cout << (cface->SurfaceArea() < 1.0e-4 ? "[Info] Ray of small face hit nothing." : "[Info] Ray hit nothing.") << " "
                          << cface->Info() << "\t" << round_double_two_digits(PS[0]) << " "
                          << round_double_two_digits(PS[1]) << " " << round_double_two_digits(PS[2]) << " | "
                          << round_double_two_digits(PE[0]) << " " << round_double_two_digits(PE[1]) << " "
                          << round_double_two_digits(PE[2]) << " | " << cface->SurfaceArea() << "\n";
            }

            for (const auto &hit: hits) {

                cFace *hit_cface = M[hit.second];
                if (cface->ID() == hit_cface->ID()) continue; // skip self-hit

                // identify parent face by checking if first level face of HIT was parent of first level face of CURRENT
                if (cface->Superface()->Parent() != nullptr && cface->Superface()->Parent() == hit_cface->Superface()) {
                    cface->SetParent(hit_cface);
                    continue;
                }

                if (cface->Superface() == hit_cface->Superface()->Parent()) continue; // first level face of CURRENT is parent face of first level face of HIT (wall hitting own child opening face)
                if (cface->Superface()->Parent() == nullptr && hit_cface->Superface()->Parent() != nullptr) continue; // wall hitting an opening
                if (cface->Superface()->Parent() != nullptr && hit_cface->Superface()->Parent() == nullptr) continue; // opening hitting a wall

                // filter corresponding if not parallel
                double angle = cface->FixedFaceNormal().Angle(hit_cface->FixedFaceNormal());
                if (angle > min_angle)
                    is_corresponding = true;

                // because map is ordered, hit_cface is the "corresponding" face
                corr = hit_cface;
                dist = hit.first - tol; // correct distance by length origin was moved
                break;
            }

            cface->SetPropertiesSpaceBoundary(corr, is_corresponding, dist, bounds_min.Z(), cface->PointOnFace());
        }
    }

    // virtual faces
    for (auto &v: virtuals) {
        if (v.second.size() < 2) {
            std::cerr << "[Warning] No virtual face pair." << std::endl;
            continue;
        }
        if (v.second.size() > 2)
            std::cerr << "[Warning] Too many faces belonging to virtual element." << std::endl;

        v.second[0]->SetCorresponding(v.second[1]);
        v.second[1]->SetCorresponding(v.second[0]);
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Perform ray tracing", std::to_string(send.size()));
}

void Kernel::perform_ray_tracing_first_level(IntersectorInterface &intersector, const std::vector<cFace *> &send, std::vector<IntersectionRay> &rays, double tol, std::unordered_map<unsigned int, cFace *> &M, gp_XYZ bounds_min) {

    auto start = std::chrono::high_resolution_clock::now();

    std::unordered_map<Product *, std::vector<cFace * >> virtuals;

    for (unsigned int i = 0; i < send.size(); i++) {

        cFace *cface = send[i];

        if (cface->IsVirtual()) {
            cface->SetPropertiesVirtualElement();
            virtuals[cface->RelProduct()].push_back(cface);
        } else {

            IntersectionRay &ray = rays[i];
            cFace *corr = nullptr;
            double dist = -1;

            std::multimap<double, unsigned int> hits = intersector.Perform(ray);

            if (hits.empty())
                std::cout << "[Info] Ray hit nothing. " << cface->Info() << "\n";

            for (const auto &hit: hits) {

                cFace *hit_cface = M[hit.second];
                if (cface->ID() == hit_cface->ID()) continue; // skip self-hit
                if (cface == hit_cface->Parent()) continue; // skip own parent face
                if (cface->Parent() == nullptr && hit_cface->Parent() != nullptr) continue; // wall hitting an opening
                if (cface->Parent() != nullptr && hit_cface->Parent() == nullptr) continue; // opening hitting a wall

                // because map is ordered, hit_cface is the "corresponding" face
                corr = hit_cface;
                dist = hit.first - tol; // correct distance by length origin was moved
                break;
            }
            // no check for angle between corresponding faces
            cface->SetPropertiesSpaceBoundary(corr, true, dist, bounds_min.Z(), cface->PointOnFace());
        }
    }

    // virtual faces
    for (auto &v: virtuals) {
        if (v.second.size() < 2) {
            std::cerr << "[Warning] No virtual face pair." << std::endl;
            continue;
        }
        if (v.second.size() > 2)
            std::cerr << "[Warning] Too many faces belonging to virtual element." << std::endl;

        v.second[0]->SetCorresponding(v.second[1]);
        v.second[1]->SetCorresponding(v.second[0]);
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Perform ray tracing", std::to_string(send.size()));
}

void Kernel::remove_holes(std::list<cFace> &cFaces, unsigned int &fid) {

    auto start = std::chrono::high_resolution_clock::now();

    unsigned int c = 0;

    std::map<cFace *, std::vector<cFace * >> inner_boundaries;

    for (auto &cface: cFaces)
        if (cface.Parent() != nullptr)
            inner_boundaries[cface.Parent()].push_back(&cface);

    for (auto &cface: cFaces) {
        if (!cface.HasHoles() || cface.IsTrash()) continue;
        const auto &i = inner_boundaries[&cface];
        cface_linked_hole_removal(cFaces, cface, fid, i, c);
        if (!cface.IsTrash())
            cface_linked_naive_triangulation(cFaces, cface, fid, i, c);

        if (!cface.IsTrash()) {
            if (inner_boundaries.empty()) std::cerr << "[Warning] Hole removal failed for face " << cface.Info() << std::endl;
            // else ... implement inner boundaries
        }
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Remove holes", std::to_string(c) + "/" + std::to_string(cFaces.size()));
}

void Kernel::decompose_concave_polygons(std::list<cFace> &cFaces, unsigned int &fid) {

    auto start = std::chrono::high_resolution_clock::now();

    unsigned int c = 0;
    std::map<cFace *, std::vector<cFace * >> inner_boundaries;

    for (auto &cface: cFaces) {
        if (cface.Parent() != nullptr) inner_boundaries[cface.Parent()].push_back(&cface);
        cface.SetWasVisited(false); // prevent recursion. the decomposer library sometimes returns wrong results
    }

    for (auto &cface: cFaces) {
        if (cface.IsConvex() || cface.IsTrash() || cface.WasVisited()) continue;
        const auto &i = inner_boundaries[&cface];

        if (!cface.HasHoles()) {
            cface_linked_concave_decomposing(cFaces, cface, fid, i, c);
            if (!cface.IsTrash())
                cface_linked_naive_triangulation(cFaces, cface, fid, i, c);
        } else
            cface_linked_naive_triangulation(cFaces, cface, fid, i, c);

        if (!cface.IsTrash()) {
            if (inner_boundaries.empty()) std::cerr << "[Warning] Convex decomposition failed for face " << cface.Info() << std::endl;
            // else ... implement inner boundaries
        }
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Decompose concave faces", std::to_string(c) + "/" + std::to_string(cFaces.size()));
}

void Kernel::decompose_concave_polygons_triangulation(std::list<cFace> &cFaces, unsigned int &fid) {

    auto start = std::chrono::high_resolution_clock::now();

    unsigned int c = 0;
    std::map<cFace *, std::vector<cFace * >> inner_boundaries;

    for (auto &cface: cFaces) {
        if (cface.Parent() != nullptr) inner_boundaries[cface.Parent()].push_back(&cface);
        cface.SetWasVisited(false); // prevent recursion
    }

    for (auto &cface: cFaces) {
        if (cface.IsConvex() || cface.IsTrash() || cface.WasVisited()) continue;
        const auto &i = inner_boundaries[&cface];
        cface_linked_naive_triangulation(cFaces, cface, fid, i, c);
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Decompose concave faces by triangulation", std::to_string(c) + "/" + std::to_string(cFaces.size()));
}

void Kernel::cface_linked_hole_removal(std::list<cFace> &cFaces, cFace &cface, unsigned int &fid, const std::vector<cFace *> &inner_boundaries, unsigned int &c) {

    if (!inner_boundaries.empty()) {
        return; // Veronika doesnt want a split in BIM2SIm tool chain
        // std::cerr << "[Warning] Face has children. Paying attention to this is not implemented in this function, yet. But changes are easy!" << std::endl;
    }

    c++;
    auto a = cface.ProjectionAxis();
    TopoDS_Wire outerWire = BRepTools::OuterWire(cface.face);

    ifc2sb_triangle_hole_remove::PointChainWithSpot wire;
    ifc2sb_triangle_hole_remove::PointChainsWithSpots holes;
    ifc2sb_triangle_hole_remove::PointChainsWithSpots children; // inner boundaries of parent face. At their position parent face shall not be split

    /*
    Triangle will segfault or at least fail, when there is a duplicate vertex (hole having vertex on outer wire). In fact, it seems that Triangle sometimes cannot work with more than two segments meeting at one vertex:
        Bailing out after 159 iterations in finddirection().
        Please report this bug to jrs@cs.berkeley.edu
        Include the message above, your input data set, and the exact
        command line you used to run Triangle.
     So it didn't work to remove holepoint and make hole segments to a normal region. (this would have been the solution. ignore "hole regions" when getting triangles and nbs)
    */
    TopoDS_ListOfShape vertices = Topo(cface.face).vertices();
    std::set<unsigned int> vertices2;
    for (auto &v: vertices)
        vertices2.insert(hash(v));
    if (vertices.Size() > 2 * vertices2.size()) {
        std::cerr << "[Warning] Face " << cface.Info() << " contains duplicate vertices. Triangulation may fail or even segfault!" << std::endl;
        return;
    }

    for (auto &W: Topo(cface.face).wires()) {

        ifc2sb_triangle_hole_remove::PointChainWithSpot chain;

        if (a == 0)
            for (const auto &vertex: Topo(W).ordered_vertices_of_wire()) {
                gp_Pnt P = BRep_Tool::Pnt(TopoDS::Vertex(vertex));
                chain.first.push_back({P.Y(), P.Z()});
            }
        else if (a == 1)
            for (const auto &vertex: Topo(W).ordered_vertices_of_wire()) {
                gp_Pnt P = BRep_Tool::Pnt(TopoDS::Vertex(vertex));
                chain.first.push_back({P.X(), P.Z()});
            }
        else
            for (const auto &vertex: Topo(W).ordered_vertices_of_wire()) {
                gp_Pnt P = BRep_Tool::Pnt(TopoDS::Vertex(vertex));
                chain.first.push_back({P.X(), P.Y()});
            }

        if (chain.first.size() < 2) {
            std::cerr << "[Warning] Wire has less than three edges. Skip!" << std::endl;
            continue;
        }

        // calculate point near random edge within polygon
        // ... center point
        ifc2sb_triangle_hole_remove::Point P1 = chain.first.back();
        ifc2sb_triangle_hole_remove::Point P2 = chain.first.front();
        ifc2sb_triangle_hole_remove::Point M = {0.5 * (P1[0] + P2[0]), 0.5 * (P1[1] + P2[1])};

        // ... vector
        gp_Dir2d n(P2[0] - P1[0], P2[1] - P1[1]);
        bool counter_clockwise_orientation = ifc2sb_triangle_hole_remove::orientation(chain.first); // if ccw, then rotate vector to the left. else to the right
        gp_Vec2d v = counter_clockwise_orientation ? gp_Vec2d(-n.Y(), n.X()) : gp_Vec2d(n.Y(), -n.X());
        v.Scale(1.0e-7);

        // ... point
        chain.second = {M[0] + v.X(), M[1] + v.Y()};

        if (W == outerWire)
            wire = chain;
        else
            holes.push_back(chain);
    }

    // TODO add children by inners map
    // for ... inners ... do the same as for all other wires
    // ...

    if (wire.first.empty()) {
        std::cerr << "[Warning] No outer wire. Hole removal failed!" << std::endl;
        return;
    }

    auto result = ifc2sb_triangle_hole_remove::run(wire, holes, children);

    if (result.size() < 2) {
        std::cerr << "[Warning] Hole removal failed!" << std::endl;
        return;
    }

    std::set<cFace *> NEW;

    for (const auto &r: result) {

        BRepBuilderAPI_MakePolygon m;
        gp_Pnt R = BRep_Tool::Pnt(TopoDS::Vertex(Topo(outerWire).vertices().First()));
        // TODO remove collinear points if quality mesh 'q' is used for triangulation creating new points

        for (const auto &P_2D: r)
            m.Add(project_inv(P_2D[0], P_2D[1], a, cface.FixedFaceNormal(), R));

        TopoDS_Face F = BRepBuilderAPI_MakeFace(m.Wire()).Face();

        // new sb
        cFace new_cface(cface);
        new_cface.SetID(fid);
        fid++;
        // set face
        if (face_normal(F).Angle(cface.FaceNormal()) > 1.57) F.Complement(); // could be removed. CheckSetFaceNormal takes care of this now
        new_cface.face = F;
        new_cface.CheckSetFaceNormal();
        // add to list and space
        cFaces.push_back(new_cface);
        cface.RelSpace()->AddSecondLvlFace(&cFaces.back());
        cface.Superface()->AppendToSubfaces(&cFaces.back());
        NEW.insert(&cFaces.back());

        // TODO analyze result and the regions (position in inners vector) and link cface (attribute parent)
        // inners[xxx]->SetParent(&cFaces.back());
    }

    cface.SetIsTrash(true);
    cface.RemoveFromSubfaces(&cface);

    // Corresponding
    split_corresponding_equivalent(cface, NEW, fid, cFaces);  // TODO analyze result and the regions (position in inners vector) and link cface (attribute parent)
}

void Kernel::cface_linked_concave_decomposing(std::list<cFace> &cFaces, cFace &cface, unsigned int &fid, const std::vector<cFace *> &inner_boundaries, unsigned int &c) {

    cface.SetWasVisited(true);

    if (!inner_boundaries.empty()) {
        return;
        // std::cerr << "[Warning] Face has children. Paying attention to this is not implemented in this function, yet. But changes are easy!" << std::endl;
    }

    c++;

    std::vector<gp_Pnt2d> L = cface.ProjectOuterWireTo2D();
    std::list<std::vector<gp_Pnt2d>> result;
    std::vector<cxdsb::Vertex> verts;
    verts.reserve(L.size());

    for (const auto &P: L)
        verts.emplace_back(cxdsb::Vec2({P.X(), P.Y()}));

    cxdsb::ConcavePolygon concavePoly(verts);
    unsigned int count = 0;
    concavePoly.convexDecomp(count); // the library can't handle (nearly) collinear and coincident (pointwise isct) situations

    // access result
    std::vector<cxdsb::ConcavePolygon> subPolygonList;
    concavePoly.returnLowestLevelPolys(subPolygonList);

    if (subPolygonList.size() < 2) {
        std::cerr << "[Warning] No solution for face " << subPolygonList.size() << "\t" << cface.Info() << std::endl;
        return;
    }

    bool fail = false;

    for (auto &subPolygon: subPolygonList)
        if (!subpolygon_check(subPolygon, result, L, cface)) {
            fail = true;
            break;
        }

    if (fail) {
        std::cerr << "[Warning] No solution for face because of an invalid polygon " << subPolygonList.size() << "\t" << cface.Info() << std::endl;
        return;
    }

    std::set<cFace *> NEW;
    auto a = cface.ProjectionAxis();
    TopoDS_Wire outerWire = BRepTools::OuterWire(cface.face);

    for (const auto &r: result) {
        BRepBuilderAPI_MakePolygon m;
        gp_Pnt R = BRep_Tool::Pnt(TopoDS::Vertex(Topo(outerWire).vertices().First()));

        for (const auto &P_2D: r)
            m.Add(project_inv(P_2D.X(), P_2D.Y(), a, cface.FixedFaceNormal(), R));
        m.Close();
        TopoDS_Face F = BRepBuilderAPI_MakeFace(m.Wire()).Face();

        // new sb
        cFace new_cface(cface);
        new_cface.SetID(fid);
        fid++;
        // set face
        if (face_normal(F).Angle(cface.FaceNormal()) > 1.57) F.Complement(); // could be removed. CheckSetFaceNormal takes care of this now
        new_cface.face = F;
        new_cface.CheckSetFaceNormal();
        // add to list and space
        cFaces.push_back(new_cface);
        cface.RelSpace()->AddSecondLvlFace(&cFaces.back());
        cface.Superface()->AppendToSubfaces(&cFaces.back());
        NEW.insert(&cFaces.back());

        // TODO analyze result and the regions (position in inners vector) and link cface (attribute parent)
        // inners[xxx]->SetParent(&cFaces.back());
    }

    cface.SetIsTrash(true);
    cface.RemoveFromSubfaces(&cface);

    // Corresponding
    split_corresponding_equivalent(cface, NEW, fid, cFaces);  // TODO analyze result and the regions (position in inners vector) and link cface (attribute parent)
}

void Kernel::cface_linked_naive_triangulation(std::list<cFace> &cFaces, cFace &cface, unsigned int &fid, const std::vector<cFace *> &inner_boundaries, unsigned int &c) {

    cface.SetWasVisited(true);

    if (!inner_boundaries.empty()) {
        return;
        // std::cerr << "[Warning] Face has children. Paying attention to this is not implemented in this function, yet. But changes are easy!" << std::endl;
    }

    c++;

    // triangulate
    TopLoc_Location Loc;
    Handle(Poly_Triangulation) T = BRep_Tool::Triangulation(cface.face, Loc);
    BRepMesh_IncrementalMesh(cface.face, 0.01, false, 0.01, true);
    std::list<std::array<gp_Pnt, 3>> triangles = triangles_from_face(cface.face);
    if (triangles.empty()) return;

    std::set<cFace *> NEW;

    for (const auto &triangle: triangles) {

        TopoDS_Wire wire = BRepBuilderAPI_MakePolygon(triangle[0], triangle[1], triangle[2], true).Wire();
        TopoDS_Face F = BRepBuilderAPI_MakeFace(wire).Face();

        // new sb
        cFace new_cface(cface);
        new_cface.SetID(fid);
        fid++;
        // set face
        if (face_normal(F).Angle(cface.FaceNormal()) > 1.57) F.Complement(); // could be removed. CheckSetFaceNormal takes care of this now
        new_cface.face = F;
        new_cface.CheckSetFaceNormal();
        // add to list and space
        cFaces.push_back(new_cface);
        cface.RelSpace()->AddSecondLvlFace(&cFaces.back());
        cface.Superface()->AppendToSubfaces(&cFaces.back());
        NEW.insert(&cFaces.back());

        // TODO analyze result and the regions (position in inners vector) and link cface (attribute parent)
        // inners[xxx]->SetParent(&cFaces.back());
    }

    cface.SetIsTrash(true);
    cface.RemoveFromSubfaces(&cface);

    // Corresponding
    split_corresponding_equivalent(cface, NEW, fid, cFaces);  // TODO analyze result and the regions (position in inners vector) and link cface (attribute parent)
}

void Kernel::split_corresponding_equivalent(const cFace &cface, const std::set<cFace *> &NEW, unsigned int &fid, std::list<cFace> &cFaces) {

    cFace *corresponding = cface.Corresponding();

    if (corresponding == nullptr) return;

    gp_Vec v(corresponding->FixedFaceNormal());
    v.Scale(50);

    bool new_faces_were_created = false;

    for (const auto &N: NEW) {

        TopoDS_Shape prism = BRepPrimAPI_MakePrism(N->face, v).Shape();
        BRepAlgoAPI_Common common_tool;

        TopoDS_ListOfShape tools, argms;
        argms.Append(corresponding->face);
        tools.Append(prism);
        common_tool.SetArguments(argms);
        common_tool.SetTools(tools);
        common_tool.SetRunParallel(true);
        common_tool.SetFuzzyValue(1e-5);
        common_tool.Build();
        TopoDS_Shape common = common_tool.Shape();
        // TopoDS_Shape common = BRepAlgoAPI_Common(corresponding->face, prism).Shape();

        if (common.IsNull()) {
            std::cerr << "[Warning] Equivalent splitting of corresponding didn't create a shape. This will create a corresponding mismatch and maybe a missing face" << cface.Info() << " " << N->Info() << " " << corresponding->Info() << ". Remove corresponding info." << std::endl;
            N->SetCorresponding(nullptr);
            N->SetSBType(SB_TYPE_2B);
            continue;
        }

        Topo top(common);
        TopoDS_ListOfShape CL = top.faces();

        if (CL.IsEmpty()) { // sometime occ gives a null result (containing only a compound) for perfectly good situations, use clipper instead

            std::cout << "[Info] Equivalent splitting of corresponding using occ didn't create a face (" << top.edges().Size() << " edges and " << top.vertices().Size() << " vertices). " << cface.Info() << " | " << N->Info() << " | " << corresponding->Info() << "\n";

            Clipper clip(corresponding->face, corresponding->FixedFaceNormal(), true);
            clip.clip(N->face, N->FixedFaceNormal());
            TopoDS_ListOfShape Clips;
            TopoDS_Face Fc;

            if (clip.success)
                for (auto &t: clip.Result)
                    for (auto &t2: Topo(BRepAlgoAPI_Common(corresponding->face, t).Shape()).faces())
                        Clips.Append(t2);

            if (Clips.Size() == 1) {
                auto faces = Topo(Clips.First()).faces();
                if (faces.Size() == 1 && !faces.First().IsNull())
                    CL.Append(faces.First());
            }

            if (CL.IsEmpty()) {
                std::cerr << "[Warning] Equivalent splitting of corresponding didn't create a face (" << top.edges().Size() << " edges and " << top.vertices().Size() << " vertices). This will create a corresponding mismatch and maybe a missing face. " << cface.Info() << " " << N->Info()
                          << " " << corresponding->Info() << ". Remove corresponding info." << std::endl;
                N->SetCorresponding(nullptr);
                N->SetSBType(SB_TYPE_2B);
                continue;
            } else
                CL.Append(Fc);
        }
        TopoDS_Face F = TopoDS::Face(CL.First());

        // new sb
        cFace new_cface(*corresponding);
        new_cface.SetID(fid);
        fid++;
        new_faces_were_created = true;

        // set face
        new_cface.face = F;
        new_cface.CheckSetFaceNormal();

        // add to list and space
        cFaces.push_back(new_cface);
        corresponding->RelSpace()->AddSecondLvlFace(&cFaces.back());
        corresponding->Superface()->AppendToSubfaces(&cFaces.back());

        // set corresponding
        cFaces.back().SetCorresponding(N);
        N->SetCorresponding(&cFaces.back());

        // TODO analyze result and the regions (position in inners vector) and link cface (attribute parent)
        // inner_boundaries[corresponding)][xxx]->SetParent(&cFaces.back());
    }

    if (new_faces_were_created) {
        corresponding->SetIsTrash(true);
        corresponding->RemoveFromSubfaces(corresponding);
    } else {
        corresponding->SetWasVisited(true);
        corresponding->SetCorresponding(*NEW.begin());
        (*NEW.begin())->SetCorresponding(corresponding);
        std::cerr << "[Warning] The corresponding face " << corresponding->Info() << " could not be split at all. Faces " << (*NEW.begin())->Info() << " and " << corresponding->Info() << " will be linked visa versa to prevent memory error." << std::endl;
    }
}

bool Kernel::subpolygon_check(const cxdsb::ConcavePolygon &subPolygon, std::list<std::vector<gp_Pnt2d>> &result, const std::vector<gp_Pnt2d> &L, const cFace &cface) {

    std::vector<gp_Pnt2d> chain;
    for (const auto &P: subPolygon.getVertices())
        chain.emplace_back(P.position.x, P.position.y);

    if (!polygon_is_convex(chain)) {
        std::cerr << "[Warning] Result for " << cface.Info() << " is not convex! Skip decomposition" << std::endl;
        return false;
    }

    if (chain.size() < 3) {
        std::cerr << "[Warning] New convex wire of " << cface.Info() << " consists of less than three points (" << chain.size() << ")! Skip decomposition" << std::endl;
        return false;
    }

    chain = remove_coincident_points_from_polygon(chain, 1.0e-5);

    if (chain.size() < 3) {
        std::cerr << "[Warning] New convex wire of " << cface.Info() << " consists of less than three points after removing coincidence (" << chain.size() << ")! Skip decomposition" << std::endl;
        return false;
    }

    chain = remove_colinear_points_from_polygon(chain, 1.0e-5);

    if (chain.size() < 3) {
        std::cerr << "[Warning] New convex wire of " << cface.Info() << " consists of less than three points after removing collinearity (" << chain.size() << ")! Skip decomposition" << std::endl;
        return false;
    }

    if (!is_point_inside_2D_polygon(L, centroid_2D_polygon(chain))) { // sometimes result of the decomposition library does not lie on the original concavee face
        std::cerr << "[Warning] Generated concave face does not lie on original face " << cface.Info() << "! Skip decomposition" << std::endl;
        return false;
    }

    result.push_back(chain);
    return true;
}

void Kernel::simplify_fenestration_faces(std::list<cFace> &cFaces, unsigned int &fid) {

    auto start = std::chrono::high_resolution_clock::now();

    int c = 0;

    std::map<cFace *, std::vector<cFace * >> inner_boundaries;

    for (auto &cface: cFaces) {
        if (cface.Parent() != nullptr) inner_boundaries[cface.Parent()].push_back(&cface);
        cface.SetWasVisited(false); // prevent recursion
    }

    for (auto &cface: cFaces) {

        if (cface.HasHoles() || !cface.IsOpening() || cface.WasVisited()) continue;

        cface.SetWasVisited(true);
        c++;

        // vertex to points
        TopoDS_ListOfShape vertices = Topo(cface.OuterWire()).ordered_vertices_of_wire();

        if (vertices.Size() < 3) {
            std::cerr << "[Warning] Outer wire consists of less than three vertices!";
            continue;
        }

        std::vector<gp_Pnt> Pnts = vertex_to_point_list(vertices);

        if (cface.Corresponding() != nullptr) Pnts = sort_points_from_polygon_by_dist(Pnts, &cface); // sort so point removal will be the same for corresponding boundaries

        double tol = 1.0e-4;
        Pnts = remove_coincident_points_from_polygon(Pnts, tol);

        if (Pnts.size() < 3) {
            std::cerr << "[Warning] Outer wire of " << cface.Info() << " consists of less than three points after removing coincident points (" << Pnts.size() << ")! Skip space boundary" << std::endl;
            continue;
        }

        Pnts = remove_colinear_points_from_polygon(Pnts, tol);

        if (Pnts.size() < 3) {
            std::cerr << "[Warning] Outer wire of " << cface.Info() << " consists of less than three points after removing collinear points (" << Pnts.size() << ")! Skip space boundary" << std::endl;
            continue;
        }

        // face creation
        BRepBuilderAPI_MakePolygon M;
        for (const auto &P: Pnts) M.Add(P);
        M.Close();
        TopoDS_Face F = BRepBuilderAPI_MakeFace(M.Wire());

        // set
        cface.face = F;
        cface.CheckSetFaceNormal();

        if (Topo(cface.OuterWire()).edges().Size() > 16) {

            TopoDS_ListOfShape vertices = Topo(cface.OuterWire()).ordered_vertices_of_wire();

            if (vertices.Size() < 3) {
                std::cerr << "[Warning] Outer wire consists of less than three vertices!";
                continue;
            }

            // define plane
            std::vector<gp_Pnt> Pnts = vertex_to_point_list(vertices);
            gp_Pnt O = Pnts[0]; // Origin of plane
            gp_Pnt T = Pnts[1];
            gp_Dir n = cface.FixedFaceNormal();
            gp_Dir u_plane, v_plane;
            calc_plane_uv_vectors(O, T, n, u_plane, v_plane);

            // 2D of center
            double uC, vC;
            calc_point_uv_parameters_on_plane(cface.Center(), O, u_plane, v_plane, uC, vC);
            gp_Pnt2d C2D(uC, vC);

            // transform to 2D
            std::vector<gp_Pnt2d> Pnts2D;
            for (const auto &P: Pnts) {
                double u, v;
                calc_point_uv_parameters_on_plane(P, O, u_plane, v_plane, u, v);
                Pnts2D.emplace_back(u, v);
            }

            // create 2D rectangle
            double s = sqrt(cface.SurfaceArea());
            double h = 0.5 * s;
            std::vector<gp_Pnt2d> Rect2D;
            Rect2D.emplace_back(uC - h, vC - h);
            Rect2D.emplace_back(uC + h, vC - h);
            Rect2D.emplace_back(uC + h, vC + h);
            Rect2D.emplace_back(uC - h, vC + h);

            // transform rectangle to 3D
            std::vector<gp_Pnt> Rect;
            for (const auto &P2D: Rect2D)
                Rect.emplace_back((gp_Vec(O.XYZ()) + P2D.X() * u_plane + P2D.Y() * v_plane).XYZ());

            // face creation
            BRepBuilderAPI_MakePolygon M;
            for (const auto &P: Rect) M.Add(P);
            M.Close();
            cface.face = BRepBuilderAPI_MakeFace(M.Wire());
            cface.CheckSetFaceNormal();

            // rotate rectangle to be perpendicular to global coordinate system
            auto axis = cface.ProjectionAxis();
            gp_Vec r(gp_Vec(Rect[1].XYZ()) - gp_Vec(Rect[0].XYZ()));
            double angle = axis == 0 ? r.Angle(gp_Dir(0, 1, 0)) : r.Angle(gp_Dir(1, 0, 0));
            rotate(cface.face, cface.Center(), n, angle);
            std::cout << "Fenestration face " << cface.Info() << "was simplified to rectangle\n";

        } else if (Topo(cface.OuterWire()).edges().Size() > 4) {

            auto &inners = inner_boundaries[&cface];
            if (!inners.empty()) continue;

            // triangulate
            TopLoc_Location Loc;
            Handle(Poly_Triangulation) T = BRep_Tool::Triangulation(cface.face, Loc);
            BRepMesh_IncrementalMesh(cface.face, 0.01, false, 0.01, true);
            std::list<std::array<gp_Pnt, 3>> triangles = triangles_from_face(cface.face);
            if (triangles.empty()) continue;

            std::set<cFace *> NEW;
            std::cout << "Fenestration face " << cface.Info() << "was simplified to triangulation\n";

            for (const auto &triangle: triangles) {

                TopoDS_Wire wire = BRepBuilderAPI_MakePolygon(triangle[0], triangle[1], triangle[2], true).Wire();
                TopoDS_Face F = BRepBuilderAPI_MakeFace(wire).Face();

                // new sb
                cFace new_cface(cface);
                new_cface.SetID(fid);
                fid++;
                // set face
                new_cface.face = F;
                new_cface.CheckSetFaceNormal();
                // add to list and space
                cFaces.push_back(new_cface);
                cface.RelSpace()->AddSecondLvlFace(&cFaces.back());
                cface.Superface()->AppendToSubfaces(&cFaces.back());
                NEW.insert(&cFaces.back());
            }

            cface.SetIsTrash(true);
            cface.RemoveFromSubfaces(&cface);

            // Corresponding
            split_corresponding_equivalent(cface, NEW, fid, cFaces);
        }
/*
        TopoDS_ListOfShape vertices = Topo(cface.OuterWire()).ordered_vertices_of_wire();

        if (vertices.Size() < 3) {
            std::cerr << "[Warning] Outer wire consists of less than three vertices!";
            continue;
        }

        // define plane
        std::vector<gp_Pnt> Pnts = vertex_to_point_list(vertices);
        gp_Pnt O = Pnts[0]; // Origin of plane
        gp_Pnt T = Pnts[1];
        gp_Dir n = cface.FixedFaceNormal();
        gp_Dir u_plane, v_plane;
        calc_plane_uv_vectors(O, T, n, u_plane, v_plane);

        std::cerr << "n: " << n.X() << " " << n.Y() << " " << n.Z() << std::endl;
        std::cerr << "u_plane: " << u_plane.X() << " " << u_plane.Y() << " " << u_plane.Z() << std::endl;
        std::cerr << "v_plane: " << v_plane.X() << " " << v_plane.Y() << " " << v_plane.Z() << std::endl;

        // 2D of center
        double uC, vC;
        calc_point_uv_parameters_on_plane(cface.Center(), O, u_plane, v_plane, uC, vC);
        gp_Pnt2d C2D(uC, vC);

        for (const auto &P: Pnts)
            std::cerr << "3Da: " << P.X() << " " << P.Y() << " " << P.Z() << std::endl;

        // transform to 2D
        std::vector<gp_Pnt2d> Pnts2D;
        std::vector<double> U;
        U.reserve(Pnts2D.size());

        for (const auto &P: Pnts) {
            double u, v;
            calc_point_uv_parameters_on_plane(P, O, u_plane, v_plane, u, v);
            Pnts2D.emplace_back(u, v);
            U.push_back(u);
        }

        // create 2D rectangle
        const auto lim = std::minmax_element(begin(U), end(U));
        double s1 = *lim.second - *lim.first;
        double s2 = cface.SurfaceArea() / s1;
        double h2 = 0.5 * s2;
        double h1 = 0.5 * s1;
        std::vector<gp_Pnt2d> Rect2D;
        Rect2D.emplace_back(uC - h2, vC - h1);
        Rect2D.emplace_back(uC + h2, vC - h1);
        Rect2D.emplace_back(uC + h2, vC + h1);
        Rect2D.emplace_back(uC - h2, vC + h1);

        // transform rectangle to 3D
        std::vector<gp_Pnt> Rect;
        for (const auto &P2D: Rect2D) {
            gp_Vec t = gp_Vec(O.XYZ()) + P2D.X() * u_plane + P2D.Y() * v_plane; // back to 3D
            gp_Pnt P(t.XYZ());
            Rect.emplace_back(P);
            std::cerr << "Rect1: " << gp_Pnt(t.XYZ()).X() << " " << gp_Pnt(t.XYZ()).Y() << " " << gp_Pnt(t.XYZ()).Z() << std::endl;
        }

        // face creation
        BRepBuilderAPI_MakePolygon M;
        for (const auto &P: Rect) M.Add(P);
        M.Close();
        TopoDS_Face F = BRepBuilderAPI_MakeFace(M.Wire());

        // rotate rectangle to be perpendicular to global coordinate system
        // double d = Rect[0].Distance(Rect[1]);
        // double angle = fabs(d - h) < 1.0e-7 ? gp_Vec(gp_Vec(Rect[2].XYZ()) - gp_Vec(Rect[1].XYZ())).Angle(n) : gp_Vec(gp_Vec(Rect[1].XYZ()) - gp_Vec(Rect[0].XYZ())).Angle(n);
        // std::cerr << "angle: " << angle * 180 / PI << std::endl;
        //  rotate(F, cface.Center(), n, -angle);

        // set
        cface.face = F;
        if (!cface.CheckFaceNormal()) cface.ComplementFace();
        */
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Simplify fenestration faces", std::to_string(c) + "/" + std::to_string(cFaces.size()));
}

void Kernel::calculate_shading(std::list<cFace> &cFaces, std::list<Space> &spaces, std::unordered_map<std::string, std::list<TopoDS_Face>> &shadings, std::unique_ptr<IfcParse::IfcFile> &model, IfcGeom::IteratorSettings settings) {

    auto start = std::chrono::high_resolution_clock::now();

    // todo std::list<TopoDS_Face> shadings_devices = calculate_shading_devices(cFaces, spaces));  ... trivial (read ifcshadingdevices from ifc model. generate shapes. maybe simplify. move to nearest ifcspace)
    std::list<TopoDS_Face> shadings_extern = calculate_shading_extern(spaces, model, settings);
    std::list<TopoDS_Face> shadings_facade = calculate_shading_facade(cFaces, spaces);

    shadings["Shading:Site"].insert(shadings["Shading:Site"].end(), shadings_extern.begin(), shadings_extern.end());
    shadings["Shading:Building"].insert(shadings["Shading:Building"].end(), shadings_facade.begin(), shadings_facade.end());

/*    std::list<viewerHelper::DisplayShapes> ds;

*//*    for (const auto &space: spaces) {
        if (!space.is_facade) continue;
        ds.emplace_back();
        ds.back().shape = space.shell;
        ds.back().transparency = 0;
        break;
    }*//*

//    for (const auto &space: spaces) {
//        if (space.is_facade) continue;
//        ds.emplace_back();
//        ds.back().shape = space.shell;
//        ds.back().transparency = 0.0;
//    }

    for (const auto &cface: cFaces) {
        if (cface.RelSpace()->is_facade) continue;
        ds.emplace_back();
        ds.back().shape = cface.face;
        ds.back().clr_by_string = true;
        ds.back().transparency = 0.0;
        if(cface.IsOpening()) {
            if(cface.Ancestor()->Opening()->declaration().is("IfcDoor")) {
                ds.back().clr_string = "BROWN";
            }
            else{
                ds.back().clr_string = "BLUE";
                ds.back().transparency = 0.8;
            }
        }else if(cface.InternalOrExternal() == SB_IE_EXTERNAL && cface.FixedFaceNormal().Z()>0.7) {
            ds.back().clr_string = "OPENSTUDIORED";
            ds.back().transparency = 0.8;
        }
        else
        ds.back().clr_string = "EGG";
    }

    for (const auto &cface: shadings) {
        ds.emplace_back();
        ds.back().shape = cface;
        ds.back().clr_by_string = true;
        ds.back().transparency = 0.0;
        ds.back().clr_string = "PURPLE";
    }
    ViewerMain::start_viewer(ds);*/

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Calculate shading faces", std::to_string(shadings.size()) + "/" + std::to_string(cFaces.size()));
}

std::list<TopoDS_Face> Kernel::calculate_shading_extern(std::list<Space> &spaces, std::unique_ptr<IfcParse::IfcFile> &model, IfcGeom::IteratorSettings settings) const {

    std::list<TopoDS_Face> shadings;

    /*
    Space *facade = nullptr;
    for (auto &space: spaces)
        if (space.is_facade) {
            facade = &space;
            break;
        }

    if (facade == nullptr)
        return shadings;
    */

    std::list<Product> products;
    std::set<std::string> guids;
    gp_XYZ bounds_min;
    gp_XYZ bounds_max;

    boost::shared_ptr<IfcEntityList> IfcEntityList = nullptr;
    try { IfcEntityList = model->instances_by_type("IfcBuildingElementProxy"); }
    catch (...) {}

    if (IfcEntityList != nullptr)
        for (auto E: *IfcEntityList) { // for the entities, check if they are child of IfcElement class, only those can have or be an void filling product

            std::string name = E->data().getArgument(2)->toString();
            std::string guid = E->data().getArgument(0)->toString();
            guid = guid.substr(1, guid.size() - 2); // remove quote signs

            if (name.find("Baum") != std::string::npos || name.find("baum") != std::string::npos)
                guids.insert(guid);
        }

    if (guids.empty())
        return shadings;

    if (!generate_shapes_from_ifc_guids(model, settings, products, guids, bounds_min, bounds_max))
        return shadings;

    if (products.empty())
        return shadings;

    for (auto &product: products) {
        std::cout << "[Info] Add product " << product.guid << " as shading.\n";
        for (auto &face: Topo(product.shape).faces())
            if (fabs(face_normal(TopoDS::Face(face)).Z()) < 0.5)
                shadings.push_back(TopoDS::Face(face));
    }

    return shadings;
}

std::list<TopoDS_Face> Kernel::calculate_shading_facade(std::list<cFace> &cFaces, std::list<Space> &spaces) {

    std::list<TopoDS_Face> shadings;

    // get facade
    Space *facade = nullptr;
    for (auto &space: spaces)
        if (space.is_facade) {
            facade = &space;
            break;
        }
    if (facade == nullptr) return shadings;

    // create some lists and maps
    std::set<cFace *> F;
    std::set<Product *> prod; // contains products that provide space boundaries related to ifcspace
    TopoDS_ListOfShape Lfacade;
    NCollection_DataMap<TopoDS_Shape, cFace *> N;

    for (auto &cface: cFaces)
        if (!cface.RelSpace()->is_facade)
            prod.insert(cface.RelProduct());
        else {
            Lfacade.Append(cface.face);
            N.Bind(Lfacade.Last(), &cface);
        }

    // prepare intersector
    TopoDS_Shape Comp = compound_from_shape_list(Lfacade);
    IntCurvesFace_ShapeIntersector intersector;
    intersector.Load(Comp, 1.0e-5);

    // find shadings
    for (auto &cface: cFaces)
        if (cface.IsPotentiallyShading(prod))
            F.insert(&cface);

    // cluster to shells and store also free edges
    TopoDS_Shape shells = sew_cfaces(F, 1.0e-4, false);
    std::list<std::pair<TopoDS_Shape, TopoDS_ListOfShape >> free;
    for (auto &shell: Topo(shells).shells())
        free.emplace_back(shell, ShapeChecker(shell).free_edges());

    /*  single faces
    TopTools_IndexedDataMapOfShapeListOfShape shellMap;
    TopExp::MapShapesAndAncestors(shells, TopAbs_FACE, TopAbs_SHELL, shellMap);
    for (auto &f: Topo(shells).faces()) // single faces having no shell
        if (shellMap.FindFromKey(f).IsEmpty()) free.emplace_back(f, ShapeChecker(f).free_edges());
    */

    // find planar open wires formed by shading faces
    for (auto &it: free) {

        TopoDS_ListOfShape used, planes;

        // find planes constructed from the open wires of the shell
        for (auto &edge: it.second) {

            TopoDS_ListOfShape vertices, C;
            vertices.Append(TopExp::FirstVertex(TopoDS::Edge(edge), Standard_True));
            C.Append(edge);

            while (true) {

                auto last_vertex = vertices.Last();
                bool found = false;

                for (auto &edge2: it.second) {
                    if (C.Contains(edge2) or used.Contains(edge2)) continue;
                    auto P1 = TopExp::FirstVertex(TopoDS::Edge(edge2));
                    auto P2 = TopExp::LastVertex(TopoDS::Edge(edge2));

                    if (last_vertex.IsSame(P1)) {
                        found = true;
                        vertices.Append(P2);
                        C.Append(edge2);
                        break;
                    }
                    if (last_vertex.IsSame(P2)) {
                        found = true;
                        vertices.Append(P1);
                        C.Append(edge2);
                        break;
                    }
                }
                if (!found) break;

            }

            if (C.Size() == 1) continue;

            BRepBuilderAPI_MakeWire mW;
            for (auto &c: C)
                mW.Add(TopoDS::Edge(c));
            mW.Build();
            if (!mW.IsDone()) continue;

            TopoDS_Wire Wire = mW.Wire();
            if (!Wire.Closed()) continue;

            for (auto &c: C)
                used.Append(c);

            BRepBuilderAPI_MakePolygon m;
            for (auto &v: vertices)
                m.Add(TopoDS::Vertex(v));
            m.Close();

            if (!m.IsDone()) continue;

            TopoDS_Face temp_plane = BRepBuilderAPI_MakeFace(m.Wire()).Face();
            if (temp_plane.IsNull()) continue;

            planes.Append(temp_plane);
        }

        if (planes.IsEmpty()) continue;

        // move to ifcspace
        for (auto &plane: planes) {

            // distance to facade must be zero
            gp_Pnt center = face_center(TopoDS::Face(plane));
            double dist = minimal_distance_between_shapes(facade->shell, BRepBuilderAPI_MakeVertex(center).Vertex());
            if (dist > 1.0e-9) continue;

            // get the nearest face on facade
            gp_Vec normal = gp_Vec(face_normal(TopoDS::Face(plane)));
            gp_Lin Line1(center, normal);
            intersector.PerformNearest(Line1, -1.0e-5, 1.0e-5);
            if (intersector.NbPnt() == 0) continue;

            const auto &hit_face = intersector.Face(1);
            cFace *hit_cface = N.Find(hit_face);

            // TODO skip if hit_cface is no opening?

            // get 2a corresponding and done
            if (hit_cface->Corresponding() == nullptr) continue;
            if (hit_cface->Corresponding()->RelSpace()->is_facade) continue;

            double dist2 = minimal_distance_between_shapes(hit_cface->face, hit_cface->Corresponding()->face);
            if (dist2 > 1.0) continue;

            gp_Vec mover(hit_cface->FixedFaceNormal().Reversed());
            mover.Scale(dist2);
            move(it.first, mover);
            for (auto &shading_face: Topo(it.first).faces())
                shadings.push_back(TopoDS::Face(shading_face));
            break;
        }
    }

    return shadings;
}

template<typename SetType>
bool Kernel::is_element_in_set(std::set<SetType> const &s, SetType const &element) { return s.find(element) != s.end(); }

void Kernel::log_spaces(const std::list<Space> &spaces) {
    std::cout << "\tID\tFacade\tFaces\tOld\tSet\tV\t\t\tZmin\tZmax\tNbs\n";
    for (const auto &space: spaces)
        space.Log();
}

void Kernel::log_zones(const std::vector<std::tuple<std::array<double, 3>, TopoDS_Shape, std::list<cFace>>

> &zones) {
    unsigned int i = 0;
    for (
        const auto &zone
            : zones) {
        auto mid = std::get<0>(zone);
        auto fuse = std::get<1>(zone);
        auto cfaces = std::get<2>(zone);
        std::cout << "\t" << std::setfill('0') << std::setw(6) <<
                  i;
        std::cout << "\t" << std::boolalpha << fuse.

                IsNull();

        std::cout << "\t" << std::setfill('0') << std::setw(6) <<
                  Topo(fuse)
                          .

                                  faces()

                          .

                                  Size();

        std::cout << "\t" << std::setfill('0') << std::setw(6) << cfaces.

                size();

        if (!cfaces.

                empty()

                ) {
            std::cout << "\t" << std::setfill('0') << std::setw(6) << cfaces.

                            front()

                    .

                            ID();

            std::cout << "\t" << std::setfill('0') << std::setw(6) << cfaces.

                            back()

                    .

                            ID();

        }
        std::cout << "\t" << std::setfill('0') << std::setw(5) << std::fixed << std::setprecision(2) << mid[0];
        std::cout << "\t" << std::setfill('0') << std::setw(5) << std::fixed << std::setprecision(2) << mid[1];
        std::cout << "\t" << std::setfill('0') << std::setw(5) << std::fixed << std::setprecision(2) << mid[2];
        std::cout << "\n";
        i++;
    }
}

void Kernel::collect_spaces_and_ifcfaces(std::list<Product> &products, ifcspaceInfoList &ifcspaces, std::list<oFace> &ifcspace_ifc_faces, std::list<cFace> &ifcspace_faces, std::list<Space> &spaces, unsigned int &space_id_counter, std::unordered_map<Space *, ifcspaceInfo *> &M,
                                         unsigned int &fid) {

    for (auto &space: ifcspaces) {

        products.emplace_back(space.product, space.guid, space.shape); // add spaces to products
        std::set<cFace *> spaceFaces;

        for (const auto &f: Topo(space.shape).faces()) {

            ifcspace_ifc_faces.emplace_back(TopoDS::Face(f.Complemented()), &products.back(), 0); // complemented to be consistent with Graph.cpp. normals pointing into space now
            ifcspace_ifc_faces.back().SetNormalStatus(FACE_NORMAL_KNOWN);

            ifcspace_faces.emplace_back(ifcspace_ifc_faces.back().face, &ifcspace_ifc_faces.back(), fid);
            // ifcspace_faces.back().ComplementFace(); // to be consistent with Graph.cpp. normals pointing into space now
            // ifcspace_faces.back().SetNormalStatus(FACE_NORMAL_KNOWN);
            // ifcspace_faces.back().SetFixedFaceNormal(ifcspace_faces.back().FaceNormal());
            fid++;
            spaceFaces.insert(&ifcspace_faces.back());
        }

        spaces.emplace_back(space_id_counter, spaceFaces, false);
        for (auto &spaceFace: spaceFaces) spaceFace->SetSpace(&(spaces.back()));
        M[&spaces.back()] = &space;
        space_id_counter++;
    }
}

void Kernel::find_clipping_pairs_all(cfaceSetMap &M_planar, cfaceSetMap &M_curved, std::list<cFace> &space_faces, std::list<cFace> &element_faces, double depth_2b) {

    auto start = std::chrono::high_resolution_clock::now();

    cface_tree3D tree_inner;  // fill rtree
    for (auto &f: element_faces) {
        Bnd_Box bnd = aabb(f.face, 0.001);
        double min[3] = {bnd.CornerMin().X(), bnd.CornerMin().Y(), bnd.CornerMin().Z()};
        double max[3] = {bnd.CornerMax().X(), bnd.CornerMax().Y(), bnd.CornerMax().Z()};
        tree_inner.Insert(min, max, &f);
    }

    for (auto &cface: space_faces) {
        gp_Vec v(cface.FixedFaceNormal().Reversed());
        v.Scale(depth_2b);
        Bnd_Box bnd = aabb(BRepPrimAPI_MakePrism(cface.face, v).Shape(), 0.001);
        double min[3] = {bnd.CornerMin().X(), bnd.CornerMin().Y(), bnd.CornerMin().Z()};
        double max[3] = {bnd.CornerMax().X(), bnd.CornerMax().Y(), bnd.CornerMax().Z()};
        std::list<cFace *> found_inner_faces;
        tree_inner.Search(min, max, [&found_inner_faces](cFace *found_inner_face) {
            found_inner_faces.push_back(found_inner_face);
            return true;
        });

        for (auto &f: found_inner_faces) {
            if (cface.IsPolygon() && f->IsPolygon()) M_planar[&cface].insert(f);
            else {
                std::cerr << "non_planar" << std::endl;
                M_curved[&cface].insert(f);
            }
        }
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Find clip face pairs", std::to_string(M_planar.size()) + ", " + std::to_string(M_curved.size()));
}

void Kernel::find_clipping_pairs_space_approach(cfaceSetMap &M_planar, cfaceSetMap &M_curved, std::list<cFace> &cFaces, double depth_2b, double min_angle) {

    auto start = std::chrono::high_resolution_clock::now();

    cface_tree3D tree_1st;  // fill rtree
    for (auto &f: cFaces) {
        gp_Vec v2(f.FixedFaceNormal().Reversed());
        v2.Scale(depth_2b);
        Bnd_Box bnd2 = aabb(BRepPrimAPI_MakePrism(f.face, v2).Shape(), 0.001);
        double min2[3] = {bnd2.CornerMin().X(), bnd2.CornerMin().Y(), bnd2.CornerMin().Z()};
        double max2[3] = {bnd2.CornerMax().X(), bnd2.CornerMax().Y(), bnd2.CornerMax().Z()};
        tree_1st.Insert(min2, max2, &f);
    }

    for (auto &cface: cFaces) {
        Bnd_Box bnd = aabb(cface.face, 0.001);
        double min[3] = {bnd.CornerMin().X(), bnd.CornerMin().Y(), bnd.CornerMin().Z()};
        double max[3] = {bnd.CornerMax().X(), bnd.CornerMax().Y(), bnd.CornerMax().Z()};
        std::list<cFace *> found_faces;
        tree_1st.Search(min, max, [&found_faces](cFace *found_inner_face) {
            found_faces.push_back(found_inner_face);
            return true;
        });

        for (auto &f: found_faces) {
            if (cface.ID() == f->ID() || cface.FixedFaceNormal().Angle(f->FixedFaceNormal()) < min_angle) continue;
            if (cface.IsPolygon() && f->IsPolygon()) M_planar[&cface].insert(f);
            else {
                std::cerr << "non_planar" << std::endl;
                M_curved[&cface].insert(f);
            }
        }
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Find clip face pairs", std::to_string(M_planar.size()) + ", " + std::to_string(M_curved.size()));
}

void Kernel::identify_sb_types_ray_tracing_space_approach(std::list<cFace> &cFaces, gp_XYZ bounds_min, double transmission_length, double min_angle) {

    // empty mesh data containers
    std::list<std::array<double, 3>> vertices;
    std::list<std::array<unsigned int, 3>> triangles;
    std::list<unsigned int> attributes;
    std::unordered_map<unsigned int, cFace *> M;

    // pre-processing
    triangulate_cfaces_for_ray_tracer(cFaces);
    create_mesh_for_ray_tracer(cFaces, vertices, triangles, attributes, M);

    // build intersector
    auto start = std::chrono::high_resolution_clock::now();

    IntersectorInterface intersector(vertices, triangles, attributes, 7, 10, false);

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Build intersector", std::to_string(triangles.size()));

    // calculate rays
    const double tol = 1.0e-5;
    std::vector<cFace *> send;
    std::vector<IntersectionRay> rays;
    create_rays_for_ray_tracer(cFaces, send, rays, tol, transmission_length);

    // perform intersections
    perform_ray_tracing_space_approach(intersector, send, rays, tol, M, min_angle, bounds_min);
}

void Kernel::perform_ray_tracing_space_approach(IntersectorInterface &intersector, const std::vector<cFace *> &send, std::vector<IntersectionRay> &rays, double tol, std::unordered_map<unsigned int, cFace *> &M, double min_angle, gp_XYZ bounds_min) {

    auto start = std::chrono::high_resolution_clock::now();

    std::unordered_map<Product *, std::vector<cFace * >> virtuals;

    for (unsigned int i = 0; i < send.size(); i++) {

        cFace *cface = send[i];

        IntersectionRay &ray = rays[i];
        cFace *corr = nullptr;
        bool is_corresponding = false;
        double dist = -1;

        std::multimap<double, unsigned int> hits = intersector.Perform(ray);

        if (hits.empty())
            std::cout << "[Info] Ray hit nothing. " << cface->Info() << "\n";

        for (const auto &hit: hits) {

            cFace *hit_cface = M[hit.second];
            if (cface->ID() == hit_cface->ID()) continue; // skip self-hit

            // identify parent face by checking if first level face of HIT was parent of first level face of CURRENT
            if (cface->Superface()->Parent() != nullptr && cface->Superface()->Parent() == hit_cface->Superface()) {
                cface->SetParent(hit_cface);
                continue;
            }

            if (cface->Superface() == hit_cface->Superface()->Parent()) continue; // first level face of CURRENT is parent face of first level face of HIT (wall hitting own child opening face)
            if (cface->Superface()->Parent() == nullptr && hit_cface->Superface()->Parent() != nullptr) continue; // wall hitting an opening
            if (cface->Superface()->Parent() != nullptr && hit_cface->Superface()->Parent() == nullptr) continue; // opening hitting a wall

            // filter corresponding if not parallel
            double angle = cface->FixedFaceNormal().Angle(hit_cface->FixedFaceNormal());
            if (angle > min_angle)
                is_corresponding = true;

            // because map is ordered, hit_cface is the "corresponding" face
            corr = hit_cface;
            dist = hit.first - tol; // correct distance by length origin was moved
            break;
        }
        cface->SetPropertiesSpaceApproach(corr, is_corresponding, dist, bounds_min.Z(), cface->PointOnFace());

    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Perform ray tracing", std::to_string(send.size()));
}

void Kernel::perform_ray_tracing_space_approach_elements(std::list<cFace> &element_faces, std::list<cFace> &cFaces_2ndLvl, double transmission_length) {

    // empty mesh data containers
    std::list<std::array<double, 3>> vertices;
    std::list<std::array<unsigned int, 3>> triangles;
    std::list<unsigned int> attributes;
    std::unordered_map<unsigned int, cFace *> M;

    // pre-processing
    triangulate_cfaces_for_ray_tracer(element_faces);
    create_mesh_for_ray_tracer(element_faces, vertices, triangles, attributes, M);

    // build intersector
    auto start = std::chrono::high_resolution_clock::now();

    IntersectorInterface intersector(vertices, triangles, attributes, 7, 10, false);
    // IntersectorInterface intersector("/home/fluid/Schreibtisch/a.stl", 7, 10, false);

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Build intersector", std::to_string(triangles.size()));

    // calculate rays
    std::vector<cFace *> send;
    std::vector<IntersectionRay> rays;
    const double tol = 1.0e-5;
    create_rays_for_ray_tracer_space_approach(cFaces_2ndLvl, send, rays, tol, transmission_length);

    // perform intersections
    perform_ray_tracing_space_approach(intersector, send, rays, tol, M);
}

void Kernel::create_rays_for_ray_tracer_space_approach(std::list<cFace> &cFaces_2ndLvl, std::vector<cFace *> &send, std::vector<IntersectionRay> &rays, double tol, double transmission_length) {

    auto start = std::chrono::high_resolution_clock::now();

    unsigned int n = cFaces_2ndLvl.size();
    rays.reserve(n);
    send.reserve(n);

    // create ray
    for (auto &cface: cFaces_2ndLvl) {
        IntersectionRay ray = cface.RayBehindSpaceApproach(tol, transmission_length);
        rays.push_back(ray);
        send.push_back(&cface);
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Create rays", std::to_string(rays.size()));
}

void Kernel::perform_ray_tracing_space_approach(IntersectorInterface &intersector, const std::vector<cFace *> &send, std::vector<IntersectionRay> &rays, double tol, std::unordered_map<unsigned int, cFace *> &M) {

    auto start = std::chrono::high_resolution_clock::now();

    for (unsigned int i = 0; i < send.size(); i++) {

        IntersectionRay &ray = rays[i];
        cFace *cface = send[i];

        cface->ClearMaterials();

        std::multimap<double, unsigned int> hits = intersector.Perform(ray);

        if (hits.empty()) continue;

        std::map<double, std::set<cFace *>> behinds;

        for (const auto &hit: hits) {
            cFace *hit_cface = M[hit.second];
            double dist = hit.first - tol; // correct distance by length origin was moved
            dist = round_double_to_n_decimal_places(dist, 5);

            if (hit_cface->IsIfcClass("IfcSite")) {
                cface->SetInternalOrExternal(SB_IE_EXTERNAL_EARTH); // TODO set this later and test if distance is same as wall thickness (biggest distance value in material list)
                continue;
            }

            cface->AppendToMaterial(hit_cface->RelProduct(), dist);
            behinds[round_double_to_n_decimal_places(dist, 3)].insert(hit_cface);

            if (dist < 1e-6) {
                cface->SetAncestor(hit_cface->Ancestor()); // ancestor is now a building element instead of an ifcspace
                if (cface->IsIfcClass("IfcOpeningElement"))
                    cface->SetPhysicalOrVirtual(SB_PV_VIRTUAL);
            }
        }

        // for 2b faces transmission length was used to find hit elements. remove elements that occur behind air gaps // TODO ray tracing with spaces to define internal/external and distance for 2b
        std::map<double, std::list<Product *>> D;
        for (auto &mat: cface->Materials())
            for (auto &d: mat.second)
                D[round_double_to_n_decimal_places(d, 3)].push_back(mat.first);

        double maximum_depth = D.begin()->first;
        std::set<Product *> odd, even;
        std::map<Product *, std::set<double>> NewMaterial;

        // if odd is empty, break, that's the air gap
        for (auto &d: D) {
            for (auto &p: d.second) {
                const bool is_in_odd = odd.find(p) != odd.end();
                const bool is_in_even = even.find(p) != even.end();
                if (!is_in_odd && !is_in_even) odd.insert(p);
                else if (is_in_odd) {
                    even.insert(p);
                    odd.erase(p);
                } else {
                    odd.insert(p);
                    even.erase(p);
                }
            }

            for (auto &p: d.second) NewMaterial[p].insert(d.first);
            maximum_depth = d.first;
            if (odd.empty()) break;
        }

        cface->SetMaterials(NewMaterial);

        // identify external faces as faces without behind space face (currently 2b) and no odd number of hits of one single product (outside element). The material layers must be consecutive, if there is a gap the face is external
        // there will be 2b faces flagged as internal despite being inside but they are ignored anyways for simulation // TODO ray tracing with ifcspaces
        if (cface->SBType() == SB_TYPE_2B && cface->InternalOrExternal() == SB_IE_INTERNAL) {
            cface->SetInternalOrExternal(SB_IE_EXTERNAL);
            for (auto &mat: cface->Materials())
                if (mat.second.size() % 2 != 0) cface->SetInternalOrExternal(SB_IE_INTERNAL);
        }

        // identify external 2a as external faces that have a corresponding face
        if (cface->InternalOrExternal() != SB_IE_INTERNAL) {
            cFace *behind_cface = *behinds[maximum_depth].begin();
            double angle = cface->FixedFaceNormal().Angle(behind_cface->FixedFaceNormal());

            if (angle > 179 * M_PI / 180 || angle < 1 * M_PI / 180) {
                cface->SetCorresponding(behind_cface);
                cface->SetSBType(SB_TYPE_2A);
            }
        }
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Perform ray tracing", std::to_string(send.size()));
}

void Kernel::create_virtual_elements_space_approach(std::list<cFace> &cFaces_2ndLvl, std::list<Product> &products, std::list<oFace> &ifcspace_ifc_faces) {

    auto start = std::chrono::high_resolution_clock::now();

    if (ifcSchema == IFC2X3)
        create_virtual_elements_space_approach_worker<Ifc2x3>(cFaces_2ndLvl, products, ifcspace_ifc_faces);
    else if (ifcSchema == IFC4)
        create_virtual_elements_space_approach_worker<Ifc4>(cFaces_2ndLvl, products, ifcspace_ifc_faces);
    else if (ifcSchema == IFC4X1)
        create_virtual_elements_space_approach_worker<Ifc4x1>(cFaces_2ndLvl, products, ifcspace_ifc_faces);
    else if (ifcSchema == IFC4X2)
        create_virtual_elements_space_approach_worker<Ifc4x2>(cFaces_2ndLvl, products, ifcspace_ifc_faces);
    else if (ifcSchema == IFC4X3_RC1)
        create_virtual_elements_space_approach_worker<Ifc4x3_rc1>(cFaces_2ndLvl, products, ifcspace_ifc_faces);
    else {
        std::cerr << "[Error] Schema not implemented yet! " << std::endl;
        return;
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Create virtual elements", std::to_string(cFaces_2ndLvl.size()));
}

template<typename Schema>
void Kernel::create_virtual_elements_space_approach_worker(std::list<cFace> &cFaces_2ndLvl, std::list<Product> &products, std::list<oFace> &ifcspace_ifc_faces) {

    for (auto &cface: cFaces_2ndLvl)
        if (cface.IsIfcClass("IfcSpace")) {

            auto IfcVirtualElement = ice::IfcVirtualElement<Schema>();
            products.emplace_back(IfcVirtualElement, IfcVirtualElement->GlobalId(), TopoDS_Shape());

            ifcspace_ifc_faces.emplace_back(cface.face, &products.back(), 0);
            cface.SetAncestor(&ifcspace_ifc_faces.back());
            cface.SetPhysicalOrVirtual(SB_PV_VIRTUAL);

            if (cface.Corresponding() != nullptr) {
                ifcspace_ifc_faces.emplace_back(cface.Corresponding()->face, &products.back(), 0);
                cface.Corresponding()->SetAncestor(&ifcspace_ifc_faces.back());
                cface.Corresponding()->SetPhysicalOrVirtual(SB_PV_VIRTUAL);
            }
        }
}

void Kernel::check_corresponding_face_pairs_space_approach(std::list<cFace> &cFaces) {

    auto start = std::chrono::high_resolution_clock::now();

    unsigned int i = 0;

    for (auto &cface: cFaces) {

        if (cface.PhysicalOrVirtual() == SB_PV_VIRTUAL && cface.Corresponding() == nullptr) std::cout << "[Info] Virtual face does not have a corresponding face. " << cface.Info() << "\n";

        if (cface.SpaceBehind() != nullptr) if (cface.SpaceBehind()->is_facade) continue;

        if (cface.Corresponding() != nullptr) {
            if (cface.Corresponding()->IsTrash()) std::cerr << "[Warning] Corresponding face " << cface.Corresponding()->Info() << " of face " << cface.Info() << " is trash." << std::endl;

            bool good = cface.CheckCorresponding(0.001, 3.05);

            if (!good) {
                i++;
                // nullify corresponding's attributes
                cface.Corresponding()->SetCorresponding(nullptr);
                cface.Corresponding()->SetSBType(SB_TYPE_2B);
                // nullify own attributes
                cface.SetCorresponding(nullptr);
                cface.SetSBType(SB_TYPE_2B);
            }
        } else if (cface.SBType() == SB_TYPE_2A) {
            std::cerr << "[Warning] No corresponding face for 2a face " << cface.Info() << "." << std::endl;
            cface.SetSBType(SB_TYPE_2B);
        }
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Check corresponding faces ", std::to_string(i) + "/" + std::to_string(cFaces.size()));
}

void Kernel::identify_duplicate_faces_unknown_normals(std::list<cFace> &cFaces, std::unordered_map<std::string, unsigned int> ranks) {

    // identify mat-mat faces by checking face normals of faces with same id (which were created during fuse, when to faces were coplanar).

    auto start = std::chrono::high_resolution_clock::now();

    // find duplicate faces
    std::unordered_map<unsigned int, std::list<cFace *>> id2pairs;

    for (auto &cface: cFaces) {
        unsigned int id = cface.FaceID();
        if (id2pairs.find(id) == id2pairs.end()) {
            std::list<cFace *> l = {&cface};
            id2pairs[id] = l;
        } else
            id2pairs[id].push_back(&cface);
    }

    // check for mat-mat faces using face normals
    for (auto &id2pair: id2pairs) {

        // skip faces that have no duplicate ids
        if (id2pair.second.size() == 1) continue;

        // get all mat faces of the duplicate faces
        std::list<cFace *> mat_faces, off_faces;

        for (auto &it_face: id2pair.second)
            it_face->IsOffset() ? off_faces.push_back(it_face) : mat_faces.push_back(it_face);

        if (mat_faces.size() == 1) { // if only one mat face paired with offset faces, then delete all offsets and keep mat
            for (auto &off_face: off_faces)
                off_face->SetIsTrash(true);

        } else if (mat_faces.size() > 1) { // if more than one mat face paired with offset faces, then ...

            // delete all offsets and ...
            for (auto &off_face: off_faces)
                off_face->SetIsTrash(true);

            // if all same normal, keep one
            std::multimap<int, cFace *> rank_to_mat_face = sort_mat_mat_faces(mat_faces, ranks);
            auto it_face = rank_to_mat_face.begin();
            for (std::advance(it_face, 1); it_face != rank_to_mat_face.end(); ++it_face) { // skip first face
                it_face->second->SetIsTrash(true);
                it_face->second->SetIsCoplanar(true);
            }

        } else { // if no mat faces and only offset faces

            // keep one
            auto it_face = off_faces.begin();
            for (std::advance(it_face, 1); it_face != off_faces.end(); ++it_face) // skip first face
                (*it_face)->SetIsTrash(true);
        }
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Identify duplicate faces", std::to_string(cFaces.size()));
}

void Kernel::ifaces_nbs(iFace &iface, used_orientation o1, std::unordered_map<cFace *, iFace *> &M) {

    oriFace null_nb(nullptr, USEDFORWARD);

    if (o1 == USEDREVERSED && iface.cface->NormalStatus() == FACE_NORMAL_KNOWN) {
        for (auto &a: iface.cface->adjacentfaces)
            iface.nb[o1][a.first] = null_nb;
    } else {

        for (auto &a: iface.cface->adjacentfaces) {

            auto edge_id = a.first;
            auto e1 = TopoDS::Edge(iface.cface->halfedges[edge_id]);
            //auto v1 = o1 == USEDFORWARD ? TopExp::FirstVertex(e1, Standard_True) : TopExp::LastVertex(e1, Standard_True);
            auto r1 = o1 == USEDFORWARD ? e1.Orientation() : (e1.Orientation() == TopAbs_FORWARD ? TopAbs_REVERSED : TopAbs_FORWARD);

            std::list<oriFace> L;

            for (auto &nb: a.second) {

                // neighbour face
                iFace *f2 = M[nb];

                // check edge orientations to find adequate orientation and add to nbs
                auto e2 = TopoDS::Edge(nb->halfedges[edge_id]);
                // auto v2 = TopExp::FirstVertex(e2, Standard_True);
                // used_orientation o2 = v1.IsSame(v2) ? USEDREVERSED : USEDFORWARD;
                used_orientation o2 = e2.Orientation() == r1 ? USEDREVERSED : USEDFORWARD;
                L.emplace_back(f2, o2);
            }

            if (L.empty()) {
                iface.nb[o1][edge_id] = null_nb;
            } else {
                oriFace nb = L.size() == 1 ? nb = L.front() : *(iface.cface->FaceWithMinimalAngleConsiderBothOrientations(edge_id, o1, L));
                iface.nb[o1][edge_id] = (nb.orient == USEDREVERSED && nb.iface->cface->NormalStatus() == FACE_NORMAL_KNOWN) ? null_nb : nb; // remove nb if it has an incompatible orientation
            }
        }
    }

}

void Kernel::ifaces_nbs(iFace &iface, std::unordered_map<cFace *, iFace *> &M) {
    ifaces_nbs(iface, USEDFORWARD, M);
    ifaces_nbs(iface, USEDREVERSED, M);
}

std::stack<oriFace> Kernel::find_spaces_normals_unknown_start_face(std::queue<oriFace> &q) {

    std::stack<oriFace> stack;

    while (!q.empty()) {

        oriFace oface = q.front();
        q.pop();

        if (oface.iface->bad) continue;

        if (oface.orient == USEDFORWARD) {
            if (!oface.iface->isUsedForward) {
                stack.push(oface);
                return stack;
            }
        } else {
            if (!oface.iface->isUsedReversed) {
                stack.push(oface);
                return stack;
            }
        }
    }

    return stack;
}

bool Kernel::find_spaces_normals_unknown_accept_solution(const std::set<oriFace> &L) {

    if (L.size() < 2)
        return false;

    std::set<Product *> P;
    for (auto &l: L)
        P.insert(l.iface->cface->RelProduct());
    if (P.size() < 3) return false;

    // closed if for every edge there is same amount of half-edges with opposite orientation
    std::unordered_map<unsigned int, std::map<unsigned int, std::list<bool>>> ids;
    for (auto &l: L)
        for (auto &edge: l.iface->cface->halfedges) {
            unsigned int v1 = hash(l.orient == USEDFORWARD ? TopExp::FirstVertex(TopoDS::Edge(edge.second), Standard_True) : TopExp::LastVertex(TopoDS::Edge(edge.second), Standard_True));
            ids[edge.first][v1].emplace_back(true);
        }

    for (auto &id: ids) {

        const auto &m = id.second;

        if (m.size() != 2) return false;
        if (m.begin()->second.size() != m.rbegin()->second.size()) return false;
    }

    return true;
}

double Kernel::find_spaces_normals_unknown_volume(const std::set<oriFace> &L) {

    BRepBuilderAPI_Sewing sew;

    for (auto l: L)
        l.orient == USEDFORWARD ? sew.Add(l.iface->cface->face) : sew.Add(l.iface->cface->face.Complemented());

    sew.Perform();
    TopoDS_Shape shape = sew.SewedShape();
    TopExp_Explorer Ex;
    TopoDS_Shell S;

    for (Ex.Init(shape, TopAbs_SHELL); Ex.More(); Ex.Next()) {
        S = TopoDS::Shell(Ex.Current());
        break;
    }

    return volume(S);
}

std::queue<oriFace> Kernel::find_spaces_normals_unknown_queue(std::list<iFace> &ifaces) {

    // usedforward means normal is same as face normal calculated from occ shape (considering topabs_orientation)

    std::queue<oriFace> q;

    // prioritize faces of known orientation and wall faces.
    for (auto &iface: ifaces)
        if (iface.cface->NormalStatus() == FACE_NORMAL_KNOWN && iface.cface->IsIfcClass("IfcWall"))
            q.emplace(&iface, USEDFORWARD);

    for (auto &iface: ifaces)
        if (iface.cface->NormalStatus() == FACE_NORMAL_KNOWN && !iface.cface->IsIfcClass("IfcWall"))
            q.emplace(&iface, USEDFORWARD);

    for (auto &iface: ifaces)
        if (iface.cface->NormalStatus() != FACE_NORMAL_KNOWN) {
            q.emplace(&iface, USEDFORWARD);
            q.emplace(&iface, USEDREVERSED);
        }

    return q;
}

void Kernel::find_spaces_normals_unknown(std::list<cFace> &cFaces, std::list<Space> &spaces, unsigned int &space_id, unsigned int &fid) {

    auto start = std::chrono::high_resolution_clock::now();

    // Rules:
    // (A) Every face could be used twice but with different orientation. Exception: Faces that can only be used once with that specific orientation: 1) Faces with known orientation 2) Building element faces (non offset) when use_non_virtual_faces_twice is false
    // (B) Accept solutions for space, when none of the following criteria applies:
    // (B.1) - Only faces of one product in space
    // (B.2) - Contains at least one faces that can not be used with necessary orientation, because of rule (A)
    // (B.3) - Prefer useful solutions to useless solutions (e.g. if all faces of a product used in spaces point inwards, the solution is not useful)
    // (C) - Within a space one face can only be used with one orientation

    std::list<iFace> ifaces;
    std::unordered_map<cFace *, iFace *> M;

    for (auto &cface: cFaces) {
        ifaces.emplace_back(&cface);
        M[&cface] = &ifaces.back();
    }

    // calculate chosen adjacent nb for each edge on an iface for both orientations
#ifdef PARALLEL_PROCESSING
    std::vector<iFace *> V;
    V.reserve(ifaces.size());
    for (auto &iface: ifaces)
        V.push_back(&iface);
    unsigned int chunk = std::ceil(0.1 * ifaces.size() / num_threads);

#pragma omp parallel for default(none) shared(V, M, chunk) schedule(dynamic, chunk) num_threads(num_threads)
    for (auto &i: V)
        ifaces_nbs(*i, M);
#else
    for (auto &iface: ifaces)
        ifaces_nbs(iface, M);
#endif

    // generate queue for finding start face. prioritize faces of known orientation.
    // TODO use octree and ray tracer (without offsets) to get spaces, fluidpoints and start faces with known orientation (especially for facade). if start face is found by ray all non-virtual faces of found component can be set to known
    // (ray tracer uses point as start for ray. hit face is start face and gets needed orientation)
    std::queue<oriFace> q = find_spaces_normals_unknown_queue(ifaces);

    // find components building closed volumes
    std::multimap<double, std::set<oriFace>> components = find_spaces_normals_unknown_recursion(ifaces, cFaces, q);

    std::cout << "Volumes: ";
    for (auto &i: components)
        std::cout << round_double_one_digit(i.first) << " ";
    std::cout << "\n";

    // detect number of occurrences of faces in all components and set duplicate attribute for ifaces
    find_spaces_normals_unknown_face_usages(components);

    // save spaces
    find_spaces_normals_unknown_save_spaces(cFaces, spaces, space_id, fid, components, M);

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << print_time(elapsed.count(), "Search components", std::to_string(spaces.size()));
}

void Kernel::find_spaces_normals_unknown_face_usages(const std::multimap<double, std::set<oriFace>> &components) {

    std::unordered_map<iFace *, std::list<used_orientation>> usages;

    for (auto &v: components)
        for (auto &oface: v.second)
            usages[oface.iface].push_back(oface.orient);

    for (auto &u: usages) {

        const auto &os = u.second;
        const auto &iface = u.first;

        if (os.size() > 1) {

            iface->duplicate = true;

            auto n1 = std::count(os.begin(), os.end(), USEDFORWARD);
            auto n2 = std::count(os.begin(), os.end(), USEDREVERSED);

            std::cout << "[Info] Face " << iface->cface->Info() << " occurs " << os.size() << " times (" << n1 << " forward, " << n2 << " reversed) in all spaces.\t" << iface->cface->IfcClass() << "\n";
            if (n1 > 1) std::cerr << "[Warning] Redundant use forward (" << n1 << ") " << iface->cface->Info() << std::endl; // bad_duplicates.insert(iface); // USEDFORWARD
            if (n2 > 1) std::cerr << "[Warning] Redundant use reversed (" << n2 << ") " << iface->cface->Info() << std::endl; // bad_duplicates.insert(iface); // USEDREVERSED
        }
    }

    //for (auto &iface: duplicates)
    //    std::cout << "[Info] Face " << iface->cface->ID() << " (" << iface->cface->IfcClass() << ", " << iface->cface->IsOffset() << ") occurs multiple times.\n";
}

void Kernel::find_spaces_normals_unknown_save_spaces_component(const std::pair<const double, std::set<oriFace>> &i, std::list<cFace> &cFaces, std::list<Space> &spaces, unsigned int &space_id, unsigned int &fid, const std::multimap<double, std::set<oriFace>> &components) {

    // TODO additional criteria: percent of fixed normal faces, percent of already used faces, (space lying within an object. ray tracer detects if space is bounded or within two object faces with normals pointing outward)

    // skip rules
    if (i.first > 0 && i != *components.rbegin()) { // volume
        std::cout << "[Info] Skip non-facade component\n";
        return; // only one space with positive volume allowed (facade)
    }

    for (auto &oface: i.second)
        if (!oface.iface->cface->IsVirtual() && !oface.iface->cface->IsOffset()) // virtual faces and offset faces can be used more than once (TODO maximum two times, once in each direction)
            if (oface.iface->used) {
                std::cout << "[Info] Already used " << oface.iface->cface->Info() << ", " << oface.iface->cface->IsOffset() << "). Skip component (" << i.first << ").\n";
                return;
            }

    std::set<cFace *> space_faces;

    for (auto &oface: i.second) {

        const auto &cface = oface.iface->cface;
        const auto &o = oface.orient;

        oface.iface->used = true;

        if (!oface.iface->duplicate) { // face only used once
            if (o == USEDREVERSED)
                cface->ComplementFace();
            space_faces.insert(cface);
        } else {
            if (o == USEDREVERSED) {
                add_face_reversed(cFaces, cface, fid);
                space_faces.insert(&cFaces.back());
            } else
                space_faces.insert(cface);
        }
    }

    spaces.emplace_back(space_id, space_faces, false);
    for (auto &sf: space_faces)
        sf->SetSpace(&(spaces.back()));
    space_id++;
}

void Kernel::find_spaces_normals_unknown_save_spaces(std::list<cFace> &cFaces, std::list<Space> &spaces, unsigned int &space_id, unsigned int &fid, const std::multimap<double, std::set<oriFace>> &components, std::unordered_map<cFace *, iFace *> &M) {

    // store facade first
    find_spaces_normals_unknown_save_spaces_component(*components.rbegin(), cFaces, spaces, space_id, fid, components);

    for (auto i = components.begin(); i != components.end(); ++i)
        if (i->first < 0) // skip facade
            find_spaces_normals_unknown_save_spaces_component(*i, cFaces, spaces, space_id, fid, components);

    for (auto &cface: cFaces)
        if (cface.RelSpace() != nullptr) {
            cface.SetNormalStatus(FACE_NORMAL_KNOWN);
        } else {
            cface.SetIsTrash(true);

            // add to inner if no counter-oriented version of this face is used in any space
            if (!M[&cface]->used)
                cface.SetIsInner(true);
        }
}

std::multimap<double, std::set<oriFace>> Kernel::find_spaces_normals_unknown_recursion(std::list<iFace> &ifaces, const std::list<cFace> &cFaces, std::queue<oriFace> &q) {

    std::multimap<double, std::set<oriFace>> components;

    while (true) {

        std::stack<oriFace> stack = find_spaces_normals_unknown_start_face(q);

        if (stack.empty())
            break;

        std::set<oriFace> L, tested;

        while (!stack.empty()) {

            oriFace oface = stack.top();
            L.insert(oface);
            stack.pop();

            for (auto &a: oface.iface->cface->adjacentfaces) {

                oriFace nb = oface.iface->nb[oface.orient][a.first];

                // skip faces under certain conditions
                if (nb.iface == nullptr) continue; // no nb on that edge
                if (nb.iface->bad) continue;
                if (tested.find(nb) != tested.end()) continue; // face was already tested
                tested.insert(nb);
                if (nb.orient == USEDREVERSED && nb.iface->cface->NormalStatus() == FACE_NORMAL_KNOWN) continue; // normal of nb face is locked (because it is known) but face had to be used reversed
                if (L.find(nb) != L.end()) continue; // face is already part of this component
                if (nb.orient == USEDFORWARD && nb.iface->isUsedForward || nb.orient == USEDREVERSED && nb.iface->isUsedReversed) continue; // face is not allowed to use. // not necessary, because multiple solutions with a face allowed

                stack.push(nb);
            }
        }

        if (!find_spaces_normals_unknown_accept_solution(L)) {
            for (auto &oface: L)
                oface.iface->bad = true;
            continue;
        }

        double V = find_spaces_normals_unknown_volume(L);
        if (fabs(V) > 0.05) components.insert(std::pair<double, std::set<oriFace>>(V, L));

        // Viewer::visualize_component_search(cFaces, L, oface);
        // Viewer::visualize_component_search(cFaces, L, *L.begin());

        for (auto &oface: L)
            if (oface.orient == USEDFORWARD) oface.iface->isUsedForward = true;
            else oface.iface->isUsedReversed = true;
    }

    // find facade space
    auto min = components.begin()->first;
    auto max = components.rbegin()->first;
    if (max > 0 && fabs(max) > fabs(min))
        std::cout << "[Info] Facade found (" << min << "\t" << max << ").\n";
    else
        std::cerr << "[Warning] No Facade found (" << min << "\t" << max << ")!\n";

    return components;
}

void Kernel::add_face_reversed(std::list<cFace> &cFaces, cFace *cface, unsigned int &fid) {

    cFaces.push_back(*cface);
    cFaces.back().SetID(fid);
    fid++;
    cFaces.back().face = TopoDS::Face(cface->face.Complemented());
    cFaces.back().UpdateHalfEdges();
    for (auto &it: cFaces.back().adjacentfaces)
        for (auto &a: it.second) // add adjacency in nb faces
            a->adjacentfaces[it.first].push_back(&cFaces.back());
}

std::string Kernel::remove_first_and_last_char(std::string s) { return s.size() < 2 ? s : s.substr(1, s.size() - 2); }