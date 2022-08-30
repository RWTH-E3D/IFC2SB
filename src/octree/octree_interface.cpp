#include "octree_interface.h"


int octree_interface::process_stl_file(const std::string &filename, const int oct_depth, const bool write_vtk) {

    //const char *infile = "/home/fluid/Schreibtisch/octree/VisGeom/AC20-FZK-Haus_with_SB.stl";
    //int oct_depth = 4;
    //bool write_vtk = true;

    const char *infile = filename.c_str();

    std::ifstream f(infile);
    if (!f.good()) {
        std::cerr << "[Error] File not found!";
        return 1;
    }

    model *my_model;    // geometric model
    octree *my_octree;    // voxel model

    char *voxfile = new char[strlen(infile) + 1];
    strcpy(voxfile, infile);
    voxfile[strlen(voxfile) - 3] = 'v';
    voxfile[strlen(voxfile) - 2] = 'o';
    voxfile[strlen(voxfile) - 1] = 'x';

    char *vtkfile(nullptr);
    if (write_vtk) {
        vtkfile = new char[strlen(infile) + 1];
        strcpy(vtkfile, infile);
        vtkfile[strlen(vtkfile) - 3] = 'v';
        vtkfile[strlen(vtkfile) - 2] = 't';
        vtkfile[strlen(vtkfile) - 1] = 'k';
    }

    std::cout << "Info: Generating octree for file " << infile << " at level depth " << oct_depth << ".\n";

    std::cout << "Info: Reading in model..." << std::flush;

    my_model = new model();
    my_model->load_model(infile, oct_depth, true);

    my_octree = new octree();
    my_octree->set_model(my_model);
    my_octree->set_filenames(voxfile, vtkfile);
    if (my_octree->gen_octree(oct_depth) != OCTREE_DONE) { std::cerr << " ERROR >> couldn't generate octree" << std::endl; }

    // clean up
    delete my_octree;
    delete my_model;
    if (vtkfile) delete[] vtkfile;
    delete[] voxfile;

    std::cout << "Info: Generation and writeout finished. Terminating execution!\n";
    return 0;
}

void octree_interface::process_mesh(std::vector<std::pair<std::array<double, 3>, std::set<std::string>>> &spaces_guids, std::vector<std::pair<std::array<double, 3>, std::set<unsigned int>>> &spaces_triangles, std::vector<std::set<unsigned int>> &zones,
                                    const std::vector<std::tuple<double, double, double>> &vertices, const std::vector<std::tuple<int, int, int>> &faces, const std::vector<std::string> &attrs, const std::set<std::array<double, 3>> &fluid_points, const int oct_depth, const bool write_vtk) {

    if (faces.size() != attrs.size() || 3 * faces.size() != vertices.size()) {
        std::cerr << "[Error] Corrupt mesh data!" << std::endl;
        std::vector<std::pair<std::array<double, 3>, std::set<std::string>>> e;
        return;
    }

    const char *infile = "/home/fluid/Schreibtisch/octree/VisGeom/Output.xyz";

    char *voxfile = new char[strlen(infile) + 1];
    strcpy(voxfile, infile);
    voxfile[strlen(voxfile) - 3] = 'v';
    voxfile[strlen(voxfile) - 2] = 'o';
    voxfile[strlen(voxfile) - 1] = 'x';

    char *vtkfile(nullptr); // if vtkfile stays nullptr, no try to write vtk file is undertaken
    if (write_vtk) {
        vtkfile = new char[strlen(infile) + 1];
        strcpy(vtkfile, infile);
        vtkfile[strlen(vtkfile) - 3] = 'v';
        vtkfile[strlen(vtkfile) - 2] = 't';
        vtkfile[strlen(vtkfile) - 1] = 'k';
    }

#ifdef PRINT
    std::cout << "\tGenerating octree for file at level depth " << oct_depth << ".\n";
    std::cout << "\tReading in model..." << std::flush;
#endif

    model *my_model;    // geometric model
    octree *my_octree;    // voxel model

    my_model = new model();
    my_model->load_model(vertices, faces, attrs, oct_depth, true);

#ifdef PRINT
    std::cout << "\tdone." << std::endl;
#endif

    my_octree = new octree();
    my_octree->set_model(my_model);
    my_octree->set_filenames(voxfile, vtkfile);
    my_octree->set_maximum_triangles_in_a_zone(6000);
    my_octree->set_fluid_points(fluid_points);
    if (my_octree->gen_octree(oct_depth) != OCTREE_DONE) { std::cerr << " ERROR >> couldn't generate octree" << std::endl; }

    // get spaces
    spaces_guids = my_octree->get_spaces_guids();
    spaces_triangles = my_octree->get_spaces_triangles();
    zones = my_octree->get_zones();

    // clean up
    delete my_octree;
    delete my_model;
    delete[] vtkfile;
    delete[] voxfile;

#ifdef PRINT
    std::cout << "Info: Generation and writeout finished. Terminating execution!\n";
#endif
}