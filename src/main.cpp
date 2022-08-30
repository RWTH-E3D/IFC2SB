// Copyright 2022 Eric Fichter
#include "headers.h"

//!\mainpage Emerged from BIM2SIM – Development of methods for the generation of simulation models using Building Information Modeling data
//! Project content: Building performance simulation (BPS) represents an efficient tool to predict energy needs from building systems to district level. Additionally, with BPS in the planning phase alternatives can be analyzed for power generation, storage, distribution, control and emission systems. As result, BPS is used to identify and quantify the optimization potential of such systems. However, there is currently a lack of automation techniques that convert digital information from a planning or as built model to a simulation model in a cost-effective manner. The intention of this project is to develop the basic methods and tools to derive directly BPS model ready for simulation, starting from a digital building information model (BIM). The resulting BPS model shall be used to improve design and to analyze optimization of building energy systems in the field. This is intended to significantly reduce the manual efforts required for modeling and to support the dissemination of various simulation applications in the building sector as well. An additional goal of this project is to link the simulation results as target data back into the model using a BIM-compliant format. Besides the energetic view, the developed methods and tools shall as well be utilized to evaluate the energy expenditure over the building life cycle and to assess the embodied energy.
//! 
//! Project term: 05/2018 – 10/2021
//! 
//! Project partners:\n 
//!     RWTH Aachen University – Institute for Energy Efficient Buildings and Indoor Climate EBC\n 
//!     RWTH Aachen University – Institute of Energy Efficiency and Sustainable Building E3D\n 
//!     Rud. Otto Meyer Technik Ltd. & Co. KG (ROM)
//! 
//!  Funding/Client: BMWi, Förderkennzeichen 03ET1562A
//! 
//!  Developer of IFC2SB:\n Eric Fichter\n RWTH Aachen University – Institute of Energy Efficiency and Sustainable Building E3D\n fichter@e3d.rwth-aachen.de


namespace po = boost::program_options;

void print_usage(bool suggest_help = true) {
    std::cout << "\nUsage: IFC2SB [options] <input.ifc> [<output>]\n"
              << "\n"
              << "Generation of 2nd Level Space Boundaries based on IFC files.\n"
              << "If no output filename given, <input>_with_SB.ifc will be used as the output file.\n";
    if (suggest_help) std::cout << "\nRun 'IFC2SB --help' for more information.";
    std::cout << std::endl;
}

void print_options(const po::options_description &options) {
    std::cout << "\n" << options;
    std::cout << std::endl;
}

void print_version() {
    std::cout << IFC2SB_NAME << " " << IFC2SB_VERSION << " (uses OCCT " << OCC_VERSION_STRING_EXT << " and IfcOpenShell " << IFCOPENSHELL_VERSION << ")\n";
}

bool file_exists(const std::string &filename) {
    std::ifstream file(IfcUtil::path::from_utf8(filename).c_str());
    return file.good();
}

std::string change_extension(const std::string &fn, const std::string &ext) {
    typename std::string::size_type dot = fn.find_last_of('.');
    if (dot != std::string::npos) {
        return fn.substr(0, dot) + ext;
    } else {
        return fn + ext;
    }
}

int main(int argc, char **argv) {

    typedef po::command_line_parser command_line_parser;

    unsigned int num_threads;
    double transmission_length;

    // Generic options
    po::options_description generic_options("Command line options");
    generic_options.add_options()
            ("help,h", "display usage information")
            ("version,v", "display version information")
            ("yes,y", "answer 'yes' automatically to possible confirmation queries (e.g. overwriting an existing output file)")
            ("threads,j", po::value<unsigned int>(&num_threads)->default_value(1), "Number of parallel processing threads for geometry interpretation");

    // File options
    po::options_description fileio_options;
    fileio_options.add_options()
            ("input-file", new po::typed_value<std::string, char>(nullptr), "input IFC file")
            ("output-file", new po::typed_value<std::string, char>(nullptr), "output geometry file");

    // IFC options
    po::options_description ifc_options("Calculation method");
    ifc_options.add_options()
            ("graph", "Use graph approach to calculate space boundaries")
            ("clip", "Use clip approach to calculate space boundaries (only basic structure implemented)");

    // Geometry options
    po::options_description geom_options("Geometry options");
    geom_options.add_options()
            ("stl", "Output an .stl file instead of .ifc (will set first_level_only to true)")
            ("first", "Generate first level space boundaries only")
            ("split", "Use original IfcSpaces to split new spaces")
            ("noholes", "Remove holes from space boundaries")
            ("convex", "Decompose concave space boundaries to convex (noholes must be activated)")
            ("simpleopening", "Simplify fenestration space boundaries")
            ("shadings", "Calculate shadings")
            ("virtualopening", "Use Ifcopeningelelements to create virtual space boundaries")
            ("virtualbyspace", "Use IfcSpaces to create virtual space boundaries")
            ("length,l", po::value<double>(&transmission_length)->default_value(0.41, "0.41"), "Length to distinguish 2a/2b SB type");

    // Command line options
    po::options_description cmdline_options;
    cmdline_options.add(generic_options).add(fileio_options).add(geom_options).add(ifc_options);

    // Positional options
    po::positional_options_description positional_options;
    positional_options.add("input-file", 1);
    positional_options.add("output-file", 1);

    po::variables_map vmap;
    try {
        po::store(command_line_parser(argc, argv).options(cmdline_options).positional(positional_options).run(), vmap);
    } catch (const po::unknown_option &e) {
        std::cerr << "[Error] Unknown option '" << e.get_option_name().c_str() << "'\n\n";
        print_usage();
        return EXIT_FAILURE;
    } catch (const po::error_with_option_name &e) {
        std::cerr << "[Error] Invalid usage of '" << e.get_option_name().c_str() << "': " << e.what() << "\n\n";
        return EXIT_FAILURE;
    } catch (const std::exception &e) {
        std::cerr << "[Error] " << e.what() << "\n\n";
        print_usage();
        return EXIT_FAILURE;
    } catch (...) {
        std::cerr << "[Error] Unknown error parsing command line options\n\n";
        print_usage();
        return EXIT_FAILURE;
    }

    po::notify(vmap);

    if (vmap.count("version")) {
        print_version();
        return EXIT_SUCCESS;
    } else if (vmap.count("help")) {
        print_usage(false);
        print_options(generic_options.add(geom_options).add(ifc_options));
        return EXIT_SUCCESS;
    } else if (!vmap.count("input-file")) {
        std::cerr << "[Error] Input file not specified!" << std::endl;
        print_usage();
        return EXIT_FAILURE;
    }

    // input filename
    const std::string input_filename = vmap["input-file"].as<std::string>();
    if (!file_exists(IfcUtil::path::to_utf8(input_filename))) {
        std::cerr << "[Error] Input file '" << input_filename << "' does not exist!" << std::endl;
        return EXIT_FAILURE;
    }

    // output format
    bool stl = false;
    std::string ending;
    if (vmap.count("stl") == 1) {
        stl = true;
        ending = ".stl";
    } else
        ending = ".ifc";

    bool first_level_only = vmap.count("first") == 1 && !stl; // first level only
    bool space_split = vmap.count("split") == 1; // split spaces using old IfcSpaces
    bool remove_holes_by_face_split = vmap.count("noholes") == 1;
    bool decompose_concave_polygons = vmap.count("convex") == 1;
    bool simplify_fenestration_faces = vmap.count("simpleopening") == 1;
    bool calculate_shadings = vmap.count("shadings") == 1;
    bool use_ifcopeningelelements_for_virtual_boundaries = vmap.count("virtualopening") == 1;
    bool use_spaces_for_virtual_boundaries = vmap.count("virtualbyspace") == 1;

    // output filename
    std::string output_filename;

    if (vmap.count("output-file") == 1) // user made input
        output_filename = change_extension(vmap["output-file"].as<std::string>(), ending);
    else
        output_filename = change_extension(input_filename, "_with_SB" + ending); // use default path and name

    if (output_filename.size() < 5) {
        std::cerr << "[Error] Invalid or unsupported output file '" << output_filename << "' given!" << std::endl;
        print_usage();
        return EXIT_FAILURE;
    }

    if (file_exists(IfcUtil::path::to_utf8(output_filename)) && !vmap.count("yes")) {
        std::string answer;
        std::cout << "A file '" << output_filename << "' already exists. Overwrite the existing file (y/n)?" << std::endl;
        std::cin >> answer;
        if (!boost::iequals(answer, "yes") && !boost::iequals(answer, "y")) {
            return EXIT_SUCCESS;
        }
    }

    // number of threads
    if (num_threads <= 0 || num_threads > std::thread::hardware_concurrency())
        num_threads = std::thread::hardware_concurrency();

    // approach
    const bool graph = vmap.count("graph") != 0;
    const bool clip = vmap.count("clip") != 0;

    // user info
    std::cout << "\n";
    std::cout << "\033[1m\033[37m";
    std::cout << "-------------------------------------------------\n";
    std::cout << "### " << IFC2SB_FULLNAME << " " << IFC2SB_NAME << " ###\n";
    std::cout << "Generation of second level space boundaries for Ifc4 files.\n";
    std::cout << "\n";
    std::cout << "#####  ######  ####   ####   ####  #####  \n";
    std::cout << "  #    #      #    # #    # #      #    # \n";
    std::cout << "  #    #####  #           #  ####  #####  \n";
    std::cout << "  #    #      #       ####       # #    # \n";
    std::cout << "  #    #      #    # #      #    # #    # \n";
    std::cout << "#####  #       ####  ######  ####  #####  \n";
    std::cout << "\n";
    std::cout << "Version:     " << IFC2SB_VERSION << " (2022)" << "\n";
    std::cout << "Developer:   " << IFC2SB_AUTHOR << "\n";
    std::cout << "Institute:   E3D - Institute of Energy Efficiency and Sustainable Building,\n";
    std::cout << "             RWTH Aachen University \n";
    std::cout << "-------------------------------------------------\n";
    std::cout << "\033[0m";
    std::cout << "\n";
    std::cout << "-------------------------------------------------\n";
    std::cout << "# Arguments\n";
    std::cout << "Input file:          " << input_filename << "\n";
    std::cout << "Output file:         " << output_filename << "\n";
    std::cout << "Threads:             " << std::to_string(num_threads) << "\n";
    if (graph) std::cout << "Graph:               " << std::boolalpha << graph << "\n";
    else std::cout << "Clip:                " << std::boolalpha << clip << "\n";
    if(!clip) std::cout << "First level only:    " << std::boolalpha << first_level_only << "\n";
    if(!clip) std::cout << "STL:                 " << std::boolalpha << stl << "\n";
    std::cout << "Transmission length: " << transmission_length << "\n";
    if(!clip) std::cout << "Split:               " << std::boolalpha << space_split << "\n";
    std::cout << "No holes:            " << std::boolalpha << remove_holes_by_face_split << "\n";
    std::cout << "Convex:              " << std::boolalpha << decompose_concave_polygons << "\n";
    std::cout << "Simple Opening:      " << std::boolalpha << simplify_fenestration_faces << "\n";
    std::cout << "Shadings:            " << std::boolalpha << calculate_shadings << "\n";
    std::cout << "Virtual Opening:     " << std::boolalpha << use_ifcopeningelelements_for_virtual_boundaries << "\n";
    std::cout << "Virtual by Spaces:   " << std::boolalpha << use_spaces_for_virtual_boundaries << "\n";
    std::cout << "-------------------------------------------------\n";

    std::cout << "\n-------------------------------------------------" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();

    bool status = true;

    if (graph) {

        Graph G(input_filename,
                output_filename,
                num_threads,
                stl,
                space_split,
                transmission_length,
                first_level_only,
                remove_holes_by_face_split,
                decompose_concave_polygons,
                simplify_fenestration_faces,
                calculate_shadings,
                use_ifcopeningelelements_for_virtual_boundaries,
                use_spaces_for_virtual_boundaries);
        status = G.run();
    } else if (clip) {

        Clip C(input_filename,
               output_filename,
               num_threads,
               transmission_length,
               remove_holes_by_face_split,
               decompose_concave_polygons,
               simplify_fenestration_faces,
               calculate_shadings,
               use_ifcopeningelelements_for_virtual_boundaries);
        status = C.run();
    } else
        std::cout << "Nothing done! Graph and partially Clip are the only implemented space boundary generators." << std::endl;

    if (status) std::cout << "Program finished successfully!" << std::endl;
    else std::cout << "Program interrupted by error!" << std::endl;

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Total Elapsed time: " << elapsed.count() << " s\n";
    std::cout << "-------------------------------------------------\n";

    return status ? EXIT_SUCCESS : EXIT_FAILURE;
}
