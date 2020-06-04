//Compile with 
//g++ -o config config.cpp -lboost_program_options


/*

g++ -g `pkg-config --cflags ibsimu-1.0.6dev` -c -o config.o config.cpp -lboost_program_options
*/

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>


#include <geometry.hpp>
#include <dxf_solid.hpp>
#include <mydxffile.hpp>
#include <epot_field.hpp>
#include <epot_efield.hpp>
#include <meshvectorfield.hpp>


#include <string>
#include <vector>

/*
*/

// https://www.boost.org/doc/libs/1_60_0/doc/html/program_options/tutorial.html
namespace bpo = boost::program_options;
 
typedef double origin_t; //Start [3]
typedef double mesh_cell_size_h_t; //h
typedef double geometry_value_t; //sizereq[3]
typedef double voltage_t;
typedef double dxf_scale_factor_t;
typedef double bfield_scale_factor_t;
typedef std::string bound_type_string_t;
typedef bound_e bound_type_t;
typedef bpo::variables_map run_parameters_t;

bpo::options_description config_file_options_m() 
{
    bpo::options_description config_file_options_o("Config file options");
    config_file_options_o.add_options()
        ("help", "print help message")
        ("geometry-mode", bpo::value<std::string>(), "Geometry Mode (MODE_3D, MODE_CYL...)")
        ("origin-x", bpo::value<origin_t>(), "geometry origin x")
        ("origin-y", bpo::value<origin_t>(), "geometry origin y")
        ("origin-z", bpo::value<origin_t>(), "geometry origin z")
        ("mesh-cell-size-h", bpo::value<mesh_cell_size_h_t>(), "mesh cell size h")
        ("geometry-start-x", bpo::value<geometry_value_t>(), "geometry start x")
        ("geometry-start-y", bpo::value<geometry_value_t>(), "geometry start y")
        ("geometry-start-z", bpo::value<geometry_value_t>(), "geometry start z")
        ("geometry-size-x", bpo::value<geometry_value_t>(), "geometry size x")
        ("geometry-size-y", bpo::value<geometry_value_t>(), "geometry size y")
        ("geometry-size-z", bpo::value<geometry_value_t>(), "geometry size z")

        ("wall-1-bound-type", bpo::value<bound_type_string_t>(), "wall bound type [BOUND_DIRICHLET, BOUND_NEUMANN]")
        ("wall-1-bound-voltage", bpo::value<voltage_t>(), "wall voltage V")
        ("wall-2-bound-type", bpo::value<bound_type_string_t>(), "wall bound type [BOUND_DIRICHLET, BOUND_NEUMANN]")
        ("wall-2-bound-voltage", bpo::value<voltage_t>(), "wall voltage V")
        ("wall-3-bound-type", bpo::value<bound_type_string_t>(), "wall bound type [BOUND_DIRICHLET, BOUND_NEUMANN]")
        ("wall-3-bound-voltage", bpo::value<voltage_t>(), "wall voltage V")
        ("wall-4-bound-type", bpo::value<bound_type_string_t>(), "wall bound type [BOUND_DIRICHLET, BOUND_NEUMANN]")
        ("wall-4-bound-voltage", bpo::value<voltage_t>(), "wall voltage V")
        ("wall-5-bound-type", bpo::value<bound_type_string_t>(), "wall bound type [BOUND_DIRICHLET, BOUND_NEUMANN]")
        ("wall-5-bound-voltage", bpo::value<voltage_t>(), "wall voltage V")
        ("wall-6-bound-type", bpo::value<bound_type_string_t>(), "wall bound type [BOUND_DIRICHLET, BOUND_NEUMANN]")
        ("wall-6-bound-voltage", bpo::value<voltage_t>(), "wall voltage V")

        ("dxf-filename", bpo::value<std::string>(), "dxf filename relative to the executable")
        ("dxfsolid-layername", bpo::value<std::vector<std::string>>(), "DSFSolid layername")
        ("dxfsolid-scalefactor", bpo::value<std::vector<dxf_scale_factor_t>>(), "DSFSolid scale factor")
        ("dxfsolid-bound-type", bpo::value<std::vector<bound_type_string_t>>(), "DSFSolid boundtype [BOUND_DIRICHLET, BOUND_NEUMANN]")
        ("dxfsolid-bound-voltage", bpo::value<std::vector<voltage_t>>(), "DSFSolid bound voltage")
        
        ("bfield-mode", bpo::value<std::string>(), "Bfield Mode (MODE_3D, MODE_CYL...)")
        ("bfield-filename", bpo::value<std::string>(), "Bfield filename, path relative to the executable")
        ("bfield-xscale", bpo::value<bfield_scale_factor_t>(), "DSFSolid scale xfactor")
        ("bfield-fscale", bpo::value<bfield_scale_factor_t>(), "DSFSolid scale ffactor")        
        ("bfield-sel-f1", bpo::value<bool>(), "bfield-sel-f1")
        ("bfield-sel-f2", bpo::value<bool>(), "bfield-sel-f2")
        ("bfield-sel-f3", bpo::value<bool>(), "bfield-sel-f3")
        ("bfield-translate-x", bpo::value<geometry_value_t>(), "bfield translate x")
        ("bfield-translate-y", bpo::value<geometry_value_t>(), "bfield translate y")
        ("bfield-translate-z", bpo::value<geometry_value_t>(), "bfield translate z")



    ;
    return config_file_options_o;

}

bpo::options_description command_line_options_m() 
{
    bpo::options_description command_line_options_o("Command line options");
    command_line_options_o.add_options()
        ("help", "print help message")
        ("config-file", bpo::value<std::string>(), "configuration file with all the simulation parameters, path relative to executable");

    return command_line_options_o;
}

geom_mode_e geometry_mode_m(std::string gm_string_o)
{
    geom_mode_e geometry_mode;
    if (gm_string_o == "MODE_3D")
        geometry_mode = MODE_3D;
    else if (gm_string_o == "MODE_CYL")
        geometry_mode = MODE_CYL;
    else if (gm_string_o == "MODE_2D")
        geometry_mode = MODE_2D;
    else if (gm_string_o == "MODE_1D")
        geometry_mode = MODE_1D;
    else
        throw new Error("Config-file error: geometry-mode not implemented");
    return geometry_mode;
}

Geometry* geometry_m(bpo::variables_map &vm_o) 
{

    

    const std::string &gm_string_o = vm_o["geometry-mode"].as<std::string>();
    geom_mode_e geometry_mode = geometry_mode_m(gm_string_o);

    const origin_t origin_x = vm_o["origin-x"].as<origin_t>();
    const origin_t origin_y = vm_o["origin-y"].as<origin_t>();
    const origin_t origin_z = vm_o["origin-z"].as<origin_t>();
    
    const mesh_cell_size_h_t mesh_cell_size_h = vm_o["mesh-cell-size-h"].as<mesh_cell_size_h_t>();
    const geometry_value_t geometry_start_x = vm_o["geometry-start-x"].as<geometry_value_t>();
    const geometry_value_t geometry_start_y = vm_o["geometry-start-y"].as<geometry_value_t>();
    const geometry_value_t geometry_start_z = vm_o["geometry-start-z"].as<geometry_value_t>();
    const geometry_value_t geometry_size_x = vm_o["geometry-size-x"].as<geometry_value_t>() - geometry_start_x;
    const geometry_value_t geometry_size_y = vm_o["geometry-size-y"].as<geometry_value_t>() - geometry_start_y;
    const geometry_value_t geometry_size_z = vm_o["geometry-size-z"].as<geometry_value_t>() - geometry_start_z;
    
    
    const Vec3D origin_o(origin_x, origin_y, origin_z );

    const Int3D size_3D_o(
        (int)floor(geometry_size_x/mesh_cell_size_h)+1,
        (int)floor(geometry_size_y/mesh_cell_size_h)+1,
        (int)floor(geometry_size_z/mesh_cell_size_h)+1
        );
    
    Geometry *geometry_op = new Geometry(geometry_mode, size_3D_o, origin_o, mesh_cell_size_h);
    
    return geometry_op;
}

bound_type_t bound_type_m(const bound_type_string_t &s) {
    if(s == "BOUND_DIRICHLET")
        return BOUND_DIRICHLET;
    else if (s == "BOUND_NEUMANN")
        return BOUND_NEUMANN;
    else
        throw new Error("Config-file error: Bound type unknown ");
     
}

void wall_bound_m(Geometry &geometry_o, bpo::variables_map &vm_o)
{
    //ITERATING through [wall-1-bound-type... wall-6-bound-type]
    
    std::string label_begin_o = "wall-";
    std::string label_type_end_o = "-bound-type";
    std::string label_voltage_end_o = "-bound-voltage";
    for(int i = 1; i <= 6; i++) {
        std::string label_type_o    = label_begin_o + std::to_string(i) + label_type_end_o;
        std::string label_voltage_o = label_begin_o + std::to_string(i) + label_voltage_end_o;
        if(vm_o.count(label_type_o)) {
            const bound_type_string_t bound_type_string_o = vm_o[label_type_o].as<bound_type_string_t>();
            const voltage_t bound_voltage = vm_o[label_voltage_o].as<voltage_t>();    
            const bound_type_t bound_type = bound_type_m(bound_type_string_o);

            geometry_o.set_boundary( i, Bound(bound_type, bound_voltage) );
        }
    }

}


void add_solid_m(Geometry &geometry_o,
                 int progressive, 
                 MyDXFFile *dxffile_op,
                 const std::string &layer_name_o, 
                 dxf_scale_factor_t dxf_scale_factor,
                 bound_type_t bound_type,
                 voltage_t  bound_V  ) 
{
    
    DXFSolid *s = new DXFSolid( dxffile_op, layer_name_o);
    s->scale( (double) dxf_scale_factor );
    geometry_o.set_solid(progressive, s);
    geometry_o.set_boundary(progressive, 
                Bound(bound_type, (double) bound_V));

}



void dxfsolids_m(Geometry &geometry_o, bpo::variables_map &vm_o) 
{

    const std::string &dxf_filename_o = vm_o["dxf-filename"].as<std::string>();

    MyDXFFile *dxffile_op = new MyDXFFile;
    dxffile_op->set_warning_level( 1 );
    dxffile_op->read( dxf_filename_o);



    const std::vector<std::string> &layer_names_o = 
                vm_o["dxfsolid-layername"].as<std::vector<std::string>>();

    const std::vector<dxf_scale_factor_t> &dxf_scale_factors_o = 
                vm_o["dxfsolid-scalefactor"].as<std::vector<dxf_scale_factor_t>>();

    const std::vector<voltage_t> &dxf_bound_voltages_o = 
                vm_o["dxfsolid-bound-voltage"].as<std::vector<voltage_t>>();

    const std::vector<std::string> &dxf_bound_type_strings_o = 
                vm_o["dxfsolid-bound-type"].as<std::vector<std::string>>();
  
    
    for(unsigned int i=0; i< layer_names_o.size(); i++) {

        const int progressive = i+ 7; // 0..6 are reserved
        const bound_type_t bound_type = bound_type_m(dxf_bound_type_strings_o[i]);

        std::cout<<i<<": "<<layer_names_o[i]<<" "
                    <<dxf_scale_factors_o[i]<<" "
                    <<bound_type<<": "
                    <<dxf_bound_voltages_o[i]<<std::endl;

        add_solid_m(geometry_o,
                    progressive, 
                    dxffile_op,
                    layer_names_o[i], 
                    dxf_scale_factors_o[i],
                    bound_type,
                    dxf_bound_voltages_o[i]) ;


    }

}



MeshVectorField* bfield_m(Geometry &geometry_o, bpo::variables_map &vm_o)
{
    const std::string &bfield_mode_string_o = vm_o["bfield-mode"].as<std::string>();
    const std::string &bfield_filename_o = vm_o["bfield-filename"].as<std::string>();
    const bfield_scale_factor_t bfield_xscale = vm_o["bfield-xscale"].as<bfield_scale_factor_t>();
    const bfield_scale_factor_t bfield_fscale = vm_o["bfield-fscale"].as<bfield_scale_factor_t>();
    const bool sel_f1 = vm_o["bfield-sel-f1"].as<bool>();
    const bool sel_f2 = vm_o["bfield-sel-f2"].as<bool>();    
    const bool sel_f3 = vm_o["bfield-sel-f3"].as<bool>();    
    const geometry_value_t bfield_translate_x = vm_o["bfield-translate-x"].as<geometry_value_t>();
    const geometry_value_t bfield_translate_y = vm_o["bfield-translate-y"].as<geometry_value_t>();
    const geometry_value_t bfield_translate_z = vm_o["bfield-translate-z"].as<geometry_value_t>();

    geom_mode_e bfield_mode = geometry_mode_m(bfield_mode_string_o);

    bool fsel[3] = {sel_f1, sel_f2, sel_f3};
    MeshVectorField *bfield = 
        new MeshVectorField(
            bfield_mode, 
            fsel, 
            (double) bfield_xscale, 
            (double) bfield_fscale, 
            bfield_filename_o);
    bfield->translate( Vec3D(
                        bfield_translate_x,
                        bfield_translate_y,
                        bfield_translate_z) );

    return bfield;
}




run_parameters_t* run_parameters_m(int argc, char *argv[]) 
{

    bpo::options_description command_line_options_o("Command line options");
    command_line_options_o.add_options()
        ("help", "print help message")
        ("config-file", bpo::value<std::string>(), "configuration file, path relative to executable");


    bpo::variables_map vm_cmdl_o;
    store(parse_command_line(argc, argv, command_line_options_o), vm_cmdl_o);

    bpo::options_description config_file_options_o = config_file_options_m();

    if (vm_cmdl_o.count("config-file")) {
        run_parameters_t* vm_op;

        const char* config_filename = vm_cmdl_o["config-file"].as<std::string>().c_str();
        store(parse_config_file(config_filename, config_file_options_o, true), *vm_op);
        bpo::notify(*vm_op);
        return vm_op;
    }

    std::cout << command_line_options_o <<std::endl;
    std::cout << config_file_options_o;
    return nullptr;

    

}