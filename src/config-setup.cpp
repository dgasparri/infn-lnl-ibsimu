//Compile with 
//g++ -o config config.cpp -lboost_program_options


/*

g++ -g `pkg-config --cflags ibsimu-1.0.6dev` -c -o config.o config.cpp -lboost_program_options
*/

#include "config.h"

namespace ic_setup = ibsimu_client::setup

bound_type_t ic_setup::bound_type_m(const bound_type_string_t &s) 
{
    if(s == "BOUND_DIRICHLET")
        return BOUND_DIRICHLET;
    if (s == "BOUND_NEUMANN")
        return BOUND_NEUMANN;
    
    throw new Error("Config-file error: Bound type unknown ");
     
}



geom_mode_e ic_setup::geometry_mode_m(std::string gm_string_o)
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

Geometry* ic_setup::geometry_m(bpo::variables_map &vm_o) 
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


void ic_setup::wall_bound_m(Geometry &geometry_o, bpo::variables_map &vm_o)
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


void ic_setup::add_solid_m(Geometry &geometry_o,
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



void ic_setup::dxfsolids_m(Geometry &geometry_o, bpo::variables_map &vm_o) 
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



MeshVectorField* ic_setup::bfield_m(Geometry &geometry_o, bpo::variables_map &vm_o)
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


int ic_setup::num_cores_m(bpo::variables_map &vm_o, int default_v)
{
    if(vm_o.count("ibsimu-cores"))
        return vm_o["ibsimu-cores"].as<int>();    
    else
        return default_v;
}

message_type_e message_threshold_m(bpo::variables_map &vm_o, message_type_e default_v)
{
    if(vm_o.count("ibsimu-message-threshold")) {
        const std::string &s = vm_o["ibsimu-message-threshold"].as<std::string>();
        if(s == "MSG_VERBOSE") 
            return MSG_VERBOSE;
    }
    return default_v;

}


physics_parameters_t* physics_parameters_m(bpo::variables_map &vm_o)
{
    double Te = vm_o["electron-temperature-Te"].as<double>();
    double Up = vm_o["plasma-potential-Up"].as<double>();
    voltage_t gndV = vm_o["ground-voltage"].as<voltage_t>();
    double plasma_init_x  = vm_o["plasma-init-x"].as<double>();
    double plasma_init_y  = vm_o["plasma-init-y"].as<double>();
    double plasma_init_z  = vm_o["plasma-init-z"].as<double>();
    physics_parameters_t* phypars = new physics_parameters_t;
    phypars->electron_temperature_Te = Te;
    phypars->plasma_potential_Up = Up;
    phypars->ground_V = gndV;
    phypars->plasma_init_x  = plasma_init_x ;
    phypars->plasma_init_y  = plasma_init_y ;
    phypars->plasma_init_z  = plasma_init_z ;
    return phypars;
}


bpo::options_description options_commandline_simulation_m()
{
    bpo::options_description command_line_options_o("Command line options");
    command_line_options_o.add_options()
        ("help", "print help message")
        ("run", bpo::value<std::string>()->default_value(""), "directory containing the config.ini of the file. Output files will be written in that directory ")
        ("config-file", bpo::value<std::string>()->default_value(""), "configuration file, path relative to executable")
        ("run-output", bpo::value<std::string>()->default_value("OUT_NORMAL"), "output files generated in the run [OUT_NORMAL (default, only final files), OUT_EVOLUTION (every 10 loops and last), OUT_BEGIN (first 3 loops and final), OUT_VERBOSE (first 3, every 10 loops and last)]")
        ("loop-output", bpo::value<std::string>()->default_value("LOOP_END"), "output files generated in the loop [LOOP_END (default, only at the end of the loop), LOOP_VERBOSE (every step of the loop)]")
        ;
    return command_line_options_o;

}

parameters_commandline_simulation_o* parameters_commandline_simulation_m(int argc, char *argv[])
{

    command_line_options_o = options_commandline_simulation_m();
    
    bpo::variables_map vm_cmdl_o;
    store(parse_command_line(argc, argv, command_line_options_o), vm_cmdl_o);

    parameters_commandline_simulation_t *options_op =
            new parameters_commandline_simulation_t; 
    
    options_op->run_o = vm_cmdl_o["run"].as<std::string>();
    options_op->config_filename_o = vm_cmdl_o["config-file"].as<std::string>();
    std::string run_output_o = vm_cmdl_o["run-output"].as<std::string>();
    if(run_output_o == "OUT_NORMAL") 
        options_op->run_output = OUT_NORMAL; 
    if(run_output_o == "OUT_EVOLUTION") 
        options_op->run_output = OUT_EVOLUTION; 
    if(run_output_o == "OUT_BEGIN") 
        options_op->run_output = OUT_BEGIN; 
    if(run_output_o == "OUT_VERBOSE") 
        options_op->run_output = OUT_VERBOSE; 
    
    std::string loop_output_o = vm_cmdl_o["loop-output"].as<std::string>();
    if(loop_output_o == "LOOP_END") 
        options_op->loop_output = LOOP_END;
    if(loop_output_o == "LOOP_VERBOSE") 
        options_op->loop_output = LOOP_VERBOSE;

    return options_op;
} 

simulation_parameters_t* simulation_parameters_m(int argc, char *argv[]) 
{

    bpo::options_description config_file_options_o = config_file_options_m();

    if (vm_cmdl_o.count("config-file")) {
        run_parameters_t* vm_op = new run_parameters_t;

        const char* config_filename = vm_cmdl_o["config-file"].as<std::string>().c_str();
        store(parse_config_file(config_filename, config_file_options_o, true), *vm_op);
        bpo::notify(*vm_op);
        return vm_op;
    }

    std::cout << command_line_options_o <<std::endl;
    std::cout << config_file_options_o;
    return nullptr;

    

}




analysis_parameters_t* analysis_parameters_m(int argc, char *argv[]) 
{

    bpo::options_description command_line_options_o("Command line options");
    command_line_options_o.add_options()
        ("help", "print help message")
        ("config-file", bpo::value<std::string>(), "configuration file, path relative to executable")
        ("epot-file", bpo::value<std::string>(), "epot file, path relative to executable")
        ("pdb-file", bpo::value<std::string>(), "pdb file, path relative to executable")
        ("bfield-file", bpo::value<std::string>(), "bfield file, path relative to executable")

    ;


    bpo::variables_map vm_cmdl_o;
    store(parse_command_line(argc, argv, command_line_options_o), vm_cmdl_o);

    bpo::options_description config_file_options_o = config_file_options_m();

    bool all_files = true;
    if (!vm_cmdl_o.count("config-file")) 
        all_files = false;
    if (!vm_cmdl_o.count("epot-file")) 
        all_files = false;
    if (!vm_cmdl_o.count("pdb-file")) 
        all_files = false;

    if (!all_files) {
        std::cout<<"Usage:"<<std::endl<<" ./analysis --config-file config.ini --epot-file epot.dat --pdb-file pdb.dat"<<std::endl;
        std::cout << command_line_options_o <<std::endl;
        std::cout << config_file_options_o;
        return nullptr;
    }
    
    analysis_parameters_t *params = new analysis_parameters_t;

    bpo::variables_map vm_o;
    std::string epot_filename_o;
    std::string pdb_filename_o;


    run_parameters_t* vm_op = new run_parameters_t;
    const char* config_filename = vm_cmdl_o["config-file"].as<std::string>().c_str();
    store(parse_config_file(config_filename, config_file_options_o, true), *vm_op);
    bpo::notify(*vm_op);
    params->vm_op = vm_op;

    params->epot_filename_o = vm_cmdl_o["epot-file"].as<std::string>();
    params->pdb_filename_o = vm_cmdl_o["pdb-file"].as<std::string>();
   
    
    return params;


    

}
