#pragma once


#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include <ibsimu.hpp>
#include <error.hpp>
#include <geometry.hpp>
#include <dxf_solid.hpp>
#include <mydxffile.hpp>
#include <epot_field.hpp>
#include <epot_efield.hpp>
#include <meshvectorfield.hpp>


#include <string>
#include <vector>
#include <unistd.h>

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




struct physics_parameters_t {
    //double electron_charge_density_rhoe;
    double electron_temperature_Te;
    double plasma_potential_Up;
    double ground_V;
    double plasma_init_x;
    double plasma_init_y;
    double plasma_init_z;
};

struct analysis_parameters_t {
    bpo::variables_map *vm_op;
    std::string epot_filename_o;
    std::string pdb_filename_o;
};

namespace ibsimu_client {
    bpo::variables_map* parameters_configfile_m(std::string configfile_o);
}


bpo::options_description command_line_options_m();
geom_mode_e geometry_mode_m(std::string gm_string_o);
Geometry* geometry_m(bpo::variables_map &vm_o) ;

bound_type_t bound_type_m(const bound_type_string_t &s);
void wall_bound_m(Geometry &geometry_o, bpo::variables_map &vm_o);
void add_solid_m(Geometry &geometry_o,
                 int progressive, 
                 MyDXFFile *dxffile_op,
                 const std::string &layer_name_o, 
                 dxf_scale_factor_t dxf_scale_factor,
                 bound_type_t bound_type,
                 voltage_t  bound_V  );


void dxfsolids_m(Geometry &geometry_o, bpo::variables_map &vm_o);
MeshVectorField* bfield_m(Geometry &geometry_o, bpo::variables_map &vm_o);

int num_cores_m(bpo::variables_map &vm_o, int default_v = 2);
message_type_e message_threshold_m(bpo::variables_map &vm_o, message_type_e default_v);

physics_parameters_t* physics_parameters_m(bpo::variables_map &vm_o);


namespace ibsimu_client::simulation {

    enum run_output_t {OUT_NORMAL, OUT_EVOLUTION, OUT_BEGIN, OUT_VERBOSE};
    enum loop_output_t {LOOP_END, LOOP_VERBOSE};

    struct parameters_commandline_t {
        std::string run_o;
        std::string config_filename_o;
        run_output_t run_output;
        loop_output_t loop_output;
    };

    
    parameters_commandline_t* parameters_commandline_m(int argc, char *argv[]);
    void show_help();

}




analysis_parameters_t* analysis_parameters_m(int argc, char *argv[]); 
