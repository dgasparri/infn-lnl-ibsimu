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

bpo::options_description config_file_options_m(); 

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

run_parameters_t* run_parameters_m(int argc, char *argv[]); 
