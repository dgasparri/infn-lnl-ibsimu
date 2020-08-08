#pragma once

#include <ibsimu.hpp>
#include <error.hpp>
#include <geometry.hpp>
#include <epot_field.hpp>
#include <epot_efield.hpp>
#include <meshvectorfield.hpp>
#include <epot_bicgstabsolver.hpp>
#include <gtkplotter.hpp>
#include <trajectorydiagnostics.hpp>

#include <boost/program_options/options_description.hpp>


#include "datatype.h"
#include "config.h"

namespace bpo = boost::program_options;


typedef std::function<void(int,const char*, EpotField&, ParticleDataBaseCyl&)> save_output_prototype_t;

namespace ibsimu_client::output {

    void output_options_m(bpo::options_description& command_line_options_o);

    //loop_number == -1 -> salva a prescindere
    void save_output_base_m(
        std::string run_directory_o, 
        ibsimu_client::run_output_t run_output, 
        ibsimu_client::loop_output_t loop_output,
        std::string geom_prefix_o,
        std::string epot_prefix_o,
        std::string pdb_prefix_o,
        Geometry* geometry_op,
        int loop_number,
        const char* stage,
        EpotField& epot_o,
        ParticleDataBaseCyl& pdb_o);



}
