/*

Compile in dir ibsimu/sandbox con altri file in ibsimu/src:
g++ -lboost_program_options -Wall -g `pkg-config --cflags ibsimu-1.0.6dev`    -c -o bfielddebug.o bfielddebug.cpp

Link:
g++ -std=c++17 -o bfielddebug bfielddebug.o ../src/config.o ../src/config-setup.o ../src/beam.o `pkg-config --libs ibsimu-1.0.6dev` -lboost_program_options

*/

#include <fstream>
#include <iostream>
#include <vector>
#include <functional>
//#include <numeric>


#include <ibsimu.hpp>
#include <error.hpp>
#include <geometry.hpp>
#include <epot_field.hpp>
#include <epot_efield.hpp>
#include <meshvectorfield.hpp>
#include <epot_bicgstabsolver.hpp>
#include <gtkplotter.hpp>
#include <trajectorydiagnostics.hpp>

#include "../src/datatype.h"
#include "../src/config.h"
#include "../src/config-setup.h"
#include "../src/beam.h"



namespace ic = ibsimu_client;
namespace ic_config = ibsimu_client::config;
namespace ic_setup = ibsimu_client::setup;
namespace ic_beam = ibsimu_client::beam;


int main(int argc, char *argv[]) 
{

    ic::parameters_commandline_t* cmdlp_op = 
            ic_config::parameters_commandline_m(argc, argv, true);
    
    const int buffer_len = 2500;
    char current_directory[buffer_len];
    getcwd(current_directory, buffer_len);


    if (!ic_config::clean_runpath_m(current_directory, cmdlp_op)) 
    {
        ic_config::show_help(true);
        return 0;
    }


    std::cout<<"Current directory: "<<current_directory<<std::endl;
    std::cout<<"Run Directory: "<<cmdlp_op->run_o<<std::endl;
    std::cout<<"Config filename: "<<cmdlp_op->config_filename_o<<std::endl;

    bpo::variables_map* params_op;
    try {
        params_op = ic_config::parameters_configfile_m(cmdlp_op->config_filename_o);
    } catch (boost::exception& e) {
        std::cout<<"Error in config file"<<std::endl;
        std::cout<<diagnostic_information(e)<<std::endl;
        std::cout<<"./simulation --help for help"<<std::endl;
        return 1;
    }
    

    try {
    	ibsimu.set_message_threshold( ic_config::message_threshold_m(*params_op, MSG_VERBOSE), 1 );
	    ibsimu.set_thread_count( ic_config::num_cores_m(*params_op, 2) );

        Geometry* geometry_op = ic_setup::geometry_m(*params_op); 

        ic_setup::wall_bound_m(*geometry_op, *params_op);
        ic_setup::dxfsolids_m(*geometry_op, *params_op, cmdlp_op->run_o);

        MeshVectorField* bfield_op = ic_setup::bfield_m(*geometry_op, *params_op, cmdlp_op->run_o);

        geometry_op->build_mesh();

        std::ofstream outfile("bfield-original.txt");
        for(double x = -0.175; x<= 0.331; x+=0.001) 
        for(double y = 0.0; y<=0.11; y+=0.001) {
            const Vec3D v(x, y);
            const Vec3D b = (*bfield_op)(v); 
            outfile<<v[0]<<","<<v[1]<<","<<v[2]<<","<<b[0]<<","<<b[1]<<","<<b[2]<<std::endl;
        }
        outfile.close();

        std::ofstream outfile2("bfield-translated.txt");
        for(double x = -0.175; x<= 0.331; x+=0.001) 
        for(double y = 0.0; y<=0.11; y+=0.001) {
            const Vec3D v(x+0.0003, y+0.0003);
            const Vec3D b = (*bfield_op)(v); 
            outfile2<<x<<","<<y<<","<<v[2]<<","<<b[0]<<","<<b[1]<<","<<b[2]<<std::endl;
        }
        outfile2.close();



        
        // std::cout<<"Vector v: x: "<<v[0]<<" y: "<<v[1]<<" z "<<v[2]<<std::endl;
        // std::cout<<"Vector b: x: "<<b[0]<<" y: "<<b[1]<<" z "<<b[2]<<std::endl;



    } catch( Error e ) {
	    e.print_error_message( ibsimu.message(0) );
    }


}
