#include <fstream>
#include <iomanip>
#include <limits>
#include "meshvectorfield.hpp"
#include "dxf_solid.hpp"
#include "mydxffile.hpp"
#include "gtkplotter.hpp"
#include "geomplotter.hpp"
#include "geometry.hpp"
#include "func_solid.hpp"
#include "epot_efield.hpp"
#include "error.hpp"
#include "ibsimu.hpp"
#include "trajectorydiagnostics.hpp"
#include "particledatabase.hpp"
#include "particlediagplotter.hpp"

#include "config.h"





void analysis(Geometry &geometry_o, std::string &epot_filename_o, std::string &pdb_filename_o, std::string bfield_filename_o)
{

    std::ifstream is_epot(epot_filename_o);
    if( !is_epot.good() )
	    throw( Error( ERROR_LOCATION, (std::string)"couldn\'t open file \'" + epot_filename_o + "\'" ) );
    
    EpotField epot( is_epot, geometry_o );
    is_epot.close();

    EpotEfield efield( epot );
    field_extrpl_e efldextrpl[6] = { FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE, 
				     FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE, 
				     FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE };
    efield.set_extrapolation( efldextrpl );

    std::ifstream is_pdb(pdb_filename_o);
    if( !is_pdb.good() )
	throw( Error( ERROR_LOCATION, (std::string)"couldn\'t open file \'" + pdb_filename_o + "\'" ) );
    ParticleDataBase3D pdb( is_pdb, geometry_o );
    is_pdb.close();

    VectorField *bfield = NULL;
    if( false ) {
    	bool fout[3] = {true, true, true};
	    MeshVectorField *mesh_bfield = new MeshVectorField( MODE_3D, fout, 1.0e-3, 1.0, bfield_filename_o );
	    field_extrpl_e bfldextrpl[6] = { FIELD_ZERO, FIELD_ZERO, 
					 FIELD_ZERO, FIELD_ZERO, 
					 FIELD_ZERO, FIELD_ZERO };
	    mesh_bfield->set_extrapolation( bfldextrpl );
	    bfield = mesh_bfield;
    }

    MeshScalarField tdens( geometry_o );
    pdb.build_trajectory_density_field( tdens );

    int temp_params = 0;
    GTKPlotter plotter( &temp_params, nullptr );
    plotter.set_geometry( &geometry_o );
    plotter.set_epot( &epot );
    plotter.set_efield( &efield );
    //if( bfield )
	//plotter.set_bfield( bfield );
    plotter.set_trajdens( &tdens );
    plotter.set_particledatabase( &pdb );
    plotter.new_geometry_plot_window();
    plotter.run();
}





int main(int argc, char *argv[]) 
{

    analysis_parameters_t *analysis_parameters_op;
    analysis_parameters_op = analysis_parameters_m(argc, argv);
    if(! analysis_parameters_op)
        return 0; //No config file with run parameters provided

    run_parameters_t *run_parameters_op = analysis_parameters_op->vm_op;

    try {
    	ibsimu.set_message_threshold( message_threshold_m(*run_parameters_op, MSG_VERBOSE), 1 );
	    ibsimu.set_thread_count( num_cores_m(*run_parameters_op, 2) );

        Geometry* geometry_op = geometry_m(*run_parameters_op); 

        wall_bound_m(*geometry_op, *run_parameters_op);
        dxfsolids_m(*geometry_op, *run_parameters_op);

        
        /*
        MeshVectorField* bfield_op = bfield_m(*geometry_op, *run_parameters_op);

        geometry_op->build_mesh();

        physics_parameters_t phy_pars = physics_parameters_m(*run_parameters_op);

        simulation(*geometry_op, *bfield_op, phy_pars);
        */

    	
    } catch( Error e ) {
	    e.print_error_message( ibsimu.message(0) );
    }

    

}
