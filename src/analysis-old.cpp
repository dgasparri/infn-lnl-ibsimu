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

int main()
{
    
}

/*


void analysis(Geometry &geometry_o,EpotField &epot_o,ParticleDataBase &pdb_o, MeshVectorField &bfield_o, physics_parameters_t &phy_pars)
{

    EpotEfield efield_o( epot_o );
    field_extrpl_e efldextrpl[6] = { FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE, 
                    FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE, 
                    FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE };
    efield_o.set_extrapolation( efldextrpl );


	// field_extrpl_e bfldextrpl[6] = { FIELD_ZERO, FIELD_ZERO, 
	// 				 FIELD_ZERO, FIELD_ZERO, 
	// 				 FIELD_ZERO, FIELD_ZERO };
	// bfield_o->set_extrapolation( bfldextrpl );

    MeshScalarField tdens_o( geometry_o );
    pdb_o.build_trajectory_density_field( tdens_o );

    int temp_params = 0;
    GTKPlotter plotter( &temp_params, nullptr );
    plotter.set_geometry( &geometry_o );
    plotter.set_epot( &epot_o );
    plotter.set_efield( &efield_o );
    
	plotter.set_bfield( &bfield_o );
    plotter.set_trajdens( &tdens_o );
    plotter.set_particledatabase( &pdb_o );
    plotter.new_geometry_plot_window();
    plotter.run();
}





int main(int argc, char *argv[]) 
{

    analysis_parameters_t *cmdl_parameters_op;
    cmdl_parameters_op = analysis_parameters_m(argc, argv);
    if(! cmdl_parameters_op)
        return 0; //No config file with run parameters provided

    run_parameters_t *run_parameters_op = cmdl_parameters_op->vm_op;

    try {
    	ibsimu.set_message_threshold( message_threshold_m(*run_parameters_op, MSG_VERBOSE), 1 );
	    ibsimu.set_thread_count( num_cores_m(*run_parameters_op, 2) );

        Geometry* geometry_op = geometry_m(*run_parameters_op); 

        wall_bound_m(*geometry_op, *run_parameters_op);
        dxfsolids_m(*geometry_op, *run_parameters_op);

        
        
        MeshVectorField* bfield_op = bfield_m(*geometry_op, *run_parameters_op);

        geometry_op->build_mesh();

        std::ifstream is_epot(cmdl_parameters_op->epot_filename_o);
        if( !is_epot.good() )
            throw( Error( ERROR_LOCATION, (std::string)"couldn\'t open file \'" + cmdl_parameters_op->epot_filename_o + "\'" ) );
        
        EpotField epot_o( is_epot, *geometry_op );
        is_epot.close();


        std::ifstream is_pdb(cmdl_parameters_op->pdb_filename_o);
        if( !is_pdb.good() )
            throw( Error( ERROR_LOCATION, (std::string)"couldn\'t open file \'" + cmdl_parameters_op->pdb_filename_o + "\'" ) );

        ParticleDataBase *pdb_op;
        if(false) {
            pdb_op = new ParticleDataBase3D( is_pdb, *geometry_op );
        } else if(true) {
            pdb_op = new ParticleDataBaseCyl( is_pdb, *geometry_op );
        }
        
        
        is_pdb.close();


        physics_parameters_t &phy_pars = *physics_parameters_m(*run_parameters_op);

        analysis(
            *geometry_op,
            epot_o,
            *pdb_op,
            *bfield_op, 
            phy_pars);
        

    	
    } catch( Error e ) {
	    e.print_error_message( ibsimu.message(0) );
    }

    

}
*/