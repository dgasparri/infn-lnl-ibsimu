#include <fstream>
#include <iostream>
#include <ibsimu.hpp>
#include <error.hpp>
#include <geometry.hpp>
#include <epot_field.hpp>
#include <epot_efield.hpp>
#include <meshvectorfield.hpp>
#include <epot_bicgstabsolver.hpp>
#include <gtkplotter.hpp>
#include <trajectorydiagnostics.hpp>

#include "config.h"




int Nrounds = 50;
double h = 0.5e-3;
double nperh = 500;
double r0 = 10e-3;
double Npart = r0/h*nperh;
double Jtotal = 30.0;
double Te = 10.0;
double E0 = 10.0;
double Tt = 1.0;
double Up = 20.0;
double sc_alpha = 0.7;


void add_beam( Geometry &geom, ParticleDataBaseCyl &pdb, double q, double m, 
	       double Jtotal, double frac )
{
    pdb.add_2d_beam_with_energy( Npart*frac, Jtotal*frac, q, m,
				 E0, 0.0, Tt,
				 geom.origo(0), 0,
				 geom.origo(0), r0 );
}

void simulation( Geometry &geometry_o, MeshVectorField &bfield_o )
{
    double q = 6;
    double m = 15;
    double B0 = 0.9;
    double r_aperture = 4.0e-3;
    double vz = sqrt(-2.0*q*CHARGE_E*Vgnd/(m*MASS_U));
    double Erms = q*CHARGE_E*B0*r_aperture*r_aperture/(8*m*MASS_U*vz);
    ibsimu.message(1) << "Erms = "<< Erms << " m rad\n";


    //exit(0);
    std::cout<<"After mesh"<<std::endl;

    EpotField epot( geometry_o );
    MeshScalarField scharge( geometry_o );
    MeshScalarField scharge_ave( geometry_o );

    EpotEfield efield( epot );
    field_extrpl_e extrapl[6] = { FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE,
				  FIELD_SYMMETRIC_POTENTIAL, FIELD_EXTRAPOLATE,
				  FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE };
    efield.set_extrapolation( extrapl );
  
    EpotBiCGSTABSolver solver( geometry_o );
    InitialPlasma init_plasma( AXIS_X, 0.5e-3 );
    solver.set_initial_plasma( Up, &init_plasma );

    ParticleDataBaseCyl pdb( geometry_o );
    bool pmirror[6] = {false, false,
		       true, false,
		       false, false};
    pdb.set_mirror( pmirror );

    for( int a = 0; a < Nrounds; a++ ) {

        ibsimu.message(1) << "Major cycle " << a << "\n";
        ibsimu.message(1) << "-----------------------\n";

        if( a == 1 ) {
            double rhoe = pdb.get_rhosum();
            solver.set_pexp_plasma( rhoe, Te, Up );
        }

        solver.solve( epot, scharge_ave );
        if( solver.get_iter() == 0 ) {
            ibsimu.message(1) << "No iterations, breaking major cycle\n";
            break;
        }
        
        efield.recalculate();

        pdb.clear();
        add_beam( geometry_o, pdb, 1, 15, Jtotal, 0.050 );
        add_beam( geometry_o, pdb, 2, 15, Jtotal, 0.100 );
        add_beam( geometry_o, pdb, 3, 15, Jtotal, 0.200 );
        add_beam( geometry_o, pdb, 4, 15, Jtotal, 0.310 );
        add_beam( geometry_o, pdb, 5, 15, Jtotal, 0.250 );
        add_beam( geometry_o, pdb, 6, 15, Jtotal, 0.085 );
        add_beam( geometry_o, pdb, 7, 15, Jtotal, 0.005 );
        //add_beam( geometry_o, pdb, 6, 15, Jtotal, 1.0 );
        pdb.iterate_trajectories( scharge, efield, bfield_o );

        
        TrajectoryDiagnosticData tdata;
        vector<trajectory_diagnostic_e> diag;
        diag.push_back( DIAG_R );
        diag.push_back( DIAG_RP );
        diag.push_back( DIAG_AP );
        diag.push_back( DIAG_CURR );
        pdb.trajectories_at_plane( tdata, AXIS_X, geom.max(0), diag );
        EmittanceConv emit( 100, 100, 
                    tdata(0).data(), tdata(1).data(), 
                    tdata(2).data(), tdata(3).data() );

        ofstream dataout( "emit.txt", ios_base::app );
        dataout << emit.alpha() << " "
            << emit.beta() << " "
            << emit.epsilon() << "\n";
        dataout.close();

        if( a == 0 ) {
            scharge_ave = scharge;
        } else {
            double sc_beta = 1.0-sc_alpha;
            uint32_t nodecount = scharge.nodecount();
            for( uint32_t b = 0; b < nodecount; b++ ) {
            scharge_ave(b) = sc_alpha*scharge(b) + sc_beta*scharge_ave(b);
            }
        }
    }

    geometry_o.save( "geom.dat" );
    epot.save( "epot.dat" );
    pdb.save( "pdb.dat" );

    MeshScalarField tdens( geometry_o );
    pdb.build_trajectory_density_field( tdens );

    GTKPlotter plotter( argc, argv );
    plotter.set_geometry( &geometry_o );
    plotter.set_epot( &epot );
    plotter.set_particledatabase( &pdb );
    plotter.set_efield( &efield );
    plotter.set_trajdens( &tdens );
    plotter.set_bfield( &bfield );
    plotter.set_scharge( &scharge );
    plotter.new_geometry_plot_window();
    plotter.run();
}



int main(int argc, char *argv[]) 
{

    remove( "emit.txt" );

    run_parameters_t *run_parameters_op;
    run_parameters_op = run_parameters_m(argc, argv);
    if(! run_parameters_op)
        return 0; //No config file with run parameters provided

    try {
    	ibsimu.set_message_threshold( message_threshold_m(*run_parameters_op, MSG_VERBOSE), 1 );
	    ibsimu.set_thread_count( num_cores_m(*run_parameters_op, 2) );

        Geometry* geometry_op = geometry_m(*run_parameters_op); 

        wall_bound_m(*geometry_op, *run_parameters_op);
        dxfsolids_m(*geometry_op, *run_parameters_op);

        MeshVectorField* bfield_op = bfield_m(*geometry_op, *run_parameters_op);

        geometry_op->build_mesh();

        simulation(*geometry_op, *bfield_op)
    	
    } catch( Error e ) {
	    e.print_error_message( ibsimu.message(0) );
    }

    

}
