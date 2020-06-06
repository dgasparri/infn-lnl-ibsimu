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

void simulation( Geometry &geometry_o, MeshVectorField &bfield_o, physics_parameters_t &phy_params_o )
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
    
    //http://ibsimu.sourceforge.net/manual_1_0_6/classInitialPlasma.html
    //Initial plasma exists at coordinates less than val in axis direction
    InitialPlasma initial_plasma_m( AXIS_X, 0.5e-3 );
    solver.set_initial_plasma( 
        phy_params_o.plasma_potential_Up, 
        &initial_plasma_m );

    ParticleDataBaseCyl pdb( geometry_o );
    bool pmirror[6] = {false, false,
		       true, false,
		       false, false};
    pdb.set_mirror( pmirror );



    for( int a = 0; a < Nrounds; a++ ) {


        geometry_o.save( std::string("geom-") + std::string(a) + std::string("-init.dat") );
        epot.save( std::string("epot-") + std::string(a) + std::string("-init.dat") );
        pdb.save( std::string("pdb-") + std::string(a) + std::string("-init.dat") );


        ibsimu.message(1) << "Major cycle " << a << "\n";
        ibsimu.message(1) << "-----------------------\n";

        if( a == 1 ) {
            double electron_charge_density_rhoe = pdb.get_rhosum();
            solver.set_pexp_plasma( 
                electron_charge_density_rhoe, 
                phy_params_o.electron_temperature_Te, 
                phy_params_o.plasma_potential_Up );
        }


        solver.solve( epot, scharge_ave );
        if( solver.get_iter() == 0 ) {
            ibsimu.message(1) << "No iterations, breaking major cycle\n";
            break;
        }

        geometry_o.save( std::string("geom-") + std::string(a) + std::string("-doposolve.dat") );
        epot.save( std::string("epot-") + std::string(a) + std::string("-doposolve.dat") );
        pdb.save( std::string("pdb-") + std::string(a) + std::string("-doposolve.dat") );

        efield.recalculate();

        geometry_o.save( std::string("geom-") + std::string(a) + std::string("-doporecalculate.dat") );
        epot.save( std::string("epot-") + std::string(a) + std::string("-doporecalculate.dat") );
        pdb.save( std::string("pdb-") + std::string(a) + std::string("-doporecalculate.dat") );

        pdb.clear();

        geometry_o.save( std::string("geom-") + std::string(a) + std::string("-dopopdbclear.dat") );
        epot.save( std::string("epot-") + std::string(a) + std::string("-dopopdbclear.dat") );
        pdb.save( std::string("pdb-") + std::string(a) + std::string("-dopopdbclear.dat") );


        add_beam( geometry_o, pdb, 1, 15, Jtotal, 0.050 );
        add_beam( geometry_o, pdb, 2, 15, Jtotal, 0.100 );
        add_beam( geometry_o, pdb, 3, 15, Jtotal, 0.200 );
        add_beam( geometry_o, pdb, 4, 15, Jtotal, 0.310 );
        add_beam( geometry_o, pdb, 5, 15, Jtotal, 0.250 );
        add_beam( geometry_o, pdb, 6, 15, Jtotal, 0.085 );
        add_beam( geometry_o, pdb, 7, 15, Jtotal, 0.005 );
        geometry_o.save( std::string("geom-") + std::string(a) + std::string("-dopoaddbean.dat") );
        epot.save( std::string("epot-") + std::string(a) + std::string("-dopoaddbean.dat") );
        pdb.save( std::string("pdb-") + std::string(a) + std::string("-dopoaddbean.dat") );

        //add_beam( geometry_o, pdb, 6, 15, Jtotal, 1.0 );
        pdb.iterate_trajectories( scharge, efield, bfield_o );

        geometry_o.save( std::string("geom-") + std::string(a) + std::string("-dopopdbtraj.dat") );
        epot.save( std::string("epot-") + std::string(a) + std::string("-dopopdbtraj.dat") );
        pdb.save( std::string("pdb-") + std::string(a) + std::string("-dopopdbtraj.dat") );

        
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
        break;
    }

    geometry_o.save( "geom.dat" );
    epot.save( "epot.dat" );
    pdb.save( "pdb.dat" );

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

        physics_parameters_t phy_pars = physics_parameters_m(*run_parameters_op);

        simulation(*geometry_op, *bfield_op, phy_pars)
    	
    } catch( Error e ) {
	    e.print_error_message( ibsimu.message(0) );
    }

    

}
