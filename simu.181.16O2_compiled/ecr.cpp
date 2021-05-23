#include <fstream>
#include <iostream>
#include <ibsimu.hpp>
#include <error.hpp>
#include <geometry.hpp>
#include <dxf_solid.hpp>
#include <mydxffile.hpp>
#include <epot_field.hpp>
#include <epot_efield.hpp>
#include <meshvectorfield.hpp>
#include <epot_bicgstabsolver.hpp>
#include <gtkplotter.hpp>
#include <trajectorydiagnostics.hpp>



using namespace std;



int Nrounds = 50;
double h = 0.2e-3;
double nperh = 200;
double r0 = 15e-3;
double Npart = r0/h*nperh;
double Jtotal = 50.0;
double Te = 10.0;
double E0 = 9.0;
double Tt = 1.0;
double Up = 24020.0;
double sc_alpha = 0.7;




void simu( int *argc, char ***argv )
{
    double q = 2;
    double m = 16;
    /*double B0 = 0.9;
    double r_aperture = 4.0e-3;
    double vz = sqrt(-2.0*q*CHARGE_E*Vgnd/(m*MASS_U));
    double Erms = q*CHARGE_E*B0*r_aperture*r_aperture/(8*m*MASS_U*vz);
    ibsimu.message(1) << "Erms = "<< Erms << " m rad\n";*/

    Vec3D origin( 4e-3, 0, 0 );
    Vec3D sizereq( 530e-3, 60e-3, 0 );
    Int3D size( floor(sizereq[0]/h)+1,
		floor(sizereq[1]/h)+1,
		1 );
    Geometry geom( MODE_CYL, size, origin, h );

    MyDXFFile *dxffile = new MyDXFFile;
    dxffile->set_warning_level( 1 );
    dxffile->read( "CNAO_MOD_Luca.dxf" );




    //extraction_line
    DXFSolid *s1 = new DXFSolid( dxffile, "7" );
    s1->scale( 1e-3 );
    geom.set_solid( 7, s1 );
    //puller_line
    DXFSolid *s2 = new DXFSolid( dxffile, "8" );
    s2->scale( 1e-3 );
    geom.set_solid( 8, s2 );
    //focus_line
    DXFSolid *s3 = new DXFSolid( dxffile, "9" );
    s3->scale( 1e-3 );
    geom.set_solid( 9, s3 );
    //ground_line
    DXFSolid *s4 = new DXFSolid( dxffile, "10" );
    s4->scale( 1e-3 );
    geom.set_solid( 10, s4 );

    //geom.set_boundary( 1, Bound(BOUND_NEUMANN, 0.0) ); // xmin
    //geom.set_boundary( 2, Bound(BOUND_NEUMANN, 0.0) ); // xmax
    //geom.set_boundary( 3, Bound(BOUND_NEUMANN, 0.0) ); // rmin
    //geom.set_boundary( 4, Bound(BOUND_NEUMANN, 0.0) ); // rmax

    //Vestrazione
    geom.set_boundary(  7, Bound(BOUND_DIRICHLET, 24.0e3) );
    geom.set_boundary(  8, Bound(BOUND_DIRICHLET, -1.0e3) );
    geom.set_boundary(  9, Bound(BOUND_DIRICHLET, 700) );
    geom.set_boundary( 10, Bound(BOUND_DIRICHLET, 0.0) );

    geom.build_mesh();

    EpotField epot( geom );
    MeshScalarField scharge( geom );
    MeshScalarField scharge_ave( geom );
    bool fsel[3] = {true, true, false};
    MeshVectorField bfield( MODE_CYL, fsel, 1.0, 1.0, "../assets/matrice_campo_10mm.16.04.2021.txt" );
    

    EpotEfield efield( epot );
    field_extrpl_e extrapl[6] = { FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE,
				  FIELD_SYMMETRIC_POTENTIAL, FIELD_EXTRAPOLATE,
				  FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE };
    efield.set_extrapolation( extrapl );
  
    EpotBiCGSTABSolver solver( geom );
    InitialPlasma init_plasma( AXIS_X, 0.5e-3 );
    solver.set_initial_plasma( Up, &init_plasma );

    ParticleDataBaseCyl pdb( geom );
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

    pdb.add_2d_beam_with_energy( 
        200000, //Npart
        Jtotal, //Current density A/m^2
        2, // charge of beam particle in multipels of e
        16, // mass (u)
        E0, //Mean energy E
        0.5, //Tp
        Tt, //Tt
        0, //x 0.5e-3
        0,
        0, // x 0.5e-3, 
        r0 );

	pdb.iterate_trajectories( scharge, efield, bfield );

    
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

    geom.save( "geom.dat" );
    epot.save( "epot.dat" );
    pdb.save( "pdb.dat" );

    MeshScalarField tdens( geom );
    pdb.build_trajectory_density_field( tdens );

    GTKPlotter plotter( argc, argv );
    plotter.set_geometry( &geom );
    plotter.set_epot( &epot );
    plotter.set_particledatabase( &pdb );
    plotter.set_efield( &efield );
    plotter.set_trajdens( &tdens );
    plotter.set_bfield( &bfield );
    plotter.set_scharge( &scharge );
    plotter.new_geometry_plot_window();
    plotter.run();
}


int main( int argc, char **argv )
{
    remove( "emit.txt" );

    try {
	ibsimu.set_message_threshold( MSG_VERBOSE, 1 );
	ibsimu.set_thread_count( 2 );
	simu( &argc, &argv );
    } catch( Error e ) {
	e.print_error_message( ibsimu.message(0) );
    }

    return( 0 );
}
