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

#include "config.cpp"




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

/*
void simu( int *argc, char ***argv )
{
    double q = 6;
    double m = 15;
    double B0 = 0.9;
    double r_aperture = 4.0e-3;
    double vz = sqrt(-2.0*q*CHARGE_E*Vgnd/(m*MASS_U));
    double Erms = q*CHARGE_E*B0*r_aperture*r_aperture/(8*m*MASS_U*vz);
    ibsimu.message(1) << "Erms = "<< Erms << " m rad\n";

    Vec3D origin( -5e-3, 0, 0 );
    Vec3D sizereq( 505e-3, 82e-3, 0 );
    Int3D size( floor(sizereq[0]/h)+1,
		floor(sizereq[1]/h)+1,
		1 );
    Geometry geom( MODE_CYL, size, origin, h );

    geom.set_boundary( 1, Bound(BOUND_DIRICHLET, Vplasma) ); // xmin
    geom.set_boundary( 2, Bound(BOUND_NEUMANN, 0.0) ); // xmax
    geom.set_boundary( 3, Bound(BOUND_NEUMANN, 0.0) ); // rmin
    geom.set_boundary( 4, Bound(BOUND_NEUMANN, 0.0) ); // rmax


    MyDXFFile *dxffile = new MyDXFFile;
    dxffile->set_warning_level( 1 );
    dxffile->read( "geom.dxf" );
    
    
    DXFSolid *s_extraction = new DXFSolid( dxffile, "extraction_line" );
    s_extraction->scale( 1e-3 );
    geom.set_solid( EXTRACTION, s_extraction );
    geom.set_boundary(  EXTRACTION, Bound(BOUND_DIRICHLET, Vextraction) );

    DXFSolid *s_puller = new DXFSolid( dxffile, "puller_line" );
    s_puller->scale( 1e-3 );
    geom.set_solid( PULLER, s_puller );
    geom.set_boundary(  PULLER, Bound(BOUND_DIRICHLET, Vpuller) );

    //Debug, da https://sourceforge.net/p/ibsimu/code/ci/master/tree/src/dxf_solid.cpp#l57


    DXFSolid *s_focus = new DXFSolid( dxffile,  "focus_line" );
    s_focus->scale( 1e-3 );
    geom.set_solid( FOCUS, s_focus );
    geom.set_boundary(  FOCUS, Bound(BOUND_DIRICHLET, Vfocus) );
   
    DXFSolid *s_ground = new DXFSolid( dxffile, "ground_line" );
    s_ground->scale( 1e-3 );
    geom.set_solid( GROUND, s_ground );
    geom.set_boundary(  GROUND, Bound(BOUND_DIRICHLET, Vgnd) );

    




    geom.build_mesh();
    //exit(0);
    std::cout<<"After mesh"<<std::endl;

    EpotField epot( geom );
    MeshScalarField scharge( geom );
    MeshScalarField scharge_ave( geom );
    bool fsel[3] = {true, true, false};
    MeshVectorField bfield( MODE_CYL, fsel, 1.0e-3, 1.0, "bfield.txt" );
    bfield.translate( Vec3D(-19.3e-3,0,0) );

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
	add_beam( geom, pdb, 1, 15, Jtotal, 0.050 );
	add_beam( geom, pdb, 2, 15, Jtotal, 0.100 );
	add_beam( geom, pdb, 3, 15, Jtotal, 0.200 );
	add_beam( geom, pdb, 4, 15, Jtotal, 0.310 );
	add_beam( geom, pdb, 5, 15, Jtotal, 0.250 );
	add_beam( geom, pdb, 6, 15, Jtotal, 0.085 );
	add_beam( geom, pdb, 7, 15, Jtotal, 0.005 );
	//add_beam( geom, pdb, 6, 15, Jtotal, 1.0 );
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
*/


int main(int argc, char *argv[]) 
{

    remove( "emit.txt" );

    run_parameters_t *run_parameters_op;
    run_parameters_op = run_parameters_m(argc, argv);
    if(! run_parameters_op)
        return 0; //No config file with run parameters provided

    try {
    	ibsimu.set_message_threshold( MSG_VERBOSE, 1 );
	    ibsimu.set_thread_count( 2 );

        Geometry* geometry_op = geometry_m(*run_parameters_op); 

        wall_bound_m(*geometry_op, *run_parameters_op);
        dxfsolids_m(*geometry_op, *run_parameters_op);

        MeshVectorField* bfield = bfield_m(*geometry_op, *run_parameters_op);

    	//simu( &argc, &argv );
    } catch( Error e ) {
	    e.print_error_message( ibsimu.message(0) );
    }

    

}
