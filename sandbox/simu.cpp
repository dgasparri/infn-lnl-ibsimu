#include <fstream>
#include <iomanip>
#include <limits>
#include "epot_bicgstabsolver.hpp"
#include "meshvectorfield.hpp"
#include "dxf_solid.hpp"
#include "stl_solid.hpp"
#include "stlfile.hpp"
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



using namespace std;


const std::string bfieldfn = "Radis_final.dat";


const double Vbias = 19.0e3;
const double Vpuller = 10.0e3;
const double Vdump = 6.0e3;
const double Veinzel = 20e3;


class ForcedPot : public CallbackFunctorB_V {

public:

    ForcedPot() {}
    ~ForcedPot() {}

    virtual bool operator()( const Vec3D &x ) const {
        return( x[2] < 0.2e-3 && x[0]*x[0]+x[1]*x[1] > 6e-3*6e-3 );
    }
};


void simu( int argc, char **argv )
{
    double start = -3.0e-3;
    double h = 0.4e-3;
    double sizereq[3] = { 50.0e-3,
                          50.0e-3, 
                          125.0e-3-start };
    Int3D meshsize( (int)floor(sizereq[0]/h)+1,
                    (int)floor(sizereq[1]/h)+1,
                    (int)floor(sizereq[2]/h)+1 );
    Vec3D origo( -25.0e-3, -25.0e-3, start );
    Geometry geom( MODE_3D, meshsize, origo, h );

    double dx = 0.00405834+0.00505786;
    double dy = -0.00808772+0.0140138;
    double angle = atan2( dx, dy ) + 0.5*M_PI;
    std::cout << "angle = " << angle << "\n";
    Transformation T;
    T.translate( Vec3D( 0, -197, 0 ) );
    T.scale( Vec3D( 1.0e-3, -1.0e-3, 1.0e-3 ) );
    T.rotate_x( 0.5*M_PI );
    T.translate( Vec3D( 0, 0, h/50.0 ) );
    T.rotate_z( angle );

    STLFile *fplasma = new STLFile( "plasma.stl" );
    STLFile *fpuller = new STLFile( "puller.stl" );
    STLFile *fedump = new STLFile( "edump.stl" );
    STLFile *fedump2 = new STLFile( "edump2.stl" );
    STLFile *fmaa1a = new STLFile( "maa1a.stl" );
    STLFile *fmaa1b = new STLFile( "maa1b.stl" );
    STLFile *fmaa2a = new STLFile( "maa2a.stl" );
    STLFile *fmaa2b = new STLFile( "maa2b.stl" );
    STLFile *feinzela = new STLFile( "einzela.stl" );
    STLFile *feinzelb = new STLFile( "einzelb.stl" );

    STLSolid *plasma = new STLSolid;
    plasma->set_transformation( T );
    plasma->add_stl_file( fplasma );
    geom.set_solid( 7, plasma );

    STLSolid *puller = new STLSolid;
    puller->set_transformation( T );
    puller->add_stl_file( fpuller );
    geom.set_solid( 8, puller );

    STLSolid *edump = new STLSolid;
    edump->set_transformation( T );
    edump->add_stl_file( fedump );
    edump->add_stl_file( fedump2 );
    geom.set_solid( 9, edump );

    STLSolid *maa1 = new STLSolid;
    maa1->set_transformation( T );
    maa1->add_stl_file( fmaa1a );
    maa1->add_stl_file( fmaa1b );
    geom.set_solid( 10, maa1 );

    STLSolid *einzel = new STLSolid;
    einzel->set_transformation( T );
    einzel->add_stl_file( feinzela );
    einzel->add_stl_file( feinzelb );
    geom.set_solid( 11, einzel );

    STLSolid *maa2 = new STLSolid;
    maa2->set_transformation( T );
    maa2->add_stl_file( fmaa2a );
    maa2->add_stl_file( fmaa2b );
    geom.set_solid( 12, maa2 );

    geom.set_boundary(  1,  Bound(BOUND_NEUMANN,     0.0) );
    geom.set_boundary(  2,  Bound(BOUND_NEUMANN,     0.0) );
    geom.set_boundary(  3,  Bound(BOUND_NEUMANN,     0.0) );
    geom.set_boundary(  4,  Bound(BOUND_NEUMANN,     0.0) );
    geom.set_boundary(  5,  Bound(BOUND_DIRICHLET,   0.0) );
    geom.set_boundary(  6,  Bound(BOUND_DIRICHLET,  Vbias) );

    geom.set_boundary(  7,  Bound(BOUND_DIRICHLET,   0.0) );
    geom.set_boundary(  8,  Bound(BOUND_DIRICHLET,  Vpuller) );
    geom.set_boundary(  9,  Bound(BOUND_DIRICHLET,  Vdump) );

    geom.set_boundary( 10,  Bound(BOUND_DIRICHLET,  Vbias) );
    geom.set_boundary( 11,  Bound(BOUND_DIRICHLET,  (Vbias+Veinzel)) );
    geom.set_boundary( 12,  Bound(BOUND_DIRICHLET,  Vbias) );
    geom.build_mesh();
    geom.build_surface();

    EpotBiCGSTABSolver solver( geom );
    InitialPlasma initp( AXIS_Z, 0.2e-3 );
    solver.set_nsimp_initial_plasma( &initp );
    ForcedPot force;
    solver.set_forced_potential_volume( 0.0, &force );
    solver.set_gnewton( true );

    EpotField epot( geom );
    MeshScalarField scharge( geom );
    MeshScalarField scharge_ave( geom );

    // Define magnetic field
    bool fout[3] = {true, true, true};
    MeshVectorField bfield( MODE_3D, fout, 1.0e-3, 1.0, bfieldfn );
    field_extrpl_e bfldextrpl[6] = { FIELD_ZERO, FIELD_ZERO, 
                                     FIELD_ZERO, FIELD_ZERO, 
                                     FIELD_ZERO, FIELD_ZERO };
    bfield.set_extrapolation( bfldextrpl );

    EpotEfield efield( epot );
    field_extrpl_e efldextrpl[6] = { FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE, 
				     FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE, 
				     FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE };
    efield.set_extrapolation( efldextrpl );

    ParticleDataBase3D pdb( geom );
    pdb.set_max_steps( 1000 );
    bool pmirror[6] = { false, false, false, false, false, false };
    pdb.set_mirror( pmirror );

    // Suppress effects of magnetic field at volume where phi<1000. This is needed
    // because of erroneous field values close to the plasma electrode magnetic
    // steel (bug in Radia-3D).
    NPlasmaBfieldSuppression psup( epot, 1000.0 );
    pdb.set_bfield_suppression( &psup );

    double rho_h, rho_tot;
    for( size_t i = 0; i < 15; i++ ) {
	
	if( i == 1 ) {
	    std::vector<double> Ei, rhoi;
	    Ei.push_back( 2.0 );
	    rhoi.push_back( 0.5*rho_h );
	    double rhop = rho_tot - rho_h*0.5;
            solver.set_nsimp_plasma( rhop, 10.0, rhoi, Ei );
        }
	
	solver.solve( epot, scharge_ave );
	efield.recalculate();

	double J = 35.0;
        pdb.clear(); 
	pdb.add_cylindrical_beam_with_energy( 50000, J, -1.0, 1.0, 
					      2.0, 0.0, 1.0, 
					      Vec3D(0,0,start),
					      Vec3D(1,0,0), 
					      Vec3D(0,1,0),
					      8e-3 );
	rho_h = pdb.get_rhosum();
	pdb.add_cylindrical_beam_with_energy( 50000, J*100.0, -1.0, 1.0/1836.15, 
					      2.0, 0.0, 1.0, 
					      Vec3D(0,0,start),
					      Vec3D(1,0,0), 
					      Vec3D(0,1,0),
					      8e-3 );
        pdb.iterate_trajectories( scharge, efield, bfield );
	rho_tot = pdb.get_rhosum();

        // Space charge averaging
	if( i == 0 ) {
	    scharge_ave = scharge;
	} else {
	    double coef = 0.3;
	    scharge *= coef;
	    scharge_ave += scharge;
	    scharge_ave *= (1.0/(1.0+coef));
	}
    }
    
    geom.save( "geom.dat" );
    epot.save( "epot.dat" );
    pdb.save( "pdb.dat" );
}


int main( int argc, char **argv )
{
    try {
	//ibsimu.set_message_output( "ibsimu.txt" );
        ibsimu.set_message_threshold( MSG_VERBOSE, 1 );
	ibsimu.set_thread_count( 4 );
	simu( argc, argv );
    } catch( Error e ) {
	e.print_error_message( ibsimu.message( 0 ) );
        exit( 1 );
    }

    return( 0 );
}
