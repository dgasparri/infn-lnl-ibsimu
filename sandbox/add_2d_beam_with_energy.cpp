void ParticleDataBaseCylImp::add_2d_beam_with_energy( uint32_t N, double J, double q, double m, 
						      double E, double Tp, double Tt, 
						      double x1, double y1, double x2, double y2 )
{
    ibsimu.message( 1 ) << "Defining a cylindrical beam\n";

    _particles.reserve( _particles.size()+N );

    // Convert input parameters
    m *= MASS_U;
    q *= CHARGE_E;
    E *= CHARGE_E;
    Tp *= CHARGE_E;
    Tt *= CHARGE_E;
    double v = energy_to_velocity(E,m);
    double dvp = sqrt(Tp/m);
    double dvt = sqrt(Tt/m);

    // 0: temperature parallel to defining line (transverse to beam)
    // 1: temperature perpendicular to line (parallel to beam)
    // 2: temperature in angular direction (skew velocity)
    Random *rng;
    if( ibsimu.get_rng_type() == RNG_SOBOL )
	rng = new QRandom( 3 );
    else
	rng = new MTRandom( 3 );
    rng->set_transformation( 0, Gaussian_Transformation() );
    rng->set_transformation( 1, Gaussian_Transformation() );
    rng->set_transformation( 2, Gaussian_Transformation() ); 

    double s = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) );
    _rhosum += J/v;
    double IQ = 2.0*M_PI*J*s/N;
    double vt[3];
    ParticlePCyl x;
    x[0] = 0.0;
    Vec3D transverse( x2-x1, y2-y1, 0.0 );
    transverse /= transverse.norm2();
    Vec3D parallel( transverse[1], -transverse[0], 0.0 );

    double Isum = 0.0;
    for( uint32_t a = 0; a < N; a++ ) {
	x[1] = x1 + (x2-x1)*(a+0.5)/((double)N);
	x[3] = y1 + (y2-y1)*(a+0.5)/((double)N);

	rng->get( vt );
	double pveld = dvp*vt[1];
        double pvel = sqrt( 2.0*E/m + pveld*pveld );
	x[2] = transverse[0]*dvt*vt[0] + parallel[0]*pvel;
	x[4] = transverse[1]*dvt*vt[0] + parallel[1]*pvel;
	if( x[3] == 0.0 )
	    x[5] = 0.0;
	else
	    x[5] = dvt*vt[2]/x[3];

	add_particle( ParticleCyl( IQ*x[3], q, m, x ) );
	Isum += IQ*x[3];
    }

    delete( rng );
    ibsimu.message( 1 ) << "  Total beam current " << Isum << " A\n";
}
