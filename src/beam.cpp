#include "beam.h"



/*
    void add_2d_beams_bare_m(
                    std::vector<beam_t> beams,
                    ParticleDataBaseCyl &pdb
                        )
    {
        for(int i= 0; i < beams; i++) {
            pdb.add_2d_beam_with_energy( 
                beams[i].n_particles, // Adds a beam consisting of N particles
                beams[i].current_density, // beam current density is J (A/m^2)
                beams[i].particle_charge, // charge of beam particle is q (in multiples of e)
                beams[i].mass, //mass is m (u)
                beams[i].mean_energy, //mean energy of the beam E (eV)
                beams[i].par_temp_Tp, //parallel temperature Tp (eV)
                beams[i].trans_temp_Tt, //transverse temperature Tt (eV)
                beams[i].x1, 
                beams[i].y1,
                beams[i].x2, 
                beams[i].y2);

        }

    }
*/

ibsimu_client::beam::add_2d_beams_mt ibsimu_client::beam::add_2d_beams_helper_m(std::vector<ibsimu_client::beam::beam_t> beams) 
{
    return [beams](ParticleDataBaseCyl &pdb) 
    {
        for(beam_t beam: beams) {        
            pdb.add_2d_beam_with_energy( 
                beam.n_particles,
                beam.current_density_Am2,
                beam.particle_charge,
                beam.mass,
                beam.mean_energy,
                beam.par_temp_Tp,
                beam.trans_temp_Tt,
                beam.x1, 
                beam.y1,
                beam.x2, 
                beam.y2);

        }
    };

}




/*
void add_beam( Geometry &geom, ParticleDataBaseCyl &pdb, double q, double m, 
	       double Jtotal, double frac )
{

/*
add_beam( geometry_o, pdb, 
double q 1
double m 15
Jtotal Jtotal  
frac 0.050 
* /

    pdb.add_2d_beam_with_energy( 
        Npart*frac, // Adds a beam consisting of N particles
        Jtotal*frac, // beam current density is J (A/m^2)
        q, // charge of beam particle is q (in multiples of e)
        m, //mass is m (u)
		E0, //mean energy of the beam E (eV)
        0.0, //parallel temperature Tp (eV)
        Tt, //transverse temperature Tt (eV)
		geom.origo(0), 
        0,
		geom.origo(0), 
        r0 );
}
*/