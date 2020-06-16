#include <fstream>
#include <iostream>
#include <vector>
#include <functional>

#include <ibsimu.hpp>
#include <error.hpp>
#include <geometry.hpp>
#include <epot_field.hpp>
#include <epot_efield.hpp>
#include <meshvectorfield.hpp>
#include <epot_bicgstabsolver.hpp>
#include <gtkplotter.hpp>
#include <trajectorydiagnostics.hpp>

#include "datatype.h"
#include "config.h"
#include "config-setup.h"
#include "beam.h"



namespace ic = ibsimu_client;
namespace ic_config = ibsimu_client::config;
namespace ic_setup = ibsimu_client::setup;
namespace ic_beam = ibsimu_client::beam;


typedef std::function<void(int,const char*, EpotField&, ParticleDataBaseCyl&)> save_output_prototype_t;




int Nrounds = 50;
double h = 0.5e-3;
double nperh = 500;
double r0 = 10e-3;
double Npart = r0/h*nperh;
double Jtotal = 3000.0;

double E0 = 10.0;
double Tt = 1.0;

double sc_alpha = 0.7;





//loop_number == -1 -> salva a prescindere
void save_output_base_m(
    std::string run_directory_o, 
    ic::run_output_t run_output, 
    ic::loop_output_t loop_output,
    std::string geom_prefix_o,
    std::string epot_prefix_o,
    std::string pdb_prefix_o,
    Geometry* geometry_op,
    int loop_number,
    const char* stage,
    EpotField& epot_o,
    ParticleDataBaseCyl& pdb_o) 
{
    bool save = false;
    std::string stage_o = stage;
    
    switch(run_output) {
        case ic::OUT_EVOLUTION:
            if(! (loop_number % 10))
                save = true;
            break;
        case ic::OUT_BEGIN:
            if(loop_number < 3)
                save = true;
            break;
        case ic::OUT_VERBOSE:
            if(loop_number < 3 
                || ! (loop_number % 10))
                save = true;
            break;
        case ic::OUT_NORMAL:
        default:
            break;
    }


    //if LOOP_ENd saves only at the beginning of the loop
    if(loop_output==ic::LOOP_END) {
        if (stage_o != "A.init")
            save = false;
    }
    
    if(loop_number == -1) {
        save = true;
    }

    /*
        geometry_o.save( geomstr(a,"init-loop") );
        epot.save( epotstr(a,"init-loop"));
        pdb.save( pdbstr(a,"init-loop") );
    */
    std::string suffix;
    if(loop_number != -1) {
        suffix = to_string(".") + to_string(loop_number) + "." + stage;
    }
    suffix += ".dat";

    if(save) {
        geometry_op->save( run_directory_o + geom_prefix_o+suffix );
        epot_o.save( run_directory_o + epot_prefix_o + suffix);
        pdb_o.save( run_directory_o + pdb_prefix_o+suffix );

    }

}


void simulation( 
    Geometry &geometry_o, 
    MeshVectorField &bfield_o, 
    ic::physics_parameters_t &phy_params_o,
    save_output_prototype_t save_output_m,
    std::ofstream& emittance_csv_stream_o,
    ic_beam::add_2d_beams_mt add_2b_beam_m 
     )
{


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
    
    InitialPlasma *initial_plasma_op;
    if(phy_params_o.plasma_init_x)
        initial_plasma_op = new InitialPlasma( AXIS_X, phy_params_o.plasma_init_x );
    if(phy_params_o.plasma_init_y)
        initial_plasma_op = new InitialPlasma( AXIS_Y, phy_params_o.plasma_init_y );
    if(phy_params_o.plasma_init_z)
        initial_plasma_op = new InitialPlasma( AXIS_Z, phy_params_o.plasma_init_z );
    
    solver.set_initial_plasma( 
        phy_params_o.plasma_potential_Up, 
        initial_plasma_op);

    ParticleDataBaseCyl pdb( geometry_o );
    bool pmirror[6] = {false, false,
		       true, false,
		       false, false};
    pdb.set_mirror( pmirror );


    for( int a = 0; a < Nrounds; a++ ) {
        save_output_m(a,"A.init", epot, pdb);

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

        save_output_m(a,"B.aftersolver", epot, pdb);

        efield.recalculate();

        save_output_m(a,"C.afterefieldrecalculate", epot, pdb);


        pdb.clear();

        save_output_m(a,"D.afterpdbclear", epot, pdb);

        add_2b_beam_m(pdb);

        save_output_m(a,"E.afteraddbeam", epot, pdb);

        pdb.iterate_trajectories( scharge, efield, bfield_o );

        save_output_m(a,"F.aftertrajectories", epot, pdb);
        
        TrajectoryDiagnosticData tdata;
        std::vector<trajectory_diagnostic_e> diag;
        diag.push_back( DIAG_R );
        diag.push_back( DIAG_RP );
        diag.push_back( DIAG_AP );
        diag.push_back( DIAG_CURR );
        pdb.trajectories_at_plane( 
                tdata, 
                AXIS_X, 
                geometry_o.max(0), 
                diag );
        
        //Emittance statistics are from original data, not gridded data
        EmittanceConv emit( 
                100, // n Grid array n x m
                100, // m
                tdata(0).data(), //DIAG_R TrajectoryDiagnosticColumn
                tdata(1).data(), //DIAG_RP
                tdata(2).data(), //DIAG_AP
                tdata(3).data()  //DIAG_CURR I
                );


        emittance_csv_stream_o << emit.alpha() << ",";
        emittance_csv_stream_o << emit.beta() << ",";
        emittance_csv_stream_o << emit.epsilon() << "";
        emittance_csv_stream_o << std::endl;
        emittance_csv_stream_o.flush();

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
    
    save_output_m(-1,"", epot, pdb);

    MeshScalarField tdens( geometry_o );
    pdb.build_trajectory_density_field( tdens );
    
    
    delete(initial_plasma_op);

    int temp_a = 1;
    GTKPlotter plotter( &temp_a, nullptr );
    plotter.set_geometry( &geometry_o );
    plotter.set_epot( &epot );
    plotter.set_efield( &efield );
	plotter.set_bfield( &bfield_o );
    plotter.set_trajdens( &tdens );
    plotter.set_particledatabase( &pdb );
    plotter.new_geometry_plot_window();
    plotter.run();


}


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

        ic::physics_parameters_t &phy_pars = *ic_setup::physics_parameters_m(*params_op);

      
        std::string& lbd_run_o = cmdlp_op->run_o;
        ibsimu_client::run_output_t lbd_run_output = cmdlp_op->run_output;
        ibsimu_client::loop_output_t lbd_loop_output = cmdlp_op->loop_output;
        const std::string& prefix_geom_o = (*params_op)["ibsimu-file-prefix-geometry"].as<std::string>();
        const std::string& prefix_epot_o = (*params_op)["ibsimu-file-prefix-epot"]    .as<std::string>();
        const std::string& prefix_pdb_o  = (*params_op)["ibsimu-file-prefix-pdb"]     .as<std::string>();

        save_output_prototype_t save_output_lambda_m
             = [
                    lbd_run_o,
                    lbd_run_output,
                    lbd_loop_output,
                    prefix_geom_o,
                    prefix_epot_o,
                    prefix_pdb_o,
                    geometry_op
            ](
                    int loop_number,
                    const char* stage,
                    EpotField& epot_o,
                    ParticleDataBaseCyl& pdb_o
            ) {
                save_output_base_m(
                    lbd_run_o,
                    lbd_run_output,
                    lbd_loop_output,
                    prefix_geom_o,
                    prefix_epot_o,
                    prefix_pdb_o,
                    geometry_op,
                    loop_number,
                    stage,
                    epot_o,
                    pdb_o
                );
            };


        const std::string& stats_filename_o = (*params_op)["ibsimu-file-emittance-statistics"].as<std::string>();
        std::string fullpath_stats_filename_o;
        if(stats_filename_o[0]=='/')
            fullpath_stats_filename_o = stats_filename_o;
        else
            fullpath_stats_filename_o = cmdlp_op->run_o + stats_filename_o;

        std::ofstream emittance_csv(
            fullpath_stats_filename_o, 
            std::ios_base::out | std::ios_base::trunc );
        
        std::vector<ic_beam::beam_t> beams = ic_setup::beams_m(*params_op);
        ic_beam::add_2d_beams_mt add_2b_beam_m = ic_beam::add_2d_beams_helper_m(beams);


        simulation(
            *geometry_op, 
            *bfield_op, 
            phy_pars,
            save_output_lambda_m,
            emittance_csv,
            add_2b_beam_m
            );
    	
        emittance_csv.close();

    } catch( Error e ) {
	    e.print_error_message( ibsimu.message(0) );
    }


}
