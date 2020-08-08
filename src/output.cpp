#include "output.h"


namespace ic_output = ibsimu_client::output;

void ic_output::output_options_m(
    bpo::options_description &command_line_options_o)
{
        command_line_options_o.add_options()
            ("run-output", bpo::value<std::string>()->default_value("OUT_NORMAL"), "output files generated in the run [OUT_NORMAL (default, only final files), OUT_EVOLUTION (every 10 loops and last), OUT_BEGIN (first 3 loops and final), OUT_VERBOSE (first 3, every 10 loops and last)]")
            ("loop-output", bpo::value<std::string>()->default_value("LOOP_END"), "output files generated in the loop [LOOP_END (default, only at the end of the loop), LOOP_VERBOSE (every step of the loop)]")
        ;

}


//loop_number == -1 -> salva a prescindere
void ic_output::save_output_base_m(
    std::string run_directory_o, 
    ibsimu_client::run_output_t run_output, 
    ibsimu_client::loop_output_t loop_output,
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
        case ibsimu_client::OUT_EVOLUTION:
            if(! (loop_number % 10))
                save = true;
            break;
        case ibsimu_client::OUT_BEGIN:
            if(loop_number < 3)
                save = true;
            break;
        case ibsimu_client::OUT_VERBOSE:
            if(loop_number < 3 
                || ! (loop_number % 10))
                save = true;
            break;
        case ibsimu_client::OUT_NORMAL:
        default:
            break;
    }


    //if LOOP_ENd saves only at the beginning of the loop
    if(loop_output==ibsimu_client::LOOP_END) {
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
