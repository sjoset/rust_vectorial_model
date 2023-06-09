// Runs vectorial models on a grid of heliocentric distances and reports the total number of
// fragments at each heliocentric distance, given a configuration of values at 1 AU
use itertools_num::linspace;
mod vectorial_model;
use crate::vectorial_model::VectorialModel;
use crate::vectorial_model::VectorialModelConfig;
use crate::vectorial_model::VectorialModelResult;
use std::io::Write;

struct Cli {
    yaml_input_filename: std::path::PathBuf,
    vectorial_model_output_filename: std::path::PathBuf,
    rh_start: f64,
    rh_end: f64,
    rh_gridsize: usize,
}

fn main() {
    let cli = parse_arguments();
    let vmc_base: VectorialModelConfig =
        vectorial_model::vectorial_model_config_from_yaml(cli.yaml_input_filename);

    let results: Vec<(f64, f64)> =
        run_models_on_rh_grid(cli.rh_start, cli.rh_end, cli.rh_gridsize, &vmc_base);

    write_results(&results, &cli.vectorial_model_output_filename);
}

fn parse_arguments() -> Cli {
    let yaml_input_fname = std::env::args().nth(1).expect("No yaml input file given");
    let yaml_output_fname = std::env::args().nth(2).expect("No yaml output file given");
    let rh_start_arg = std::env::args()
        .nth(3)
        .expect("No starting heliocentric distance (AU) given");
    let rh_end_arg = std::env::args()
        .nth(4)
        .expect("No ending heliocentric distance (AU) given");
    let rh_gridsize_arg = std::env::args()
        .nth(5)
        .expect("No grid size for heliocentric distance given");
    let args = Cli {
        yaml_input_filename: std::path::PathBuf::from(yaml_input_fname),
        vectorial_model_output_filename: std::path::PathBuf::from(yaml_output_fname),
        rh_start: rh_start_arg.parse().unwrap(),
        rh_end: rh_end_arg.parse().unwrap(),
        rh_gridsize: rh_gridsize_arg.parse().unwrap(),
    };
    args
}

fn run_models_on_rh_grid(
    rh_start: f64,
    rh_end: f64,
    rh_gridsize: usize,
    vmc_base: &VectorialModelConfig,
) -> Vec<(f64, f64)> {
    // vector of tuples: (rh, total_fragment_number)
    let mut results: Vec<(f64, f64)> = Vec::new();
    // grid of r_h at which to run these vmcs
    let rhs: Vec<f64> = linspace::<f64>(rh_start, rh_end, rh_gridsize).collect();

    for (i, rh) in rhs.iter().enumerate() {
        let mut new_vmc: VectorialModelConfig = vmc_base.clone();
        let r_squared: f64 = rh.powi(2);

        //adjust lifetimes and parent outflow velocity for heliocentric distance using empirical
        //relation for outflow
        new_vmc.p_tau_d *= r_squared;
        new_vmc.p_tau_t *= r_squared;
        new_vmc.f_tau_t *= r_squared;
        new_vmc.v_outflow = 0.85 / rh.sqrt();

        print!(
            "Computing model {:05} at r_h = {:02.4} AU : v_outflow = {:02.4} ...",
            i, rh, new_vmc.v_outflow
        );

        let (_, vmr): (VectorialModel, VectorialModelResult) =
            vectorial_model::run_vectorial_model(&new_vmc);

        println!("total fragments: {:02.25e}", vmr.total_fragment_number);

        results.push((rh.clone(), vmr.total_fragment_number))
    }

    results
}

fn write_results(results: &Vec<(f64, f64)>, fname: &std::path::PathBuf) {
    let f = std::fs::OpenOptions::new()
        .write(true)
        .truncate(true)
        .create(true)
        .open(fname)
        .unwrap();
    let mut writer = std::io::BufWriter::new(&f);

    write_results_fragment_density(&mut writer, &results);
}

fn write_results_fragment_density(
    w: &mut std::io::BufWriter<&std::fs::File>,
    results: &Vec<(f64, f64)>,
) {
    write!(w, "r (AU),total_fragment_number\n").unwrap();
    for (rh, total_fragment_number) in results {
        write!(w, "{:02.25e},{:02.25e}\n", rh, total_fragment_number).unwrap();
    }
}
