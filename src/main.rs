mod vectorial_model;
use crate::vectorial_model::VectorialModel;
use crate::vectorial_model::VectorialModelConfig;
use crate::vectorial_model::VectorialModelResult;
use std::io::Write;

struct Cli {
    yaml_input_filename: std::path::PathBuf,
    vectorial_model_output_filename: std::path::PathBuf,
}

fn main() {
    let cli = parse_arguments();
    let vmc = vectorial_model::vectorial_model_config_from_yaml(cli.yaml_input_filename);

    let (coma, vmr): (VectorialModel, VectorialModelResult) =
        vectorial_model::run_vectorial_model(&vmc);

    write_results(&vmc, &coma, &vmr, &cli.vectorial_model_output_filename);
}

fn parse_arguments() -> Cli {
    let yaml_input_fname = std::env::args().nth(1).expect("No yaml input file given");
    let yaml_output_fname = std::env::args().nth(2).expect("No output filename given");
    let args = Cli {
        yaml_input_filename: std::path::PathBuf::from(yaml_input_fname),
        vectorial_model_output_filename: std::path::PathBuf::from(yaml_output_fname),
    };
    args
}

fn write_results(
    vmc: &VectorialModelConfig,
    coma: &VectorialModel,
    vmr: &VectorialModelResult,
    fname: &std::path::PathBuf,
) {
    let f = std::fs::OpenOptions::new()
        .write(true)
        .truncate(true)
        .create(true)
        .open(fname)
        .unwrap();
    let mut writer = std::io::BufWriter::new(&f);

    write_results_config(&mut writer, &vmc);
    write_results_model_info(&mut writer, &coma);
    write_results_fragment_density(&mut writer, &vmr);
    write_results_fragment_sputter(&mut writer, &coma, &vmr);
}

fn write_results_config(w: &mut std::io::BufWriter<&std::fs::File>, vmc: &VectorialModelConfig) {
    write!(
        w,
        "Radial grid size: {}\tAngular grid size: {}\tRadial substeps: {}\n\n",
        vmc.radial_points, vmc.angular_points, vmc.radial_substeps,
    )
    .unwrap();
    write!(
        w,
        "Q: {:02.5e} mol/s\tParent dissociative lifetime: {:02.5e} s\tParent total lifetime: {:02.5e} s\tParent outflow velocity {:02.5e} m/s\n",
        vmc.base_q,
        vmc.p_tau_d,
        vmc.p_tau_t,
        vmc.v_outflow
    )
    .unwrap();
    write!(
        w,
        "Parent cross sectional area: {:02.5e} m^2\tFragment total lifetime: {:02.5e} s\tFragment velocity: {:02.5e} m/s\n",
        vmc.sigma,
        vmc.f_tau_t,
        vmc.v_photo,
        ).unwrap();
    write!(
        w,
        "Parent destruction level: {:03.1}%\tFragment destruction level {:03.1}%\n\n",
        vmc.parent_destruction_level * 100.0,
        vmc.fragment_destruction_level * 100.0,
    )
    .unwrap();
}

fn write_results_model_info(w: &mut std::io::BufWriter<&std::fs::File>, coma: &VectorialModel) {
    write!(
        w,
        "Collision sphere radius: {:02.5e} km\ncoma radius: {:02.5e} km\nmax grid extent: {:02.5e} km\n\n",
        coma.collision_sphere_radius / 1000.0,
        coma.coma_radius / 1000.0,
        coma.max_grid_radius / 1000.0,
    ).unwrap();
}

fn write_results_fragment_density(
    w: &mut std::io::BufWriter<&std::fs::File>,
    vmr: &VectorialModelResult,
) {
    write!(
        w,
        "Total number of fragments: {:02.25e}\n\n",
        vmr.total_fragment_number
    )
    .unwrap();

    let results = vmr
        .volume_density_grid
        .iter()
        .zip(vmr.volume_density.iter());
    write!(w, "r (meters)\t\t\tvolume density (1/m^3)\n").unwrap();
    for (r, n) in results {
        write!(w, "{:02.25e},\t{:02.25e}\n", r, n).unwrap();
    }
}

fn write_results_fragment_sputter(
    w: &mut std::io::BufWriter<&std::fs::File>,
    coma: &VectorialModel,
    vmr: &VectorialModelResult,
) {
    write!(
        w,
        "r (meters),\t\t\ttheta (radians),\t\tfragment density (1/m^3)\n"
    )
    .unwrap();

    for i in 0..coma.radial_grid.len() - 20 {
        for j in 0..coma.angular_grid.len() - 20 {
            write!(
                w,
                "{:02.25e},\t{:02.25e},\t{:02.25e}\n",
                coma.radial_grid[i] / 1000.0,
                coma.angular_grid[j],
                vmr.fragment_sputter_grid[i][j]
            )
            .unwrap();
        }
    }
}
