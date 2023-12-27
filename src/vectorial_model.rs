use itertools_num::linspace;
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VectorialModelConfig {
    pub base_q: f64,    // parent production, molecules/second
    pub p_tau_d: f64,   // dissociative lifetime of parent, seconds
    pub p_tau_t: f64,   // total lifetime of parent, seconds
    pub v_outflow: f64, // meters/second
    pub sigma: f64,     // cross section of parent molecule, m^2
    pub f_tau_t: f64,   // total lifetime of fragment, seconds
    pub v_photo: f64,   // photodissociative velocity kick, meters/second
    pub radial_points: usize,
    pub angular_points: usize,
    pub radial_substeps: usize,

    pub parent_destruction_level: f64,
    pub fragment_destruction_level: f64,
    // pub backflow_mode: bool,
}

#[derive(Clone)]
pub struct VectorialModel {
    pub radial_grid: Vec<f64>,
    pub angular_grid: Vec<f64>,
    pub collision_sphere_radius: f64, // meters
    pub coma_radius: f64,             //meters
    pub max_grid_radius: f64,         // meters
    pub epsilon_max: f64,             // radians
}

#[derive(Clone)]
pub struct VectorialModelResult {
    pub volume_density_grid: Vec<f64>,
    pub volume_density: Vec<f64>,

    pub total_fragment_number: f64,
    pub fragment_sputter_grid: Vec<Vec<f64>>,
}

pub fn construct_vectorial_model(vmc: &VectorialModelConfig) -> VectorialModel {
    // collision sphere
    let csr: f64 = (vmc.sigma * vmc.base_q) / (4.0 * vmc.v_outflow);

    // how far parent, fragments travel before destruction_level percent are destroyed
    let coma_r: f64 =
        -1.0 * vmc.v_outflow * vmc.p_tau_t * ((1.0 - vmc.parent_destruction_level).log2());
    let fragment_beta_r: f64 = -1.0 * (1.0 - vmc.fragment_destruction_level).log2();
    let fragment_travel_dist: f64 = fragment_beta_r * vmc.f_tau_t * (vmc.v_outflow + vmc.v_photo);

    // extend the radial grid out far enough to catch the desired destruction_level of fragments
    let grid_max: f64 = coma_r + fragment_travel_dist;

    // if backflow_mode = false, we compute within the collision sphere for better results
    let min_r: f64 = 1000.0;
    // if vmc.backflow_mode {
    //     min_r = 2.0 * csr
    // };

    // make the coordinate grids
    let r_grid = geomspace(min_r, grid_max, vmc.radial_points);
    let ang_grid_prelim: Vec<f64> = linspace::<f64>(0.0, PI, vmc.angular_points + 1).collect();
    // bump the angles up to the midpoints of the prelim grid by adding half the difference of
    // adjacent elements
    // this pushes the angular grid away from the outflow axis at theta = 0
    let ang_grid: Vec<f64> = ang_grid_prelim
        .windows(2)
        .map(|x| x[0] + (x[1] - x[0]) / 2.0)
        .collect();

    // maximum ejection angle in comet frame
    let e_max: f64;
    if vmc.v_photo < vmc.v_outflow {
        e_max = (vmc.v_photo / vmc.v_outflow).asin();
    } else {
        e_max = PI;
    }

    VectorialModel {
        radial_grid: r_grid,
        angular_grid: ang_grid,
        collision_sphere_radius: csr,
        coma_radius: coma_r,
        max_grid_radius: grid_max,
        epsilon_max: e_max,
    }
}

pub fn run_vectorial_model(vmc: &VectorialModelConfig) -> (VectorialModel, VectorialModelResult) {
    let coma: VectorialModel = construct_vectorial_model(vmc);

    // initialize fragment sputter grid
    let mut fsg: Vec<Vec<f64>> = vec![vec![0.0; vmc.angular_points]; vmc.radial_points];

    let d_alpha: f64 = coma.epsilon_max / (vmc.angular_points as f64);
    let integration_factor: f64 = (1.0 / (4.0 * PI * vmc.p_tau_d)) * d_alpha / (4.0 * PI);

    // compute fragment sputter
    for i in 0..vmc.radial_points {
        for j in 0..vmc.angular_points {
            fsg[i][j] =
                fragment_sputter_at_point(coma.radial_grid[i], coma.angular_grid[j], &vmc, &coma)
                    * integration_factor;
        }
    }

    // compute fragment density as function of r
    let vdg: Vec<f64> = coma.radial_grid.clone();
    let mut vd: Vec<f64> = vec![0.0; vdg.len()];
    for i in 0..vmc.radial_points {
        for j in 0..vmc.angular_points {
            vd[i] += 2.0 * PI * fsg[i][j] * coma.angular_grid[j].sin();
        }
    }

    let total_number: f64 = calculate_total_fragment_count(&coma.radial_grid, &vd);

    let vmr = VectorialModelResult {
        volume_density_grid: vdg,
        volume_density: vd,
        total_fragment_number: total_number,
        fragment_sputter_grid: fsg,
    };

    (coma, vmr)
}

pub fn fragment_sputter_at_point(
    r: f64,
    theta: f64,
    vmc: &VectorialModelConfig,
    coma: &VectorialModel,
) -> f64 {
    let mut sputter: f64 = 0.0;
    let vp: f64 = vmc.v_outflow;
    let vf: f64 = vmc.v_photo;
    let x: f64 = r * theta.sin();
    let y: f64 = r * theta.cos();
    let (outflow_axis_points, drs) = outflow_axis_sample(
        x,
        y,
        theta,
        vmc.radial_substeps,
        coma.max_grid_radius,
        coma.epsilon_max,
    );

    let outflow_samples = outflow_axis_points.iter().zip(drs.iter());

    for (slice_r, dr) in outflow_samples {
        let sep_dist: f64 = (x.powi(2) + (slice_r - y).powi(2)).sqrt();
        let cos_eject: f64 = (y - slice_r) / sep_dist;
        let sin_eject: f64 = x / sep_dist;

        let p_extinction: f64 = ((-1.0 * slice_r) / (vmc.p_tau_t * vp)).exp();
        let v_factor: f64 = (vf.powi(2) - vp.powi(2) * sin_eject.powi(2)).sqrt();

        let v_one: f64 = vp * cos_eject + v_factor;
        let v_two: f64 = vp * cos_eject - v_factor;

        let velocities: Vec<f64>;
        if vf > vp {
            velocities = vec![v_one];
        } else {
            velocities = vec![v_one, v_two];
        }

        for v in &velocities {
            let t_frag: f64 = sep_dist / v;
            // let t_total: f64 = (slice_r / vp) + t_frag;
            // let q = production_at_time(t_total) / vp;
            let q = vmc.base_q / vp;
            let q_r_eps = (v.powi(2) * q) / (vf * (v - vp * cos_eject).abs());
            let f_extinction = (-1.0 * t_frag / vmc.f_tau_t).exp();
            let n_r = (p_extinction * f_extinction * q_r_eps) / (sep_dist.powi(2) * v);
            sputter = sputter + n_r * dr;
        }
    }

    sputter
}

pub fn outflow_axis_sample(
    x: f64,
    y: f64,
    theta: f64,
    radial_substeps: usize,
    max_grid_radius: f64,
    epsilon_max: f64,
) -> (Vec<f64>, Vec<f64>) {
    let grid_edge_angle: f64 = x.atan2(y - max_grid_radius);
    let max_subangle: f64;
    if grid_edge_angle < epsilon_max {
        max_subangle = grid_edge_angle;
    } else {
        max_subangle = epsilon_max;
    }

    // divide the space into (radial_substeps + 1) slices and take the midpoints
    let subangles: Vec<f64> = linspace::<f64>(theta, max_subangle, radial_substeps + 1).collect();
    let ejection_grid: Vec<f64> = subangles.iter().map(|t| (y - (x / t.tan()))).collect();
    let ejection_sites: Vec<f64> = ejection_grid
        .windows(2)
        .map(|x| x[0] + (x[1] - x[0]) / 2.0)
        .collect();
    // size of the slices along the axis
    let drs: Vec<f64> = ejection_grid.windows(2).map(|x| x[1] - x[0]).collect();

    (ejection_sites, drs)
}

pub fn vectorial_model_config_from_yaml(yaml_filename: std::path::PathBuf) -> VectorialModelConfig {
    let f = std::fs::File::open(yaml_filename).expect("Could not read yaml config file!");
    let vmc: VectorialModelConfig =
        serde_yaml::from_reader(f).expect("Format error in yaml config file!");
    vmc
}

fn calculate_total_fragment_count(
    volume_density_grid: &Vec<f64>,
    volume_density: &Vec<f64>,
) -> f64 {
    // difference between radial density grid points
    let drs: Vec<f64> = volume_density_grid
        .windows(2)
        .map(|x| x[1] - x[0])
        .collect();
    // take average density of left and right grid point
    let density_at_midpoints: Vec<f64> = volume_density
        .windows(2)
        .map(|x| (x[1] + x[0]) / 2.0)
        .collect();
    // calculate the mid-points between our radial density grid points
    let rs_at_midpoints: Vec<f64> = volume_density_grid
        .windows(2)
        .map(|x| (x[1] + x[0]) / 2.0)
        .collect();
    // radial part of spherical integral over the coma
    let total_number: f64 = rs_at_midpoints
        .iter()
        .zip(density_at_midpoints.iter().zip(drs.iter()))
        .map(|(r, (nr, dr))| r * r * nr * dr)
        .sum();

    4.0 * PI * total_number
}

pub fn logspace(power_start: f64, power_end: f64, num: usize) -> Vec<f64> {
    let powers = linspace::<f64>(power_start, power_end, num);
    let points: Vec<f64> = powers.into_iter().map(|p| 10.0_f64.powf(p)).collect();

    points
}

pub fn geomspace(value_start: f64, value_end: f64, num: usize) -> Vec<f64> {
    logspace(value_start.log10(), value_end.log10(), num)
}
