use crate::particle::Particle;
use rand_xoshiro::Xoshiro256StarStar;
use rand_distr::{Uniform, Distribution, Normal};
use std::fmt::{Display, Formatter};
use std::fmt;
use crate::asc::{Asc, save_asc_from_opt};
use crate::schedule::{Schedule, write_sweep_log};
use crate::PI;
use std::convert::TryInto;
use crate::common_util::{apply_pbc, global_to_relative2, relative_to_global2};

// https://stackoverflow.com/questions/26958178/how-do-i-automatically-implement-comparison-for-structs-with-floats-in-rust
#[derive(Debug, Clone)]
pub struct Disk {
    rel_pos: [f64; 2], // The relative position is the "source of truth"
    global_pos: [f64; 2], // We cache the global position to avoid calculating unecessarily
    radius: f64, // Radius in global distances
}

impl Disk {
    pub fn make_shape(r: f64) -> Self {
        Disk { rel_pos: [0.0, 0.0], global_pos: [0.0, 0.0], radius: r }
    }
}

impl Display for Disk {
    // From rust docs
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(f, "{} {} {}",
               self.rel_pos[0], self.rel_pos[1],  self.radius
        )
    }
}

impl Particle for Disk {
    const TYPE: &'static str = "Disk";

    fn parse(line: &str, unit_cell: &[f64]) -> Self {
        let params: Vec<f64> = line.split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect();
        let rel_pos = [params[0], params[1]];
        let global_pos = relative_to_global2(unit_cell, &rel_pos);
        Disk { rel_pos, global_pos, radius: params[2] } 
    }

    fn check_overlap(&self, other: &Self, offset: &[f64]) -> bool {
        let image_x = other.global_pos[0] + offset[0];
        let image_y = other.global_pos[1] + offset[1];
        (self.global_pos[0] - image_x).powi(2) 
            + (self.global_pos[1] - image_y).powi(2)
            <= (self.radius + other.radius).powi(2)
    }
    
    fn copy_shape_random_coord(&self,
                               cell: &[f64],
                               rng: &mut Xoshiro256StarStar
    ) -> Self
    {
        let uni_dist = Uniform::new(0.0, 1.0);
        let lat_x = uni_dist.sample(rng);
        let lat_y = uni_dist.sample(rng);
        let rel_pos = [lat_x, lat_y];
        let global_pos = relative_to_global2(cell, &rel_pos);
        // Turn lattice coords into euclidean coords
        Disk { rel_pos, global_pos, radius: self.radius }
    }

    fn perturb(&mut self,
               cell: &[f64],
               param: &[f64],
               rng: &mut Xoshiro256StarStar,
    ) -> (Self, usize)
    {
        let old_disk = self.clone();
        let normal = Normal::new(0.0, param[0]).unwrap();
        self.global_pos[0] += normal.sample(rng);
        self.global_pos[1] += normal.sample(rng);
        let uncorrected_rel = global_to_relative2(cell, &self.global_pos);
        // Handle pbc
        self.rel_pos = apply_pbc(&uncorrected_rel).as_slice().try_into().unwrap();
        // Recalculate global because relative is primary
        self.global_pos = relative_to_global2(cell, &self.rel_pos);
        (old_disk, 0)
    }

    // From S. Torquato and Y. Jiao PRE 80, 041104 (2009)
    fn apply_strain(&mut self, new_cell: &[f64]) {
        // Just need to change global coords
        self.global_pos = relative_to_global2(new_cell, &self.rel_pos);
    }

    fn init_obs() -> Vec<f64> {
        vec![0.0, 0.0] // [samples, sum of volume]
    }

    fn sample_obs_sweep<C: Asc<Self>>(schedule: &mut Schedule<Self>, config: &C) {
        let vol = schedule.running_obs[1]/schedule.running_obs[0];
        schedule.running_obs = vec![0.0,0.0];
        println!("Cell volume over sweep: {}", vol);
        let phi = (config.n_particles() as f64)
            *PI*config.first_particle().radius.powi(2)
            /vol;
        schedule.phi = phi;
        println!("Phi over sweep: {}", phi);
        let logline = format!("{} {} {}", schedule.current_sweep, vol, phi);
        write_sweep_log(&logline);
        save_asc_from_opt(config, &format!("sweep_{}", schedule.current_sweep));
    }

    fn sample_obs_failed_move<C: Asc<Self>>(
        schedule: &mut Schedule<Self>,
        config: &C
    )
    {
        schedule.running_obs[0] += 1.0;
        schedule.running_obs[1] += config.cell_volume();
    }

    fn sample_obs_accepted_pmove<C: Asc<Self>>(
        schedule: &mut Schedule<Self>,
        config: &C,
        _changed_idx: usize,
        _old_p: &Self
    )
    {
        schedule.running_obs[0] += 1.0;
        schedule.running_obs[1] += config.cell_volume();
    }
    
    fn sample_obs_accepted_cmove<C: Asc<Self>>(
        schedule: &mut Schedule<Self>,
        config: &C,
        _old_c: &[f64]
    )
    {
        schedule.running_obs[0] += 1.0;
        schedule.running_obs[1] += config.cell_volume();
    }

    fn hint_lower(&self) -> f64 {
        self.radius
    }

    fn hint_upper(&self) -> f64 {
        self.radius
    }

    fn lat_coord(&self) -> Vec<f64> {
        self.rel_pos.to_vec()
    }
}
