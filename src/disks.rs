use crate::particle::Particle;
use rand_xoshiro::Xoshiro256StarStar;
use rand_distr::{Uniform, Distribution, Normal};
use std::fmt::{Display, Formatter};
use std::fmt;
use crate::asc::{Asc, save_asc_from_opt};
use crate::schedule::{ObservableTracker, Schedule, write_sweep_log};
// https://stackoverflow.com/questions/31208465/pi-constant-is-ambiguous
use std::f64::consts::PI;
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

    fn init_obs<C: Asc<Self>>(_config: &C) -> ObservableTracker {
        ObservableTracker::new()
    }

    fn sample_obs_sweep<C: Asc<Self>>(schedule: &mut Schedule<Self, C>, config: &C) {
        let vol = schedule.running_obs.sum_of_vol/schedule.running_obs.n_samples;
        schedule.running_obs = Self::init_obs(config);
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
        schedule: &mut Schedule<Self, C>,
        config: &C
    )
    {
        schedule.running_obs.n_samples += 1.0;
        schedule.running_obs.sum_of_vol += config.cell_volume();
    }

    fn sample_obs_accepted_pmove<C: Asc<Self>>(
        schedule: &mut Schedule<Self, C>,
        config: &C,
        _changed_idx: usize,
        _old_p: &Self
    )
    {
        schedule.running_obs.n_samples += 1.0;
        schedule.running_obs.sum_of_vol += config.cell_volume();
    }
    
    fn sample_obs_accepted_cmove<C: Asc<Self>>(
        schedule: &mut Schedule<Self, C>,
        config: &C,
        _old_c: &[f64]
    )
    {
        schedule.running_obs.n_samples += 1.0;
        schedule.running_obs.sum_of_vol += config.cell_volume();
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

#[cfg(test)]
mod tests {
    use super::*;
    use rand_xoshiro::rand_core::SeedableRng;

    #[test]
    fn zero_inversion_perturb() {
        // Cell selected by changing by hand until test passes,
        // since this is just a regression test
        // Need to change example because I changed
        // matrix inversion method to nalgebra default
        let mut rng = Xoshiro256StarStar::seed_from_u64(0);
        let cell = [2100.0232, -2000.21983, -230.0, 5.2];
        let mut disk: Disk = Particle::parse("0.0 0.3 1.0", &cell);
        disk.perturb(&cell, &[0.0], &mut rng);
        let local = global_to_relative2(&cell, &disk.global_pos);
        assert!(local[0] < 0.0);
        assert!(disk.rel_pos[0] >= 0.0 && disk.rel_pos[0] < 1.0);
        assert!(disk.rel_pos[1] >= 0.0 && disk.rel_pos[1] < 1.0);
        let local_pbc = apply_pbc(&local);
        assert!(local_pbc[0] >= 0.0 && local_pbc[0] < 1.0);
        assert!(local_pbc[1] >= 0.0 && local_pbc[1] < 1.0);
    }
    
    #[test]
    fn zero_inversion_strain() {
        // Cell selected by changing by hand until test passes,
        // since this is just a regression test
        let mut rng = Xoshiro256StarStar::seed_from_u64(0);
        let old_cell = [1.0, 0.0, 0.0, 1.0];
        let cell = [2100.0232, -2000.21983, -230.0, 5.2];
        let mut disk: Disk = Particle::parse("0.0 0.3 1.0", &old_cell);
        disk.apply_strain(&cell);
        let local = global_to_relative2(&cell, &disk.global_pos);
        assert!(local[0] < 0.0);
        assert!(disk.rel_pos[0] >= 0.0 && disk.rel_pos[0] < 1.0);
        assert!(disk.rel_pos[1] >= 0.0 && disk.rel_pos[1] < 1.0);
        let local_pbc = apply_pbc(&local);
        assert!(local_pbc[0] >= 0.0 && local_pbc[0] < 1.0);
        assert!(local_pbc[1] >= 0.0 && local_pbc[1] < 1.0);
    }

}
