use crate::particle::Particle;
use rand_xoshiro::Xoshiro256StarStar;
use rand_distr::{Uniform, Distribution, Normal};
use std::fmt::{Display, Formatter};
use std::fmt;
use std::convert::TryInto;
use crate::asc::{Asc, save_asc_from_opt};
use crate::schedule::{Schedule, write_sweep_log, write_data_file};
use crate::PI;
use crate::common_util::{apply_pbc, relative_to_global3, global_to_relative3};

// https://stackoverflow.com/questions/26958178/how-do-i-automatically-implement-comparison-for-structs-with-floats-in-rust
#[derive(Debug, Clone)]
pub struct Sphere {
    pub rel_pos: [f64; 3],
    pub global_pos: [f64; 3],
    pub radius: f64,
}

impl Sphere {
    pub fn make_shape(r: f64) -> Self {
        Sphere { rel_pos: [0.0, 0.0, 0.0], global_pos: [0.0, 0.0, 0.0], radius: r }
    }
}

impl Display for Sphere {
    // From rust docs
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(f, "{} {} {} {}",
               self.rel_pos[0], self.rel_pos[1], self.rel_pos[2], self.radius
        )
    }
}

impl Particle for Sphere {
    const TYPE: &'static str = "Sphere";

    fn parse(line: &str, unit_cell: &[f64]) -> Self {
        let params: Vec<f64> = line.split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect();
        let rel_pos = [params[0], params[1], params[2]];
        let global_pos = relative_to_global3(unit_cell, &rel_pos);
        Sphere { rel_pos, global_pos, radius: params[3] } 
    }

    fn check_overlap(&self, other: &Self, offset: &[f64]) -> bool {
        let image_x = other.global_pos[0] + offset[0];
        let image_y = other.global_pos[1] + offset[1];
        let image_z = other.global_pos[2] + offset[2];
        (self.global_pos[0] - image_x).powi(2) 
            + (self.global_pos[1] - image_y).powi(2)
            + (self.global_pos[2] - image_z).powi(2)
            <= (self.radius + other.radius).powi(2)
    }
    
    fn overlap_scale(&self, other: &Self, offset: &[f64]) -> f64 {
        let image_x = other.global_pos[0] + offset[0];
        let image_y = other.global_pos[1] + offset[1];
        let image_z = other.global_pos[2] + offset[2];
        let displacement = ((self.global_pos[0] - image_x).powi(2) 
            + (self.global_pos[1] - image_y).powi(2)
            + (self.global_pos[2] - image_z).powi(2)).sqrt();
        let extent = self.radius + other.radius;
        (extent/displacement).powi(3) 
    }

    fn copy_shape_random_coord(&self,
                               cell: &[f64],
                               rng: &mut Xoshiro256StarStar
    ) -> Self
    {
        let uni_dist = Uniform::new(0.0, 1.0);
        let lat_x = uni_dist.sample(rng);
        let lat_y = uni_dist.sample(rng);
        let lat_z = uni_dist.sample(rng);
        let rel_pos = [lat_x, lat_y, lat_z];
        let global_pos = relative_to_global3(cell, &rel_pos);
        Sphere { rel_pos, global_pos, radius: self.radius }
    }

    fn perturb(&mut self,
               cell: &[f64],
               param: &[f64],
               rng: &mut Xoshiro256StarStar,
    ) -> (Self, usize)
    {
        let old_sphere = self.clone();
        let normal = Normal::new(0.0, param[0]).unwrap();
        self.global_pos[0] += normal.sample(rng);
        self.global_pos[1] += normal.sample(rng);
        self.global_pos[2] += normal.sample(rng);
        let uncorrected_rel = global_to_relative3(cell, &self.global_pos);
        // Handle pbc
        self.rel_pos = apply_pbc(&uncorrected_rel).as_slice().try_into().unwrap();
        // Recalculate global
        self.global_pos = relative_to_global3(cell, &self.rel_pos);
        (old_sphere, 0)
    }

    // From S. Torquato and Y. Jiao PRE 80, 041104 (2009)
    // Also need nalgebra
    fn apply_strain(&mut self, new_cell: &[f64]) {
        self.global_pos = relative_to_global3(new_cell, &self.rel_pos);
    }

    fn init_obs<C: Asc<Self>>(config: &C) -> Vec<f64> {
        let (_, min_nn_gap, _, _) = config.nn_gap_distribution();
        vec![0.0, 0.0, min_nn_gap, 0.0] // [samples, sum of volume, current min_nn_gap, sum of min_nn_gap]
    }

    fn sample_obs_sweep<C: Asc<Self>>(schedule: &mut Schedule<Self, C>, config: &C) {
        let vol = schedule.running_obs[1]/schedule.running_obs[0];
        let avg_min_gap = schedule.running_obs[3]/schedule.running_obs[0];
        schedule.running_obs = Self::init_obs(config);
        println!("Cell volume averaged over sweep: {}", vol);
        let phi = config.packing_fraction();
        schedule.phi = phi;
        println!("Phi at end of sweep: {}", phi);
        let (avg_nn_gap, min_nn_gap, max_nn_gap, gap_distr) = config.nn_gap_distribution();
        schedule.avg_gap = avg_nn_gap;
        let gap_string: String = gap_distr.iter()
                                          .map(|x| format!("{} {} {}\n", x[0], x[1], x[2]))
                                          .fold(String::new(), |acc, x| acc + &x);
        let fname = format!("nn_gap_{}", schedule.current_sweep);
        write_data_file(&gap_string, &fname);
        let (pressure, unc_pressure, chisq) = config.instantaneous_pressure();
        let logline = format!("{} {} {} {} {} {} {} {} {}", schedule.current_sweep,
                                                vol,
                                                phi,
                                                pressure,
                                                unc_pressure,
                                                avg_min_gap,
                                                avg_nn_gap,
                                                min_nn_gap,
                                                max_nn_gap);
        write_sweep_log(&logline);
        save_asc_from_opt(config, &format!("sweep_{}", schedule.current_sweep));
    }

    fn sample_obs_failed_move<C: Asc<Self>>(
        schedule: &mut Schedule<Self, C>,
        config: &C
    )
    {
        schedule.running_obs[0] += 1.0;
        schedule.running_obs[1] += config.cell_volume();
        // The new min_nn_gap is the same as the old
        schedule.running_obs[3] += schedule.running_obs[2];
    }

    fn sample_obs_accepted_pmove<C: Asc<Self>>(
        schedule: &mut Schedule<Self, C>,
        config: &C,
        changed_idx: usize,
        _old_p: &Self
    )
    {
        schedule.running_obs[0] += 1.0;
        schedule.running_obs[1] += config.cell_volume();
        let possible_gap = config.particle_gap(changed_idx, &config.particle_slice()[changed_idx]);
        if possible_gap < schedule.running_obs[2] {
            schedule.running_obs[2] = possible_gap;
        }
        schedule.running_obs[3] += schedule.running_obs[2];
    }
    
    fn sample_obs_accepted_cmove<C: Asc<Self>>(
        schedule: &mut Schedule<Self, C>,
        config: &C,
        _old_c: &[f64]
    )
    {
        schedule.running_obs[0] += 1.0;
        schedule.running_obs[1] += config.cell_volume();
        let (_, possible_gap, _, _) = config.nn_gap_distribution();
        if possible_gap < schedule.running_obs[2] {
            schedule.running_obs[2] = possible_gap;
        }
        schedule.running_obs[3] += schedule.running_obs[2];
    }

    fn hint_lower(&self) -> f64 {
        self.radius
    }

    fn hint_upper(&self) -> f64 {
        self.radius
    }

    fn vol(&self) -> f64 {
        4.0*PI*self.radius.powi(3)/3.0
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
        // Needs to be different from the one below
        let mut rng = Xoshiro256StarStar::seed_from_u64(0);
        let cell = [219.0232, -4200.21983, -2300.123123, 8443.5, 2213.4, -231.15, 5200.2, -2002.1, 231.24];
        let mut sphere: Sphere = Particle::parse("0.0 0.3 0.22 1.0", &cell);
        println!("{}", sphere);
        sphere.perturb(&cell, &[0.0], &mut rng);
        println!("{}", sphere);
        let local = global_to_relative3(&cell, &sphere.global_pos);
        assert!(local[0] < 0.0);
        assert!(sphere.rel_pos[0] >= 0.0 && sphere.rel_pos[0] < 1.0);
        assert!(sphere.rel_pos[1] >= 0.0 && sphere.rel_pos[1] < 1.0);
        assert!(sphere.rel_pos[2] >= 0.0 && sphere.rel_pos[2] < 1.0);
        let local_pbc = apply_pbc(&local);
        assert!(local_pbc[0] >= 0.0 && local_pbc[0] < 1.0);
        assert!(local_pbc[1] >= 0.0 && local_pbc[1] < 1.0);
        assert!(local_pbc[2] >= 0.0 && local_pbc[2] < 1.0);
    }
    
    #[test]
    fn zero_inversion_strain() {
        // Cell selected by changing by hand until test passes,
        // since this is just a regression test
        let mut rng = Xoshiro256StarStar::seed_from_u64(0);
        let old_cell = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let cell = [2100.0232, -2000.21983, -230.0, 833.5, 2.4, -231.15, 5.2, -200.1, 21.24];
        let mut sphere: Sphere = Particle::parse("0.0 0.3 0.22 1.0", &old_cell);
        sphere.apply_strain(&cell);
        let local = global_to_relative3(&cell, &sphere.global_pos);
        assert!(local[0] < 0.0);
        assert!(sphere.rel_pos[0] >= 0.0 && sphere.rel_pos[0] < 1.0);
        assert!(sphere.rel_pos[1] >= 0.0 && sphere.rel_pos[1] < 1.0);
        assert!(sphere.rel_pos[2] >= 0.0 && sphere.rel_pos[2] < 1.0);
        let local_pbc = apply_pbc(&local);
        println!("{:?}", local_pbc);
        assert!(local_pbc[0] >= 0.0 && local_pbc[0] < 1.0);
        assert!(local_pbc[1] >= 0.0 && local_pbc[1] < 1.0);
        assert!(local_pbc[2] >= 0.0 && local_pbc[2] < 1.0);
    }
}
