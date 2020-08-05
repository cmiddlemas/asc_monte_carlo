use crate::particle::Particle;
use rand_xoshiro::Xoshiro256StarStar;
use rand_distr::{Uniform, Distribution, Normal};
use std::fmt::{Display, Formatter};
use std::fmt;
use std::convert::TryInto;
use nalgebra::{Matrix3, Vector3};
use crate::asc::{Asc, save_asc_from_opt};
use crate::schedule::{Schedule, write_sweep_log};
use crate::PI;
use crate::common_util::apply_pbc;

// https://stackoverflow.com/questions/26958178/how-do-i-automatically-implement-comparison-for-structs-with-floats-in-rust
#[derive(Debug, Clone)]
pub struct Sphere {
    pos: [f64; 3],
    radius: f64,
}

impl Sphere {
    pub fn make_shape(r: f64) -> Self {
        Sphere { pos: [0.0, 0.0, 0.0], radius: r }
    }
}

impl Display for Sphere {
    // From rust docs
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(f, "{} {} {} {}",
               self.pos[0], self.pos[1], self.pos[2], self.radius
        )
    }
}

impl Particle for Sphere {
    const TYPE: &'static str = "Sphere";

    fn parse(line: &str, unit_cell: &[f64]) -> Self {
        let params: Vec<f64> = line.split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect();
        Sphere { pos: [params[0], params[1], params[2]], radius: params[3] } 
    }

    fn check_overlap(&self, other: &Self, offset: &[f64]) -> bool {
        //if *self == *other && offset.iter().all(|&x| x == 0.0) {
        //    return false; 
        //}
        let image_x = other.pos[0] + offset[0];
        let image_y = other.pos[1] + offset[1];
        let image_z = other.pos[2] + offset[2];
        (self.pos[0] - image_x).powi(2) 
            + (self.pos[1] - image_y).powi(2)
            + (self.pos[2] - image_z).powi(2)
            <= (self.radius + other.radius).powi(2)
    }
    
    fn copy_shape_random_coord(&self,
                               cell: &[f64],
                               rng: &mut Xoshiro256StarStar
    ) -> Self
    {
        const DIM: usize = 3;
        let uni_dist = Uniform::new(0.0, 1.0);
        let lat_x = uni_dist.sample(rng);
        let lat_y = uni_dist.sample(rng);
        let lat_z = uni_dist.sample(rng);
        // Turn lattice coords into euclidean coords
        Sphere { pos: 
            [ lat_x*cell[DIM*0 + 0] + lat_y*cell[DIM*1 + 0] + lat_z*cell[DIM*2 + 0],
              lat_x*cell[DIM*0 + 1] + lat_y*cell[DIM*1 + 1] + lat_z*cell[DIM*2 + 1],
              lat_x*cell[DIM*0 + 2] + lat_y*cell[DIM*1 + 2] + lat_z*cell[DIM*2 + 2]],
            radius: self.radius }
    }

    fn perturb(&mut self,
               cell: &[f64],
               param: &[f64],
               rng: &mut Xoshiro256StarStar,
    ) -> (Self, usize)
    {
        let old_sphere = self.clone();
        let normal = Normal::new(0.0, param[0]).unwrap();
        self.pos[0] += normal.sample(rng);
        self.pos[1] += normal.sample(rng);
        self.pos[2] += normal.sample(rng);
        // Handle pbc
        self.pos = apply_pbc(&self.pos).as_slice().try_into().unwrap();
        (old_sphere, 0)
    }

    // From S. Torquato and Y. Jiao PRE 80, 041104 (2009)
    // Also need nalgebra
    fn apply_strain(&mut self, new_cell: &[f64]) {
    }

    fn init_obs() -> Vec<f64> {
        vec![0.0, 0.0] // [samples, sum of volume]
    }

    fn sample_obs_sweep<C: Asc<Self>>(schedule: &mut Schedule<Self>, config: &C) {
        let vol = schedule.running_obs[1]/schedule.running_obs[0];
        schedule.running_obs = vec![0.0,0.0];
        println!("Cell volume over sweep: {}", vol);
        let phi = (config.n_particles() as f64)
            *4.0*PI*config.first_particle().radius.powi(3)
            /(3.0*vol);
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
        unimplemented!()
    }
}
