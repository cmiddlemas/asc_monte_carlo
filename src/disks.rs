use crate::particle::Particle;
use rand_xoshiro::Xoshiro256StarStar;
use rand_distr::{Uniform, Distribution, Normal};
use std::fmt::{Display, Formatter};
use std::fmt;
use crate::asc::{Asc, save_asc_from_opt};
use crate::schedule::{Schedule, write_sweep_log};
use crate::PI;
use nalgebra::{Matrix2, Vector2};
use std::convert::TryInto;

// https://stackoverflow.com/questions/26958178/how-do-i-automatically-implement-comparison-for-structs-with-floats-in-rust
#[derive(Debug, Clone)]
pub struct Disk {
    pos: [f64; 2],
    radius: f64,
}

impl Disk {
    pub fn make_shape(r: f64) -> Self {
        Disk { pos: [0.0, 0.0], radius: r }
    }

    // c is the unit cell, given as (u_rc
    // u_00 u_01 u_10 u_11 
    fn apply_pbc(&mut self, c: &[f64]) {
        let u = Matrix2::from_column_slice(c);
        let u_inv = u.lu().try_inverse().expect("unit cell matrix must be invertible");
        let r = Vector2::from_column_slice(&self.pos);
        // convert to lattice coords
        let mut lat_c = u_inv*r;
        // put lattice coords back in unit square
        lat_c.apply(|x| x - x.floor());
        lat_c.apply(|x| x - x.floor());
        // convert back to euclidean coords
        // https://stackoverflow.com/questions/25428920/how-to-get-a-slice-as-an-array-in-rust
        self.pos = (u*lat_c).as_slice().try_into().unwrap();
    }
}

impl Display for Disk {
    // From rust docs
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(f, "{} {} {}",
               self.pos[0], self.pos[1], self.radius
        )
    }
}

impl Particle for Disk {
    const TYPE: &'static str = "Disk";

    fn parse(line: &str) -> Self {
        let params: Vec<f64> = line.split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect();
        Disk { pos: [params[0], params[1]], radius: params[2] } 
    }

    fn check_overlap(&self, other: &Self, offset: &[f64]) -> bool {
        //if *self == *other && offset.iter().all(|&x| x == 0.0) {
        //    return false; 
        //}
        let image_x = other.pos[0] + offset[0];
        let image_y = other.pos[1] + offset[1];
        (self.pos[0] - image_x).powi(2) 
            + (self.pos[1] - image_y).powi(2)
            <= (self.radius + other.radius).powi(2)
    }
    
    fn copy_shape_random_coord(&self,
                               cell: &[f64],
                               rng: &mut Xoshiro256StarStar
    ) -> Self
    {
        const DIM: usize = 2;
        let uni_dist = Uniform::new(0.0, 1.0);
        let lat_x = uni_dist.sample(rng);
        let lat_y = uni_dist.sample(rng);
        // Turn lattice coords into euclidean coords
        Disk { pos: 
            [ lat_x*cell[DIM*0 + 0] + lat_y*cell[DIM*1 + 0],
              lat_x*cell[DIM*0 + 1] + lat_y*cell[DIM*1 + 1]],
            radius: self.radius }
    }

    fn perturb(&mut self,
               cell: &[f64],
               param: &[f64],
               rng: &mut Xoshiro256StarStar,
    ) -> (Self, usize)
    {
        let old_disk = self.clone();
        let normal = Normal::new(0.0, param[0]).unwrap();
        self.pos[0] += normal.sample(rng);
        self.pos[1] += normal.sample(rng);
        // Handle pbc
        self.apply_pbc(cell);
        (old_disk, 0)
    }

    // From S. Torquato and Y. Jiao PRE 80, 041104 (2009)
    fn apply_strain(&mut self, old_cell: &[f64], new_cell: &[f64]) {
        // get old lattice coords
        let u_old_inv = Matrix2::from_column_slice(old_cell)
            .lu()
            .try_inverse()
            .expect("Unit cells must be invertible");
        let lat_c = u_old_inv*Vector2::from_column_slice(&self.pos);
        // set to new global coords
        let u_new = Matrix2::from_column_slice(new_cell);
        self.pos = (u_new*lat_c).as_slice().try_into().unwrap();
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

    fn lat_coord(&self, cell: &[f64]) -> Vec<f64> {
        let u = Matrix2::from_column_slice(cell);
        let u_inv = u.lu().try_inverse().expect("unit cell matrix must be invertible");
        let r = Vector2::from_column_slice(&self.pos);
        // convert to lattice coords
        let lat_c = u_inv*r;
        lat_c.as_slice().try_into().unwrap()
    }
}
