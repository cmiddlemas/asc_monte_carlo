use crate::asc::{Particle};
use rand_xoshiro::Xoshiro256StarStar;
use rand_distr::{Uniform, Distribution, Normal};
use std::fmt::{Display, Formatter};
use std::fmt;
use std::cmp::PartialEq;
use crate::asc::{Asc, save_asc_from_opt};
use crate::schedule::Schedule;
use crate::OPT;

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
        const DIM: usize = 2;
        // convert to lattice coords
        // https://www.mathsisfun.com/algebra/matrix-inverse.html
        let det = c[DIM*0 + 0]*c[DIM*1 + 1] 
                - c[DIM*1 + 0]*c[DIM*0 + 1];
        let mut lat_x = (c[DIM*1 + 1]*self.pos[0] 
               - c[DIM*0 + 1]*self.pos[1])
                /det;
        let mut lat_y = (-c[DIM*1 + 0]*self.pos[0]
                 +c[DIM*0 + 0]*self.pos[1])
                /det;
        // put lattice coords back in unit square
        lat_x = lat_x - lat_x.floor();
        lat_y = lat_y - lat_y.floor();
        // convert back to euclidean coords
        self.pos[0] = lat_x*c[DIM*0 + 0] + lat_y*c[DIM*0 + 1];
        self.pos[1] = lat_x*c[DIM*1 + 0] + lat_y*c[DIM*1 + 1];
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
            [ lat_x*cell[DIM*0 + 0] + lat_y*cell[DIM*0 + 1],
              lat_x*cell[DIM*1 + 0] + lat_y*cell[DIM*1 + 1]],   
            radius: self.radius }
    }

    fn perturb(&mut self,
               cell: &[f64],
               param: &[f64],
               rng: &mut Xoshiro256StarStar,
    ) -> Self
    {
        let old_disk = self.clone();
        let normal = Normal::new(0.0, param[0]).unwrap();
        self.pos[0] += normal.sample(rng);
        self.pos[1] += normal.sample(rng);
        // Handle pbc
        self.apply_pbc(cell);
        old_disk
    }

    // From S. Torquato and Y. Jiao PRE 80, 041104 (2009)
    fn apply_strain(&mut self, old_cell: &[f64], new_cell: &[f64]) {
        const DIM: usize = 2;
        // convert to lattice coords
        // https://www.mathsisfun.com/algebra/matrix-inverse.html
        let det = old_cell[DIM*0 + 0]*old_cell[DIM*1 + 1] 
                - old_cell[DIM*1 + 0]*old_cell[DIM*0 + 1];
        let lat_x = (old_cell[DIM*1 + 1]*self.pos[0] 
               - old_cell[DIM*0 + 1]*self.pos[1])
                /det;
        let lat_y = (-old_cell[DIM*1 + 0]*self.pos[0]
                 +old_cell[DIM*0 + 0]*self.pos[1])
                /det;
        // set to new cell
        self.pos[0] = lat_x*new_cell[DIM*0 + 0] + lat_y*new_cell[DIM*0 + 1];
        self.pos[1] = lat_x*new_cell[DIM*1 + 0] + lat_y*new_cell[DIM*1 + 1];
    }

    fn init_obs() -> Vec<f64> {
        vec![0.0, 0.0] // [samples, sum of volume]
    }

    fn sample_obs_sweep(schedule: &mut Schedule<Self>, config: &Asc<Self>) {
        let vol = schedule.running_obs[1]/schedule.running_obs[0];
        schedule.running_obs = vec![0.0,0.0];
        println!("Cell volume over sweep: {}", vol);
        save_asc_from_opt(config, &format!("sweep_{}", schedule.current_sweep));
    }

    fn sample_obs_failed_move(
        schedule: &mut Schedule<Self>,
        config: &Asc<Self>
    )
    {
        schedule.running_obs[0] += 1.0;
        schedule.running_obs[1] += config.cell_volume();
    }

    fn sample_obs_accepted_pmove(
        schedule: &mut Schedule<Self>,
        config: &Asc<Self>,
        changed_idx: usize,
        old_p: &Self
    )
    {
        schedule.running_obs[0] += 1.0;
        schedule.running_obs[1] += config.cell_volume();
    }
    
    fn sample_obs_accepted_cmove(
        schedule: &mut Schedule<Self>,
        config: &Asc<Self>,
        old_c: &[f64]
    )
    {
        schedule.running_obs[0] += 1.0;
        schedule.running_obs[1] += config.cell_volume();
    }
}
