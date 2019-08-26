use std::io::prelude::*;
use std::fmt::{Debug, Display, Formatter};
use std::fmt;
use rand_xoshiro::Xoshiro256StarStar;
use rand_distr::{Uniform, Distribution};
use rand::SeedableRng;

mod asc;
mod spheres;

use asc::{Asc, Particle};
use spheres::{Disk};

// 1D point struct for extremely simple debugging
#[derive(Debug)]
struct Point1D {
    pos: [f64; 1],
}

impl Display for Point1D {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(f, "{}",
               self.pos[0]
        )
    }
}

impl Particle for Point1D {
    fn check_overlap(&self, other: &Self, offset: &[f64]) -> bool {
        false
    }
    fn copy_shape_random_coord(&self,
                               cell: &[f64],
                               rng: &mut Xoshiro256StarStar)
    -> Self
    {
        let x = Uniform::new(0.0, cell[0]).sample(rng);
        Point1D { pos: [x] }
    }
}

fn main() {
    let shape = Disk::make_shape(0.1);
    let mut init_cell = vec![10.0, 0.0, 0.0, 10.0];
    let mut rng = Xoshiro256StarStar::seed_from_u64(0);
    let mut a = Asc::make_rsa(1200, &shape, 2, init_cell, &mut rng);
    a.print_asc();
}
