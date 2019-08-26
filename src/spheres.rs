use crate::asc::{Particle};
use rand_xoshiro::Xoshiro256StarStar;
use rand_distr::{Uniform, Distribution};
use std::fmt::{Display, Formatter};
use std::fmt;

#[derive(Debug)]
pub struct Disk {
    pos: [f64; 2],
    radius: f64,
}

impl Disk {
    pub fn make_shape(r: f64) -> Self {
        Disk { pos: [0.0, 0.0], radius: r }
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
    
    fn check_overlap(&self, other: &Self, offset: &[f64]) -> bool {
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

}
