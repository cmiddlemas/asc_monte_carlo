// asc.rs
// Author: Timothy Middlemas
// Implements ASC algorithm described in
// Jiao et. al. Nature 2009
// Generic over particle shape and simulation dimension
use std::fmt::{Debug, Display};
use rand_xoshiro::Xoshiro256StarStar;
use itertools::{Itertools, Position};

// Free helper functions

fn calc_offset(n: usize, dim: usize, cell: &[f64]) -> Vec<f64>
{
    let mut offset = Vec::with_capacity(dim);
    for i in 0..dim {
        let mut coord = 0.0;
        let mut k = n as i32;
        for j in 0..dim {
            let scalar = (k % 3 - 1) as f64;
            k /= 3;
            coord += scalar*cell[dim*i + j];
        }
        offset.push(coord);
    }
    offset
}

// cell stored as (dim*row + column)
pub struct Asc<P> {
    dim: usize, // Dimension of configuration
    overbox: usize, // # of overboxes, needed for small unit cell
    cell: Vec<f64>, // Unit Cell
    p_vec: Vec<P>, // List of particles
}

impl<P: Particle + Debug + Display> Asc<P> {
    
    // Makes a trivial, generic config for testing
    pub fn make() -> Asc<P> { 
        Asc { dim: 0, overbox: 0, cell: Vec::new(), p_vec: Vec::new()}
    }

    // Make an rsa config
    pub fn make_rsa(n: usize, // Number of particles to insert
                    shape: &P, // Contains the shape information
                    dim: usize, // dimension of config
                    init_cell: Vec<f64>, // Initial cell shape
                    rng: &mut Xoshiro256StarStar // Reproducible Rng
    ) -> Asc<P> {
        let mut new_asc = Asc 
            { dim: dim, overbox: 1, cell: init_cell.clone(), p_vec: Vec::new() };
        while new_asc.p_vec.len() < n {
            if new_asc.try_add_particle(
                shape.copy_shape_random_coord(
                    &init_cell, rng
                )
            ) {
                println!("Added particle {}", new_asc.p_vec.len());
            }
        }
        new_asc
    }
    
    // Dumps all the info contained in Asc
    pub fn debug(&self) {
        println!("dim: {}", self.dim);
        println!("{:?}", &self.cell);
        for p in &self.p_vec{
            println!("{:?}", p);
        }
    }

    // Prints the config in a nice, ascii format
    pub fn print_asc(&self) {
        println!("{}", self.dim);
        for entry in self.cell.iter().with_position() {
            match entry {
                Position::Last(x) => println!("{}", x),
                Position::Middle(x) => print!("{} ", x),
                Position::First(x) => print!("{} ", x),
                Position::Only(x) => println!("{}", x),
            }
        }
        for p in &self.p_vec {
            println!("{}", p);
        }
    }

    // returns number of overlaps
    pub fn check_particle(&self, fixed: &P) -> usize {
        let mut count = 0;
        let n_offset = (2*self.overbox + 1).pow(self.dim as u32);
        let offset_vec: Vec<Vec<f64>>  = (0..n_offset)
            .map(|x| calc_offset(x, self.dim, &self.cell))
            .collect();
        for imaged in &self.p_vec {
            for i in 0..n_offset {
                // https://stackoverflow.com/questions/55461617/how-do-i-convert-a-boolean-to-an-integer-in-rust
                count += fixed.check_overlap(
                    imaged,
                    &offset_vec[i]
                ) as usize;
            }
        }
        //println!("{}", count);
        count
    }

    // returns true if add is successful
    pub fn try_add_particle(&mut self, p: P) -> bool {
        if self.check_particle(&p) > 0 {
            return false;
        }
        self.p_vec.push(p);
        true
    }
}

pub trait Particle {
    // offset is a vector to add to translational coord of 
    // other
    fn check_overlap(&self, other: &Self, offset: &[f64]) -> bool;
    fn copy_shape_random_coord(&self,
                               cell: &[f64],
                               rng: &mut Xoshiro256StarStar) 
    -> Self;
}
