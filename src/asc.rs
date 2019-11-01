// asc.rs
// Author: Timothy Middlemas
// Implements ASC algorithm described in
// Jiao et. al. Nature 2009
// Generic over particle shape and simulation dimension
use std::fmt::{Debug, Display};
use rand_xoshiro::Xoshiro256StarStar;
use rand::Rng;
use rand::seq::SliceRandom;
use itertools::{Itertools, Position};
use std::path::Path;
use std::fs::{File, OpenOptions};
use std::io::{Write, BufWriter};
use crate::OPT;
use rayon::prelude::*;
use std::ops::Range;
use rand_distr::{Uniform, Normal, Distribution};
use nalgebra::{Matrix3};
use crate::schedule::Schedule;

// Free helper functions

fn calc_offset(n: usize, dim: usize, overbox: usize, cell: &[f64]) -> Vec<f64>
{
    let mut offset = Vec::with_capacity(dim);
    let copies = 2*(overbox as i32) + 1;
    for i in 0..dim {
        let mut coord = 0.0;
        let mut k = n as i32;
        for j in 0..dim {
            let scalar = (k % copies - 1) as f64;
            k /= copies;
            coord += scalar*cell[dim*i + j];
        }
        offset.push(coord);
    }
    offset
}

// cell stored as (dim*row + column)
// and columns are interpreted as the
// lattice vectors
#[derive(Clone)]
pub struct Asc<P> {
    dim: usize, // Dimension of configuration
    overbox: usize, // # of overboxes, needed for small unit cell
    pub cell: Vec<f64>, // Unit Cell
    pub p_vec: Vec<P>, // List of particles
}

// https://old.reddit.com/r/rust/comments/98jldy/question_about_rayon_and_intopariterator_not/
// https://github.com/rust-lang/rust/issues/61768
impl<P: Particle + Debug + Display + Send + Sync + Clone> Asc<P> {
    
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
            ) {} 
            /* Some code for printing out RSA data, unlikely to want to use this
             * in a real MC run
            {
                println!("Added particle {}", new_asc.p_vec.len());
                if let Some(path) = &OPT.savefiles {
                    if path.is_dir() {
                        println!("Root filename cannot be empty/ you specified a dir. Skipping save.");
                        continue;
                    }
                    let mut full_path = path.clone();
                    // https://users.rust-lang.org/t/what-is-right-ways-to-concat-strings/3780/4
                    full_path.set_file_name(
                        format!("{}_rsa_{}.dat", 
                                path.file_name().expect("Must give a valid root filename.")
                                    .to_str().expect("Must give valid UTF-8 str."), 
                                new_asc.p_vec.len()
                    ));
                    new_asc.save_asc(&full_path);
                }
            }
            */
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
        println!("{} {} {}", self.dim, self.overbox, P::TYPE);
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

    // Saves the config in a nice, ascii format given path
    // Panics on any error
    // Tries to wait until data hits disk
    // https://doc.rust-lang.org/std/io/struct.BufWriter.html
    pub fn save_asc(&self, path: &Path) {
        let mut file = BufWriter::new(File::create(path).expect("Must specify valid path to save to."));
        writeln!(&mut file, "{} {} {}", self.dim, self.overbox, P::TYPE).expect("Failed write during save.");
        for entry in self.cell.iter().with_position() {
            match entry {
                Position::Last(x) => writeln!(&mut file, "{}", x),
                Position::Middle(x) => write!(&mut file, "{} ", x),
                Position::First(x) => write!(&mut file, "{} ", x),
                Position::Only(x) => writeln!(&mut file, "{}", x),
            }.expect("Failed write during save.");
        }
        for p in &self.p_vec {
            writeln!(&mut file, "{}", p).expect("Failed write during save.");
        }
        let mut f = file.into_inner().expect("Failed to unwrap buffer during save");
        f.flush().expect("Failed to flush file writer during save");
        f.sync_all().expect("Failed to sync during save.");
        let mut dir = OpenOptions::new()
            .read(true)
            .open(path.parent().expect("Must have parent directory"))
            .expect("Failed to open parent directory");
        dir.sync_all().expect("Failed to sync directory during save.");
    }

    // returns number of overlaps
    pub fn check_particle(&self, fixed: &P) -> usize {
        let n_offset = (2*self.overbox + 1).pow(self.dim as u32);
        
        (0..n_offset)
            .map(|x| calc_offset(x, self.dim, self.overbox, &self.cell))
            .map(|offset|
                self.p_vec.iter()
                    .map(|imaged| fixed.check_overlap(imaged, &offset) as usize)
                    .sum::<usize>()
                )
            .sum()

            // For future reference, to parallelize range:
            // https://users.rust-lang.org/t/rayon-parallel-sum-from-range/6367

            // Obsolete comment, was for when single particle check was parallelized
            // For some reason faster than into_par_iter(), also see
            // https://users.rust-lang.org/t/for-loops-in-rust/8217/4

    /* Old code, reversed inner and outer loops
        self.p_vec.par_iter().map(|imaged|
            // https://users.rust-lang.org/t/auto-vectorization-in-rust/24379/2
            offset_vec.iter().map(|offset|
                // https://stackoverflow.com/questions/55461617/how-do-i-convert-a-boolean-to-an-integer-in-rust
                fixed.check_overlap(imaged, offset) as usize
            // https://stackoverflow.com/questions/51283403/cannot-infer-type-for-b-for-filter-map-sum
            ).sum::<usize>()
        ).sum()
    */
    }

    // returns true if add is successful
    pub fn try_add_particle(&mut self, p: P) -> bool {
        if self.check_particle(&p) > 0 {
            return false;
        }
        self.p_vec.push(p);
        true
    }

    // TODO: Maybe should replace by appropriate nalgebra calls?
    pub fn cell_volume(&self) -> f64 {
        match self.dim {
            2 => { // Simple formula for determinant
                (self.cell[0]*self.cell[3]
                    - self.cell[1]*self.cell[2]).abs()
            }
            3 => {
                let u = Matrix3::from_row_slice(&self.cell);
                (u.row(0).dot(&u.row(1).cross(&u.row(2)))).abs()
            }
            _ => unimplemented!(),
        }
    }

    fn is_valid(&self) -> bool {
        self.p_vec.par_iter()
            .map(|p| self.check_particle(p))
            .all(|x| x <= 1)
    }

    pub fn try_cell_move(&mut self,
                         schedule: &Schedule<P>,
                         rng: &mut Xoshiro256StarStar
    ) -> bool
    {
        let iso_dist = Normal::new(0.0, schedule.cell_param[0])
            .unwrap();
        let shear_dist = Normal::new(0.0, schedule.cell_param[1])
            .unwrap();
        let axi_dist = Normal::new(0.0, schedule.cell_param[2])
            .unwrap();
        let uni_dist = Uniform::new(0.0, 1.0); // for probabilities

        let old_asc = self.clone();
        
        match self.dim {
            2 => {
                // Choose strain
                let iso = iso_dist.sample(rng);
                let shear = shear_dist.sample(rng);
                let axi = axi_dist.sample(rng);
                let strain = [iso + axi, shear, shear, iso - axi];
                // Change unit cell
                for (e, s) in self.cell.iter_mut().zip(strain.iter()) {
                    *e += *s;
                }
            }
            3 => {
                // Choose strain
                let iso = iso_dist.sample(rng);
                let shear1 = shear_dist.sample(rng);
                let shear2 = shear_dist.sample(rng);
                let shear3 = shear_dist.sample(rng);
                let axi1 = axi_dist.sample(rng);
                let axi2 = axi_dist.sample(rng);
                
                let axis_choice: usize = rng.gen_range(0,3);
                let d1; let d2; let d3;
                match axis_choice {
                    0 => {d1 = iso + axi1; d2 = iso + axi2; d3 = iso - axi1 - axi2},
                    1 => {d1 = iso - axi1 - axi2; d2 = iso + axi1; d3 = iso + axi2},
                    2 => {d1 = iso + axi2; d2 = iso - axi1 - axi2; d3 = iso + axi1},
                    _ => unreachable!(),
                }
                
                let strain = [d1, shear1, shear2,
                              shear1, d2, shear3,
                              shear2, shear3, d3];
                // Change unit cell
                for (e, s) in self.cell.iter_mut().zip(strain.iter()) {
                    *e += *s;
                }
            }
            _ => unimplemented!(),
        }
        // Kinda weird, need to do ref outside of closure
        // https://stackoverflow.com/questions/48717833/how-to-use-struct-self-in-member-method-closure
        let new_cell = &self.cell;
        // Change particles
        self.p_vec.par_iter_mut().for_each(|p| {
            p.apply_strain(&old_asc.cell, new_cell);
        });
        if self.is_valid() { // Accept probabilistically
            // http://www.pages.drexel.edu/~cfa22/msim/node31.html
            let new_vol = self.cell_volume();
            let old_vol = old_asc.cell_volume();
            let n_particles = self.p_vec.len() as f64;
            let vol_factor = 
                (-schedule.beta*schedule.pressure*(new_vol - old_vol)
                 +n_particles*(new_vol/old_vol).ln()).exp();
            if uni_dist.sample(rng) < vol_factor { //Keep config
                true
            } else { //Reset config
                *self = old_asc;
                false
            }
        } else { // Reset config
            *self = old_asc;
            false
        }
    }

    pub fn try_particle_move(&mut self,
                             schedule: &mut Schedule<P>,
                             rng: &mut Xoshiro256StarStar
    ) -> bool
    {
        let r_idx: usize = Uniform::new(0, self.p_vec.len())
            .sample(rng);
        let old_p = self.p_vec[r_idx]
            .perturb(&self.cell, &schedule.particle_param, rng);
        if self.check_particle(&self.p_vec[r_idx]) > 1 { //reject
            self.p_vec[r_idx] = old_p; //roll back
            P::sample_obs_failed_move(
                schedule,
                self
            );
            false
        } else { //accept
            P::sample_obs_accepted_pmove(
                schedule,
                self,
                r_idx,
                &old_p
            );
            true
        }
    }
}

pub trait Particle 
where Self: std::clone::Clone + std::marker::Sized {
    const TYPE: &'static str;

    // offset is a vector to add to translational coord of 
    // other
    fn check_overlap(&self, other: &Self, offset: &[f64]) -> bool;
    
    fn copy_shape_random_coord(&self,
                               cell: &[f64],
                               rng: &mut Xoshiro256StarStar) 
    -> Self;
    
    // Returns a copy of the original coordinates,
    // changes the given object
    fn perturb(&mut self,
               cell: &[f64],
               param: &[f64],
               rng: &mut Xoshiro256StarStar
    ) -> Self;

    fn apply_strain(&mut self, old_cell: &[f64], new_cell: &[f64]);
   
    fn init_obs() -> Vec<f64>;

    fn sample_obs_sweep(
        schedule: &mut Schedule<Self>,
        config: &Asc<Self>
    );
    
    fn sample_obs_failed_move(
        schedule: &mut Schedule<Self>,
        config: &Asc<Self>
    );

    fn sample_obs_accepted_pmove(
        schedule: &mut Schedule<Self>,
        config: &Asc<Self>,
        changed_idx: usize,
        old_p: &Self
    );

    fn sample_obs_accepted_cmove(
        schedule: &mut Schedule<Self>,
        config: &Asc<Self>,
        old_c: &[f64]
    );
}
