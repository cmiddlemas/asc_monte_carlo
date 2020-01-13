// asc.rs
// Author: Timothy Middlemas
// Defines Asc trait, to be generic
// over list of particles implementation
// and Particle, to be generic over shape
use std::fmt::{Debug, Display};
use rand_xoshiro::Xoshiro256StarStar;
use std::path::Path;
use crate::OPT;
use crate::schedule::Schedule;
use nalgebra::Matrix3;
use std::fs::{OpenOptions, rename};

pub fn save_asc_from_opt<C, P: Particle + Debug + Display + Send + Sync + Clone>
    (config: &C, annotation: &str)
    where C: Asc<P>
{
    if let Some(path) = &OPT.savefiles {
        if path.is_dir() {
            println!("Root filename cannot be empty/ you specified a dir. Skipping save.");
        } else {
            let mut full_path = path.clone();
            
            if (!OPT.save_trajectory) && (annotation != "initial" && annotation != "final") {
                full_path.set_file_name(
                    format!("{}_backup_temp.dat",
                        path.file_name().expect("Must give a valid root filename.")
                            .to_str().expect("Must give valid UTF-8 str.")
                ));
                config.save_asc(&full_path, Some(annotation));
                let mut final_path = path.clone();
                final_path.set_file_name(
                    format!("{}_backup.dat",
                        path.file_name().expect("Must give a valid root filename.")
                            .to_str().expect("Must give valid UTF-8 str.")
                ));
                // This code tries and get a little bit of safety
                // https://unix.stackexchange.com/questions/297632/is-it-broken-to-replace-an-existing-file-without-fsync
                // https://lwn.net/Articles/457667/
                // https://stackoverflow.com/questions/3764822/how-to-durably-rename-a-file-in-posix
                rename(&full_path, &final_path).expect("Rename should be successful");
                let canonical = final_path.canonicalize().expect("Must be able to canonicalize");
                let dir = OpenOptions::new()
                    .read(true)
                    .open(canonical.parent().expect("Canonical form should have parent"))
                    .expect("Failed to open parent directory");
                dir.sync_all().expect("Failed to sync directory during save.");
            } else {
            
                // https://users.rust-lang.org/t/what-is-right-ways-to-concat-strings/3780/4
                full_path.set_file_name(
                    format!("{}_{}.dat",
                        path.file_name().expect("Must give a valid root filename.")
                            .to_str().expect("Must give valid UTF-8 str."),
                        annotation
                ));
                config.save_asc(&full_path, None);
            }
        }
    }
}

// TODO: maybe make 2d with nalgebra too?
pub fn volume(dim: usize, unit_cell: &[f64]) -> f64 {
    match dim {
        2 => { // Simple formula for determinant
            (unit_cell[0]*unit_cell[3]
                - unit_cell[1]*unit_cell[2]).abs()
        }
        3 => {
            let u = Matrix3::from_row_slice(unit_cell);
            (u.row(0).dot(&u.row(1).cross(&u.row(2)))).abs()
        }
        _ => unimplemented!(),
    }
}

pub trait Asc<P> 
where P: Particle + Debug + Display + Send + Sync + Clone 
{
    // Print configuration onto stdout
    fn print_asc(&self);

    // Save configuration to file
    fn save_asc(&self, path: &Path, annotation: Option<&str>);

    // Check particle over overlaps in Asc, return number of overlaps
    fn check_particle(&self, fixed: &P) -> usize;

    // Return cell volume
    fn cell_volume(&self) -> f64;

    // True if no overlaps
    fn is_valid(&self) -> bool;

    // Try to change the cell by straining
    fn try_cell_move(&mut self, schedule: &Schedule<P>, rng: &mut Xoshiro256StarStar) -> bool;

    // Try to move a particle
    fn try_particle_move(&mut self, schedule: &mut Schedule<P>, rng: &mut Xoshiro256StarStar) -> bool;

    // Return number of particles in Asc
    fn n_particles(&self) -> usize;

    // Return a reference to the first particle stored in Asc
    fn first_particle(&self) -> &P;
}

pub trait Particle 
where Self: Clone + Send + Sync + Debug + Display + std::marker::Sized {
    const TYPE: &'static str;

    // Parse from line in file
    fn parse(line: &str) -> Self;
    
    // offset is a vector to add to translational coord of 
    // other
    fn check_overlap(&self, other: &Self, offset: &[f64]) -> bool;
    
    fn copy_shape_random_coord(&self,
                               cell: &[f64],
                               rng: &mut Xoshiro256StarStar) 
    -> Self;
    
    // Returns a copy of the original coordinates,
    // changes the given object
    // Also returns the move_type, which is
    // a usize corresponding to which particle param
    // the move corresponded to, or always zero if
    // OPT.combined_move is true
    fn perturb(&mut self,
               cell: &[f64],
               param: &[f64],
               rng: &mut Xoshiro256StarStar
    ) -> (Self, usize);

    fn apply_strain(&mut self, old_cell: &[f64], new_cell: &[f64]);
   
    fn init_obs() -> Vec<f64>;

    fn sample_obs_sweep<C: Asc<Self>>(
        schedule: &mut Schedule<Self>,
        config: &C
    );
    
    fn sample_obs_failed_move<C: Asc<Self>>(
        schedule: &mut Schedule<Self>,
        config: &C
    );

    fn sample_obs_accepted_pmove<C: Asc<Self>>(
        schedule: &mut Schedule<Self>,
        config: &C,
        changed_idx: usize,
        old_p: &Self
    );

    fn sample_obs_accepted_cmove<C: Asc<Self>>(
        schedule: &mut Schedule<Self>,
        config: &C,
        old_c: &[f64]
    );

    // Lower bound on dimension of particle
    fn hint_lower(&self) -> f64;

    // Upper bound on dimension of particle
    fn hint_upper(&self) -> f64;

    // Returns relative translational coordinates
    // Type signature is a bit of a hack, will waste space
    // in 2d and cause problems if we try to implement
    // d > 3
    // TODO: proper type level integers
    fn lat_coord(&self, cell: &[f64]) -> Vec<f64>;
}
