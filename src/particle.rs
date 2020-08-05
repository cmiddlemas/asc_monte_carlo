use crate::schedule::Schedule;
use crate::asc::Asc;
use std::fmt::{Debug, Display};
use rand_xoshiro::Xoshiro256StarStar;

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
