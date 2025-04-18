use crate::schedule::{ObservableTracker, Schedule};
use crate::asc::Asc;
use std::fmt::{Debug, Display};
use rand_xoshiro::Xoshiro256StarStar;
use kiss3d::window::Window;

pub trait Particle 
where Self: Clone + Send + Sync + Debug + Display + std::marker::Sized {
    const TYPE: &'static str;

    // Parse from line in file
    fn parse(line: &str, unit_cell: &[f64]) -> Self;
    
    // offset is a vector to add to translational coord of 
    // other
    fn check_overlap(&self, other: &Self, offset: &[f64]) -> bool;

    // Gives the scaling factor V'/V needed to bring particles into contact
    fn overlap_scale(&self, other: &Self, offset: &[f64]) -> f64 { unimplemented!() }
    
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

    fn apply_strain(&mut self, new_cell: &[f64]);
   
    fn init_obs<C: Asc<Self>>(config: &C) -> ObservableTracker;

    fn sample_obs_sweep<C: Asc<Self>>(
        schedule: &mut Schedule<Self, C>,
        config: &C
    );
    
    fn sample_obs_failed_move<C: Asc<Self>>(
        schedule: &mut Schedule<Self, C>,
        config: &C
    );

    fn sample_obs_accepted_pmove<C: Asc<Self>>(
        schedule: &mut Schedule<Self, C>,
        config: &C,
        changed_idx: usize,
        old_p: &Self
    );

    fn sample_obs_accepted_cmove<C: Asc<Self>>(
        schedule: &mut Schedule<Self, C>,
        config: &C,
        old_c: &[f64]
    );

    // Return code is just the propagation of window.render()
    fn render_packing<C: Asc<Self>>(window: &mut Window, config: &C) -> bool { unimplemented!() }

    // Lower bound on dimension of particle
    fn hint_lower(&self) -> f64;

    // Upper bound on dimension of particle
    fn hint_upper(&self) -> f64;

    // Returns the particle volume
    fn vol(&self) -> f64 { unimplemented!() }

    // Returns relative translational coordinates
    // Type signature is a bit of a hack, will waste space
    // in 2d and cause problems if we try to implement
    // d > 3
    // TODO: proper type level integers
    fn lat_coord(&self) -> Vec<f64>;
}
