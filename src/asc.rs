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
use crate::particle::Particle;
use crate::common_util::{min_float_slice, max_float_slice, linear_fit};
use std::fs::{OpenOptions, rename};
use rand_distr::{Uniform, Distribution};

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

pub trait Asc<P> 
where 
    P: Particle + Debug + Display + Send + Sync + Clone,
    Self: std::marker::Sized + Clone
{
    // Print configuration onto stdout
    fn print_asc(&self);

    // Save configuration to file
    fn save_asc(&self, path: &Path, annotation: Option<&str>);

    // Check particle over overlaps in Asc, return number of overlaps
    fn check_particle(&self, fixed: &P) -> usize;

    // Give nearest particle gap in terms of delta_phi
    // Need to give the unique index of the particle as well for
    // exclusion purposes
    fn particle_gap(&self, exclude_idx: usize, particle: &P) -> f64 { unimplemented!() }

    // Return cell volume
    fn cell_volume(&self) -> f64;

    // True if no overlaps
    fn is_valid(&self) -> bool;

    // Applies a random strain to configuration. Does not check overlap, but will rebuild
    // acceleration structures if necessary
    // Returns an option with (Asc object, trace_strain), if None, it means that the
    // strain attempt failed to produce a valid object, and the given object is consumed.
    fn apply_random_strain(self, schedule: &Schedule<P>, rng: &mut Xoshiro256StarStar) -> Option<(Self, f64)>;

    // Try to change the cell by straining
    fn try_cell_move(&mut self, schedule: &mut Schedule<P>, rng: &mut Xoshiro256StarStar) -> bool {
        let new_asc = self.clone();

        // This implementation switched because cannot drop through an &mut
        // https://doc.rust-lang.org/book/ch15-03-drop.html
        if let Some((new_asc, trace_strain)) = new_asc.apply_random_strain(schedule, rng) {
            let uni_dist = Uniform::new(0.0, 1.0); // for probabilities
            
            if new_asc.is_valid() { // Accept probabilistically
                // http://www.pages.drexel.edu/~cfa22/msim/node31.html
                let new_vol = new_asc.cell_volume();
                let old_vol = self.cell_volume();
                let n_particles = new_asc.n_particles() as f64;
                let vol_factor = if OPT.log_volume_step {
                    // This taken from Frenkel and Smit UMS 2nd Ed. (2002)
                    // among other places
                    (-schedule.beta*schedule.pressure*(new_vol - old_vol)
                        +(n_particles + 1.0)*((new_vol/old_vol).ln())).exp()
                } else {
                    if OPT.linear_acceptance {
                        (-schedule.beta*schedule.pressure*old_vol*trace_strain
                        +n_particles*((1.0 + trace_strain).ln())).exp()
                    } else {
                        (-schedule.beta*schedule.pressure*(new_vol - old_vol)
                        +n_particles*((new_vol/old_vol).ln())).exp()
                    }
                };
                if uni_dist.sample(rng) < vol_factor { //Keep config
                    //eprintln!("Old config! {}", self.cell_volume());
                    //eprintln!("Look out, a config: {}", new_asc.cell_volume());
                    P::sample_obs_accepted_cmove(
                        schedule,
                        &new_asc,
                        self.unit_cell()
                    );
                    // Set self to successful configuration
                    *self = new_asc;
                    true
                } else {
                    P::sample_obs_failed_move(
                        schedule,
                        self
                    );
                    false
                }
            } else {
                P::sample_obs_failed_move(
                    schedule,
                    self
                );
                false
            }
        } else {
            // Assume that the failure of apply_random_strain is due to
            // a configuration that would have been invalid anyway, if
            // it was possible for the underlying data structure to keep
            // working
            // Log a warning, since this is unexpected, but not necessarily incorrect
            // behavior, it is likely due to setting cell changes too large
            eprintln!("Warning: failure to strain cell because data structure became invalid");
            false
        }
    }

    // Try to move a particle
    fn try_particle_move(&mut self, schedule: &mut Schedule<P>, rng: &mut Xoshiro256StarStar) -> bool;

    // Return number of particles in Asc
    fn n_particles(&self) -> usize;

    // Return phi
    fn packing_fraction(&self) -> f64 {
        let sum_particle_vol: f64 = self.particle_slice().iter()
                                   .map(|p| p.vol())
                                   .sum();
        sum_particle_vol/self.cell_volume()
    }

    // Return a slice that gives the unit cell
    fn unit_cell(&self) -> &[f64];

    // Return a reference to the first particle stored in Asc
    fn first_particle(&self) -> &P;

    // Return a slice of particles
    fn particle_slice(&self) -> &[P] { unimplemented!() }

    // Returns nearest-neighbor gap distribution function
    // in form [delta phi, P(delta phi), unc P(delta phi)]
    // Working on this using the ASC code Duyu gave me
    // Also returns average nn gap in first variable, smallest
    // gap in second, and largest in third
    fn nn_gap_distribution(&self) -> (f64, f64, f64, Vec<[f64; 3]>) {
        let delta_phi: Vec<f64> = self.particle_slice().iter().enumerate()
                       .map(|(i, p)| self.particle_gap(i, p))
                       .collect();
        
        let n_obs = delta_phi.len() as f64;
        // https://stackoverflow.com/questions/56872714/iter-sum-stops-working-as-soon-as-i-add-further-operation
        let avg_nn_gap: f64 = delta_phi.iter().sum::<f64>()/n_obs;
        let max_delta_phi = max_float_slice(&delta_phi);
        let min_delta_phi = min_float_slice(&delta_phi);
        let delta_phi_step = (max_delta_phi-min_delta_phi)/(OPT.n_bins_gap as f64);
        
        let midpoint: Vec<f64> = (0..OPT.n_bins_gap).map(|x|
                                                          min_delta_phi
                                                            + (x as f64)*delta_phi_step 
                                                            + delta_phi_step/2.0)
                                                    .collect();
        let mut histogram = vec![0.0; OPT.n_bins_gap];
        
        for obs in delta_phi {
            // Safely cast float to int
            // https://github.com/rust-lang/rust/issues/10184
            let proposed_bin_f64 = (obs - min_delta_phi)/delta_phi_step;
            // eprintln!("{} {}", min_delta_phi, max_delta_phi);
            assert!(proposed_bin_f64.is_finite());
            assert!(proposed_bin_f64 >= 0.0 && proposed_bin_f64 <= 10000000.0);
            let proposed_bin = proposed_bin_f64 as usize;
            if proposed_bin < OPT.n_bins_gap {
                histogram[proposed_bin] += 1.0;
            }
        }

        let distribution = midpoint.iter().zip(histogram.iter())
           .map(|(x, y)| [*x, (*y)/n_obs, (*y).sqrt()/n_obs])
           .collect();
        
        (avg_nn_gap, min_delta_phi, max_delta_phi, distribution)
    }

    // Returns the instantaneous pressure, its uncertainty, chisq,
    // and R^2 of fit
    fn instantaneous_pressure(&self) -> (f64, f64, f64) {
        let (_, _, _, gap_distr) = self.nn_gap_distribution();
        // Need to fit ln P1(delta phi) as indicated in
        // Eppenga and Frenkel, Mol. Phys. 52 (1984)
        // and Duyu's paper on bulk equilibrium and MRJ polyhedra
        let ln_gap: Vec<[f64; 3]> = gap_distr.iter().map(|x| [x[0], x[1].ln(), x[2]/x[1]])
                                     .collect();
        let trimmed_data = &ln_gap[0..OPT.n_bins_fit];
        let (_, b, _, unc_b, chisq) = linear_fit(trimmed_data);
        
        let phi = self.packing_fraction();
        let alpha = -b;
        let unc_alpha = unc_b;
        let p = 1.0 + phi*alpha/2.0;
        let unc_p = phi*unc_alpha/2.0;
        (p, unc_p, chisq)
    }
}
