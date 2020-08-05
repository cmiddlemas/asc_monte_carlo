use crate::particle::Particle;
use rand_xoshiro::Xoshiro256StarStar;
use rand_distr::{Uniform, Distribution, Normal};
use std::fmt::{Display, Formatter};
use std::fmt;
use std::convert::TryInto;
use nalgebra::{Matrix2, Vector2};
use rgsl::{Minimizer, MinimizerType, Value, minimizer};
use crate::asc::{Asc, save_asc_from_opt};
use crate::schedule::{Schedule, write_sweep_log};
use crate::{PI,OPT};
use crate::common_util::apply_pbc;

// https://stackoverflow.com/questions/26958178/how-do-i-automatically-implement-comparison-for-structs-with-floats-in-rust
#[derive(Debug, Clone)]
pub struct Ellipse {
    pos: [f64; 2], // [x, y]
    theta: f64,
    semi_axes: [f64; 2], // for no rot, [x, y]
}

impl Ellipse {
    pub fn make_shape(a: f64, b: f64) -> Self {
        Ellipse { pos: [0.0, 0.0], theta: 0.0, semi_axes: [a, b] }
    }

    fn minor_semi_axis(&self) -> f64 {
        *self.semi_axes.iter()
            .min_by(|a, b| a.partial_cmp(b).expect("No NaNs"))
            .expect("Not empty")
    }

    fn major_semi_axis(&self) -> f64 {
        *self.semi_axes.iter()
            .max_by(|a, b| a.partial_cmp(b).expect("No NaNs"))
            .expect("Not empty")
    }
}

impl Display for Ellipse {
    // From rust docs
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(f, "{} {} {} {} {}",
               self.pos[0], self.pos[1], self.theta,
               self.semi_axes[0], self.semi_axes[1]
        )
    }
}

// Overlap function defined in Donev. et. al. J. Comp. Phys. 202 (2005).
// for params, we have (X_A^-1, X_B^-1, r_AB)
// https://stackoverflow.com/questions/25272392/wrong-number-of-type-arguments-expected-1-but-found-0
fn pw_overlap(lambda: f64, params: &mut (Matrix2<f64>, Matrix2<f64>, Vector2<f64>)) -> f64 {
    let y_mat = lambda*params.0 + (1.0f64 - lambda)*params.1;
    let y_inv = y_mat.lu().try_inverse().expect("Must have invertible Y matrix");
    -(lambda*(1.0f64 - lambda)*params.2.dot(&(y_inv*params.2)) - 1.0f64)
}

impl Particle for Ellipse {
    
    const TYPE: &'static str = "Ellipse";
 
    fn parse(line: &str, unit_cell: &[f64]) -> Self {
        let params: Vec<f64> = line.split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect();
        Ellipse { pos: [params[0], params[1]], 
                    theta: params[2],
                    semi_axes: [params[3], params[4]] } 
    }

    fn check_overlap(&self, other: &Self, offset: &[f64]) -> bool {
        let image_x = other.pos[0] + offset[0];
        let image_y = other.pos[1] + offset[1];
        let disp_x = self.pos[0] - image_x;
        let disp_y = self.pos[1] - image_y;
        let disp2 = disp_x.powi(2) + disp_y.powi(2);
        // Pre-check, as suggested in Hard Convex Body Fluids
        // Try and avoid evaluating PW potential by checking inner and outer
        // spheres
        if disp2 <= (self.minor_semi_axis() + other.minor_semi_axis()).powi(2) {
            return true;
        }
        if disp2 > (self.major_semi_axis() + other.major_semi_axis()).powi(2) {
            return false;
        }
        // Do brent's algo
        // based on example
        // https://github.com/GuillaumeGomez/rust-GSL/blob/master/examples/minimization.rs
        let mut iter = 0usize;
        let mut lambda;
        let mut upper = 1.0f64;
        let mut lower = 0.0f64;
        let mut status = Value::Continue;

        let min_type: MinimizerType<(Matrix2<f64>, Matrix2<f64>, Vector2<f64>)> = MinimizerType::brent();
        let mut min_instance: Minimizer<(Matrix2<f64>, Matrix2<f64>, Vector2<f64>)> = Minimizer::new(&min_type)
            .expect("Couldn't alloc minimizer");
        
        // Create X_A^-1 X_B^-1 and r_AB
        let self_rot = Matrix2::from_row_slice(
            &[self.theta.cos(), -self.theta.sin(),
            self.theta.sin(), self.theta.cos()]);
        let other_rot = Matrix2::from_row_slice(
            &[other.theta.cos(), -other.theta.sin(),
            other.theta.sin(), other.theta.cos()]);
        let self_diagonal = Vector2::new(self.semi_axes[0].powi(2),
                                         self.semi_axes[1].powi(2));
        let self_shape = Matrix2::from_diagonal(&self_diagonal);
        let other_diagonal = Vector2::new(other.semi_axes[0].powi(2),
                                          other.semi_axes[1].powi(2));
        let other_shape = Matrix2::from_diagonal(&other_diagonal);
        let x_a_inv = self_rot*self_shape*(self_rot.transpose());
        let x_b_inv = other_rot*other_shape*(other_rot.transpose());
        let r_ab = Vector2::new(disp_x, disp_y);
        let mut param_tuple = (x_a_inv, x_b_inv, r_ab);
        let mut param_tuple2 = param_tuple.clone();
        if pw_overlap(0.0, &mut param_tuple2) < pw_overlap(1.0, &mut param_tuple2) {
            lambda = 0.0 + OPT.brent_abs_tol;
        } else {
            lambda = 1.0 - OPT.brent_abs_tol;
        }
        min_instance.set(pw_overlap, &mut param_tuple, lambda, lower, upper);
        
        while status == Value::Continue && iter < OPT.brent_max_iter {
            iter += 1;
            status = min_instance.iterate();
            if status != Value::Success {
                panic!("Error occurred in brent method: {:?}", status);
            }
            lower = min_instance.x_lower();
            upper = min_instance.x_upper();
            status = minimizer::test_interval(lower, upper, OPT.brent_abs_tol, OPT.brent_rel_tol);
            if status == Value::Success || status == Value::Continue {
                lambda = min_instance.x_minimum();
                if pw_overlap(lambda, &mut param_tuple2) < 0.0 {
                    return false;
                } else {
                    if status == Value::Success {
                        return true;
                    }
                }
            }
        }
        panic!("Error occurred in brent method or reached max iteration: {:?}", status);
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
        let theta_dist = Uniform::new(0.0, 2.0*PI);
        let theta = theta_dist.sample(rng);
        // Turn lattice coords into euclidean coords
        Ellipse { pos: 
            [ lat_x*cell[DIM*0 + 0] + lat_y*cell[DIM*1 + 0],
              lat_x*cell[DIM*0 + 1] + lat_y*cell[DIM*1 + 1]],
            theta,
            semi_axes: self.semi_axes }
    }

    fn perturb(&mut self,
               cell: &[f64],
               param: &[f64],
               rng: &mut Xoshiro256StarStar,
    ) -> (Self, usize)
    {
        let old_ell = self.clone();
        let normal_trans = Normal::new(0.0, param[0]).expect("Need valid translation sigma");
        let normal_rot = Normal::new(0.0, param[1]).expect("Need valid rotation sigma");
        let uni_dist = Uniform::new(0.0, 1.0);
        // Equally likely to be a rotation or translation
        // 0 -> translation, 1 -> rotation
        let move_type = uni_dist.sample(rng) <= 0.5;
        // Apply translation
        if OPT.combined_move || !move_type {
            for x in &mut self.pos {
                *x += normal_trans.sample(rng);
            }
            self.pos = apply_pbc(&self.pos).as_slice().try_into().unwrap();
        }
        // Apply rotation
        if OPT.combined_move || move_type {
            self.theta += normal_rot.sample(rng);
            self.theta -= 2.0*PI*(self.theta/(2.0*PI)).floor();
        }

        if OPT.combined_move {
            (old_ell, 0)
        } else {
            (old_ell, move_type as usize)
        }
    }

    // From S. Torquato and Y. Jiao PRE 80, 041104 (2009)
    fn apply_strain(&mut self, new_cell: &[f64]) {
    }

    fn init_obs() -> Vec<f64> {
        vec![0.0, 0.0] // [samples, sum of volume]
    }

    fn sample_obs_sweep<C: Asc<Self>>(schedule: &mut Schedule<Self>, config: &C) {
        let vol = schedule.running_obs[1]/schedule.running_obs[0];
        schedule.running_obs = vec![0.0,0.0];
        println!("Cell volume over sweep: {}", vol);
        let semi_prod: f64 = config.first_particle().semi_axes.iter().product();
        let phi = (config.n_particles() as f64)*PI*semi_prod/(vol);
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
        self.minor_semi_axis()
    }

    fn hint_upper(&self) -> f64 {
        self.major_semi_axis()
    }

    fn lat_coord(&self) -> Vec<f64> { 
        unimplemented!()
    }
}
