use crate::asc::Particle;
use rand_xoshiro::Xoshiro256StarStar;
use rand_distr::{Uniform, Distribution, Normal};
use std::fmt::{Display, Formatter};
use std::fmt;
use std::convert::TryInto;
use std::cmp::PartialEq;
use nalgebra::{Matrix3, Vector3, Quaternion, UnitQuaternion};
use rgsl::{Minimizer, MinimizerType, Value, minimizer};
use crate::asc::{Asc, save_asc_from_opt};
use crate::schedule::{Schedule, write_sweep_log};
use crate::{PI, OPT};

// https://stackoverflow.com/questions/26958178/how-do-i-automatically-implement-comparison-for-structs-with-floats-in-rust
#[derive(Debug, Clone, PartialEq)]
pub struct Ellipsoid {
    pos: [f64; 3], // [x, y, z]
    quat: [f64; 4], // a + bi + cj + dk -> [a, b, c, d]
    semi_axes: [f64; 3], // for no rot, [x, y, z]
}

impl Ellipsoid {
    pub fn make_shape(a: f64, b: f64, c: f64) -> Self {
        Ellipsoid { pos: [0.0, 0.0, 0.0], quat: [1.0, 0.0, 0.0, 0.0], semi_axes: [a, b, c] }
    }

    // c is the unit cell, given as (u_rc)
    // u_00 u_01 u_10 u_11 in 2d, etc.
    // Need to use nalgebra here
    fn apply_pbc(&mut self, c: &[f64]) {
        let u = Matrix3::from_column_slice(c);
        let u_inv = u.lu().try_inverse().expect("unit cell matrix must be invertible");
        let r = Vector3::from_column_slice(&self.pos);
        // convert to lattice coords
        let mut lat_c = u_inv*r;
        // put lattice coords back in unit square
        lat_c.apply(|x| x - x.floor());
        // convert back to euclidean coords
        // https://stackoverflow.com/questions/25428920/how-to-get-a-slice-as-an-array-in-rust
        self.pos = (u*lat_c).as_slice().try_into().unwrap();
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

impl Display for Ellipsoid {
    // From rust docs
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(f, "{} {} {} {} {} {} {} {} {} {}",
               self.pos[0], self.pos[1], self.pos[2],
               self.quat[0], self.quat[1], self.quat[2],
               self.quat[3], self.semi_axes[0], self.semi_axes[1],
               self.semi_axes[2]
        )
    }
}

// Overlap function defined in Donev. et. al. J. Comp. Phys. 202 (2005).
// for params, we have (X_A^-1, X_B^-1, r_AB)
// https://stackoverflow.com/questions/25272392/wrong-number-of-type-arguments-expected-1-but-found-0
fn pw_overlap(lambda: f64, params: &mut (Matrix3<f64>, Matrix3<f64>, Vector3<f64>)) -> f64 {
    let y_mat = lambda*params.0 + (1.0f64 - lambda)*params.1;
    let y_inv = y_mat.lu().try_inverse().expect("Must have invertible Y matrix");
    -(lambda*(1.0f64 - lambda)*params.2.dot(&(y_inv*params.2)) - 1.0f64)
}

impl Particle for Ellipsoid {
    
    const TYPE: &'static str = "Ellipsoid";
 
    fn parse(line: &str) -> Self {
        let params: Vec<f64> = line.split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect();
        Ellipsoid { pos: [params[0], params[1], params[2]], 
                    quat: [params[3], params[4], params[5], params[6]],
                    semi_axes: [params[7], params[8], params[9]] } 
    }

    fn check_overlap(&self, other: &Self, offset: &[f64]) -> bool {
        let image_x = other.pos[0] + offset[0];
        let image_y = other.pos[1] + offset[1];
        let image_z = other.pos[2] + offset[2];
        let disp_x = self.pos[0] - image_x;
        let disp_y = self.pos[1] - image_y;
        let disp_z = self.pos[2] - image_z;
        let disp2 = disp_x.powi(2) + disp_y.powi(2) + disp_z.powi(2);
        // Pre-check, as suggested in Hard Convex Body Fluids
        // Try and avoid evaluating PW potential by checking inner and outer
        // spheres
        if disp2 <= (self.minor_semi_axis() + other.minor_semi_axis()).powi(2) {
            if OPT.check_overlap {
                if *self != *other {
                    println!("Trivial overlap:\n\
                        self: {}\n\
                        other: {}",
                        self, other
                    );
                }
            }
            return true;
        }

        if disp2 > (self.major_semi_axis() + other.major_semi_axis()).powi(2) {
            return false;
        }
        // Do brent's algo
        // based on example
        // https://github.com/GuillaumeGomez/rust-GSL/blob/master/examples/minimization.rs
        let mut iter = 0usize;
        //let mut lambda;
        let mut lambda = 0.5;
        let mut upper = 1.0f64;
        let mut lower = 0.0f64;
        let mut status = Value::Continue;

        let min_type: MinimizerType<(Matrix3<f64>, Matrix3<f64>, Vector3<f64>)> = MinimizerType::brent();
        let mut min_instance: Minimizer<(Matrix3<f64>, Matrix3<f64>, Vector3<f64>)> = Minimizer::new(&min_type)
            .expect("Couldn't alloc minimizer");
        
        // Create X_A^-1 X_B^-1 and r_AB
        let self_quat = UnitQuaternion::from_quaternion(
            Quaternion::new(self.quat[0],
                            self.quat[1],
                            self.quat[2],
                            self.quat[3]));
        let other_quat = UnitQuaternion::from_quaternion(
            Quaternion::new(other.quat[0],
                            other.quat[1],
                            other.quat[2],
                            other.quat[3]));
        let self_rot = self_quat.to_rotation_matrix();
        let other_rot = other_quat.to_rotation_matrix();
        let self_diagonal = Vector3::new(self.semi_axes[0].powi(2),
                                         self.semi_axes[1].powi(2),
                                         self.semi_axes[2].powi(2));
        let self_shape = Matrix3::from_diagonal(&self_diagonal);
        let other_diagonal = Vector3::new(other.semi_axes[0].powi(2),
                                          other.semi_axes[1].powi(2),
                                          other.semi_axes[2].powi(2));
        let other_shape = Matrix3::from_diagonal(&other_diagonal);
        let x_a_inv = self_rot.matrix()*self_shape*(self_rot.matrix().transpose());
        let x_b_inv = other_rot.matrix()*other_shape*(other_rot.matrix().transpose());
        let r_ab = Vector3::new(disp_x, disp_y, disp_z);
        let mut param_tuple = (x_a_inv, x_b_inv, r_ab);
        let mut param_tuple2 = param_tuple.clone();
        //if pw_overlap(0.0, &mut param_tuple2) < pw_overlap(1.0, &mut param_tuple2) {
        //    lambda = 0.0 + OPT.brent_abs_tol;
        //} else {
        //    lambda = 1.0 - OPT.brent_abs_tol;
        //}
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
            if OPT.check_overlap {
                if status == Value::Success {
                    lambda = min_instance.x_minimum();
                    let pw_val = pw_overlap(lambda, &mut param_tuple2);
                    if pw_val < 0.0 {
                        if -pw_val <= OPT.near_overlap_tol {
                            // From stack overflow 29483365
                            // https://stackoverflow.com/questions/29483365/what-is-the-syntax-for-a-multiline-string-literal
                            println!("Near overlap from brent\n\
                                self: {}\n\
                                other: {}\n\
                                lambda: {}\n\
                                pw: {}",
                                self, other, lambda, pw_val
                            );
                        }
                        return false;
                    } else {
                        println!("Non-trivial overlap\n\
                            self: {}\n\
                            other: {}\n\
                            lambda: {}\n\
                            pw: {}",
                            self, other, lambda, pw_val
                        );
                        return true;
                    }
                }
            } else {
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
        }
        panic!("Error occurred in brent method or reached max iteration: {:?}", status);
    }
    
    fn copy_shape_random_coord(&self,
                               cell: &[f64],
                               rng: &mut Xoshiro256StarStar
    ) -> Self
    {
        const DIM: usize = 3;
        let uni_dist = Uniform::new(0.0, 1.0);
        let lat_x = uni_dist.sample(rng);
        let lat_y = uni_dist.sample(rng);
        let lat_z = uni_dist.sample(rng);
        // Sampling strategy from, among other refs, not everything I read is listed
        // https://math.stackexchange.com/questions/2178462/proof-of-generating-random-rotations-using-quaternions
        // https://math.stackexchange.com/questions/442418/random-generation-of-rotation-matrices
        // https://en.wikipedia.org/wiki/Rotation_matrix#Uniform_random_rotation_matrices
        // https://math.stackexchange.com/questions/184086/uniform-distributions-on-the-space-of-rotations-in-3d/184685
        // https://github.com/RobotLocomotion/drake/issues/10319
        // https://stackoverflow.com/questions/31600717/how-to-generate-a-random-quaternion-quickly
        // https://stats.stackexchange.com/questions/143903/how-to-generate-uniformly-random-orthogonal-matrices-of-positive-determinant/143970#143970
        // https://math.stackexchange.com/questions/2178462/proof-of-generating-random-rotations-using-quaternions
        let normal_dist = Normal::new(0.0f64, 1.0f64).unwrap();
        let mut q0 = normal_dist.sample(rng);
        let mut q1 = normal_dist.sample(rng);
        let mut q2 = normal_dist.sample(rng);
        let mut q3 = normal_dist.sample(rng);
        let qnorm = (q0*q0 + q1*q1 + q2*q2 + q3*q3).sqrt();
        q0 /= qnorm; q1 /= qnorm; q2 /= qnorm; q3 /= qnorm;
        // Turn lattice coords into euclidean coords
        Ellipsoid { pos: 
            [ lat_x*cell[DIM*0 + 0] + lat_y*cell[DIM*0 + 1] + lat_z*cell[DIM*0 + 2],
              lat_x*cell[DIM*1 + 0] + lat_y*cell[DIM*1 + 1] + lat_z*cell[DIM*1 + 2],
              lat_x*cell[DIM*2 + 0] + lat_y*cell[DIM*2 + 1] + lat_z*cell[DIM*2 + 2]],
            quat: [ q0, q1, q2, q3 ],
            semi_axes: self.semi_axes }
    }

    fn perturb(&mut self,
               cell: &[f64],
               param: &[f64],
               rng: &mut Xoshiro256StarStar,
    ) -> Self
    {
        let old_ell = self.clone();
        let normal_trans = Normal::new(0.0, param[0]).expect("Need valid translation sigma");
        let normal_rot = Normal::new(0.0, param[1]).expect("Need valid rotation sigma");
        // Apply translation
        for x in &mut self.pos {
            *x += normal_trans.sample(rng);
        }
        // Handle pbc
        self.apply_pbc(cell);
        // Apply rotation, strategy adapted from Frenkel and Mulder he of revolution paper
        let mut qnorm = 0.0;
        for q in &mut self.quat {
            *q += normal_rot.sample(rng);
            qnorm += (*q)*(*q);
        }
        // Normalize
        qnorm = qnorm.sqrt();
        for q in &mut self.quat {
            *q /= qnorm;
        }
        old_ell
    }

    // From S. Torquato and Y. Jiao PRE 80, 041104 (2009)
    // Also need nalgebra
    fn apply_strain(&mut self, old_cell: &[f64], new_cell: &[f64]) {
        // get old lattice coords
        let u_old_inv = Matrix3::from_column_slice(old_cell)
            .lu()
            .try_inverse()
            .expect("Unit cells must be invertible");
        let lat_c = u_old_inv*Vector3::from_column_slice(&self.pos);
        // set to new global coords
        let u_new = Matrix3::from_column_slice(new_cell);
        self.pos = (u_new*lat_c).as_slice().try_into().unwrap();
    }

    fn init_obs() -> Vec<f64> {
        vec![0.0, 0.0] // [samples, sum of volume]
    }

    fn sample_obs_sweep(schedule: &mut Schedule<Self>, config: &Asc<Self>) {
        let vol = schedule.running_obs[1]/schedule.running_obs[0];
        schedule.running_obs = vec![0.0,0.0];
        println!("Cell volume over sweep: {}", vol);
        let semi_prod: f64 = config.p_vec[0].semi_axes.iter().product();
        let phi = (config.p_vec.len() as f64)*4.0*PI*semi_prod/(3.0*vol);
        println!("Phi over sweep: {}", phi);
        let logline = format!("{} {} {}", schedule.current_sweep, vol, phi);
        write_sweep_log(&logline);
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
        _changed_idx: usize,
        _old_p: &Self
    )
    {
        schedule.running_obs[0] += 1.0;
        schedule.running_obs[1] += config.cell_volume();
    }
    
    fn sample_obs_accepted_cmove(
        schedule: &mut Schedule<Self>,
        config: &Asc<Self>,
        _old_c: &[f64]
    )
    {
        schedule.running_obs[0] += 1.0;
        schedule.running_obs[1] += config.cell_volume();
    }
}
