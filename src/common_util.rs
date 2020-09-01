use nalgebra::{Matrix2, Vector2, Matrix3, Vector3};
use std::convert::TryInto;
use rand_distr::{Normal, Distribution};
use rand::Rng;
use rand_xoshiro::Xoshiro256StarStar;
use crate::schedule::Schedule;
use crate::OPT;

// For FFI
use libc::{c_double, size_t, c_int};

#[link(name = "gsl")]
extern "C" {
    fn gsl_fit_wlinear(x: *const c_double,
                       xstride: size_t,
                       w: *const c_double,
                       wstride: size_t,
                       y: *const c_double,
                       ystride: size_t,
                       n: size_t,
                       c0: *mut c_double,
                       c1: *mut c_double,
                       cov00: *mut c_double,
                       cov01: *mut c_double,
                       cov11: *mut c_double,
                       chisq: *mut c_double) 
        -> c_int;      
}

// TODO: Make return type non-dynamic for abstraction over dimension
// Applies periodic boundary conditions to a set of lattice coordinates
// See Ge's code, need to do this
// to put domain in [0.0, 1.0)
// Consider very small negative f64 values
pub fn apply_pbc(lat_coord: &[f64]) -> Vec<f64> {
    lat_coord.iter()
             .map(|x| x - x.floor())
             .map(|x| x - x.floor())
             .collect()
}

pub fn global_to_relative2(unit_cell: &[f64], global_coord: &[f64]) -> [f64; 2] {
    let u_inv = Matrix2::from_column_slice(unit_cell)
        .try_inverse()
        .unwrap();
    let g_r = Vector2::from_column_slice(global_coord);
    (u_inv*g_r).as_slice().try_into().unwrap()
}

pub fn relative_to_global2(unit_cell: &[f64], rel_coord: &[f64]) -> [f64; 2] {
    let u = Matrix2::from_column_slice(unit_cell);
    let l_r = Vector2::from_column_slice(rel_coord);
    (u*l_r).as_slice().try_into().unwrap()
}

pub fn global_to_relative3(unit_cell: &[f64], global_coord: &[f64]) -> [f64; 3] {
    let u_inv = Matrix3::from_column_slice(unit_cell)
        .try_inverse()
        .unwrap();
    let g_r = Vector3::from_column_slice(global_coord);
    (u_inv*g_r).as_slice().try_into().unwrap()
}

pub fn relative_to_global3(unit_cell: &[f64], rel_coord: &[f64]) -> [f64; 3] {
    let u = Matrix3::from_column_slice(unit_cell);
    let l_r = Vector3::from_column_slice(rel_coord);
    (u*l_r).as_slice().try_into().unwrap()
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

// Returns (Tr(strain), new_cell)
pub fn gen_random_strain<P>(dim: usize,
                         old_cell: &[f64],
                         schedule: &Schedule<P>,
                         rng: &mut Xoshiro256StarStar)
-> (f64, Vec<f64>) 
{
    // Set up random distributions
    let iso_dist = Normal::new(0.0, schedule.cell_param[0])
        .unwrap();
    let shear_dist = Normal::new(0.0, schedule.cell_param[1])
        .unwrap();
    let axi_dist = Normal::new(0.0, schedule.cell_param[2])
        .unwrap();
    
    let dimf = dim as f64;
    if OPT.exact_volume_step {
        let old_vol = volume(dim, old_cell);
        let v_step = iso_dist.sample(rng);
        let rel_change = ((v_step + old_vol)/old_vol).powf(1.0/dimf);
        let new_cell = old_cell.iter()
                               .map(|&x| x*rel_change)
                               .collect();
        (rel_change - 1.0, new_cell)
    } else if OPT.log_volume_step {
        let log_step = dimf*iso_dist.sample(rng);
        let new_cell = old_cell.iter()
                               .map(|&x| x*(log_step.exp().powf(1.0/dimf)))
                               .collect();
        (log_step, new_cell)
    } else {
        // Choose random strain
        match dim {
            2 => {
                // Choose strain
                let iso = iso_dist.sample(rng);
                let shear = shear_dist.sample(rng);
                let axi = axi_dist.sample(rng);
                let strain = Matrix2::from_row_slice(
                    &[iso + axi, shear, shear, iso - axi]
                );
                let trace_strain = 2.0*iso;
                // Change unit cell
                // Implicit transposition in read order
                // Noticed that this is necessary due to older
                // group code, courtesy of Duyu, Steve, and Yang
                let old_cell_mat = Matrix2::from_column_slice(&old_cell);
                let new_cell_mat = old_cell_mat + strain*old_cell_mat;

                (trace_strain, new_cell_mat.as_slice().to_vec())
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
                
                let strain = Matrix3::from_row_slice(&[d1, shear1, shear2,
                                                       shear1, d2, shear3,
                                                       shear2, shear3, d3]);
                let trace_strain = 3.0*iso;
                // Change unit cell
                // Implicit transposition in read order
                // Noticed that this is necessary due to older
                // group code, courtesy of Duyu, Steve, and Yang
                let old_cell_mat = Matrix3::from_column_slice(&old_cell);
                let new_cell_mat = old_cell_mat + strain*old_cell_mat;

                (trace_strain, new_cell_mat.as_slice().to_vec())
            }
            _ => unimplemented!(),
        }
    }
}

pub fn min_float_slice(slice: &[f64]) -> f64 {
    let mut current_min = f64::INFINITY;
    for f in slice {
        if *f < current_min {
            current_min = *f;
        }
    }
    current_min
}

pub fn max_float_slice(slice: &[f64]) -> f64 {
    let mut current_max = f64::NEG_INFINITY;
    for f in slice {
        if *f > current_max {
            current_max = *f;
        }
    }
    current_max
}

// Makes a linear fit y = a + bx to the given data
// Goes by least squares
// return
// (a, b, unc_a, unc_b, chisq)
// Currently implemented by calling out to GSL, may change
// later
// used many references regarding the unsafe portion
// but these were the most useful
// Also remember talking to Alex Ozdemir about something similar
// https://doc.rust-lang.org/nomicon/ffi.html
// https://doc.rust-lang.org/book/ch19-01-unsafe-rust.html
// https://s3.amazonaws.com/temp.michaelfbryan.com/arrays/index.html
pub fn linear_fit(data: &[[f64; 3]]) -> (f64, f64, f64, f64, f64) {
    let x: Vec<f64> = data.iter().map(|v| v[0]).collect();
    let y: Vec<f64> = data.iter().map(|v| v[1]).collect();
    let w: Vec<f64> = data.iter().map(|v| 1.0/(v[2]*v[2])).collect();
    let n = x.len();
    let mut c0 = 0.0;
    let mut c1 = 0.0;
    let mut cov00 = 0.0;
    let mut cov01 = 0.0;
    let mut cov11 = 0.0;
    let mut chisq = 0.0;
    unsafe { gsl_fit_wlinear(x.as_ptr(), 1, w.as_ptr(), 1, y.as_ptr(), 1, n,
        &mut c0, &mut c1, &mut cov00, &mut cov01, &mut cov11, &mut chisq) };
    (c0, c1, cov00.sqrt(), cov11.sqrt(), chisq)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn map_back_negative() {
        let x = [-1e-50, 2.0];
        let y = apply_pbc(&x);
        assert_eq!(y[0], 0.0);
        assert_eq!(y[1], 0.0);
    }

    #[test]
    fn inverse() {
        let cell = [1.4, 0.2, -3.8, 2.0];
        let l_r = [0.2, 0.3];
        let g_r = relative_to_global2(&cell, &l_r);
        let inverse = global_to_relative2(&cell, &g_r);
        let rel_err: Vec<f64> = l_r.iter().zip(inverse.iter())
            .map(|(x,y)| ((x-y)/x).abs() )
            .collect();
        assert!(rel_err[0] < 1e-14, "Error: {}", rel_err[0]);
        assert!(rel_err[1] < 1e-14, "Error: {}", rel_err[1]);
    }
    
    #[test]
    fn inverse2() {
        let cell = [1.4, 0.2, -3.8, 2.0];
        let g_r = [17.2, -100.3];
        let l_r = global_to_relative2(&cell, &g_r);
        let inverse = relative_to_global2(&cell, &l_r);
        let rel_err: Vec<f64> = g_r.iter().zip(inverse.iter())
            .map(|(x,y)| ((x-y)/x).abs() )
            .collect();
        assert!(rel_err[0] < 1e-14, "Error: {}", rel_err[0]);
        assert!(rel_err[1] < 1e-14, "Error: {}", rel_err[1]);
    }

    #[test]
    fn inverse_3d() {
        let cell = [1.4, 0.2, -4.5, -3.8, 2.0, 0.7, -20.0, 0.0, 30.5];
        let l_r = [0.2, 0.3, -0.7];
        let g_r = relative_to_global3(&cell, &l_r);
        let inverse = global_to_relative3(&cell, &g_r);
        let rel_err: Vec<f64> = l_r.iter().zip(inverse.iter())
            .map(|(x,y)| ((x-y)/x).abs() )
            .collect();
        assert!(rel_err[0] < 1e-14, "Error: {}", rel_err[0]);
        assert!(rel_err[1] < 1e-14, "Error: {}", rel_err[1]);
        assert!(rel_err[2] < 1e-14, "Error: {}", rel_err[2]);
    }
    
    // Note: using 0.2 and 2.0 as third value in g_r failed to
    // converge at the 1e-14 level. This is probably a true
    // reflection of the numerical stability. I've just decided
    // to go with the default algorithms provided by nalgebra.
    // Thus, these two tests act a regression test for numerical
    // accuracy provided by using nalgebra defaults.
    // Using Mathematica 12.0.0.0 confirms this suspicion,
    // but of course the exact numerics are not reproduced
    #[test]
    fn inverse2_3d() {
        let cell = [1.4, 0.2, -4.5, -3.8, 2.0, 0.7, -20.0, 0.0, 30.5];
        let g_r = [17.2, -100.3, 20.0];
        let l_r = global_to_relative3(&cell, &g_r);
        let inverse = relative_to_global3(&cell, &l_r);
        let rel_err: Vec<f64> = g_r.iter().zip(inverse.iter())
            .map(|(x,y)| ((x-y)/x).abs() )
            .collect();
        assert!(rel_err[0] < 1e-14, "Error: {}", rel_err[0]);
        assert!(rel_err[1] < 1e-14, "Error: {}", rel_err[1]);
        assert!(rel_err[2] < 1e-14, 
            "Error, val1, val2: {} {} {}",
            rel_err[2], g_r[2], inverse[2]
        );
    }

    // Version of above test with lighter tolerance
    #[test]
    fn inverse3_3d() {
        let cell = [1.4, 0.2, -4.5, -3.8, 2.0, 0.7, -20.0, 0.0, 30.5];
        let g_r = [17.2, -100.3, 2.0];
        let l_r = global_to_relative3(&cell, &g_r);
        let inverse = relative_to_global3(&cell, &l_r);
        let rel_err: Vec<f64> = g_r.iter().zip(inverse.iter())
            .map(|(x,y)| ((x-y)/x).abs() )
            .collect();
        assert!(rel_err[0] < 1e-14, "Error: {}", rel_err[0]);
        assert!(rel_err[1] < 1e-14, "Error: {}", rel_err[1]);
        assert!(rel_err[2] >= 1e-14, 
            "Error, val1, val2: {} {} {}",
            rel_err[2], g_r[2], inverse[2]
        );
        assert!(rel_err[2] < 1e-13, 
            "Error, val1, val2: {} {} {}",
            rel_err[2], g_r[2], inverse[2]
        );
    }


    #[test]
    fn volume2() {
        let cell = [2.0, 0.0, 0.0, 2.0];
        assert!(volume(2, &cell) == 4.0);
    }

    #[test]
    fn volume3() {
        let cell = [2.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0];
        assert!(volume(3, &cell) == 8.0);
    }
    
    #[test]
    fn slice_max() {
        let array = [0.2, 0.5, -0.2, 10.0, 2.0];
        assert!(max_float_slice(&array) == 10.0);
    }
    
    #[test]
    fn slice_min() {
        let array = [0.2, 0.5, -0.2, 10.0, 2.0];
        assert!(min_float_slice(&array) == -0.2);
    }

}
