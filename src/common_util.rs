use nalgebra::{Matrix2, Vector2, Matrix3, Vector3};
use std::convert::TryInto;

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
}
