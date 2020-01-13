// cell_list.rs
// Author: Timothy Middlemas
// Implements the Asc trait using an underlying
// cell list representation. Performs automatic resizing
// Assumes monodispersity for now
// Is probably similar to Ge's cell list implementation,
// which I have access to, but this is in native rust
use crate::asc::{Asc, Particle, volume};
use crate::overbox_list::OverboxList;
use std::fmt::{Debug, Display};
use nalgebra::{Matrix2, Matrix3, Vector2, Vector3};
use crate::OPT;
use crate::schedule::Schedule;
use rand_xoshiro::Xoshiro256StarStar;
use std::path::Path;
use itertools::{Itertools, Position};
use std::fs::{File, OpenOptions};
use std::io::{Write, BufWriter};
use std::convert::TryInto;
use rayon::prelude::*;
use rand_distr::{Uniform, Normal, Distribution};
use rand::Rng;

#[derive(Clone)]
pub struct CellList<P> {
    dim: usize,
    pub unit_cell: Vec<f64>, // Unit cell, stored as (dim*row + column), each row is a lattice vector
    pub p_list: Vec<P>,
    n_particles: usize,
    lin_subdiv: usize,
    c_assoc_list: Vec<Vec<usize>>, // Association lists given as [cell1(p1...), cell2(p1...)...]
    p_assoc_list: Vec<usize>, // Gives first (cell identity) index into c_assoc_list, same order as p_list
}

// Returns the maximum valid value of lin_subdiv given a dimension,
// the current unit cell, and an upper bound on the radius of the 
// circumscribing sphere of the largest particle in the packing
// References:
// https://scicomp.stackexchange.com/questions/3107/minimum-image-convention-for-triclinic-unit-cell
// http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.57.1696&rep=rep1&type=pdf
// https://hoomd-blue.readthedocs.io/en/stable/nlist.html
fn max_subdiv(dim: usize, unit_cell: &[f64], max_radius: f64) -> usize {
    let vol = volume(dim, unit_cell);
    let max_diam = 2.0*max_radius;
    match dim {
        2 => {
            let u = Matrix2::from_column_slice(unit_cell);
            let base0 = u.column(0).norm();
            let base1 = u.column(1).norm();
            let max_base = base0.max(base1);
            let min_width = vol/max_base;
            // https://stackoverflow.com/questions/37506672/convert-float-to-integer-in-rust
            (min_width/max_diam).floor() as usize
        }
        3 => {
            let u = Matrix3::from_column_slice(unit_cell);
            let base0 = u.column(0).cross(&u.column(1)).norm();
            let base1 = u.column(1).cross(&u.column(2)).norm();
            let base2 = u.column(2).cross(&u.column(0)).norm();
            let max_base = base0.max(base1.max(base2));
            let min_width = vol/max_base;
            (min_width/max_diam).floor() as usize
        }
        _ => unimplemented!(),
    }
}

// Cell layout:
// 7 8 9
// 4 5 6
// 1 2 3
fn assign_cell<P: Particle>
    (dim: usize,
    lin_subdiv: usize,
    unit_cell: &[f64],
    p: &P) 
-> usize 
{
    let ls_f64 = lin_subdiv as f64;
    match dim {
        2 => {
            let lat_c = p.lat_coord(unit_cell);
            let mut lat_idx = lat_c.iter().map(|x| (x*ls_f64).floor() as usize);
            let idx1 = lat_idx.next().unwrap();
            let idx2 = lat_idx.next().unwrap();
            idx1 + lin_subdiv*idx2
        }
        3 => {
            let lat_c = p.lat_coord(unit_cell);
            let mut lat_idx = lat_c.iter().map(|x| (x*ls_f64).floor() as usize);
            let idx1 = lat_idx.next().unwrap();
            let idx2 = lat_idx.next().unwrap();
            let idx3 = lat_idx.next().unwrap();
            idx1 + lin_subdiv*idx2 + lin_subdiv*lin_subdiv*idx3
        }
        _ => unimplemented!(),
    }
}

impl<P: Particle + Debug + Display + Send + Sync + Clone> CellList<P> {
    // So we don't have to duplicate code for initialization
    pub fn from_overbox_list(o_list: OverboxList<P>) -> CellList<P> {
        let unit_cell = o_list.cell;
        let dim = o_list.dim;
        let n_particles = o_list.p_vec.len();
        let p_list = o_list.p_vec;

        let mut list = CellList { dim, unit_cell, p_list, n_particles, lin_subdiv: 0, c_assoc_list: Vec::new(), p_assoc_list: Vec::new() };
        list.build_cell_list();
        list
    }

    pub fn build_cell_list(&mut self) {
        // Note: this choice assumes monodispersity!
        let max_radius = self.first_particle().hint_upper();
        
        let max_subdiv = max_subdiv(self.dim, &self.unit_cell, max_radius);
        self.lin_subdiv = if max_subdiv > OPT.subdiv_offset {
            max_subdiv - OPT.subdiv_offset
        } else {
            panic!("subdiv_offset exceeds max_subdiv");
        };

        // TODO: check truncating behavior
        self.c_assoc_list = vec![Vec::new(); self.lin_subdiv.pow(self.dim as u32)];
        self.p_assoc_list = vec![0; self.n_particles];
        for (i, p) in self.p_list.iter().enumerate() {
            let idx = assign_cell(self.dim, self.lin_subdiv, &self.unit_cell, p);
            self.c_assoc_list[idx].push(i);
            self.p_assoc_list[i] = idx;
        }

    }

    fn rebin_particle(&mut self, p_idx: usize) {
        let new_cell = assign_cell(self.dim, self.lin_subdiv, &self.unit_cell, &self.p_list[p_idx]);
        let old_cell = self.p_assoc_list[p_idx];
        if old_cell == new_cell {
            return;
        } else {
            let idx_in_old_cell = self.c_assoc_list[old_cell].iter()
                .position(|&p| p == p_idx)
                .unwrap();
            self.c_assoc_list[old_cell].remove(idx_in_old_cell);
            self.c_assoc_list[new_cell].push(p_idx);
            self.p_assoc_list[p_idx] = new_cell;
        }
    }
    
    // Dumps info for debugging purposes
    #[allow(dead_code)]
    pub fn debug(&self) {
        println!("dim: {}", self.dim);
        println!("unit_cell: {:?}", &self.unit_cell);
        println!("n_particles: {}", self.n_particles);
        println!("lin_subdiv: {}", self.lin_subdiv);
        for (i, cell) in self.c_assoc_list.iter().enumerate() {
            println!("cell {}:", i);
            for p in cell {
                println!("{:?}", &self.p_list[*p]);
            }
        }
    }
}

// Returns a list of the the cells surrounding the given cell
// Includes the given cell
// i.e.
// 1  2  3  4
// 5  6  7  8
// 9  10 11 12
// 13 14 15 16
// Querying 6 would give [1, 2, 3, 5, 6, 7, 9, 10, 11]
// in some order
// Also gives the pbc offset for cells on border, i.e. if given is 12,
// 9 would come wth a +u0 offset
fn cells_around(cell: usize, dim: usize, lin_subdiv: usize, unit_cell: &[f64])
    -> Vec<(usize, Vec<f64>)> 
{ 
    // Be careful, this cast is only correct for
    // small enough integers!
    let c = cell as isize;
    let l_s = lin_subdiv as isize;

    let mut surrounding_list = Vec::new();
    
    match dim {
        2 => {
            let u = Matrix2::from_column_slice(unit_cell);
            let u0 = u.column(0);
            let u1 = u.column(1);
            
            let given_x_idx = c % l_s;
            let given_y_idx = c / l_s;
            
            for i in -1..2 {
                let out_x_idx = given_x_idx + i;
                for j in -1..2 {
                    let out_y_idx = given_y_idx + j;
                    
                    let mut offset = Vector2::new(0.0, 0.0);
                    if out_x_idx < 0 {
                        offset -= u0;
                    } else if out_x_idx >= l_s {
                        offset += u0;
                    }
                    if out_y_idx < 0 {
                        offset -= u1;
                    } else if out_y_idx >= l_s {
                        offset += u1;
                    }

                    let pbc_x_idx = ((out_x_idx%l_s + l_s)%l_s) as usize;
                    let pbc_y_idx = ((out_y_idx%l_s + l_s)%l_s) as usize;

                    let out_idx = pbc_x_idx + lin_subdiv*pbc_y_idx;

                    surrounding_list.push((out_idx, offset.as_slice().try_into().unwrap()));

                }
            }
        }
        
        3 => {
            let u = Matrix3::from_column_slice(unit_cell);
            let u0 = u.column(0);
            let u1 = u.column(1);
            let u2 = u.column(2);
            
            let given_x_idx = c % l_s;
            let given_y_idx = (c / l_s) % l_s;
            let given_z_idx = c / (l_s * l_s);
            
            for i in -1..2 {
                let out_x_idx = given_x_idx + i;
                for j in -1..2 {
                    let out_y_idx = given_y_idx + j;
                    for k in -1..2 {
                        let out_z_idx = given_z_idx + k;
                        
                        let mut offset = Vector3::new(0.0, 0.0, 0.0);
                        if out_x_idx < 0 {
                            offset -= u0;
                        } else if out_x_idx >= l_s {
                            offset += u0;
                        }
                        if out_y_idx < 0 {
                            offset -= u1;
                        } else if out_y_idx >= l_s {
                            offset += u1;
                        }
                        if out_z_idx < 0 {
                            offset -= u2;
                        } else if out_z_idx >= l_s {
                            offset += u2;
                        }

                        let pbc_x_idx = ((out_x_idx%l_s + l_s)%l_s) as usize;
                        let pbc_y_idx = ((out_y_idx%l_s + l_s)%l_s) as usize;
                        let pbc_z_idx = ((out_z_idx%l_s + l_s)%l_s) as usize;

                        let out_idx = pbc_x_idx + lin_subdiv*pbc_y_idx + lin_subdiv*lin_subdiv*pbc_z_idx;

                        surrounding_list.push((out_idx, offset.as_slice().try_into().unwrap()));
                    }
                }
            }
        }
        
        _ => unimplemented!(),
    
    };
    
    surrounding_list
}

impl<P: Particle + Debug + Display + Send + Sync + Clone> Asc<P> for CellList<P> {

    // Print configuration onto stdout
    fn print_asc(&self) {
        println!("{} 1 {}", self.dim, P::TYPE);
        for entry in self.unit_cell.iter().with_position() {
            match entry {
                Position::Last(x) => println!("{}", x),
                Position::Middle(x) => print!("{} ", x),
                Position::First(x) => print!("{} ", x),
                Position::Only(x) => println!("{}", x),
            }
        }
        for p in &self.p_list {
            println!("{}", p);
        }
    }

    // Save configuration to file
    fn save_asc(&self, path: &Path, annotation: Option<&str>) {
        // See implementation in overbox_list.rs for references
        let mut file = BufWriter::new(File::create(path).expect("Must specify valid path to save to."));
        if let Some(annotation) = annotation {
            writeln!(&mut file, "{}", annotation).expect("Failed write during save.");
        }
        writeln!(&mut file, "{} 1 {}", self.dim, P::TYPE).expect("Failed write during save.");
        for entry in self.unit_cell.iter().with_position() {
            match entry {
                Position::Last(x) => writeln!(&mut file, "{}", x),
                Position::Middle(x) => write!(&mut file, "{} ", x),
                Position::First(x) => write!(&mut file, "{} ", x),
                Position::Only(x) => writeln!(&mut file, "{}", x),
            }.expect("Failed write during save.");
        }
        for p in &self.p_list {
            writeln!(&mut file, "{}", p).expect("Failed write during save.");
        }
        let mut f = file.into_inner().expect("Failed to unwrap buffer during save");
        f.flush().expect("Failed to flush file writer during save");
        f.sync_all().expect("Failed to sync during save.");
        let canonical = path.canonicalize().expect("Must be able to canonicalize");
        let dir = OpenOptions::new()
            .read(true)
            .open(canonical.parent().expect("Canonical form should have parent"))
            .expect("Failed to open parent directory");
        dir.sync_all().expect("Failed to sync directory during save.");
    }

    // Check particle over overlaps in Asc, return number of overlaps
    fn check_particle(&self, fixed: &P) -> usize {
        // Can't just look it up because particle may have been perturbed
        let current_cell = assign_cell(self.dim, self.lin_subdiv, &self.unit_cell, fixed);
        let mut n_overlaps = 0usize;

        for (cell, offset) in cells_around(current_cell,
                                           self.dim,
                                           self.lin_subdiv,
                                           &self.unit_cell) 
        {
            for p_idx in &self.c_assoc_list[cell] {
                n_overlaps += fixed.check_overlap(&self.p_list[*p_idx], &offset) as usize;
            }
        }
        
        n_overlaps
    }

    // Return cell volume
    fn cell_volume(&self) -> f64 { volume(self.dim, &self.unit_cell) }

    // True if no overlaps
    fn is_valid(&self) -> bool {
        self.p_list.par_iter()
            .map(|p| self.check_particle(p))
            .all(|x| x <= 1)
    }

    // Try to change the cell by straining
    fn try_cell_move(&mut self, schedule: &Schedule<P>, rng: &mut Xoshiro256StarStar) -> bool {
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
                let strain = Matrix2::from_row_slice(
                    &[iso + axi, shear, shear, iso - axi]
                );
                // Change unit cell
                // Implicit transposition in read order
                // Noticed that this is necessary due to older
                // group code, courtesy of Duyu, Steve, and Yang
                let current_cell = Matrix2::from_column_slice(&self.unit_cell);
                let new_cell = current_cell + strain*current_cell;
                self.unit_cell = new_cell.as_slice().to_vec();
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
                // Change unit cell
                // Implicit transposition in read order
                let current_cell = Matrix3::from_column_slice(&self.unit_cell);
                let new_cell = current_cell + strain*current_cell;
                self.unit_cell = new_cell.as_slice().to_vec();
            }
            _ => unimplemented!(),
        }
        // Kinda weird, need to do ref outside of closure
        // https://stackoverflow.com/questions/48717833/how-to-use-struct-self-in-member-method-closure
        let new_cell = &self.unit_cell;
        // Change particles
        self.p_list.par_iter_mut().for_each(|p| {
            p.apply_strain(&old_asc.unit_cell, new_cell);
        });

        // check to see if we need a cell list rebuild
        let max_radius = self.first_particle().hint_upper();
        if self.lin_subdiv 
            > max_subdiv(self.dim, &self.unit_cell, max_radius)
        {
            self.build_cell_list();
        }


        if self.is_valid() { // Accept probabilistically
            // http://www.pages.drexel.edu/~cfa22/msim/node31.html
            let new_vol = self.cell_volume();
            let old_vol = old_asc.cell_volume();
            let n_particles = self.n_particles as f64;
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

    // Try to move a particle
    fn try_particle_move(&mut self, schedule: &mut Schedule<P>, rng: &mut Xoshiro256StarStar) -> bool 
    {
        let r_idx: usize = Uniform::new(0, self.n_particles)
            .sample(rng);
        let (old_p, move_type) = self.p_list[r_idx]
            .perturb(&self.unit_cell, &schedule.particle_param, rng);
        self.rebin_particle(r_idx);
        schedule.particle_tries[move_type] += 1;
        if self.check_particle(&self.p_list[r_idx]) > 1 { //reject
            self.p_list[r_idx] = old_p; //roll back
            self.rebin_particle(r_idx);
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
            schedule.particle_accepts[move_type] += 1;
            true
        }
    }

    // Return number of particles in Asc
    fn n_particles(&self) -> usize { self.n_particles }

    // Return a reference to the first particle stored in Asc
    fn first_particle(&self) -> &P
    {
        &self.p_list[0]
    }

}


