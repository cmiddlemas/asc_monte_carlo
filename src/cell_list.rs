// cell_list.rs
// Author: Timothy Middlemas
// Implements the Asc trait using an underlying
// cell list representation. Performs automatic resizing
// Assumes monodispersity for now
// Is probably similar to Ge's cell list implementation,
// which I have access to, but this is in native rust
// Also looked at the M. Skoge LS code and
// http://www.acclab.helsinki.fi/~knordlun/moldyn/lecture03.pdf
// and
// http://cacs.usc.edu/education/cs596/01-1LinkedListCell.pdf
use crate::asc::Asc;
use crate::common_util::volume;
use crate::particle::Particle;
use crate::overbox_list::OverboxList;
use std::fmt::{Debug, Display};
use nalgebra::{Matrix2, Matrix3, Vector2, Vector3};
use crate::OPT;
use crate::schedule::Schedule;
use rand_xoshiro::Xoshiro256StarStar;
use rand_distr::{Uniform, Distribution};
use std::path::Path;
use itertools::{Itertools, Position};
use std::fs::{File, OpenOptions};
use std::io::{Write, BufWriter};
use rayon::prelude::*;
use crate::common_util::{min_width, gen_random_strain};

#[derive(Clone, Debug)]
pub struct CellList<P> {
    pub dim: usize,
    pub unit_cell: Vec<f64>, // Unit cell, stored as (dim*row + column), each row is a lattice vector
    pub p_list: Vec<P>,
    pub n_particles: usize,
    pub lin_subdiv: usize,
    pub c_assoc_list: Vec<Vec<usize>>, // Association lists given as [cell1(p1...), cell2(p1...)...]
    pub p_assoc_list: Vec<usize>, // Gives first (cell identity) index into c_assoc_list, same order as p_list
}

impl<P> CellList<P> {
}

// Returns the maximum valid value of lin_subdiv given a dimension,
// the current unit cell, and an upper bound on the radius of the 
// circumscribing sphere of the largest particle in the packing
// References:
// https://scicomp.stackexchange.com/questions/3107/minimum-image-convention-for-triclinic-unit-cell
// http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.57.1696&rep=rep1&type=pdf
// https://hoomd-blue.readthedocs.io/en/stable/nlist.html
fn max_subdiv(dim: usize, unit_cell: &[f64], max_radius: f64) -> usize {
    let max_diam = 2.0*max_radius;
    let min_width = min_width(dim, unit_cell);    
    // https://stackoverflow.com/questions/37506672/convert-float-to-integer-in-rust
    // Safely convert between float and usize
    let answer_f64 = (min_width/max_diam).floor();
    assert!(answer_f64.is_finite());
    assert!(answer_f64 >= 0.0 && answer_f64 <= 100000000.0);
    answer_f64 as usize
}

// Cell layout:
// 7 8 9
// 4 5 6
// 1 2 3
fn assign_cell<P: Particle>
    (dim: usize,
    lin_subdiv: usize,
    p: &P) 
-> usize 
{
    let ls_f64 = lin_subdiv as f64;
    match dim {
        2 => {
            let lat_c = p.lat_coord();
            // Safely cast between float and usize
            let mut lat_idx = lat_c.iter().map(|x| 
                { let ans = (x*ls_f64).floor();
                assert!(ans.is_finite());
                assert!(ans >= 0.0 && ans <= 100000000.0);
                ans as usize });
            let idx1 = lat_idx.next().unwrap();
            let idx2 = lat_idx.next().unwrap();
            idx1 + lin_subdiv*idx2
        }
        3 => {
            let lat_c = p.lat_coord();
            let mut lat_idx = lat_c.iter().map(|x| 
                { let ans = (x*ls_f64).floor();
                assert!(ans.is_finite());
                assert!(ans >= 0.0 && ans <= 100000000.0);
                ans as usize });
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
            eprintln!("max_subdiv = {}", max_subdiv);
            eprintln!("cell = {:?}", &self.unit_cell);
            eprintln!("vol = {:?}", self.cell_volume());
            panic!("subdiv_offset exceeds max_subdiv");
        };

        // TODO: check truncating behavior
        self.c_assoc_list = vec![Vec::new(); self.lin_subdiv.pow(self.dim as u32)];
        self.p_assoc_list = vec![0; self.n_particles];
        for (i, p) in self.p_list.iter().enumerate() {
            let idx = assign_cell(self.dim, self.lin_subdiv, p);
            self.c_assoc_list[idx].push(i);
            self.p_assoc_list[i] = idx;
        }

    }

    fn rebin_particle(&mut self, p_idx: usize) {
        let new_cell = assign_cell(self.dim, self.lin_subdiv, &self.p_list[p_idx]);
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

    // prevent resizing for up to 3 dimensions
    let mut surrounding_list = Vec::with_capacity(30);
    
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

                    // https://doc.rust-lang.org/std/convert/trait.From.html
                    // https://doc.rust-lang.org/std/vec/struct.Vec.html#impl-From%3C%26%27_%20%5BT%5D%3E
                    surrounding_list.push((out_idx, Vec::<f64>::from(offset.as_slice())));

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

                        surrounding_list.push((out_idx, Vec::<f64>::from((offset).as_slice())));
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
        let current_cell = assign_cell(self.dim, self.lin_subdiv, fixed);
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
    
    fn particle_gap(&self, exclude_idx: usize, particle: &P) -> f64 {
        let current_cell = assign_cell(self.dim, self.lin_subdiv, particle);
        let mut smallest_gap: f64 = f64::INFINITY;
        let phi = self.packing_fraction();

        for (cell, offset) in cells_around(current_cell,
                                           self.dim,
                                           self.lin_subdiv,
                                           &self.unit_cell) 
        {
            for p_idx in &self.c_assoc_list[cell] {
                if *p_idx != exclude_idx { // don't compute gap for particle with itself
                    let scale = particle.overlap_scale(&self.p_list[*p_idx], &offset);
                    //assert!(scale <= 1.0, "{} {} {} {:?} {:?} {:?}",
                    //        scale,
                    //        particle.check_overlap(&self.p_list[*p_idx], &offset),
                    //        self.is_valid(),
                    //        self.unit_cell(),
                    //        &particle,
                    //        &self.p_list[*p_idx]);
                    let gap = (phi/scale) - phi;
                    if gap < smallest_gap {
                        smallest_gap = gap;
                    }
                }
            }
        }

        // If the smallest gap is still inf, then
        // we need to resort to a brute force strategy
        // Shouldn't have much effect on performance, since
        // should be triggered rarely in the course of the simulation,
        // at least after we leave the almost ideal gas state
        if smallest_gap == f64::INFINITY {
            let overbox_list = OverboxList::from_cell_list(self.clone());
            smallest_gap = overbox_list.particle_gap(exclude_idx, particle);
        }
        
        smallest_gap
    }

    // Return cell volume
    fn cell_volume(&self) -> f64 { volume(self.dim, &self.unit_cell) }

    // True if no overlaps
    fn is_valid(&self) -> bool {
        if OPT.no_rayon {
            self.p_list.iter()
                .map(|p| self.check_particle(p))
                .all(|x| x <= 1)
        } else {
            self.p_list.par_iter()
                .map(|p| self.check_particle(p))
                .all(|x| x <= 1)
        }
    }
    
    fn apply_random_strain(mut self, schedule: &Schedule<P, Self>, rng: &mut Xoshiro256StarStar)
        -> Option<(Self, f64)>
    {
        // Choose random strain
        let (trace_strain, new_cell) = gen_random_strain(self.dim, &self.unit_cell, schedule, rng);
        // Apply strain to cell
        self.unit_cell = new_cell;
                        
        // Kinda weird, need to do ref outside of closure
        // https://stackoverflow.com/questions/48717833/how-to-use-struct-self-in-member-method-closure
        let new_cell = &self.unit_cell;
        // Change particles
        if OPT.no_rayon {
            self.p_list.iter_mut().for_each(|p| {
                p.apply_strain(new_cell);
            });
        } else {
            self.p_list.par_iter_mut().for_each(|p| {
                p.apply_strain(new_cell);
            });
        }
    
        // check to see if we need a cell list rebuild
        let max_radius = self.first_particle().hint_upper();
        let max_subdiv = max_subdiv(self.dim, &self.unit_cell, max_radius);
        if max_subdiv < 1 { 
            return None;
        } else if self.lin_subdiv > max_subdiv {
            self.build_cell_list();
        }

        Some((self, trace_strain))
    }

    // Try to move a particle
    fn try_particle_move(&mut self, schedule: &mut Schedule<P, Self>, rng: &mut Xoshiro256StarStar) -> bool 
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
            //assert!(self.is_valid());
            false
        } else { //accept
            //eprintln!("accept pmove");
            P::sample_obs_accepted_pmove(
                schedule,
                self,
                r_idx,
                &old_p
            );
            schedule.particle_accepts[move_type] += 1;
            //assert!(self.is_valid(), "{:?}, {:?}", &self.p_list[r_idx], &self.unit_cell);
            true
        }
    }

    // Return number of particles in Asc
    fn n_particles(&self) -> usize { self.n_particles }

    fn unit_cell(&self) -> &[f64] { &self.unit_cell }

    fn dim(&self) -> usize { self.dim }

    // Return a reference to the first particle stored in Asc
    fn first_particle(&self) -> &P
    {
        &self.p_list[0]
    }

    fn particle_slice(&self) -> &[P] {
        &self.p_list
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::overbox_list::OverboxList;
    use crate::spheres::Sphere; 

    #[test]
    fn particle_gap1() {
        // Find nearest gap in neighbor cell
        let overbox_list = OverboxList {
            dim: 3,
            overbox: 1,
            cell: vec![1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
            p_vec: vec![
                Sphere {
                    rel_pos: [0.201, 0.201, 0.201],
                    global_pos: [0.201, 0.201, 0.201],
                    radius: 0.1
                },
                Sphere {
                    rel_pos: [0.201, 0.402, 0.201],
                    global_pos: [0.201, 0.402, 0.201],
                    radius: 0.1
                }
            ]
        };
        let cell_list = CellList::from_overbox_list(overbox_list);
        assert!(cell_list.lin_subdiv == 5);
        let gap = cell_list.particle_gap(0, cell_list.first_particle());
        // Computed with Mathematica 12.0.0.0 Student Edition
        let true_gap = 0.0001262930718718569;
        eprintln!("gap: {}, true_gap: {}", gap, true_gap);
        assert!((gap - true_gap).abs()/true_gap <= 1e-10);
    }

    #[test]
    fn particle_gap2() {
        // Find nearest gap in other cell
        let overbox_list = OverboxList {
            dim: 3,
            overbox: 1,
            cell: vec![1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
            p_vec: vec![
                Sphere {
                    rel_pos: [0.201, 0.201, 0.201],
                    global_pos: [0.201, 0.201, 0.201],
                    radius: 0.1
                },
                Sphere {
                    rel_pos: [0.201, 0.601, 0.201],
                    global_pos: [0.201, 0.601, 0.201],
                    radius: 0.1
                }
            ]
        };
        let cell_list = CellList::from_overbox_list(overbox_list);
        assert!(cell_list.lin_subdiv == 5);
        let gap = cell_list.particle_gap(0, cell_list.first_particle());
        // Computed with Mathematica 12.0.0.0 Student Edition
        let true_gap = 0.05864306286700947;
        eprintln!("gap: {}, true_gap: {}", gap, true_gap);
        assert!((gap - true_gap).abs()/true_gap <= 1e-10);
    }

    #[test]
    fn particle_gap3() {
        // Find nearest gap in same cell
        let overbox_list = OverboxList {
            dim: 3,
            overbox: 1,
            cell: vec![1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
            p_vec: vec![
                Sphere {
                    rel_pos: [0.201, 0.201, 0.201],
                    global_pos: [0.201, 0.201, 0.201],
                    radius: 0.1
                },
                Sphere {
                    rel_pos: [0.399, 0.399, 0.399],
                    global_pos: [0.399, 0.399, 0.399],
                    radius: 0.1
                }
            ]
        };
        let cell_list = CellList::from_overbox_list(overbox_list);
        assert!(cell_list.lin_subdiv == 5);
        let gap = cell_list.particle_gap(0, cell_list.first_particle());
        // Computed with Mathematica 12.0.0.0 Student Edition
        let true_gap = 0.03386068461403755;
        eprintln!("gap: {}, true_gap: {}", gap, true_gap);
        assert!((gap - true_gap).abs()/true_gap <= 1e-10);
    }

    #[test]
    fn symmetric_invalidation1() {
        // Show two particles, which have the property that
        // overlap(1,2) = 1 and overlap(2,1) = 2 when computed
        // with a non-symmetric or non-decision making method
        // for computing overlaps
        use crate::common_util::{global_to_relative3, relative_to_global3};
        use std::convert::TryInto;
        let unit_cell = vec![3.982197460629651, 0.13715627693204338, 0.49349869179284556, 0.08310728005514077, 3.3947201254916064, -0.6848127113190684, 0.6571811393005731, -0.7374698335070697, 4.5158482079131845];
        let rel_pos1 = [0.11061712357162513, 0.2454929018265659, 0.190921542191982];
        let global_pos1 = relative_to_global3(&unit_cell, &rel_pos1);
        let rel_pos2 = [0.039411461262676944, 0.22175768105691646, 0.7727097479639695];
        let global_pos2 = relative_to_global3(&unit_cell, &rel_pos2);
        let overbox_list = OverboxList {
            dim: 3,
            overbox: 1,
            cell: unit_cell,
            p_vec: vec![
                Sphere {
                    rel_pos: rel_pos1,
                    global_pos: global_pos1,
                    radius: 1.0
                },
                Sphere {
                    rel_pos: rel_pos2,
                    global_pos: global_pos2,
                    radius: 1.0
                }
            ]
        };
        let cell_list = CellList::from_overbox_list(overbox_list);
        let p1 = &cell_list.particle_slice()[0];
        let p2 = &cell_list.particle_slice()[1];
        let check1 = cell_list.check_particle(p1);
        let check2 = cell_list.check_particle(p2);
        assert!(check1 == check2, "{} {} {} {}",
                check1,
                check2,
                cell_list.particle_gap(0, p1),
                cell_list.particle_gap(1, p2));
        assert!(cell_list.particle_gap(0, p1) >= 0.0, "{}", cell_list.particle_gap(0, p1));
        assert!(cell_list.particle_gap(1, p2) >= 0.0, "{}", cell_list.particle_gap(1, p2));
    }
    
    #[test]
    fn symmetric_invalidation2() {
        // In the initial fix of the above test, there was a leftover case not fixed.
        // This test was created to fix that, and the details of why this example is
        // interesting is documented in an internal test in spheres.rs
        use crate::common_util::{global_to_relative3, relative_to_global3};
        use std::convert::TryInto;
        let unit_cell = vec![5.07019429309184, -0.47515173025982416, 0.18018763060534648, -0.5656991499933727, 3.0853092570916774, -0.2935140522456714, 0.07612697602838228, -0.021482269817008064, 3.9939399781072584];
        let rel_pos1 = [0.5029715172380667, 0.23835146318049266, 0.08529553399114896];
        let global_pos1 = relative_to_global3(&unit_cell, &rel_pos1);
        let rel_pos2 = [0.8993901682271355, 0.42715227859966864, 0.9539818845758558];
        let global_pos2 = relative_to_global3(&unit_cell, &rel_pos2);
        let overbox_list = OverboxList {
            dim: 3,
            overbox: 1,
            cell: unit_cell,
            p_vec: vec![
                Sphere {
                    rel_pos: rel_pos1,
                    global_pos: global_pos1,
                    radius: 1.0
                },
                Sphere {
                    rel_pos: rel_pos2,
                    global_pos: global_pos2,
                    radius: 1.0
                }
            ]
        };
        let cell_list = CellList::from_overbox_list(overbox_list);
        let p1 = &cell_list.particle_slice()[0];
        let p2 = &cell_list.particle_slice()[1];
        let check1 = cell_list.check_particle(p1);
        let check2 = cell_list.check_particle(p2);
        assert!(check1 == 2);
        assert!(check1 == check2, "{} {} {} {}",
                check1,
                check2,
                cell_list.particle_gap(0, p1),
                cell_list.particle_gap(1, p2));
    }

}
