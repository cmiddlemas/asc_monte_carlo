// overbox_list.rs
// Author: Timothy Middlemas
// Implements ASC algorithm described in
// Jiao et. al. Nature 2009
// Generic over particle shape and simulation dimension
// Uses a naive "overbox" implementation to achieve pbc
use std::fmt::{Debug, Display};
use rand_xoshiro::Xoshiro256StarStar;
use itertools::{Itertools, Position};
use std::path::Path;
use std::fs::{File, OpenOptions};
use std::io::{Write, BufWriter, BufReader, BufRead};
use std::path::PathBuf;
use rayon::prelude::*;
use rand_distr::{Uniform, Distribution};
use crate::schedule::Schedule;
use crate::OPT;
use crate::asc::Asc;
use crate::common_util::{volume, gen_random_strain, min_float_slice};
use crate::particle::Particle;
use crate::cell_list::CellList;

// Free helper functions

fn calc_offset(n: usize, dim: usize, overbox: usize, cell: &[f64]) -> Vec<f64>
{
    let mut offset = Vec::with_capacity(dim);
    let copies = 2*(overbox as i32) + 1;
    for i in 0..dim {
        let mut coord = 0.0;
        let mut k = n as i32;
        for j in 0..dim {
            let scalar = (k % copies - overbox as i32) as f64;
            k /= copies;
            coord += scalar*cell[dim*j + i];
        }
        offset.push(coord);
    }
    offset
}

// https://stackoverflow.com/questions/34572784/why-can-i-iterate-over-a-slice-twice-but-not-a-vector
// https://doc.rust-lang.org/std/slice/
fn is_zero(v: &[f64]) -> bool {
    for val in v {
        if *val != 0.0 {
            return false
        }
    }
    true
}

// cell stored as (dim*row + column)
// and each row is interpreted as the
// coordinates vx vy vz ... of a lattice
// vector
#[derive(Clone)]
pub struct OverboxList<P> {
    pub dim: usize, // Dimension of configuration
    pub overbox: usize, // # of overboxes, needed for small unit cell
    pub cell: Vec<f64>, // Unit Cell
    pub p_vec: Vec<P>, // List of particles
}

// https://old.reddit.com/r/rust/comments/98jldy/question_about_rayon_and_intopariterator_not/
// https://github.com/rust-lang/rust/issues/61768
impl<P: Particle + Debug + Display + Send + Sync + Clone> OverboxList<P> {
    
    // Makes a trivial, generic config for testing
    #[allow(dead_code)]
    pub fn make() -> OverboxList<P> { 
        OverboxList { dim: 0, overbox: 0, cell: Vec::new(), p_vec: Vec::new()}
    }

    pub fn from_file(path: &PathBuf) -> OverboxList<P> {
        // Made extensive use of the rust docs for various
        // std components when writing this function
        let mut infile = BufReader::new(File::open(path).expect("Input file must be valid"));
        let mut buf = String::new();
        // 1st line, dim overbox type
        infile.read_line(&mut buf).expect("Valid utf-8");
        let mut line_one = buf.split_whitespace();
        let dim: usize = line_one.next().unwrap()
            .parse().expect("Valid dimension");
        let overbox: usize = line_one.next().unwrap()
            .parse().expect("Valid overbox");
        // 2nd line, cell
        buf.clear();
        infile.read_line(&mut buf).expect("Valid utf-8");
        let cell: Vec<f64> = buf.split_whitespace()
            .map(|x| x.parse().expect("valid unit cell"))
            .collect();
        // 3rd line and on, particle
        let p_vec: Vec<P> = infile.lines()
            .map(|x| P::parse(&x.expect("Valid utf-8"), &cell))
            .collect();
        let o_list = OverboxList { dim, overbox, cell, p_vec };
        o_list
    }

    pub fn from_cell_list(clist: CellList<P>) -> OverboxList<P> {
        OverboxList { dim: clist.dim, overbox: 1, cell: clist.unit_cell, p_vec: clist.p_list }
    }

    // Make an rsa config
    pub fn make_rsa(n: usize, // Number of particles to insert
                    shape: &P, // Contains the shape information
                    dim: usize, // dimension of config
                    init_cell: Vec<f64>, // Initial cell shape
                    rng: &mut Xoshiro256StarStar // Reproducible Rng
    ) -> OverboxList<P> {
        let mut new_asc = OverboxList 
            { dim: dim, overbox: 1, cell: init_cell.clone(), p_vec: Vec::new() };
        while new_asc.p_vec.len() < n {
            if new_asc.try_add_particle(
                shape.copy_shape_random_coord(
                    &init_cell, rng
                )
            ) {} 
            /* Some code for printing out RSA data, unlikely to want to use this
             * in a real MC run
            {
                println!("Added particle {}", new_asc.p_vec.len());
                if let Some(path) = &OPT.savefiles {
                    if path.is_dir() {
                        println!("Root filename cannot be empty/ you specified a dir. Skipping save.");
                        continue;
                    }
                    let mut full_path = path.clone();
                    // https://users.rust-lang.org/t/what-is-right-ways-to-concat-strings/3780/4
                    full_path.set_file_name(
                        format!("{}_rsa_{}.dat", 
                                path.file_name().expect("Must give a valid root filename.")
                                    .to_str().expect("Must give valid UTF-8 str."), 
                                new_asc.p_vec.len()
                    ));
                    new_asc.save_asc(&full_path);
                }
            }
            */
        }
        new_asc
    }
    
    // Dumps all the info contained in OverboxList
    #[allow(dead_code)]
    pub fn debug(&self) {
        println!("dim: {}", self.dim);
        println!("{:?}", &self.cell);
        for p in &self.p_vec{
            println!("{:?}", p);
        }
    }

    // returns true if add is successful
    pub fn try_add_particle(&mut self, p: P) -> bool {
        if self.check_particle(&p) > 0 {
            return false;
        }
        self.p_vec.push(p);
        true
    }

}

impl<P: Particle + Debug + Display + Clone + Send + Sync> Asc<P> for OverboxList<P> {
    // Prints the config in a nice, ascii format
    fn print_asc(&self) {
        println!("{} {} {}", self.dim, self.overbox, P::TYPE);
        for entry in self.cell.iter().with_position() {
            match entry {
                Position::Last(x) => println!("{}", x),
                Position::Middle(x) => print!("{} ", x),
                Position::First(x) => print!("{} ", x),
                Position::Only(x) => println!("{}", x),
            }
        }
        for p in &self.p_vec {
            println!("{}", p);
        }
    }

    // Saves the config in a nice, ascii format given path
    // Panics on any error
    // Tries to wait until data hits disk
    // https://doc.rust-lang.org/std/io/struct.BufWriter.html
    fn save_asc(&self, path: &Path, annotation: Option<&str>) {
        let mut file = BufWriter::new(File::create(path).expect("Must specify valid path to save to."));
        if let Some(annotation) = annotation {
            writeln!(&mut file, "{}", annotation).expect("Failed write during save.");
        }
        writeln!(&mut file, "{} {} {}", self.dim, self.overbox, P::TYPE).expect("Failed write during save.");
        for entry in self.cell.iter().with_position() {
            match entry {
                Position::Last(x) => writeln!(&mut file, "{}", x),
                Position::Middle(x) => write!(&mut file, "{} ", x),
                Position::First(x) => write!(&mut file, "{} ", x),
                Position::Only(x) => writeln!(&mut file, "{}", x),
            }.expect("Failed write during save.");
        }
        for p in &self.p_vec {
            writeln!(&mut file, "{}", p).expect("Failed write during save.");
        }
        let mut f = file.into_inner().expect("Failed to unwrap buffer during save");
        f.flush().expect("Failed to flush file writer during save");
        f.sync_all().expect("Failed to sync during save.");
        // https://doc.rust-lang.org/std/path/struct.Path.html#method.canonicalize
        let canonical = path.canonicalize().expect("Must be able to canonicalize");
        // some info about fsyncing directory
        // http://blog.httrack.com/blog/2013/11/15/everything-you-always-wanted-to-know-about-fsync/
        // https://old.reddit.com/r/node/comments/4r8k11/how_to_call_fsync_on_a_directory/
        // https://github.com/google/renameio/issues/11
        // https://stackoverflow.com/questions/51100698/how-to-fsync-a-directory-under-linux-in-c
        // https://stackoverflow.com/questions/20687611/is-there-fsync-but-with-path-parameter
        // https://stackoverflow.com/questions/49060587/is-there-a-way-to-call-fsync-flush-the-parent-of-a-fd
        // http://manpages.ubuntu.com/manpages/bionic/man2/fsync.2.html
        let dir = OpenOptions::new()
            .read(true)
            .open(canonical.parent().expect("Canonical form should have parent"))
            .expect("Failed to open parent directory");
        dir.sync_all().expect("Failed to sync directory during save.");
    }

    // returns number of overlaps
    fn check_particle(&self, fixed: &P) -> usize {
        let n_offset = (2*self.overbox + 1).pow(self.dim as u32);
        if OPT.parallelize_inner {
            (0..n_offset).into_par_iter()
                .map(|x| calc_offset(x, self.dim, self.overbox, &self.cell))
                .map(|offset|
                    self.p_vec.iter()
                        .map(|imaged| fixed.check_overlap(imaged, &offset) as usize)
                        .sum::<usize>()
                    )
                .sum()
            /* Chunking strategy, probably slower, probably heard about trying this
             * from somewhere on the internet
            // https://stackoverflow.com/questions/37033700/how-do-i-process-a-range-in-slices-in-rust
            (0..n_offset).collect::<Vec<usize>>()
                .par_chunks(CHUNK_SIZE)
                .map(|x| x.iter().map(|&e| calc_offset(e, self.dim, self.overbox, &self.cell)))
                .map(|offset_chunk|
                    offset_chunk.map(|offset| 
                        self.p_vec.iter()
                            .map(|imaged| fixed.check_overlap(imaged, &offset) as usize)
                            .sum::<usize>()
                    ))
                .map(|result_chunk| result_chunk.sum::<usize>())
                .sum()
            */
        } else {
            (0..n_offset)
                .map(|x| calc_offset(x, self.dim, self.overbox, &self.cell))
                .map(|offset|
                    self.p_vec.iter()
                        .map(|imaged| fixed.check_overlap(imaged, &offset) as usize)
                        .sum::<usize>()
                    )
                .sum()
        }
            // For future reference, to parallelize range:
            // https://users.rust-lang.org/t/rayon-parallel-sum-from-range/6367

            // Obsolete comment, was for when single particle check was parallelized
            // For some reason faster than into_par_iter(), also see
            // https://users.rust-lang.org/t/for-loops-in-rust/8217/4

    /* Old code, reversed inner and outer loops
        self.p_vec.par_iter().map(|imaged|
            // https://users.rust-lang.org/t/auto-vectorization-in-rust/24379/2
            offset_vec.iter().map(|offset|
                // https://stackoverflow.com/questions/55461617/how-do-i-convert-a-boolean-to-an-integer-in-rust
                fixed.check_overlap(imaged, offset) as usize
            // https://stackoverflow.com/questions/51283403/cannot-infer-type-for-b-for-filter-map-sum
            ).sum::<usize>()
        ).sum()
    */
    }

    fn particle_gap(&self, exclude_idx: usize, particle: &P) -> f64 {
        let n_offset = (2*self.overbox + 1).pow(self.dim as u32);
        let phi = self.packing_fraction();
        let minimal_list: Vec<f64> = (0..n_offset)
                .map(|x| calc_offset(x, self.dim, self.overbox, &self.cell))
                .map(|offset|
                    self.p_vec.iter().enumerate()
                        .map(|(i, imaged)|
                            if i == exclude_idx && is_zero(&offset) {
                                f64::INFINITY
                            } else {
                                (phi/particle.overlap_scale(imaged, &offset)) - phi
                            }
                        )
                        .collect::<Vec<f64>>()
                )
                .flatten()
                .collect();

        min_float_slice(&minimal_list)
    }

    fn cell_volume(&self) -> f64 {
        volume(self.dim, &self.cell)
    }

    fn is_valid(&self) -> bool {
        if OPT.parallelize_inner {
            self.p_vec.iter()
                .map(|p| self.check_particle(p))
                .all(|x| x <= 1)
        } else {
            self.p_vec.par_iter()
                .map(|p| self.check_particle(p))
                .all(|x| x <= 1)
        }
    }

    fn apply_random_strain(mut self, schedule: &Schedule<P>, rng: &mut Xoshiro256StarStar) -> Option<(Self, f64)> {
    // Choose random strain
        let (trace_strain, new_cell) = gen_random_strain(self.dim, &self.cell, schedule, rng);
        // Apply strain to cell
        self.cell = new_cell;
                        
        // Kinda weird, need to do ref outside of closure
        // https://stackoverflow.com/questions/48717833/how-to-use-struct-self-in-member-method-closure
        let new_cell = &self.cell;
        // Change particles
        if OPT.no_rayon {
            self.p_vec.iter_mut().for_each(|p| {
                p.apply_strain(new_cell);
            });
        } else {
            self.p_vec.par_iter_mut().for_each(|p| {
                p.apply_strain(new_cell);
            });
        }
     
        Some((self, trace_strain))
    }   

    fn try_particle_move(&mut self,
                             schedule: &mut Schedule<P>,
                             rng: &mut Xoshiro256StarStar
    ) -> bool
    {
        let r_idx: usize = Uniform::new(0, self.p_vec.len())
            .sample(rng);
        let (old_p, move_type) = self.p_vec[r_idx]
            .perturb(&self.cell, &schedule.particle_param, rng);
        schedule.particle_tries[move_type] += 1;
        if self.check_particle(&self.p_vec[r_idx]) > 1 { //reject
            self.p_vec[r_idx] = old_p; //roll back
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
    
    fn n_particles(&self) -> usize { 
        self.p_vec.len()
    }
    
    fn unit_cell(&self) -> &[f64] { &self.cell }

    fn first_particle(&self) -> &P {
        &self.p_vec[0]
    }

    fn particle_slice(&self) -> &[P] {
        &self.p_vec
    }
}
