use std::io::prelude::*;
use std::fmt::{Debug, Display, Formatter};
use std::fmt;
use rand_xoshiro::Xoshiro256StarStar;
use rand_distr::{Uniform, Distribution};
use rand::{SeedableRng, RngCore};
use rand::rngs::OsRng;
use std::path::PathBuf;
use std::fs::File;
use structopt::StructOpt;
use lazy_static::lazy_static;
use chrono::prelude::*;

mod asc;
mod disks;
mod schedule;
mod spheres;
mod ellipsoids;

use asc::{Asc, Particle, save_asc_from_opt};
use disks::{Disk};
use spheres::{Sphere};
use ellipsoids::{Ellipsoid};
use schedule::Schedule;

// Computed with Mathematica Student 11.2.0.0
const PI: f64 = 3.141592653589793;

// Struct for command line
#[derive(StructOpt, Debug)]
#[structopt(name = "asc_monte_carlo")]
struct Opt {
    /// Dimension
    #[structopt(short = "d", long, default_value = "3")]
    dimension: usize,

    /// Number of particles
    #[structopt(short = "n", long, default_value = "10")]
    number: usize,

    /// Work with ellipsoids?
    #[structopt(long)]
    ellipsoid: bool,

    /// x-aligned semi-axis
    #[structopt(short = "a", default_value = "1.0")]
    a_semi: f64,

    /// y-aligned semi-axis
    #[structopt(short = "b", default_value = "1.0")]
    b_semi: f64,
    
    /// x-aligned semi-axis
    #[structopt(short = "c", default_value = "1.0")]
    c_semi: f64,

    /// Side length of initial cubic cell
    #[structopt(short = "s", long, default_value = "10.0")]
    side: f64,

    /// Root filename to save configs as
    /// Default to only give initial and final
    /// config on stdout, will give all sweeps
    /// if specified
    #[structopt(short = "o", long)]
    savefiles: Option<PathBuf>,

    /// Root filename to save pretty-printed logs
    #[structopt(long)]
    logfiles: Option<PathBuf>,

    /// Number of sweeps to do
    #[structopt(long, default_value = "100")]
    sweeps: usize,

    /// Number of moves to do per sweep
    /// Recommended at least 10x particle number
    #[structopt(long, default_value = "1000")]
    moves: usize,

    /// Probability of doing a cell move
    #[structopt(long, default_value = "0.001")]
    pcell: f64,

    /// Std dev of gaussian distribution
    /// for translational moves, as fraction
    /// of particle radius/minor semiaxis
    #[structopt(long, default_value = "0.3")]
    dtrans: f64,

    /// Std dev of gaussian distribution for
    /// rotational moves
    #[structopt(long, default_value = "0.01")]
    drot: f64,

    /// Radius of spherical particles
    #[structopt(long, default_value = "1.0")]
    radius: f64,

    /// Initial pressure of system
    #[structopt(long, default_value = "1.0")]
    ipressure: f64,

    /// Initial shear step width
    #[structopt(long, default_value = "0.01")]
    shear: f64,

    /// Initial axial compression step width
    #[structopt(long, default_value = "0.01")]
    axial: f64,

    /// Initial isotropic volume change step width
    #[structopt(long, default_value = "0.01")]
    isotropic: f64,

    /// Should the program adjust the step sizes
    /// to try and keep acceptance ratio between
    /// adjust-upper and adjust-lower?
    #[structopt(long)]
    adjust: bool,

    #[structopt(long, default_value = "0.4")]
    adjust_lower: f64,

    #[structopt(long, default_value = "0.6")]
    adjust_upper: f64,

    /// Optional file holding seed in 32 u8 format
    #[structopt(long)]
    seedfile: Option<PathBuf>,

    /// Try and parallelize inner loops? Not always worth it.
    #[structopt(long)]
    parallelize_inner: bool,
}

lazy_static! {
    static ref OPT: Opt = Opt::from_args();
}

fn make_and_run_schedule<P: Particle + Clone + Debug + Display + Send + Sync>
        (mut config: Asc<P>, rng: &mut Xoshiro256StarStar) {
    println!("Initial Configuration:");
    config.print_asc();
    save_asc_from_opt(&config, "initial");
    let mut schedule = Schedule::make();
    println!("Running schedule:");
    println!("{:?}", schedule);
    schedule.run(&mut config, rng);
    println!("Ended schedule:");
    println!("{:?}", schedule);
    println!("Final Configuration:");
    config.print_asc();
    save_asc_from_opt(&config, "final");
    println!("Ending program at: {}", Utc::now());
}

fn consume<T>(arg: T) {}

fn main() {
// Immediately read command line
    consume(&*OPT);

// Report runtime parameters ---------------------------------------
    println!("Starting program at: {}", Utc::now());
    println!("Given command line options:");
    println!("{:?}", *OPT);
    
// Initialize rng --------------------------------------------------
    let mut rng = if let Some(path) = &OPT.seedfile {
        // https://stackoverflow.com/questions/31192956/whats-the-de-facto-way-of-reading-and-writing-files-in-rust-1-x
        // Found above after the fact, but previously accessed (probably), also used docs
        // https://doc.rust-lang.org/std/fs/struct.File.html
        let mut file = File::open(path).expect("Seedfile must be valid");
        let mut seed_string = String::new();
        file.read_to_string(&mut seed_string).expect("Need a successful read on seedfile");
        let mut seed = [0u8; 32];
        for (i, token) in seed_string.trim().split(' ').enumerate() {
            // https://rust-lang-nursery.github.io/rust-cookbook/text/string_parsing.html
            // https://doc.rust-lang.org/std/str/trait.FromStr.html
            seed[i] = token.parse().expect("Gave invalid u8 in seedfile");
        }
        print!("Using seed: ");
        for i in 0..32 {
            if i < 31 {
                print!("{} ", seed[i]);
            } else {
                println!("{}", seed[i]);
            }
        }
        Xoshiro256StarStar::from_seed(seed)
    } else {
        // https://docs.rs/rand/0.7.2/rand/rngs/struct.OsRng.html
        let mut seed = [0u8; 32];
        OsRng.fill_bytes(&mut seed);
        print!("Using seed: ");
        for i in 0..32 {
            if i < 31 {
                print!("{} ", seed[i]);
            } else {
                println!("{}", seed[i]);
            }
        }
        Xoshiro256StarStar::from_seed(seed)
    };
    
// Make and run the correct simulation -----------------------------------
    match OPT.dimension {
        2 => {
            if OPT.ellipsoid {
                panic!("2d ellipse packing not implemented yet!");
            } else {
                let shape = Disk::make_shape(OPT.radius);
                let mut init_cell = vec![OPT.side, 0.0, 0.0, OPT.side];
                let config = Asc::make_rsa(OPT.number, &shape, 2, init_cell, &mut rng);
                make_and_run_schedule(config, &mut rng);
            }
        }
        3 => {
            if OPT.ellipsoid {
                let shape = Ellipsoid::make_shape(OPT.a_semi, OPT.b_semi, OPT.c_semi);
                let mut init_cell = vec![OPT.side, 0.0, 0.0,
                                    0.0, OPT.side, 0.0,
                                    0.0, 0.0, OPT.side];
                let config = Asc::make_rsa(OPT.number, &shape, 3, init_cell, &mut rng);
                make_and_run_schedule(config, &mut rng);
            } else {
                let shape = Sphere::make_shape(OPT.radius);
                let mut init_cell = vec![OPT.side, 0.0, 0.0,
                                    0.0, OPT.side, 0.0,
                                    0.0, 0.0, OPT.side];
                let config = Asc::make_rsa(OPT.number, &shape, 3, init_cell, &mut rng);
                make_and_run_schedule(config, &mut rng);
            }
        }
        _ => panic!("Can't handle that dimension yet!"),
    }
}

