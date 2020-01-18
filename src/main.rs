use std::io::prelude::*;
use std::io::{BufReader, BufRead};
use std::fmt::{Debug, Display};
use rand_xoshiro::Xoshiro256StarStar;
use rand::{SeedableRng, RngCore};
use rand::rngs::OsRng;
use std::path::PathBuf;
use std::fs::File;
use structopt::StructOpt;
use lazy_static::lazy_static;
use chrono::prelude::*;

mod asc;
mod disks;
mod ellipses;
mod schedule;
mod spheres;
mod ellipsoids;
mod overbox_list;
mod cell_list;

use asc::{Asc, Particle, save_asc_from_opt};
use disks::Disk;
use ellipses::Ellipse;
use spheres::Sphere;
use ellipsoids::Ellipsoid;
use schedule::Schedule;
use overbox_list::OverboxList;
use cell_list::CellList;

// Struct for command line
#[derive(StructOpt, Debug)]
#[structopt(name = "asc_monte_carlo")]
/// A program to make packings with the adaptive shrinking
/// cell monte carlo technique. Can either take an input file
/// using --initfile or start with a low-density RSA configuration using
/// -n, -d, --ellipsoid, -s, -a, -b, -c, and --radius.
/// Can configure number of threads to run on with RAYON_NUM_THREADS
/// environment variable.
struct Opt {
    // If specified, stop after reaching this phi
    // over one sweep
    #[structopt(long)]
    max_phi: Option<f64>,

    /// If specified, use a cell list implementation
    #[structopt(long)]
    cell_list: bool,

    /// When positive, use larger cells than the smallest possible
    /// by width test in cell list
    #[structopt(long, default_value = "0")]
    subdiv_offset: usize,

    /// Check overlap, should probably specify
    /// --initfile or you'll just check an RSA
    /// configuration
    #[structopt(long)]
    check_overlap: bool,
    
    /// Near overlap tol, report overlap if closer than this
    /// Only used with --check_overlap
    /// Default value corresponds to 1% overlap
    #[structopt(long, default_value = "0.0201")]
    near_overlap_tol: f64,

    /// Dimension, only used with RSA generator
    #[structopt(short = "d", long, default_value = "3")]
    dimension: usize,

    /// Number of particles, only used with RSA generator
    #[structopt(short = "n", long, default_value = "10")]
    number: usize,

    /// Work with ellipsoids? only used with RSA generator
    #[structopt(long)]
    ellipsoid: bool,

    /// X-aligned semi-axis, only used with RSA generator
    #[structopt(short = "a", default_value = "1.0")]
    a_semi: f64,

    /// Y-aligned semi-axis, only used with RSA generator
    #[structopt(short = "b", default_value = "1.0")]
    b_semi: f64,
    
    /// Z-aligned semi-axis, only used with RSA generator
    #[structopt(short = "c", default_value = "1.0")]
    c_semi: f64,

    /// Side length of initial cubic cell, only used with RSA generator
    #[structopt(short = "s", long, default_value = "10.0")]
    side: f64,

    /// Root filename to save configs to.
    /// If unspecified, only give initial and final
    /// config on stdout. Will write all sweeps
    /// to file if specified and initial/final configuration
    /// if specified.
    #[structopt(short = "o", long)]
    savefiles: Option<PathBuf>,

    /// Option to turn on saving of entire trajectory
    /// If not, will save a backup file instead
    #[structopt(long)]
    save_trajectory: bool,

    /// Root filename to save log tables to.
    #[structopt(long)]
    logfiles: Option<PathBuf>,

    /// Number of sweeps to do
    #[structopt(long, default_value = "100")]
    sweeps: usize,

    /// Number of moves to do per sweep, would
    /// recommend at least 10x particle number
    #[structopt(long, default_value = "1000")]
    moves: usize,

    /// Decides whether to use combined or separate
    /// translations/rotations, only applicable to
    /// ellipsoids
    #[structopt(long)]
    combined_move: bool,

    /// Probability of doing a cell move
    #[structopt(long, default_value = "0.001")]
    pcell: f64,

    /// Std dev of gaussian distribution
    /// for translational moves, as fraction
    /// of particle radius/minor semiaxis
    #[structopt(long, default_value = "0.3")]
    dtrans: f64,

    /// Std dev of gaussian distribution for
    /// rotational moves, using quaternion
    /// addition algorithm in 3d
    #[structopt(long, default_value = "0.1")]
    drot: f64,

    /// Rotation cap to prevent unnecessary runaway rotations
    /// in "uniform" (dilute fluid) regime,
    /// only does anything if --combined_move is *inactive*
    #[structopt(long, default_value = "100.0")]
    drot_cap: f64,

    /// Radius of spherical particles, only used with RSA generator
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
    /// adjust_upper and adjust_lower?
    #[structopt(long)]
    adjust: bool,

    /// Gives the lower bound to use with --adjust
    #[structopt(long, default_value = "0.2")]
    adjust_lower: f64,

    /// Gives the upper bound to use with --adjust
    #[structopt(long, default_value = "0.3")]
    adjust_upper: f64,

    /// Optional file holding seed in 32 by u8 format
    #[structopt(long)]
    seedfile: Option<PathBuf>,

    /// Try and parallelize inner loops? Not always worth it.
    #[structopt(long)]
    parallelize_inner: bool,

    /// Optional file holding initial configuration
    #[structopt(long)]
    initfile: Option<PathBuf>,

    /// Turns off clamping of strain parameters, warning, may
    /// violate equilibrium conditions
    #[structopt(long)]
    no_clamp: bool,

    /// Brent absolute tolerance, used in ellipsoid overlap
    #[structopt(long, default_value = "1e-7")]
    brent_abs_tol: f64,

    /// Brent relative tolerance, used in ellipsoid overlap
    #[structopt(long, default_value = "0.0")]
    brent_rel_tol: f64,

    /// Brent max iterations, used in ellipsoid overlap
    #[structopt(long, default_value = "100")]
    brent_max_iter: usize,

    /// Multiplier for decreaasing step size
    /// when using --adjust. Multiplier for
    /// increasing step is just 1/de
    #[structopt(long, default_value = "0.9")]
    adjust_dec_mult: f64,
}

// Globals -------------------------------------------------------------------

// Computed with Mathematica Student 11.2.0.0
const PI: f64 = 3.141592653589793;

lazy_static! {
    static ref OPT: Opt = Opt::from_args();
}

// Helper functions for command line -----------------------------------------

fn make_and_run_schedule<C, P: Particle + Clone + Debug + Display + Send + Sync>
        (mut config: C, rng: &mut Xoshiro256StarStar)
    where C: Asc<P>
{
    if OPT.check_overlap {
        let validity = config.is_valid();
        println!("Found configuration is {}", validity);
    } else {
        if OPT.savefiles.is_none() {
            println!("Initial Configuration:");
            config.print_asc();
        }
        save_asc_from_opt(&config, "initial");
        let mut schedule = Schedule::make(config.first_particle());
        println!("Running schedule:");
        println!("{:?}", schedule);
        schedule.run(&mut config, rng);
        println!("Ended schedule:");
        println!("{:?}", schedule);
        if OPT.savefiles.is_none() {
            println!("Final Configuration:");
            config.print_asc();
        }
        save_asc_from_opt(&config, "final");
        println!("Ending program at: {}", Utc::now());
    }
}

fn consume<T>(_arg: T) {}

fn initialize_rng() -> Xoshiro256StarStar {
    if let Some(path) = &OPT.seedfile {
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
    }
}

// From 
// https://stackoverflow.com/questions/39204908/how-to-check-release-debug-builds-using-cfg-in-rust
// and
// https://vallentin.io/2019/06/06/versioning
#[cfg(debug_assertions)]
const BUILD_MODE: &str = "debug";
#[cfg(not(debug_assertions))]
const BUILD_MODE: &str = "release";

#[cfg(feature = "using_make")]
const USING_MAKE: &str = "true";
#[cfg(not(feature = "using_make"))]
const USING_MAKE: &str = "false";

fn main() {
// Print out build time information
    println!("Build info for asc_monte_carlo:");
    println!("Using make: {}", USING_MAKE);
    println!("Cargo version: {}", env!("C_VER"));
    println!("Commit SHA: {}", env!("VERGEN_SHA"));
    println!("Commit date: {}", env!("VERGEN_COMMIT_DATE"));
    println!("Version: {}", env!("VERGEN_SEMVER"));
    println!("Build: {}", BUILD_MODE);
    println!("Build time: {}", env!("VERGEN_BUILD_TIMESTAMP"));
    println!("Compiler version: {}", env!("V_RUSTC"));
    println!("RUSTFLAGS: {}", env!("S_RUSTFLAGS"));
    println!("Target: {}", env!("VERGEN_TARGET_TRIPLE"));
    println!("Clean working directory for build: {}", env!("WD_IS_CLEAN"));
    println!("Start program:\n");
// Immediately read command line
    consume(&*OPT);

// Report runtime parameters ---------------------------------------
    println!("Starting program at: {}", Utc::now());
    println!("Given command line options:");
    println!("{:?}", *OPT);
    
// Make and run the correct simulation -----------------------------------
    
    let mut rng = initialize_rng();
    
    // Read from file
    if let Some(path) = &OPT.initfile {
        // Figure out what type of particle from line one
        let mut infile = BufReader::new(File::open(path).expect("Input file must be valid"));
        let mut buf = String::new();
        infile.read_line(&mut buf).expect("Valid utf-8");
        let line_one = buf.split_whitespace();
        let shape: &str = line_one.last().expect("Valid shape");
        match shape {
            "Ellipse" => {
                let config: OverboxList<Ellipse> = OverboxList::from_file(path);
                if OPT.cell_list {
                    let cell_list = CellList::from_overbox_list(config);
                    make_and_run_schedule(cell_list, &mut rng);
                } else {
                    make_and_run_schedule(config, &mut rng);
                }
            }
            "Disk" => {
                let config: OverboxList<Disk> = OverboxList::from_file(path);
                if OPT.cell_list {
                    let cell_list = CellList::from_overbox_list(config);
                    make_and_run_schedule(cell_list, &mut rng);
                } else {
                    make_and_run_schedule(config, &mut rng);
                }
            }
            "Ellipsoid" => {
                let config: OverboxList<Ellipsoid> = OverboxList::from_file(path);
                if OPT.cell_list {
                    let cell_list = CellList::from_overbox_list(config);
                    make_and_run_schedule(cell_list, &mut rng);
                } else {
                    make_and_run_schedule(config, &mut rng);
                }
            }
            "Sphere" => {
                let config: OverboxList<Sphere> = OverboxList::from_file(path);
                if OPT.cell_list {
                    let cell_list = CellList::from_overbox_list(config);
                    make_and_run_schedule(cell_list, &mut rng);
                } else {
                    make_and_run_schedule(config, &mut rng);
                }
            }
            _ => panic!("Invalid shape/shape not implemented"),
        }
    } else {
        // Make initial config through RSA
        println!("Running rsa...");
        match OPT.dimension {
            2 => {
                if OPT.ellipsoid {
                    let shape = Ellipse::make_shape(OPT.a_semi, OPT.b_semi);
                    let init_cell = vec![OPT.side, 0.0, 0.0, OPT.side];
                    let config = OverboxList::make_rsa(OPT.number, &shape, 2, init_cell, &mut rng);
                    if OPT.cell_list {
                        let cell_list = CellList::from_overbox_list(config);
                        make_and_run_schedule(cell_list, &mut rng);
                    } else {
                        make_and_run_schedule(config, &mut rng);
                    }
                } else {
                    let shape = Disk::make_shape(OPT.radius);
                    let init_cell = vec![OPT.side, 0.0, 0.0, OPT.side];
                    let config = OverboxList::make_rsa(OPT.number, &shape, 2, init_cell, &mut rng);
                    if OPT.cell_list {
                        let cell_list = CellList::from_overbox_list(config);
                        make_and_run_schedule(cell_list, &mut rng);
                    } else {   
                        make_and_run_schedule(config, &mut rng);
                    }
                }
            }
            3 => {
                if OPT.ellipsoid {
                    let shape = Ellipsoid::make_shape(OPT.a_semi, OPT.b_semi, OPT.c_semi);
                    let init_cell = vec![OPT.side, 0.0, 0.0,
                                            0.0, OPT.side, 0.0,
                                            0.0, 0.0, OPT.side];
                    let config = OverboxList::make_rsa(OPT.number, &shape, 3, init_cell, &mut rng);
                    if OPT.cell_list {
                        let cell_list = CellList::from_overbox_list(config);
                        make_and_run_schedule(cell_list, &mut rng);
                    } else {
                        make_and_run_schedule(config, &mut rng);
                    }
                } else {
                    let shape = Sphere::make_shape(OPT.radius);
                    let init_cell = vec![OPT.side, 0.0, 0.0,
                                            0.0, OPT.side, 0.0,
                                            0.0, 0.0, OPT.side];
                    let config = OverboxList::make_rsa(OPT.number, &shape, 3, init_cell, &mut rng);
                    if OPT.cell_list {
                        let cell_list = CellList::from_overbox_list(config);
                        make_and_run_schedule(cell_list, &mut rng);
                    } else {
                        make_and_run_schedule(config, &mut rng);
                    }
                }
            }
            _ => panic!("Can't handle that dimension yet!"),
        }
    }
}

