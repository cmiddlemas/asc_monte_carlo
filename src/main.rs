use std::io::prelude::*;
use std::fmt::{Debug, Display, Formatter};
use std::fmt;
use rand_xoshiro::Xoshiro256StarStar;
use rand_distr::{Uniform, Distribution};
use rand::SeedableRng;
use std::path::PathBuf;
use structopt::StructOpt;
use lazy_static::lazy_static;
use chrono::prelude::*;

mod asc;
mod disks;
mod schedule;

use asc::{Asc, Particle};
use disks::{Disk};
use schedule::Schedule;

// Struct for command line
#[derive(StructOpt, Debug)]
#[structopt(name = "asc_monte_carlo")]
struct Opt {
    /// Number of particles
    #[structopt(short = "n", long, default_value = "10")]
    number: usize,

    /// Side length of initial cubic cell
    #[structopt(short = "s", long, default_value = "10.0")]
    side: f64,

    /// Root filename to save configs as
    /// Default to only give initial and final
    /// config on stdout, will give all sweeps
    /// if specified
    #[structopt(short = "o", long)]
    savefiles: Option<PathBuf>,

    /// Root filename to save logs as
    /// All logs default to go on stdout
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
    /// of particle radius
    #[structopt(long, default_value = "0.3")]
    dtrans: f64,

    /// Radius of spherical particles
    #[structopt(long, default_value = "1.0")]
    radius: f64,

    /// Initial pressure of system
    #[structopt(long, default_value = "1.0")]
    ipressure: f64,

    /// Initial shear step width
    #[structopt(long, default_value = "0.1")]
    shear: f64,

    /// Initial axial compression step width
    #[structopt(long, default_value = "0.1")]
    axial: f64,

    /// Initial isotropic volume change step width
    #[structopt(long, default_value = "0.1")]
    isotropic: f64,

    /// Should the program adjust the step sizes
    /// to try and keep acceptance ratio between
    /// 0.4 and 0.6?
    #[structopt(long)]
    adjust: bool,
}

lazy_static! {
    static ref OPT: Opt = Opt::from_args();
}

fn main() {
    println!("Starting program at: {}", Utc::now());
    println!("Given command line options:");
    println!("{:?}", *OPT);
    let shape = Disk::make_shape(OPT.radius);
    let mut init_cell = vec![OPT.side, 0.0, 0.0, OPT.side];
    let mut rng = Xoshiro256StarStar::seed_from_u64(0);
    let mut config = Asc::make_rsa(OPT.number, &shape, 2, init_cell, &mut rng);
    println!("Initial Configuration:");
    config.print_asc();
    let mut schedule = Schedule::make();
    println!("Running schedule:");
    println!("{:?}", schedule);
    schedule.run(&mut config, &mut rng);
    println!("Ended schedule:");
    println!("{:?}", schedule);
    println!("Final Configuration:");
    config.print_asc();
    println!("Ending program at: {}", Utc::now());
}

