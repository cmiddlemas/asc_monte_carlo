use std::io::prelude::*;
use std::fmt::{Debug, Display, Formatter};
use std::fmt;
use rand_xoshiro::Xoshiro256StarStar;
use rand_distr::{Uniform, Distribution};
use rand::SeedableRng;
use std::path::PathBuf;
use structopt::StructOpt;
use lazy_static::lazy_static;

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
    #[structopt(long, default_value = "0.0")]
    pcell: f64,

    /// Std dev of gaussian distribution
    /// for translational moves, as fraction
    /// of particle radius
    #[structopt(long, default_value = "1.0")]
    dtrans: f64,

    /// Radius of spherical particles
    #[structopt(long, default_value = "1.0")]
    radius: f64,
}

lazy_static! {
    static ref OPT: Opt = Opt::from_args();
}

fn main() {
    println!("{:?}", *OPT);
    let shape = Disk::make_shape(OPT.radius);
    let mut init_cell = vec![OPT.side, 0.0, 0.0, OPT.side];
    let mut rng = Xoshiro256StarStar::seed_from_u64(0);
    let mut config = Asc::make_rsa(OPT.number, &shape, 2, init_cell, &mut rng);
    //config.print_asc();
    let mut schedule = Schedule::make();
    schedule.run(&mut config, &mut rng);
    println!("{:?}", schedule);
    //config.print_asc();
}
