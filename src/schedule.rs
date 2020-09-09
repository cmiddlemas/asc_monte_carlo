use crate::asc::Asc;
use rand_xoshiro::Xoshiro256StarStar;
use rand_distr::{Uniform, Distribution};
use std::marker::PhantomData;
use std::fmt::{Debug, Display};
use std::path::PathBuf;
use std::fs::{File, OpenOptions};
use std::io::{Read, Write};
use crate::OPT;
use crate::particle::Particle;
use serde::{Serialize, Deserialize};
use kiss3d::window::Window;

pub fn write_data_file(contents: &str, suffix: &str) {
    if OPT.verbosity >= 1 {
        if let Some(logroot) = &OPT.logfiles {
            let mut logpath = logroot.clone();
            logpath.set_file_name(
                format!("{}_{}.dat",
                    logroot.file_name().expect("Must give valid root filename for logfile")
                        .to_str().expect("Valid utf-8"),
                    suffix
                ));
            // Referenced rust docs for std
            let mut logfile = OpenOptions::new()
                .append(true)
                .create(true)
                .open(logpath)
                .expect("Must be able to write to data file");
            write!(&mut logfile, "{}", contents)
                .expect("Failed write to data file");
        }
    }
}

pub fn write_sweep_log(logline: &str) {
    if let Some(logroot) = &OPT.logfiles {
        let mut logpath = logroot.clone();
        logpath.set_file_name(
            format!("{}_sweep_log.dat",
                logroot.file_name().expect("Must give valid root filename for logfile")
                    .to_str().expect("Valid utf-8")
            ));
        // Referenced rust docs for std
        let mut logfile = OpenOptions::new()
            .append(true)
            .create(true)
            .open(logpath)
            .expect("Must be able to write to logfile");
        writeln!(&mut logfile, "{}", logline)
            .expect("Failed write to logfile");
    }
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct ObservableTracker {
    pub n_samples: f64,
    pub sum_of_vol: f64,
    pub sum_of_ar: f64, // Aspect ratio of unit cell
    pub min_nn_gap: f64,
    pub idx_min_nn_gap: usize,
    pub next_to_min_nn_gap: f64,
    pub idx_next_to_min_nn_gap: usize,
    pub sum_of_min_nn_gap: f64,
}

impl ObservableTracker {
    // Basic observable tracker for those particle
    // types which only use trivial fields like
    // volume and aspect ratio
    pub fn new() -> ObservableTracker {
        ObservableTracker {
            n_samples: 0.0,
            sum_of_vol: 0.0,
            sum_of_ar: 0.0,
            min_nn_gap: 0.0,
            idx_min_nn_gap: 0,
            next_to_min_nn_gap: 0.0,
            idx_next_to_min_nn_gap: 0,
            sum_of_min_nn_gap: 0.0
        }
    }
}

// https://github.com/serde-rs/serde
#[derive(Serialize, Deserialize, Debug)]
pub struct Schedule<P,C> {
    pub current_sweep: usize,
    pub particle_accepts: Vec<usize>,
    pub particle_tries: Vec<usize>,
    cell_accepts: usize,
    cell_tries: usize,
    n_sweeps: usize,
    n_moves: usize,
    pub pressure: f64,
    p_cell_move: f64,
    pub cell_param: Vec<f64>,
    pub particle_param: Vec<f64>,
    pub running_obs: ObservableTracker,
    pub beta: f64,
    pub phi: f64,
    pub avg_gap: f64,
    _phantom1: PhantomData<P>,
    _phantom2: PhantomData<C>,
}

// https://stackoverflow.com/questions/42613974/why-cant-i-add-a-blanket-impl-on-a-trait-with-a-type-parameter
impl<P: Particle + Debug + Display + Send + Sync + Clone, C: Asc<P>> Schedule<P,C> {
    pub fn make(particle: &P, config: &C) -> Schedule<P,C> {
        // Decide on particle move parameters
        let mut p_param = Vec::new();
        let size_hint = particle.hint_lower();
        p_param.push(OPT.dtrans*size_hint);
        p_param.push(OPT.drot);
        
        // Make initial schedule
        let mut schedule = Schedule {
            current_sweep: 0,
            particle_accepts: vec![0,0],
            particle_tries: vec![0,0],
            cell_accepts: 0,
            cell_tries: 0,
            n_sweeps: OPT.sweeps,
            n_moves: OPT.moves,
            pressure: OPT.pressure,
            p_cell_move: OPT.pcell,
            cell_param: vec![OPT.isotropic, OPT.shear, OPT.axial],
            particle_param: p_param,
            running_obs: P::init_obs(config),
            beta: 1.0, // TODO: decide if this should be 
                       // accessible to user
            phi: 0.0,
            avg_gap: 0.0,
            _phantom1: PhantomData,
            _phantom2: PhantomData,
        };
        
        // Clamp strain parameters if needed
        if !OPT.no_clamp {
            for val in schedule.cell_param.iter_mut() {
                if *val > 0.01 {
                    println!("Clamping a strain parameter to 0.01!");
                    *val = 0.01;
                }
            }
        }

        // Cap rotations if need be
        if !OPT.combined_move {
            if schedule.particle_param[1] > OPT.drot_cap {
                println!("Capping rotation");
                schedule.particle_param[1] = OPT.drot_cap;
            }
        }

        schedule
    }

    // Using serde_yaml because json can't represent inf
    // properly
    // https://github.com/serde-rs/json/issues/202
    // https://stackoverflow.com/questions/1423081/json-left-out-infinity-and-nan-json-status-in-ecmascript
    pub fn save_in_log(&self) {
        if let Some(logroot) = &OPT.logfiles {
            let mut logpath = logroot.clone();
            logpath.set_file_name(
                format!("{}_schedule.yaml",
                    logroot.file_name().expect("Must give valid root filename for logfile")
                        .to_str().expect("Valid utf-8")
                ));
            // Referenced rust docs for std
            let mut logfile = OpenOptions::new()
                .append(true)
                .create(true)
                .open(logpath)
                .expect("Must be able to write to logfile");
            // https://github.com/serde-rs/serde
            let serial = serde_yaml::to_string(self).unwrap();
            write!(&mut logfile, "{}", serial)
                .expect("Failed write to logfile");
        }
    }

    // Sets a schedule using a yaml schedule. Overwrites
    // only one field using the input from the command line, which is n_sweeps,
    // since that just determines the length of the simulation,
    // and has no physical effect
    pub fn from_file(path: &PathBuf, config: &C) -> Schedule<P,C> {
        // Made extensive use of the rust docs for various
        // std components when writing this function
        let mut infile = File::open(path).expect("Schedule file must be valid");
        let mut buf = String::new();
        infile.read_to_string(&mut buf).expect("Must be able to read Schedule file");
        // https://github.com/serde-rs/serde
        let mut schedule: Schedule<P,C> = serde_yaml::from_str(&buf).unwrap();
        schedule.n_sweeps = OPT.sweeps;
        schedule.n_moves = OPT.moves;
        schedule.pressure = OPT.pressure;
        schedule.p_cell_move = OPT.pcell;

        // Need to reset some counters, since some runs
        // leave the counters unreset at the end so that
        // things like accept ratios can be computed after
        // the simulation
        schedule.particle_accepts = vec![0, 0];
        schedule.particle_tries = vec![0, 0];
        schedule.cell_accepts = 0;
        schedule.cell_tries = 0;
        schedule.running_obs = P::init_obs(config);
        schedule
    }

    fn should_terminate(&self) -> bool {
        if let Some(threshold) = OPT.max_phi {
            if self.phi >= threshold {
                return true;
            }
        } else if let Some(threshold) = OPT.gap_threshold {
            if self.avg_gap/self.phi <= threshold {
                return true;
            }
        }
        false
    }

    pub fn run(&mut self,
               config: &mut C,
               rng: &mut Xoshiro256StarStar,
               window: &mut Option<Window>)
    {
        let u_dist = Uniform::new(0.0, 1.0);
        for _i in 0..self.n_sweeps {
            for _j in 0..self.n_moves {
                if u_dist.sample(rng) < self.p_cell_move {
                    self.cell_tries += 1;
                    if config.try_cell_move(self, rng) {
                        self.cell_accepts += 1;
                    }
                } else {
                    if config.try_particle_move(self, rng) {
                    }
                }
            }
            P::sample_obs_sweep(self, config);
            // https://doc.rust-lang.org/stable/rust-by-example/scope/borrow/ref.html
            // https://www.reddit.com/r/rust/comments/bn1e5o/what_does_mut_in_mut_some_mean/
            if let Some(ref mut w) = *window {
                P::render_packing(w, config);
            }
            self.post_sweep();
            if self.should_terminate() {
                return;
            }
        }
    }

    fn post_sweep(&mut self) {
        //eprintln!("Sweep the roads!");
        self.current_sweep += 1;
        // Adjust MC move parameters based on acceptance ratio
        // This simple algorithm keeps ratio between parameters the same
        if OPT.adjust {
            let p_acc_ratio_0 = self.particle_accepts[0] as f64/(self.particle_tries[0] as f64);
            let p_acc_ratio_1 = self.particle_accepts[1] as f64/(self.particle_tries[1] as f64);
            let c_acc_ratio = self.cell_accepts as f64/(self.cell_tries as f64);
            println!("Translation (combined), rotation, cell accept ratio: {}, {}, {}",
                p_acc_ratio_0, p_acc_ratio_1, c_acc_ratio
            );
            // For next sweep
            self.particle_accepts = vec![0,0]; self.particle_tries = vec![0,0];
            self.cell_accepts = 0; self.cell_tries = 0;

            let decrease_mult = OPT.adjust_dec_mult;
            let increase_mult = 1.0/decrease_mult;

            if OPT.combined_move { // only use acc_ratio_0, for combined move
                if p_acc_ratio_0 < OPT.adjust_lower {
                    for val in self.particle_param.iter_mut() {
                        *val *= decrease_mult;
                }
                } else if p_acc_ratio_0 > OPT.adjust_upper {
                    for val in self.particle_param.iter_mut() {
                        *val *= increase_mult;
                    }
                }
            } else {
                if p_acc_ratio_0 < OPT.adjust_lower {
                    self.particle_param[0] *= decrease_mult;
                } else if p_acc_ratio_0 > OPT.adjust_upper {
                    self.particle_param[0] *= increase_mult;
                }

                if p_acc_ratio_1 < OPT.adjust_lower {
                    self.particle_param[1] *= decrease_mult;
                } else if p_acc_ratio_1 > OPT.adjust_upper {
                    self.particle_param[1] *= increase_mult;
                }
                if self.particle_param[1] > OPT.drot_cap {
                    println!("Capping rotations");
                    self.particle_param[1] = OPT.drot_cap;
                }
            }

            if c_acc_ratio < OPT.adjust_lower {
                for val in self.cell_param.iter_mut() {
                    *val *= decrease_mult;
                }
            } else if c_acc_ratio > OPT.adjust_upper {
                for val in self.cell_param.iter_mut() {
                    *val *= increase_mult;
                }
            }
            
            if !OPT.no_clamp {
                for val in self.cell_param.iter_mut() {
                    if *val > 0.01 {
                        println!("Clamping a strain parameter to 0.01!");
                        *val = 0.01;
                    }
                }
            }
            println!("Params are (p,c): {:?}, {:?}", self.particle_param, self.cell_param);
        }
    }
}
