use crate::asc::{Asc, Particle};
use rand_xoshiro::Xoshiro256StarStar;
use rand_distr::{Uniform, Distribution};
use std::marker::PhantomData;
use std::fmt::{Debug, Display};
use std::fs::OpenOptions;
use std::io::Write;
use crate::OPT;

const DECREASE_MOD: f64 = 0.9;
// Value calculated from Mathematica Student 11.2.0.0
const INCREASE_MOD: f64 = 1.111111111111111;

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

#[derive(Debug)]
pub struct Schedule<P> {
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
    pub running_obs: Vec<f64>,
    pub beta: f64,
    pub phi: f64,
    _phantom: PhantomData<P>,
}

// https://stackoverflow.com/questions/42613974/why-cant-i-add-a-blanket-impl-on-a-trait-with-a-type-parameter
impl<P: Particle + Debug + Display + Send + Sync + Clone> Schedule<P> {
    pub fn make(particle: &P) -> Schedule<P> {
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
            pressure: OPT.ipressure,
            p_cell_move: OPT.pcell,
            cell_param: vec![OPT.isotropic, OPT.shear, OPT.axial],
            particle_param: p_param,
            running_obs: P::init_obs(),
            beta: 1.0, // TODO: decide if this should be 
                       // accessible to user
            phi: 0.0,
            _phantom: PhantomData,
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

    fn should_terminate(&self) -> bool {
        if let Some(threshold) = OPT.max_phi {
            if self.phi >= threshold {
                return true;
            }
        }
        false
    }

    pub fn run<C: Asc<P>>(&mut self,
               config: &mut C,
               rng: &mut Xoshiro256StarStar)
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
            self.post_sweep();
            if self.should_terminate() {
                return;
            }
        }
    }

    fn post_sweep(&mut self) {
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

            if OPT.combined_move { // only use acc_ratio_0, for combined move
                if p_acc_ratio_0 < OPT.adjust_lower {
                    for val in self.particle_param.iter_mut() {
                        *val *= DECREASE_MOD;
                }
                } else if p_acc_ratio_0 > OPT.adjust_upper {
                    for val in self.particle_param.iter_mut() {
                        *val *= INCREASE_MOD;
                    }
                }
            } else {
                if p_acc_ratio_0 < OPT.adjust_lower {
                    self.particle_param[0] *= DECREASE_MOD;
                } else if p_acc_ratio_0 > OPT.adjust_upper {
                    self.particle_param[0] *= INCREASE_MOD;
                }

                if p_acc_ratio_1 < OPT.adjust_lower {
                    self.particle_param[1] *= DECREASE_MOD;
                } else if p_acc_ratio_1 > OPT.adjust_upper {
                    self.particle_param[1] *= INCREASE_MOD;
                }
                if self.particle_param[1] > OPT.drot_cap {
                    println!("Capping rotations");
                    self.particle_param[1] = OPT.drot_cap;
                }
            }

            if c_acc_ratio < OPT.adjust_lower {
                for val in self.cell_param.iter_mut() {
                    *val *= DECREASE_MOD;
                }
            } else if c_acc_ratio > OPT.adjust_upper {
                for val in self.cell_param.iter_mut() {
                    *val *= INCREASE_MOD;
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
