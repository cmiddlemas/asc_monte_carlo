use crate::asc::{Asc, Particle};
use rand_xoshiro::Xoshiro256StarStar;
use rand_distr::{Uniform, Distribution};
use std::marker::PhantomData;
use std::fmt::{Debug, Display};
use crate::OPT;

const ACC_LOWER_LIM: f64 = 0.4;
const ACC_UPPER_LIM: f64 = 0.6;
const DECREASE_MOD: f64 = 0.9;
// Value calculated from Mathematica Student 11.2.0.0
const INCREASE_MOD: f64 = 1.111111111111111;

#[derive(Debug)]
pub struct Schedule<P> {
    particle_accepts: usize,
    particle_tries: usize,
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
    adjust_params: bool,
    _phantom: PhantomData<P>,
}

// https://stackoverflow.com/questions/42613974/why-cant-i-add-a-blanket-impl-on-a-trait-with-a-type-parameter
impl<P: Particle + Debug + Display + Send + Sync + Clone> Schedule<P> {
    pub fn make() -> Schedule<P> {
        Schedule {
            particle_accepts: 0,
            particle_tries: 0,
            cell_accepts: 0,
            cell_tries: 0,
            n_sweeps: OPT.sweeps,
            n_moves: OPT.moves,
            pressure: OPT.ipressure,
            p_cell_move: OPT.pcell,
            cell_param: vec![OPT.isotropic, OPT.shear, OPT.axial],
            particle_param: vec![OPT.dtrans*OPT.radius],
            running_obs: P::init_obs(),
            beta: 1.0, // TODO: decide if this should be 
                       // accessible to user
            adjust_params: OPT.adjust,
            _phantom: PhantomData,
        }
    }

    pub fn run(&mut self,
               config: &mut Asc<P>,
               rng: &mut Xoshiro256StarStar)
    {
        let u_dist = Uniform::new(0.0, 1.0);
        for i in 0..self.n_sweeps {
            for j in 0..self.n_moves {
                if u_dist.sample(rng) < self.p_cell_move {
                    self.cell_tries += 1;
                    if config.try_cell_move(self, rng) {
                        self.cell_accepts += 1;
                    }
                } else {
                    self.particle_tries += 1;
                    if config.try_particle_move(self, rng) {
                        self.particle_accepts += 1;
                    }
                }
            }
            P::sample_obs_sweep(self, &config);
            self.post_sweep();
        }
    }

    fn post_sweep(&mut self) { 
        // Adjust MC move parameters based on acceptance ratio
        // This simple algorithm keeps ratio between parameters the same
        if self.adjust_params {
            let p_acc_ratio = self.particle_accepts as f64/(self.particle_tries as f64);
            let c_acc_ratio = self.cell_accepts as f64/(self.cell_tries as f64);
            println!("Particle, cell accept ratio: {}, {}", p_acc_ratio, c_acc_ratio);
            // For next sweep
            self.particle_accepts = 0; self.particle_tries = 0;
            self.cell_accepts = 0; self.cell_tries = 0;

            if p_acc_ratio < ACC_LOWER_LIM {
                for val in self.particle_param.iter_mut() {
                    *val *= DECREASE_MOD;
                }
            } else if p_acc_ratio > ACC_UPPER_LIM {
                for val in self.particle_param.iter_mut() {
                    *val *= INCREASE_MOD;
                }
            }

            if c_acc_ratio < ACC_LOWER_LIM {
                for val in self.cell_param.iter_mut() {
                    *val *= DECREASE_MOD;
                }
            } else if c_acc_ratio > ACC_UPPER_LIM {
                for val in self.cell_param.iter_mut() {
                    *val *= INCREASE_MOD;
                }
            }

            for val in self.cell_param.iter_mut() {
                if *val > 0.01 {
                    println!("Clamping a strain parameter to 0.01!");
                    *val = 0.01;
                }
            }
            println!("Params are (p,c): {:?}, {:?}", self.particle_param, self.cell_param);
        }
    }
}
