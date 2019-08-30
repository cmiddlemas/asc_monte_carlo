use crate::asc::{Asc, Particle};
use rand_xoshiro::Xoshiro256StarStar;
use rand_distr::{Uniform, Distribution};
use std::marker::PhantomData;
use std::fmt::{Debug, Display};
use crate::OPT;

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
    _phantom: PhantomData<P>,
}

// https://stackoverflow.com/questions/42613974/why-cant-i-add-a-blanket-impl-on-a-trait-with-a-type-parameter
impl<P: Particle + Debug + Display + Send + Sync> Schedule<P> {
    pub fn make() -> Schedule<P> {
        Schedule {
            particle_accepts: 0,
            particle_tries: 0,
            cell_accepts: 0,
            cell_tries: 0,
            n_sweeps: OPT.sweeps,
            n_moves: OPT.moves,
            pressure: 1.0,
            p_cell_move: OPT.pcell,
            cell_param: Vec::new(),
            particle_param: vec![OPT.dtrans*OPT.radius],
            running_obs: Vec::new(),
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
        P::sample_obs_sweep(self, &config);
    }

    fn post_sweep(&mut self) {
        return;
    }
}
