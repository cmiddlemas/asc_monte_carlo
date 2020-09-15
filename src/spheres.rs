use crate::particle::Particle;
use rand_xoshiro::Xoshiro256StarStar;
use rand_distr::{Uniform, Distribution, Normal};
use std::fmt::{Display, Formatter};
use std::fmt;
use std::convert::TryInto;
use crate::asc::{Asc, save_asc_from_opt};
use crate::schedule::{ObservableTracker, Schedule, write_sweep_log, write_data_file};
use std::f64::consts::PI;
use crate::common_util::{negate, apply_pbc, relative_to_global3, global_to_relative3};
use kiss3d::window::Window;
// https://github.com/sebcrozet/kiss3d/issues/66
// sebcrozet evidently reversed decision in that thread,
// so we don't need to manually resolve these dependencies ourselves!
use kiss3d::nalgebra::{Translation3, Point3};
use kiss3d::scene::SceneNode;
use lazy_static::lazy_static;

lazy_static! {
    static ref RED: Point3<f32> = Point3::new(1.0, 0.0, 0.0);
}

lazy_static! {
    static ref GREEN: Point3<f32> = Point3::new(0.0, 1.0, 0.0);
}

lazy_static! {
    static ref BLUE: Point3<f32> = Point3::new(0.0, 0.0, 1.0);
}

// https://stackoverflow.com/questions/26958178/how-do-i-automatically-implement-comparison-for-structs-with-floats-in-rust
#[derive(Debug, Clone)]
pub struct Sphere {
    pub rel_pos: [f64; 3],
    pub global_pos: [f64; 3],
    pub radius: f64,
}

impl Sphere {
    
    pub fn make_shape(r: f64) -> Self {
        Sphere { rel_pos: [0.0, 0.0, 0.0], global_pos: [0.0, 0.0, 0.0], radius: r }
    }
    
    fn squared_scaling(&self, other: &Self, offset: &[f64]) -> f64 {
        let image_x = other.global_pos[0] + offset[0];
        let image_y = other.global_pos[1] + offset[1];
        let image_z = other.global_pos[2] + offset[2];
        let displacement2 = (self.global_pos[0] - image_x).powi(2) 
                + (self.global_pos[1] - image_y).powi(2)
                + (self.global_pos[2] - image_z).powi(2);
        let extent2 = (self.radius + other.radius).powi(2);
        extent2/displacement2
    }

    // Must be symmetric to avoid precision problems
    // No roots for efficiency
    fn squared_symmetric_squared_scaling(&self, other: &Self, offset: &[f64]) -> f64 {
        self.squared_scaling(other, offset)*other.squared_scaling(self, &negate(offset))
    }

}

impl Display for Sphere {
    // From rust docs
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(f, "{} {} {} {}",
               self.rel_pos[0], self.rel_pos[1], self.rel_pos[2], self.radius
        )
    }
}

impl Particle for Sphere {
    const TYPE: &'static str = "Sphere";

    fn parse(line: &str, unit_cell: &[f64]) -> Self {
        let params: Vec<f64> = line.split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect();
        let rel_pos = [params[0], params[1], params[2]];
        let global_pos = relative_to_global3(unit_cell, &rel_pos);
        Sphere { rel_pos, global_pos, radius: params[3] } 
    }

    // We recommend using same numerical function for the next two,
    // since they always need to give consistent results with each other.
    fn check_overlap(&self, other: &Self, offset: &[f64]) -> bool {
        // This would probably work too
        //(self.squared_scaling(other, offset)) > 1.0 || (other.squared_scaling(self, &negate(offset)) > 1.0)
        self.squared_symmetric_squared_scaling(other, offset) > 1.0
    }
    
    fn overlap_scale(&self, other: &Self, offset: &[f64]) -> f64 {
        self.squared_symmetric_squared_scaling(other, offset).powf(0.75)
    }

    fn copy_shape_random_coord(&self,
                               cell: &[f64],
                               rng: &mut Xoshiro256StarStar
    ) -> Self
    {
        let uni_dist = Uniform::new(0.0, 1.0);
        let lat_x = uni_dist.sample(rng);
        let lat_y = uni_dist.sample(rng);
        let lat_z = uni_dist.sample(rng);
        let rel_pos = [lat_x, lat_y, lat_z];
        let global_pos = relative_to_global3(cell, &rel_pos);
        Sphere { rel_pos, global_pos, radius: self.radius }
    }

    fn perturb(&mut self,
               cell: &[f64],
               param: &[f64],
               rng: &mut Xoshiro256StarStar,
    ) -> (Self, usize)
    {
        let old_sphere = self.clone();
        let normal = Normal::new(0.0, param[0]).unwrap();
        self.global_pos[0] += normal.sample(rng);
        self.global_pos[1] += normal.sample(rng);
        self.global_pos[2] += normal.sample(rng);
        let uncorrected_rel = global_to_relative3(cell, &self.global_pos);
        // Handle pbc
        self.rel_pos = apply_pbc(&uncorrected_rel).as_slice().try_into().unwrap();
        // Recalculate global
        self.global_pos = relative_to_global3(cell, &self.rel_pos);
        (old_sphere, 0)
    }

    // From S. Torquato and Y. Jiao PRE 80, 041104 (2009)
    // Also need nalgebra
    fn apply_strain(&mut self, new_cell: &[f64]) {
        self.global_pos = relative_to_global3(new_cell, &self.rel_pos);
    }

    fn init_obs<C: Asc<Self>>(config: &C) -> ObservableTracker {
        let (_, min_nn_gap, idx_min_nn_gap, _, _) = config.nn_gap_distribution();
        ObservableTracker {
            n_samples: 0.0,
            sum_of_vol: 0.0,
            sum_of_ar: 0.0,
            min_nn_gap,
            idx_min_nn_gap,
            sum_of_min_nn_gap: 0.0
        }
    }

    fn sample_obs_sweep<C: Asc<Self>>(schedule: &mut Schedule<Self, C>, config: &C) {
        let vol = schedule.running_obs.sum_of_vol/schedule.running_obs.n_samples;
        let avg_min_gap = schedule.running_obs.sum_of_min_nn_gap/schedule.running_obs.n_samples;
        let aspect_ratio = schedule.running_obs.sum_of_ar/schedule.running_obs.n_samples;
        schedule.running_obs = Self::init_obs(config);
        println!("Cell volume averaged over sweep: {}", vol);
        println!("Aspect ratio averaged over sweep: {}", aspect_ratio);
        let phi = config.packing_fraction();
        schedule.phi = phi;
        println!("Phi at end of sweep: {}", phi);
        let (avg_nn_gap,
             min_nn_gap,
             _idx_min_nn_gap,
             max_nn_gap,
             gap_distr) = config.nn_gap_distribution();
        schedule.avg_min_gap = avg_min_gap;
        println!("Average min gap, mean gap: {} {}", avg_min_gap, avg_nn_gap);
        let gap_string: String = gap_distr.iter()
                                          .map(|x| format!("{} {} {}\n", x[0], x[1], x[2]))
                                          .fold(String::new(), |acc, x| acc + &x);
        let fname = format!("nn_gap_{}", schedule.current_sweep);
        write_data_file(&gap_string, &fname);
        let (pressure, unc_pressure, chisq) = config.instantaneous_pressure();
        let logline = format!("{} {} {} {} {} {} {} {} {} {} {}", schedule.current_sweep,
                                                vol,
                                                phi,
                                                pressure,
                                                unc_pressure,
                                                chisq,
                                                avg_min_gap,
                                                avg_nn_gap,
                                                min_nn_gap,
                                                max_nn_gap,
                                                aspect_ratio);
        write_sweep_log(&logline);
        save_asc_from_opt(config, &format!("sweep_{}", schedule.current_sweep));
    }

    fn sample_obs_failed_move<C: Asc<Self>>(
        schedule: &mut Schedule<Self, C>,
        config: &C
    )
    {
        schedule.running_obs.n_samples += 1.0;
        schedule.running_obs.sum_of_vol += config.cell_volume();
        // The new min_nn_gap is the same as the old
        schedule.running_obs.sum_of_min_nn_gap += schedule.running_obs.min_nn_gap;
        schedule.running_obs.sum_of_ar += config.aspect_ratio();
    }

    fn sample_obs_accepted_pmove<C: Asc<Self>>(
        schedule: &mut Schedule<Self, C>,
        config: &C,
        changed_idx: usize,
        _old_p: &Self
    )
    {
        schedule.running_obs.n_samples += 1.0;
        schedule.running_obs.sum_of_vol += config.cell_volume();
        
        let possible_gap = config.particle_gap(changed_idx, &config.particle_slice()[changed_idx]);
        
        if possible_gap <= schedule.running_obs.min_nn_gap {
            schedule.running_obs.min_nn_gap = possible_gap;
            schedule.running_obs.idx_min_nn_gap = changed_idx;
        // Invalidated current minimum gap
        } else if changed_idx == schedule.running_obs.idx_min_nn_gap {
            let (_, min_nn_gap, idx_min_nn_gap, _, _) = config.nn_gap_distribution();
            schedule.running_obs.min_nn_gap = min_nn_gap;
            schedule.running_obs.idx_min_nn_gap = idx_min_nn_gap;
        }
        
        schedule.running_obs.sum_of_min_nn_gap += schedule.running_obs.min_nn_gap;
        
        schedule.running_obs.sum_of_ar += config.aspect_ratio();
    }
    
    fn sample_obs_accepted_cmove<C: Asc<Self>>(
        schedule: &mut Schedule<Self, C>,
        config: &C,
        _old_c: &[f64]
    )
    {
        schedule.running_obs.n_samples += 1.0;
        schedule.running_obs.sum_of_vol += config.cell_volume();
        
        let (_, min_nn_gap, idx_min_nn_gap, _, _) = config.nn_gap_distribution();
        schedule.running_obs.min_nn_gap = min_nn_gap;
        schedule.running_obs.idx_min_nn_gap = idx_min_nn_gap;

        schedule.running_obs.sum_of_min_nn_gap += schedule.running_obs.min_nn_gap;
        
        schedule.running_obs.sum_of_ar += config.aspect_ratio();
    }

    fn render_packing<C: Asc<Self>>(window: &mut Window, config: &C) -> bool {
        // Could use the solution in
        // https://stackoverflow.com/questions/55293051/how-do-i-clear-the-scene-in-kiss3d
        // but since our needs are simple, we're just going to keep track of all
        // of the added objects and unlink them after drawing them.
        let mut node_vec: Vec<SceneNode> = Vec::with_capacity(config.n_particles());
        
        // Place particles
        for p in config.particle_slice() {
            let translation = Translation3::new(p.global_pos[0] as f32,
                                                p.global_pos[1] as f32,
                                                p.global_pos[2] as f32);
            let mut sphere = window.add_sphere(p.radius as f32);
            sphere.set_local_translation(translation);
            node_vec.push(sphere);
        }

        // Place unit cell
        let cell = config.unit_cell();
        let origin = Point3::new(0.0f32, 0.0f32, 0.0f32);
        let u1 = Point3::new(cell[0] as f32,
                             cell[1] as f32,
                             cell[2] as f32);
        window.draw_line(&origin, &u1, &*RED);
        let u2 = Point3::new(cell[3] as f32,
                             cell[4] as f32,
                             cell[5] as f32);
        window.draw_line(&origin, &u2, &*GREEN);
        let u3 = Point3::new(cell[6] as f32,
                             cell[7] as f32,
                             cell[8] as f32);
        window.draw_line(&origin, &u3, &*BLUE);

        let status = window.render();
        
        // Clean up window
        for mut node in node_vec {
            node.unlink();
        }
        
        status
    }

    fn hint_lower(&self) -> f64 {
        self.radius
    }

    fn hint_upper(&self) -> f64 {
        self.radius
    }

    fn vol(&self) -> f64 {
        4.0*PI*self.radius.powi(3)/3.0
    }

    fn lat_coord(&self) -> Vec<f64> {
        self.rel_pos.to_vec()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand_xoshiro::rand_core::SeedableRng;

    #[test]
    fn zero_inversion_perturb() {
        // Cell selected by changing by hand until test passes,
        // since this is just a regression test
        // Needs to be different from the one below
        let mut rng = Xoshiro256StarStar::seed_from_u64(0);
        let cell = [219.0232, -4200.21983, -2300.123123, 8443.5, 2213.4, -231.15, 5200.2, -2002.1, 231.24];
        let mut sphere: Sphere = Particle::parse("0.0 0.3 0.22 1.0", &cell);
        println!("{}", sphere);
        sphere.perturb(&cell, &[0.0], &mut rng);
        println!("{}", sphere);
        let local = global_to_relative3(&cell, &sphere.global_pos);
        assert!(local[0] < 0.0);
        assert!(sphere.rel_pos[0] >= 0.0 && sphere.rel_pos[0] < 1.0);
        assert!(sphere.rel_pos[1] >= 0.0 && sphere.rel_pos[1] < 1.0);
        assert!(sphere.rel_pos[2] >= 0.0 && sphere.rel_pos[2] < 1.0);
        let local_pbc = apply_pbc(&local);
        assert!(local_pbc[0] >= 0.0 && local_pbc[0] < 1.0);
        assert!(local_pbc[1] >= 0.0 && local_pbc[1] < 1.0);
        assert!(local_pbc[2] >= 0.0 && local_pbc[2] < 1.0);
    }
    
    #[test]
    fn zero_inversion_strain() {
        // Cell selected by changing by hand until test passes,
        // since this is just a regression test
        let mut rng = Xoshiro256StarStar::seed_from_u64(0);
        let old_cell = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let cell = [2100.0232, -2000.21983, -230.0, 833.5, 2.4, -231.15, 5.2, -200.1, 21.24];
        let mut sphere: Sphere = Particle::parse("0.0 0.3 0.22 1.0", &old_cell);
        sphere.apply_strain(&cell);
        let local = global_to_relative3(&cell, &sphere.global_pos);
        assert!(local[0] < 0.0);
        assert!(sphere.rel_pos[0] >= 0.0 && sphere.rel_pos[0] < 1.0);
        assert!(sphere.rel_pos[1] >= 0.0 && sphere.rel_pos[1] < 1.0);
        assert!(sphere.rel_pos[2] >= 0.0 && sphere.rel_pos[2] < 1.0);
        let local_pbc = apply_pbc(&local);
        println!("{:?}", local_pbc);
        assert!(local_pbc[0] >= 0.0 && local_pbc[0] < 1.0);
        assert!(local_pbc[1] >= 0.0 && local_pbc[1] < 1.0);
        assert!(local_pbc[2] >= 0.0 && local_pbc[2] < 1.0);
    }
    
    #[test]
    fn internal_symmetric_invalidation2() {
        // This test shows why symmetric_invalidation2 in cell_list.rs
        // was added in the first place, as it caught a bug in my first implementation
        // aimed at fixing symmetric_invalidation1
        use std::convert::TryInto;
        use crate::overbox_list::OverboxList;
        use crate::cell_list::CellList;
        let unit_cell = vec![5.07019429309184, -0.47515173025982416, 0.18018763060534648, -0.5656991499933727, 3.0853092570916774, -0.2935140522456714, 0.07612697602838228, -0.021482269817008064, 3.9939399781072584];
        let forward_offset: [f64; 3] = [-unit_cell[6], -unit_cell[7], -unit_cell[8]];
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
        let p1_scale = p1.squared_scaling(p2, &forward_offset);
        let p2_scale = p2.squared_scaling(p1, &negate(&forward_offset));
        println!("{} {}", p1_scale, p2_scale);
        assert!(p1_scale <= 1.0);
        assert!(p2_scale > 1.0);
        assert!(p2_scale*p1_scale > 1.0);
    }

}
