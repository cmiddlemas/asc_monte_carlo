#!/bin/bash

cargo build --release

export RAYON_NUM_THREADS=3
export RUST_BACKTRACE=full

TEST_DIR=test_eq_sphere_melt

# https://www.cyberciti.biz/faq/howto-check-if-a-directory-exists-in-a-bash-shellscript/
# https://stackoverflow.com/questions/13553173/whats-the-meaning-of-the-operator-in-linux-shell
# Also ubuntu default crontab
[ -d "$TEST_DIR" ] || mkdir "$TEST_DIR"

(
[ -d "$TEST_DIR/true_0_5" ] || mkdir $TEST_DIR/true_0_5
# Set up true equilibrium system, using auto adjust to get in the vicinity
# Then do a measurement run, with fixed step sizes determined by last sweep
cargo run --release -- --initfile sphere_fcc.dat --log-volume-step --no-clamp --cell-list --pressure 0.5 --pcell 0.002 --moves 50000 --sweeps 5000 --adjust --logfiles $TEST_DIR/true_0_5/equilibration -o $TEST_DIR/true_0_5/equilibration 2>$TEST_DIR/true_0_5/eq_std_err.txt >$TEST_DIR/true_0_5/eq_std_out.txt
cargo run --release -- --save-trajectory --log-volume-step --no-clamp --cell-list --initfile $TEST_DIR/true_0_5/equilibration_final.dat --sweeps 1000 --schedulefile $TEST_DIR/true_0_5/equilibration_schedule.json --logfiles $TEST_DIR/true_0_5/measure -o $TEST_DIR/true_0_5/measure 2>$TEST_DIR/true_0_5/me_std_err.txt >$TEST_DIR/true_0_5/me_std_out.txt

# Now do the same, but use the Torquato ASC method, which involves a linearization of deformations

[ -d "$TEST_DIR/linear_0_5" ] || mkdir $TEST_DIR/linear_0_5
# Set up linear equilibrium system, using auto adjust to get in the vicinity
# Then do a measurement run, with fixed step sizes determined by last sweep
cargo run --release -- --initfile sphere_fcc.dat --cell-list --pressure 0.5 --pcell 0.002 --moves 50000 --sweeps 5000 --adjust --logfiles $TEST_DIR/linear_0_5/equilibration -o $TEST_DIR/linear_0_5/equilibration 2>$TEST_DIR/linear_0_5/eq_std_err.txt >$TEST_DIR/linear_0_5/eq_std_out.txt
cargo run --release -- --save-trajectory --cell-list --initfile $TEST_DIR/linear_0_5/equilibration_final.dat --sweeps 1000 --schedulefile $TEST_DIR/linear_0_5/equilibration_schedule.json --logfiles $TEST_DIR/linear_0_5/measure -o $TEST_DIR/linear_0_5/measure 2>$TEST_DIR/linear_0_5/me_std_err.txt >$TEST_DIR/linear_0_5/me_std_out.txt

# Do one final set, which uses the linear steps, but turns off shear and axial deformations
[ -d "$TEST_DIR/isolinear_0_5" ] || mkdir $TEST_DIR/isolinear_0_5
# Set up isolinear equilibrium system, using auto adjust to get in the vicinity
# Then do a measurement run, with fixed step sizes determined by last sweep
cargo run --release -- --initfile sphere_fcc.dat --axial 0.0 --shear 0.0 --cell-list --pressure 0.5 --pcell 0.002 --moves 50000 --sweeps 5000 --adjust --logfiles $TEST_DIR/isolinear_0_5/equilibration -o $TEST_DIR/isolinear_0_5/equilibration 2>$TEST_DIR/isolinear_0_5/eq_std_err.txt >$TEST_DIR/isolinear_0_5/eq_std_out.txt
cargo run --release -- --save-trajectory --cell-list --initfile $TEST_DIR/isolinear_0_5/equilibration_final.dat --sweeps 1000 --schedulefile $TEST_DIR/isolinear_0_5/equilibration_schedule.json --logfiles $TEST_DIR/isolinear_0_5/measure -o $TEST_DIR/isolinear_0_5/measure 2>$TEST_DIR/isolinear_0_5/me_std_err.txt >$TEST_DIR/isolinear_0_5/me_std_out.txt
) &

# Do the same, for 3 more pressures. Chosen with help of Erpenbeck and Wood (1984)

(
# p = 1.0
[ -d "$TEST_DIR/true_1" ] || mkdir $TEST_DIR/true_1
cargo run --release -- --initfile sphere_fcc.dat --log-volume-step --no-clamp --cell-list --pressure 1 --pcell 0.002 --moves 50000 --sweeps 5000 --adjust --logfiles $TEST_DIR/true_1/equilibration -o $TEST_DIR/true_1/equilibration 2>$TEST_DIR/true_1/eq_std_err.txt >$TEST_DIR/true_1/eq_std_out.txt
cargo run --release -- --save-trajectory --log-volume-step --no-clamp --cell-list --initfile $TEST_DIR/true_1/equilibration_final.dat --sweeps 1000 --schedulefile $TEST_DIR/true_1/equilibration_schedule.json --logfiles $TEST_DIR/true_1/measure -o $TEST_DIR/true_1/measure 2>$TEST_DIR/true_1/me_std_err.txt >$TEST_DIR/true_1/me_std_out.txt

[ -d "$TEST_DIR/linear_1" ] || mkdir $TEST_DIR/linear_1
cargo run --release -- --initfile sphere_fcc.dat --cell-list --pressure 1 --pcell 0.002 --moves 50000 --sweeps 5000 --adjust --logfiles $TEST_DIR/linear_1/equilibration -o $TEST_DIR/linear_1/equilibration 2>$TEST_DIR/linear_1/eq_std_err.txt >$TEST_DIR/linear_1/eq_std_out.txt
cargo run --release -- --save-trajectory --cell-list --initfile $TEST_DIR/linear_1/equilibration_final.dat --sweeps 1000 --schedulefile $TEST_DIR/linear_1/equilibration_schedule.json --logfiles $TEST_DIR/linear_1/measure -o $TEST_DIR/linear_1/measure 2>$TEST_DIR/linear_1/me_std_err.txt >$TEST_DIR/linear_1/me_std_out.txt

[ -d "$TEST_DIR/isolinear_1" ] || mkdir $TEST_DIR/isolinear_1
cargo run --release -- --initfile sphere_fcc.dat --axial 0.0 --shear 0.0 --cell-list --pressure 1 --pcell 0.002 --moves 50000 --sweeps 5000 --adjust --logfiles $TEST_DIR/isolinear_1/equilibration -o $TEST_DIR/isolinear_1/equilibration 2>$TEST_DIR/isolinear_1/eq_std_err.txt >$TEST_DIR/isolinear_1/eq_std_out.txt
cargo run --release -- --save-trajectory --cell-list --initfile $TEST_DIR/isolinear_1/equilibration_final.dat --sweeps 1000 --schedulefile $TEST_DIR/isolinear_1/equilibration_schedule.json --logfiles $TEST_DIR/isolinear_1/measure -o $TEST_DIR/isolinear_1/measure 2>$TEST_DIR/isolinear_1/me_std_err.txt >$TEST_DIR/isolinear_1/me_std_out.txt
) &

(
# p = 1.5
[ -d "$TEST_DIR/true_1_5" ] || mkdir $TEST_DIR/true_1_5
cargo run --release -- --initfile sphere_fcc.dat --log-volume-step --no-clamp --cell-list --pressure 1.5 --pcell 0.002 --moves 50000 --sweeps 5000 --adjust --logfiles $TEST_DIR/true_1_5/equilibration -o $TEST_DIR/true_1_5/equilibration 2>$TEST_DIR/true_1_5/eq_std_err.txt >$TEST_DIR/true_1_5/eq_std_out.txt
cargo run --release -- --save-trajectory --log-volume-step --no-clamp --cell-list --initfile $TEST_DIR/true_1_5/equilibration_final.dat --sweeps 1000 --schedulefile $TEST_DIR/true_1_5/equilibration_schedule.json --logfiles $TEST_DIR/true_1_5/measure -o $TEST_DIR/true_1_5/measure 2>$TEST_DIR/true_1_5/me_std_err.txt >$TEST_DIR/true_1_5/me_std_out.txt

[ -d "$TEST_DIR/linear_1_5" ] || mkdir $TEST_DIR/linear_1_5
cargo run --release -- --initfile sphere_fcc.dat --cell-list --pressure 1.5 --pcell 0.002 --moves 50000 --sweeps 5000 --adjust --logfiles $TEST_DIR/linear_1_5/equilibration -o $TEST_DIR/linear_1_5/equilibration 2>$TEST_DIR/linear_1_5/eq_std_err.txt >$TEST_DIR/linear_1_5/eq_std_out.txt
cargo run --release -- --save-trajectory --cell-list --initfile $TEST_DIR/linear_1_5/equilibration_final.dat --sweeps 1000 --schedulefile $TEST_DIR/linear_1_5/equilibration_schedule.json --logfiles $TEST_DIR/linear_1_5/measure -o $TEST_DIR/linear_1_5/measure 2>$TEST_DIR/linear_1_5/me_std_err.txt >$TEST_DIR/linear_1_5/me_std_out.txt

[ -d "$TEST_DIR/isolinear_1_5" ] || mkdir $TEST_DIR/isolinear_1_5
cargo run --release -- --initfile sphere_fcc.dat --axial 0.0 --shear 0.0 --cell-list --pressure 1.5 --pcell 0.002 --moves 50000 --sweeps 5000 --adjust --logfiles $TEST_DIR/isolinear_1_5/equilibration -o $TEST_DIR/isolinear_1_5/equilibration 2>$TEST_DIR/isolinear_1_5/eq_std_err.txt >$TEST_DIR/isolinear_1_5/eq_std_out.txt
cargo run --release -- --save-trajectory --cell-list --initfile $TEST_DIR/isolinear_1_5/equilibration_final.dat --sweeps 1000 --schedulefile $TEST_DIR/isolinear_1_5/equilibration_schedule.json --logfiles $TEST_DIR/isolinear_1_5/measure -o $TEST_DIR/isolinear_1_5/measure 2>$TEST_DIR/isolinear_1_5/me_std_err.txt >$TEST_DIR/isolinear_1_5/me_std_out.txt
) &

(
# p = 2.0
[ -d "$TEST_DIR/true_2" ] || mkdir $TEST_DIR/true_2
cargo run --release -- --initfile sphere_fcc.dat --log-volume-step --no-clamp --cell-list --pressure 2 --pcell 0.002 --moves 50000 --sweeps 5000 --adjust --logfiles $TEST_DIR/true_2/equilibration -o $TEST_DIR/true_2/equilibration 2>$TEST_DIR/true_2/eq_std_err.txt >$TEST_DIR/true_2/eq_std_out.txt
cargo run --release -- --save-trajectory --log-volume-step --no-clamp --cell-list --initfile $TEST_DIR/true_2/equilibration_final.dat --sweeps 1000 --schedulefile $TEST_DIR/true_2/equilibration_schedule.json --logfiles $TEST_DIR/true_2/measure -o $TEST_DIR/true_2/measure 2>$TEST_DIR/true_2/me_std_err.txt >$TEST_DIR/true_2/me_std_out.txt

[ -d "$TEST_DIR/linear_2" ] || mkdir $TEST_DIR/linear_2
cargo run --release -- --initfile sphere_fcc.dat --cell-list --pressure 2 --pcell 0.002 --moves 50000 --sweeps 5000 --adjust --logfiles $TEST_DIR/linear_2/equilibration -o $TEST_DIR/linear_2/equilibration 2>$TEST_DIR/linear_2/eq_std_err.txt >$TEST_DIR/linear_2/eq_std_out.txt
cargo run --release -- --save-trajectory --cell-list --initfile $TEST_DIR/linear_2/equilibration_final.dat --sweeps 1000 --schedulefile $TEST_DIR/linear_2/equilibration_schedule.json --logfiles $TEST_DIR/linear_2/measure -o $TEST_DIR/linear_2/measure 2>$TEST_DIR/linear_2/me_std_err.txt >$TEST_DIR/linear_2/me_std_out.txt

[ -d "$TEST_DIR/isolinear_2" ] || mkdir $TEST_DIR/isolinear_2
cargo run --release -- --initfile sphere_fcc.dat --axial 0.0 --shear 0.0 --cell-list --pressure 2 --pcell 0.002 --moves 50000 --sweeps 5000 --adjust --logfiles $TEST_DIR/isolinear_2/equilibration -o $TEST_DIR/isolinear_2/equilibration 2>$TEST_DIR/isolinear_2/eq_std_err.txt >$TEST_DIR/isolinear_2/eq_std_out.txt
cargo run --release -- --save-trajectory --cell-list --initfile $TEST_DIR/isolinear_2/equilibration_final.dat --sweeps 1000 --schedulefile $TEST_DIR/isolinear_2/equilibration_schedule.json --logfiles $TEST_DIR/isolinear_2/measure -o $TEST_DIR/isolinear_2/measure 2>$TEST_DIR/isolinear_2/me_std_err.txt >$TEST_DIR/isolinear_2/me_std_out.txt
) &
