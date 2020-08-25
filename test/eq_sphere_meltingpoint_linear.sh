#!/bin/bash

cargo build --release

export RAYON_NUM_THREADS=3
export RUST_BACKTRACE=full

TEST_DIR=eq_sphere_meltingpoint_linear

# https://www.cyberciti.biz/faq/howto-check-if-a-directory-exists-in-a-bash-shellscript/
# https://stackoverflow.com/questions/13553173/whats-the-meaning-of-the-operator-in-linux-shell
# Also ubuntu default crontab
[ -d "$TEST_DIR" ] || mkdir "$TEST_DIR"

# See previous testing scripts for comments

(
# p = 1.12
[ -d "$TEST_DIR/true_1_12" ] || mkdir $TEST_DIR/true_1_12
cargo run --release -- --initfile sphere_fcc.dat --cell-list --pressure 1.12 --pcell 0.002 --moves 50000 --sweeps 25000 --adjust --logfiles $TEST_DIR/true_1_12/equilibration -o $TEST_DIR/true_1_12/equilibration 2>$TEST_DIR/true_1_12/eq_std_err.txt >$TEST_DIR/true_1_12/eq_std_out.txt
cargo run --release -- --save-trajectory --cell-list --initfile $TEST_DIR/true_1_12/equilibration_final.dat --sweeps 25000 --schedulefile $TEST_DIR/true_1_12/equilibration_schedule.json --logfiles $TEST_DIR/true_1_12/measure -o $TEST_DIR/true_1_12/measure 2>$TEST_DIR/true_1_12/me_std_err.txt >$TEST_DIR/true_1_12/me_std_out.txt
) &
(
# p = 1.14
[ -d "$TEST_DIR/true_1_14" ] || mkdir $TEST_DIR/true_1_14
cargo run --release -- --initfile sphere_fcc.dat --cell-list --pressure 1.14 --pcell 0.004 --moves 50000 --sweeps 25000 --adjust --logfiles $TEST_DIR/true_1_14/equilibration -o $TEST_DIR/true_1_14/equilibration 2>$TEST_DIR/true_1_14/eq_std_err.txt >$TEST_DIR/true_1_14/eq_std_out.txt
cargo run --release -- --save-trajectory --cell-list --initfile $TEST_DIR/true_1_14/equilibration_final.dat --sweeps 25000 --schedulefile $TEST_DIR/true_1_14/equilibration_schedule.json --logfiles $TEST_DIR/true_1_14/measure -o $TEST_DIR/true_1_14/measure 2>$TEST_DIR/true_1_14/me_std_err.txt >$TEST_DIR/true_1_14/me_std_out.txt
) &
(
# p = 1.16
[ -d "$TEST_DIR/true_1_16" ] || mkdir $TEST_DIR/true_1_16
cargo run --release -- --initfile sphere_fcc.dat --cell-list --pressure 1.16 --pcell 0.006 --moves 50000 --sweeps 25000 --adjust --logfiles $TEST_DIR/true_1_16/equilibration -o $TEST_DIR/true_1_16/equilibration 2>$TEST_DIR/true_1_16/eq_std_err.txt >$TEST_DIR/true_1_16/eq_std_out.txt
cargo run --release -- --save-trajectory --cell-list --initfile $TEST_DIR/true_1_16/equilibration_final.dat --sweeps 25000 --schedulefile $TEST_DIR/true_1_16/equilibration_schedule.json --logfiles $TEST_DIR/true_1_16/measure -o $TEST_DIR/true_1_16/measure 2>$TEST_DIR/true_1_16/me_std_err.txt >$TEST_DIR/true_1_16/me_std_out.txt
) &
(
# p = 1.18
[ -d "$TEST_DIR/true_1_18" ] || mkdir $TEST_DIR/true_1_18
cargo run --release -- --initfile sphere_fcc.dat --cell-list --pressure 1.18 --pcell 0.008 --moves 50000 --sweeps 25000 --adjust --logfiles $TEST_DIR/true_1_18/equilibration -o $TEST_DIR/true_1_18/equilibration 2>$TEST_DIR/true_1_18/eq_std_err.txt >$TEST_DIR/true_1_18/eq_std_out.txt
cargo run --release -- --save-trajectory --cell-list --initfile $TEST_DIR/true_1_18/equilibration_final.dat --sweeps 25000 --schedulefile $TEST_DIR/true_1_18/equilibration_schedule.json --logfiles $TEST_DIR/true_1_18/measure -o $TEST_DIR/true_1_18/measure 2>$TEST_DIR/true_1_18/me_std_err.txt >$TEST_DIR/true_1_18/me_std_out.txt
) &
(
# p = 1.2
[ -d "$TEST_DIR/true_1_2" ] || mkdir $TEST_DIR/true_1_2
cargo run --release -- --initfile sphere_fcc.dat --cell-list --pressure 1.2 --pcell 0.008 --moves 50000 --sweeps 25000 --adjust --logfiles $TEST_DIR/true_1_2/equilibration -o $TEST_DIR/true_1_2/equilibration 2>$TEST_DIR/true_1_2/eq_std_err.txt >$TEST_DIR/true_1_2/eq_std_out.txt
cargo run --release -- --save-trajectory --cell-list --initfile $TEST_DIR/true_1_2/equilibration_final.dat --sweeps 25000 --schedulefile $TEST_DIR/true_1_2/equilibration_schedule.json --logfiles $TEST_DIR/true_1_2/measure -o $TEST_DIR/true_1_2/measure 2>$TEST_DIR/true_1_2/me_std_err.txt >$TEST_DIR/true_1_2/me_std_out.txt
) &
(
# p = 1.22
[ -d "$TEST_DIR/true_1_22" ] || mkdir $TEST_DIR/true_1_22
cargo run --release -- --initfile sphere_fcc.dat --cell-list --pressure 1.22 --pcell 0.008 --moves 50000 --sweeps 25000 --adjust --logfiles $TEST_DIR/true_1_22/equilibration -o $TEST_DIR/true_1_22/equilibration 2>$TEST_DIR/true_1_22/eq_std_err.txt >$TEST_DIR/true_1_22/eq_std_out.txt
cargo run --release -- --save-trajectory --cell-list --initfile $TEST_DIR/true_1_22/equilibration_final.dat --sweeps 25000 --schedulefile $TEST_DIR/true_1_22/equilibration_schedule.json --logfiles $TEST_DIR/true_1_22/measure -o $TEST_DIR/true_1_22/measure 2>$TEST_DIR/true_1_22/me_std_err.txt >$TEST_DIR/true_1_22/me_std_out.txt
) &
(
# p = 1.24
[ -d "$TEST_DIR/true_1_24" ] || mkdir $TEST_DIR/true_1_24
cargo run --release -- --initfile sphere_fcc.dat --cell-list --pressure 1.24 --pcell 0.008 --moves 50000 --sweeps 25000 --adjust --logfiles $TEST_DIR/true_1_24/equilibration -o $TEST_DIR/true_1_24/equilibration 2>$TEST_DIR/true_1_24/eq_std_err.txt >$TEST_DIR/true_1_24/eq_std_out.txt
cargo run --release -- --save-trajectory --cell-list --initfile $TEST_DIR/true_1_24/equilibration_final.dat --sweeps 25000 --schedulefile $TEST_DIR/true_1_24/equilibration_schedule.json --logfiles $TEST_DIR/true_1_24/measure -o $TEST_DIR/true_1_24/measure 2>$TEST_DIR/true_1_24/me_std_err.txt >$TEST_DIR/true_1_24/me_std_out.txt
) &
