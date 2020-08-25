#!/bin/bash

cargo build --release

export RAYON_NUM_THREADS=3
export RUST_BACKTRACE=full

TEST_DIR=eq_sphere_meltingpoint2

# https://www.cyberciti.biz/faq/howto-check-if-a-directory-exists-in-a-bash-shellscript/
# https://stackoverflow.com/questions/13553173/whats-the-meaning-of-the-operator-in-linux-shell
# Also ubuntu default crontab
[ -d "$TEST_DIR" ] || mkdir "$TEST_DIR"

# See previous testing scripts for comments

(
# p = 1.02
[ -d "$TEST_DIR/true_1_02" ] || mkdir $TEST_DIR/true_1_02
cargo run --release -- --initfile sphere_fcc.dat --log-volume-step --no-clamp --cell-list --pressure 1.02 --pcell 0.002 --moves 50000 --sweeps 5000 --adjust --logfiles $TEST_DIR/true_1_02/equilibration -o $TEST_DIR/true_1_02/equilibration 2>$TEST_DIR/true_1_02/eq_std_err.txt >$TEST_DIR/true_1_02/eq_std_out.txt
cargo run --release -- --save-trajectory --log-volume-step --no-clamp --cell-list --initfile $TEST_DIR/true_1_02/equilibration_final.dat --sweeps 1000 --schedulefile $TEST_DIR/true_1_02/equilibration_schedule.json --logfiles $TEST_DIR/true_1_02/measure -o $TEST_DIR/true_1_02/measure 2>$TEST_DIR/true_1_02/me_std_err.txt >$TEST_DIR/true_1_02/me_std_out.txt
) &
(
# p = 1.04
[ -d "$TEST_DIR/true_1_04" ] || mkdir $TEST_DIR/true_1_04
cargo run --release -- --initfile sphere_fcc.dat --log-volume-step --no-clamp --cell-list --pressure 1.04 --pcell 0.004 --moves 50000 --sweeps 5000 --adjust --logfiles $TEST_DIR/true_1_04/equilibration -o $TEST_DIR/true_1_04/equilibration 2>$TEST_DIR/true_1_04/eq_std_err.txt >$TEST_DIR/true_1_04/eq_std_out.txt
cargo run --release -- --save-trajectory --log-volume-step --no-clamp --cell-list --initfile $TEST_DIR/true_1_04/equilibration_final.dat --sweeps 1000 --schedulefile $TEST_DIR/true_1_04/equilibration_schedule.json --logfiles $TEST_DIR/true_1_04/measure -o $TEST_DIR/true_1_04/measure 2>$TEST_DIR/true_1_04/me_std_err.txt >$TEST_DIR/true_1_04/me_std_out.txt
) &
(
# p = 1.06
[ -d "$TEST_DIR/true_1_06" ] || mkdir $TEST_DIR/true_1_06
cargo run --release -- --initfile sphere_fcc.dat --log-volume-step --no-clamp --cell-list --pressure 1.06 --pcell 0.006 --moves 50000 --sweeps 5000 --adjust --logfiles $TEST_DIR/true_1_06/equilibration -o $TEST_DIR/true_1_06/equilibration 2>$TEST_DIR/true_1_06/eq_std_err.txt >$TEST_DIR/true_1_06/eq_std_out.txt
cargo run --release -- --save-trajectory --log-volume-step --no-clamp --cell-list --initfile $TEST_DIR/true_1_06/equilibration_final.dat --sweeps 1000 --schedulefile $TEST_DIR/true_1_06/equilibration_schedule.json --logfiles $TEST_DIR/true_1_06/measure -o $TEST_DIR/true_1_06/measure 2>$TEST_DIR/true_1_06/me_std_err.txt >$TEST_DIR/true_1_06/me_std_out.txt
) &
(
# p = 1.08
[ -d "$TEST_DIR/true_1_08" ] || mkdir $TEST_DIR/true_1_08
cargo run --release -- --initfile sphere_fcc.dat --log-volume-step --no-clamp --cell-list --pressure 1.08 --pcell 0.008 --moves 50000 --sweeps 5000 --adjust --logfiles $TEST_DIR/true_1_08/equilibration -o $TEST_DIR/true_1_08/equilibration 2>$TEST_DIR/true_1_08/eq_std_err.txt >$TEST_DIR/true_1_08/eq_std_out.txt
cargo run --release -- --save-trajectory --log-volume-step --no-clamp --cell-list --initfile $TEST_DIR/true_1_08/equilibration_final.dat --sweeps 1000 --schedulefile $TEST_DIR/true_1_08/equilibration_schedule.json --logfiles $TEST_DIR/true_1_08/measure -o $TEST_DIR/true_1_08/measure 2>$TEST_DIR/true_1_08/me_std_err.txt >$TEST_DIR/true_1_08/me_std_out.txt
) &
