#!/bin/bash

cargo build --release

export RAYON_NUM_THREADS=3
export RUST_BACKTRACE=full

TEST_DIR=eq_sphere_meltingpoint

# https://www.cyberciti.biz/faq/howto-check-if-a-directory-exists-in-a-bash-shellscript/
# https://stackoverflow.com/questions/13553173/whats-the-meaning-of-the-operator-in-linux-shell
# Also ubuntu default crontab
[ -d "$TEST_DIR" ] || mkdir "$TEST_DIR"

# See previous testing scripts for comments

(
# p = 1.1
[ -d "$TEST_DIR/true_1_1" ] || mkdir $TEST_DIR/true_1_1
cargo run --release -- --initfile sphere_fcc.dat --log-volume-step --no-clamp --cell-list --pressure 1.1 --pcell 0.002 --moves 50000 --sweeps 5000 --adjust --logfiles $TEST_DIR/true_1_1/equilibration -o $TEST_DIR/true_1_1/equilibration 2>$TEST_DIR/true_1_1/eq_std_err.txt >$TEST_DIR/true_1_1/eq_std_out.txt
cargo run --release -- --save-trajectory --log-volume-step --no-clamp --cell-list --initfile $TEST_DIR/true_1_1/equilibration_final.dat --sweeps 1000 --schedulefile $TEST_DIR/true_1_1/equilibration_schedule.json --logfiles $TEST_DIR/true_1_1/measure -o $TEST_DIR/true_1_1/measure 2>$TEST_DIR/true_1_1/me_std_err.txt >$TEST_DIR/true_1_1/me_std_out.txt
) &

(
# p = 1.2
[ -d "$TEST_DIR/true_1_2" ] || mkdir $TEST_DIR/true_1_2
cargo run --release -- --initfile sphere_fcc.dat --log-volume-step --no-clamp --cell-list --pressure 1.2 --pcell 0.002 --moves 50000 --sweeps 5000 --adjust --logfiles $TEST_DIR/true_1_2/equilibration -o $TEST_DIR/true_1_2/equilibration 2>$TEST_DIR/true_1_2/eq_std_err.txt >$TEST_DIR/true_1_2/eq_std_out.txt
cargo run --release -- --save-trajectory --log-volume-step --no-clamp --cell-list --initfile $TEST_DIR/true_1_2/equilibration_final.dat --sweeps 1000 --schedulefile $TEST_DIR/true_1_2/equilibration_schedule.json --logfiles $TEST_DIR/true_1_2/measure -o $TEST_DIR/true_1_2/measure 2>$TEST_DIR/true_1_2/me_std_err.txt >$TEST_DIR/true_1_2/me_std_out.txt
) &

(
# p = 1.3
[ -d "$TEST_DIR/true_1_3" ] || mkdir $TEST_DIR/true_1_3
cargo run --release -- --initfile sphere_fcc.dat --log-volume-step --no-clamp --cell-list --pressure 1.3 --pcell 0.002 --moves 50000 --sweeps 5000 --adjust --logfiles $TEST_DIR/true_1_3/equilibration -o $TEST_DIR/true_1_3/equilibration 2>$TEST_DIR/true_1_3/eq_std_err.txt >$TEST_DIR/true_1_3/eq_std_out.txt
cargo run --release -- --save-trajectory --log-volume-step --no-clamp --cell-list --initfile $TEST_DIR/true_1_3/equilibration_final.dat --sweeps 1000 --schedulefile $TEST_DIR/true_1_3/equilibration_schedule.json --logfiles $TEST_DIR/true_1_3/measure -o $TEST_DIR/true_1_3/measure 2>$TEST_DIR/true_1_3/me_std_err.txt >$TEST_DIR/true_1_3/me_std_out.txt
) &

(
# p = 1.4
[ -d "$TEST_DIR/true_1_4" ] || mkdir $TEST_DIR/true_1_4
cargo run --release -- --initfile sphere_fcc.dat --log-volume-step --no-clamp --cell-list --pressure 1.4 --pcell 0.002 --moves 50000 --sweeps 5000 --adjust --logfiles $TEST_DIR/true_1_4/equilibration -o $TEST_DIR/true_1_4/equilibration 2>$TEST_DIR/true_1_4/eq_std_err.txt >$TEST_DIR/true_1_4/eq_std_out.txt
cargo run --release -- --save-trajectory --log-volume-step --no-clamp --cell-list --initfile $TEST_DIR/true_1_4/equilibration_final.dat --sweeps 1000 --schedulefile $TEST_DIR/true_1_4/equilibration_schedule.json --logfiles $TEST_DIR/true_1_4/measure -o $TEST_DIR/true_1_4/measure 2>$TEST_DIR/true_1_4/me_std_err.txt >$TEST_DIR/true_1_4/me_std_out.txt
) &
