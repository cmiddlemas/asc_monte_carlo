#!/bin/bash

cargo build --release

export RAYON_NUM_THREADS=2
export RUST_BACKTRACE=full

TEST_DIR=mrj_sphere_twostage_quick

# https://www.cyberciti.biz/faq/howto-check-if-a-directory-exists-in-a-bash-shellscript/
# https://stackoverflow.com/questions/13553173/whats-the-meaning-of-the-operator-in-linux-shell
# Also ubuntu default crontab
[ -d "$TEST_DIR" ] || mkdir "$TEST_DIR"

# Going to run at an initial compression rate of 1e-0,
# up to some threshold in nn_gap, then finish compressing at
# 1e-1. 1e-2, and 1e-3. Run 5 each.

# https://stackoverflow.com/questions/22768533/d-command-not-found
(
RUN_NAME=final_1e-1_1
[ -d "$TEST_DIR/$RUN_NAME" ] || mkdir $TEST_DIR/$RUN_NAME

cargo run --release -- -d 3 -n 200 -s 16 --cell-list --pressure 1e100 --pcell 0.005 --moves 40000 --sweeps 50000 --adjust --gap-threshold 1e-8 --logfiles $TEST_DIR/$RUN_NAME/initial_log -o $TEST_DIR/$RUN_NAME/initial_config 2>$TEST_DIR/$RUN_NAME/initial.stderr >$TEST_DIR/$RUN_NAME/initial.stdout

cargo run --release -- --initfile $TEST_DIR/$RUN_NAME/initial_config_final.dat --schedulefile $TEST_DIR/$RUN_NAME/initial_log_schedule.json --cell-list --pressure 1e100 --pcell 0.0005 --moves 40000 --sweeps 5000 --adjust --logfiles $TEST_DIR/$RUN_NAME/final_log -o $TEST_DIR/$RUN_NAME/final_config 2>$TEST_DIR/$RUN_NAME/final.stderr >$TEST_DIR/$RUN_NAME/final.stdout
) &

(
RUN_NAME=final_1e-1_2
[ -d "$TEST_DIR/$RUN_NAME" ] || mkdir $TEST_DIR/$RUN_NAME

cargo run --release -- -d 3 -n 200 -s 16 --cell-list --pressure 1e100 --pcell 0.005 --moves 40000 --sweeps 50000 --adjust --gap-threshold 1e-8 --logfiles $TEST_DIR/$RUN_NAME/initial_log -o $TEST_DIR/$RUN_NAME/initial_config 2>$TEST_DIR/$RUN_NAME/initial.stderr >$TEST_DIR/$RUN_NAME/initial.stdout

cargo run --release -- --initfile $TEST_DIR/$RUN_NAME/initial_config_final.dat --schedulefile $TEST_DIR/$RUN_NAME/initial_log_schedule.json --cell-list --pressure 1e100 --pcell 0.0005 --moves 40000 --sweeps 5000 --adjust --logfiles $TEST_DIR/$RUN_NAME/final_log -o $TEST_DIR/$RUN_NAME/final_config 2>$TEST_DIR/$RUN_NAME/final.stderr >$TEST_DIR/$RUN_NAME/final.stdout
) &

(
RUN_NAME=final_1e-1_3
[ -d "$TEST_DIR/$RUN_NAME" ] || mkdir $TEST_DIR/$RUN_NAME

cargo run --release -- -d 3 -n 200 -s 16 --cell-list --pressure 1e100 --pcell 0.005 --moves 40000 --sweeps 50000 --adjust --gap-threshold 1e-8 --logfiles $TEST_DIR/$RUN_NAME/initial_log -o $TEST_DIR/$RUN_NAME/initial_config 2>$TEST_DIR/$RUN_NAME/initial.stderr >$TEST_DIR/$RUN_NAME/initial.stdout

cargo run --release -- --initfile $TEST_DIR/$RUN_NAME/initial_config_final.dat --schedulefile $TEST_DIR/$RUN_NAME/initial_log_schedule.json --cell-list --pressure 1e100 --pcell 0.0005 --moves 40000 --sweeps 5000 --adjust --logfiles $TEST_DIR/$RUN_NAME/final_log -o $TEST_DIR/$RUN_NAME/final_config 2>$TEST_DIR/$RUN_NAME/final.stderr >$TEST_DIR/$RUN_NAME/final.stdout
) &

(
RUN_NAME=final_1e-1_4
[ -d "$TEST_DIR/$RUN_NAME" ] || mkdir $TEST_DIR/$RUN_NAME

cargo run --release -- -d 3 -n 200 -s 16 --cell-list --pressure 1e100 --pcell 0.005 --moves 40000 --sweeps 50000 --adjust --gap-threshold 1e-8 --logfiles $TEST_DIR/$RUN_NAME/initial_log -o $TEST_DIR/$RUN_NAME/initial_config 2>$TEST_DIR/$RUN_NAME/initial.stderr >$TEST_DIR/$RUN_NAME/initial.stdout

cargo run --release -- --initfile $TEST_DIR/$RUN_NAME/initial_config_final.dat --schedulefile $TEST_DIR/$RUN_NAME/initial_log_schedule.json --cell-list --pressure 1e100 --pcell 0.0005 --moves 40000 --sweeps 5000 --adjust --logfiles $TEST_DIR/$RUN_NAME/final_log -o $TEST_DIR/$RUN_NAME/final_config 2>$TEST_DIR/$RUN_NAME/final.stderr >$TEST_DIR/$RUN_NAME/final.stdout
) &

(
RUN_NAME=final_1e-1_5
[ -d "$TEST_DIR/$RUN_NAME" ] || mkdir $TEST_DIR/$RUN_NAME

cargo run --release -- -d 3 -n 200 -s 16 --cell-list --pressure 1e100 --pcell 0.005 --moves 40000 --sweeps 50000 --adjust --gap-threshold 1e-8 --logfiles $TEST_DIR/$RUN_NAME/initial_log -o $TEST_DIR/$RUN_NAME/initial_config 2>$TEST_DIR/$RUN_NAME/initial.stderr >$TEST_DIR/$RUN_NAME/initial.stdout

cargo run --release -- --initfile $TEST_DIR/$RUN_NAME/initial_config_final.dat --schedulefile $TEST_DIR/$RUN_NAME/initial_log_schedule.json --cell-list --pressure 1e100 --pcell 0.0005 --moves 40000 --sweeps 5000 --adjust --logfiles $TEST_DIR/$RUN_NAME/final_log -o $TEST_DIR/$RUN_NAME/final_config 2>$TEST_DIR/$RUN_NAME/final.stderr >$TEST_DIR/$RUN_NAME/final.stdout
) &
