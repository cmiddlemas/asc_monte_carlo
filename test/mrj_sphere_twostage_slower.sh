#!/bin/bash

cargo build --release

export RAYON_NUM_THREADS=2
export RUST_BACKTRACE=full

TEST_DIR=mrj_sphere_twostage_slower

# https://www.cyberciti.biz/faq/howto-check-if-a-directory-exists-in-a-bash-shellscript/
# https://stackoverflow.com/questions/13553173/whats-the-meaning-of-the-operator-in-linux-shell
# Also ubuntu default crontab
[ -d "$TEST_DIR" ] || mkdir "$TEST_DIR"

(
RUN_NAME=final_1e-4_1
[-d "$TEST_DIR/$RUN_NAME" ] || mkdir $TEST_DIR/$RUN_NAME

cargo run --release -- --no-rayon -d 3 -n 200 -s 16 --cell-list --pressure inf --pcell 0.005 --moves 4000 --sweeps 50000 --adjust --gap-threshold 0.003 --logfiles $TEST_DIR/$RUN_NAME/initial_log -o $TEST_DIR/$RUN_NAME/initial_config 2>$TEST_DIR/$RUN_NAME/initial.stderr >$TEST_DIR/$RUN_NAME/initial.stdout

cargo run --release -- --no-rayon --initfile $TEST_DIR/$RUN_NAME/initial_config_final.dat --schedulefile $TEST_DIR/$RUN_NAME/initial_log_schedule.yaml --cell-list --pressure inf --pcell 0.0000005 --moves 40000000 --sweeps 5000 --adjust --logfiles $TEST_DIR/$RUN_NAME/final_log -o $TEST_DIR/$RUN_NAME/final_config 2>$TEST_DIR/$RUN_NAME/final.stderr >$TEST_DIR/$RUN_NAME/final.stdout
) &

(
RUN_NAME=final_1e-4_2
[-d "$TEST_DIR/$RUN_NAME" ] || mkdir $TEST_DIR/$RUN_NAME

cargo run --release -- --no-rayon -d 3 -n 200 -s 16 --cell-list --pressure inf --pcell 0.005 --moves 4000 --sweeps 50000 --adjust --gap-threshold 0.003 --logfiles $TEST_DIR/$RUN_NAME/initial_log -o $TEST_DIR/$RUN_NAME/initial_config 2>$TEST_DIR/$RUN_NAME/initial.stderr >$TEST_DIR/$RUN_NAME/initial.stdout

cargo run --release -- --no-rayon --initfile $TEST_DIR/$RUN_NAME/initial_config_final.dat --schedulefile $TEST_DIR/$RUN_NAME/initial_log_schedule.yaml --cell-list --pressure inf --pcell 0.0000005 --moves 40000000 --sweeps 5000 --adjust --logfiles $TEST_DIR/$RUN_NAME/final_log -o $TEST_DIR/$RUN_NAME/final_config 2>$TEST_DIR/$RUN_NAME/final.stderr >$TEST_DIR/$RUN_NAME/final.stdout
) &

(
RUN_NAME=final_1e-4_3
[-d "$TEST_DIR/$RUN_NAME" ] || mkdir $TEST_DIR/$RUN_NAME

cargo run --release -- --no-rayon -d 3 -n 200 -s 16 --cell-list --pressure inf --pcell 0.005 --moves 4000 --sweeps 50000 --adjust --gap-threshold 0.003 --logfiles $TEST_DIR/$RUN_NAME/initial_log -o $TEST_DIR/$RUN_NAME/initial_config 2>$TEST_DIR/$RUN_NAME/initial.stderr >$TEST_DIR/$RUN_NAME/initial.stdout

cargo run --release -- --no-rayon --initfile $TEST_DIR/$RUN_NAME/initial_config_final.dat --schedulefile $TEST_DIR/$RUN_NAME/initial_log_schedule.yaml --cell-list --pressure inf --pcell 0.0000005 --moves 40000000 --sweeps 5000 --adjust --logfiles $TEST_DIR/$RUN_NAME/final_log -o $TEST_DIR/$RUN_NAME/final_config 2>$TEST_DIR/$RUN_NAME/final.stderr >$TEST_DIR/$RUN_NAME/final.stdout
) &

(
RUN_NAME=final_1e-4_4
[-d "$TEST_DIR/$RUN_NAME" ] || mkdir $TEST_DIR/$RUN_NAME

cargo run --release -- --no-rayon -d 3 -n 200 -s 16 --cell-list --pressure inf --pcell 0.005 --moves 4000 --sweeps 50000 --adjust --gap-threshold 0.003 --logfiles $TEST_DIR/$RUN_NAME/initial_log -o $TEST_DIR/$RUN_NAME/initial_config 2>$TEST_DIR/$RUN_NAME/initial.stderr >$TEST_DIR/$RUN_NAME/initial.stdout

cargo run --release -- --no-rayon --initfile $TEST_DIR/$RUN_NAME/initial_config_final.dat --schedulefile $TEST_DIR/$RUN_NAME/initial_log_schedule.yaml --cell-list --pressure inf --pcell 0.0000005 --moves 40000000 --sweeps 5000 --adjust --logfiles $TEST_DIR/$RUN_NAME/final_log -o $TEST_DIR/$RUN_NAME/final_config 2>$TEST_DIR/$RUN_NAME/final.stderr >$TEST_DIR/$RUN_NAME/final.stdout
) &

(
RUN_NAME=final_1e-4_5
[-d "$TEST_DIR/$RUN_NAME" ] || mkdir $TEST_DIR/$RUN_NAME

cargo run --release -- --no-rayon -d 3 -n 200 -s 16 --cell-list --pressure inf --pcell 0.005 --moves 4000 --sweeps 50000 --adjust --gap-threshold 0.003 --logfiles $TEST_DIR/$RUN_NAME/initial_log -o $TEST_DIR/$RUN_NAME/initial_config 2>$TEST_DIR/$RUN_NAME/initial.stderr >$TEST_DIR/$RUN_NAME/initial.stdout

cargo run --release -- --no-rayon --initfile $TEST_DIR/$RUN_NAME/initial_config_final.dat --schedulefile $TEST_DIR/$RUN_NAME/initial_log_schedule.yaml --cell-list --pressure inf --pcell 0.0000005 --moves 40000000 --sweeps 5000 --adjust --logfiles $TEST_DIR/$RUN_NAME/final_log -o $TEST_DIR/$RUN_NAME/final_config 2>$TEST_DIR/$RUN_NAME/final.stderr >$TEST_DIR/$RUN_NAME/final.stdout
) &
