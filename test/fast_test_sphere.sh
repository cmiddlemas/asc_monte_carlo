#!/bin/bash

cargo build --release

export RAYON_NUM_THREADS=3
export RUST_BACKTRACE=full

TEST_DIR=fast_test_sphere

# https://www.cyberciti.biz/faq/howto-check-if-a-directory-exists-in-a-bash-shellscript/
# https://stackoverflow.com/questions/13553173/whats-the-meaning-of-the-operator-in-linux-shell
# Also ubuntu default crontab
[ -d "$TEST_DIR" ] || mkdir "$TEST_DIR"

# Going to run at a variety of compression rates

(
cargo run --release -- --adjust --log-volume-step --n-bins-fit 5 --no-clamp -v --cell-list -n 500 -s 22.0 --pressure 0.5 --pcell 0.002 --moves 50000 --sweeps 1000 -o $TEST_DIR/config --logfiles $TEST_DIR/log 2>$TEST_DIR/stderr.txt >$TEST_DIR/stdout.txt
) &

(
cargo run --release -- --adjust --n-bins-fit 5 --cell-list -n 500 -s 22.0 --pressure inf --pcell 0.002 --moves 50000 --sweeps 5000 -o $TEST_DIR/mrj_config --logfiles $TEST_DIR/mrj_log 2>$TEST_DIR/mrj_stderr.txt >$TEST_DIR/mrj_stdout.txt
) &
