#!/bin/bash

cargo build --release

export RAYON_NUM_THREADS=3
export RUST_BACKTRACE=full

TEST_DIR=mrj_sphere

# https://www.cyberciti.biz/faq/howto-check-if-a-directory-exists-in-a-bash-shellscript/
# https://stackoverflow.com/questions/13553173/whats-the-meaning-of-the-operator-in-linux-shell
# Also ubuntu default crontab
[ -d "$TEST_DIR" ] || mkdir "$TEST_DIR"

# Going to run at a variety of compression rates

(
[ -d "$TEST_DIR/compress_1e--1" ] || mkdir $TEST_DIR/compress_1e--1
cargo run --release -- -d 3 -n 500 -s 28 --cell-list --pressure inf --pcell 0.02 --moves 1000 --sweeps 50000 --adjust --logfiles $TEST_DIR/compress_1e--1/run1 -o $TEST_DIR/compress_1e--1/run1 2>$TEST_DIR/compress_1e--1/run1_std_err.txt >$TEST_DIR/compress_1e--1/run1_std_out.txt
cargo run --release -- -d 3 -n 500 -s 28 --cell-list --pressure inf --pcell 0.02 --moves 1000 --sweeps 50000 --adjust --logfiles $TEST_DIR/compress_1e--1/run2 -o $TEST_DIR/compress_1e--1/run2 2>$TEST_DIR/compress_1e--1/run2.txt >$TEST_DIR/compress_1e--1/run2.txt
cargo run --release -- -d 3 -n 500 -s 28 --cell-list --pressure inf --pcell 0.02 --moves 1000 --sweeps 50000 --adjust --logfiles $TEST_DIR/compress_1e--1/run3 -o $TEST_DIR/compress_1e--1/run3 2>$TEST_DIR/compress_1e--1/run3.txt >$TEST_DIR/compress_1e--1/run3.txt
cargo run --release -- -d 3 -n 500 -s 28 --cell-list --pressure inf --pcell 0.02 --moves 1000 --sweeps 50000 --adjust --logfiles $TEST_DIR/compress_1e--1/run4 -o $TEST_DIR/compress_1e--1/run4 2>$TEST_DIR/compress_1e--1/run4.txt >$TEST_DIR/compress_1e--1/run4.txt
cargo run --release -- -d 3 -n 500 -s 28 --cell-list --pressure inf --pcell 0.02 --moves 1000 --sweeps 50000 --adjust --logfiles $TEST_DIR/compress_1e--1/run5 -o $TEST_DIR/compress_1e--1/run5 2>$TEST_DIR/compress_1e--1/run5.txt >$TEST_DIR/compress_1e--1/run5.txt
) &

(
[ -d "$TEST_DIR/compress_1e--2" ] || mkdir $TEST_DIR/compress_1e--2
cargo run --release -- -d 3 -n 500 -s 28 --cell-list --pressure inf --pcell 0.2 --moves 100 --sweeps 50000 --adjust --logfiles $TEST_DIR/compress_1e--2/run1 -o $TEST_DIR/compress_1e--2/run1 2>$TEST_DIR/compress_1e--2/run1_std_err.txt >$TEST_DIR/compress_1e--2/run1_std_out.txt
cargo run --release -- -d 3 -n 500 -s 28 --cell-list --pressure inf --pcell 0.2 --moves 100 --sweeps 50000 --adjust --logfiles $TEST_DIR/compress_1e--2/run2 -o $TEST_DIR/compress_1e--2/run2 2>$TEST_DIR/compress_1e--2/run2.txt >$TEST_DIR/compress_1e--2/run2.txt
cargo run --release -- -d 3 -n 500 -s 28 --cell-list --pressure inf --pcell 0.2 --moves 100 --sweeps 50000 --adjust --logfiles $TEST_DIR/compress_1e--2/run3 -o $TEST_DIR/compress_1e--2/run3 2>$TEST_DIR/compress_1e--2/run3.txt >$TEST_DIR/compress_1e--2/run3.txt
cargo run --release -- -d 3 -n 500 -s 28 --cell-list --pressure inf --pcell 0.2 --moves 100 --sweeps 50000 --adjust --logfiles $TEST_DIR/compress_1e--2/run4 -o $TEST_DIR/compress_1e--2/run4 2>$TEST_DIR/compress_1e--2/run4.txt >$TEST_DIR/compress_1e--2/run4.txt
cargo run --release -- -d 3 -n 500 -s 28 --cell-list --pressure inf --pcell 0.2 --moves 100 --sweeps 50000 --adjust --logfiles $TEST_DIR/compress_1e--2/run5 -o $TEST_DIR/compress_1e--2/run5 2>$TEST_DIR/compress_1e--2/run5.txt >$TEST_DIR/compress_1e--2/run5.txt
) &
