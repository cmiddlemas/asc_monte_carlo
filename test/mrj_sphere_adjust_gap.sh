#!/bin/bash

cargo build --release

export RUST_BACKTRACE=full

TEST_DIR=mrj_sphere_adjust_gap

# https://www.cyberciti.biz/faq/howto-check-if-a-directory-exists-in-a-bash-shellscript/
# https://stackoverflow.com/questions/13553173/whats-the-meaning-of-the-operator-in-linux-shell
# Also ubuntu default crontab
[ -d "$TEST_DIR" ] || mkdir "$TEST_DIR"

# Going to run at a variety of compression rates

(
[ -d "$TEST_DIR/compress_1" ] || mkdir $TEST_DIR/compress_1
cargo run --release -- --no-rayon -d 3 -n 500 -s 25 --cell-list --pressure inf --pcell 0.002 --axial 0.0 --shear 0.0 --moves 50000 --sweeps 5000 --adjust --adjust-gap 1.0 --uniform-moves --logfiles $TEST_DIR/compress_1/run1 -o $TEST_DIR/compress_1/run1 2>$TEST_DIR/compress_1/run1_std_err.txt >$TEST_DIR/compress_1/run1_std_out.txt

cargo run --release -- --no-rayon -d 3 -n 500 -s 25 --cell-list --pressure inf --pcell 0.002 --axial 0.0 --shear 0.0 --moves 50000 --sweeps 5000 --adjust --adjust-gap 1.0 --uniform-moves --logfiles $TEST_DIR/compress_1/run2 -o $TEST_DIR/compress_1/run2 2>$TEST_DIR/compress_1/run2_std_err.txt >$TEST_DIR/compress_1/run2_std_out.txt

cargo run --release -- --no-rayon -d 3 -n 500 -s 25 --cell-list --pressure inf --pcell 0.002 --axial 0.0 --shear 0.0 --moves 50000 --sweeps 5000 --adjust --adjust-gap 1.0 --uniform-moves --logfiles $TEST_DIR/compress_1/run3 -o $TEST_DIR/compress_1/run3 2>$TEST_DIR/compress_1/run3_std_err.txt >$TEST_DIR/compress_1/run3_std_out.txt

cargo run --release -- --no-rayon -d 3 -n 500 -s 25 --cell-list --pressure inf --pcell 0.002 --axial 0.0 --shear 0.0 --moves 50000 --sweeps 5000 --adjust --adjust-gap 1.0 --uniform-moves --logfiles $TEST_DIR/compress_1/run4 -o $TEST_DIR/compress_1/run4 2>$TEST_DIR/compress_1/run4_std_err.txt >$TEST_DIR/compress_1/run4_std_out.txt

cargo run --release -- --no-rayon -d 3 -n 500 -s 25 --cell-list --pressure inf --pcell 0.002 --axial 0.0 --shear 0.0 --moves 50000 --sweeps 5000 --adjust --adjust-gap 1.0 --uniform-moves --logfiles $TEST_DIR/compress_1/run5 -o $TEST_DIR/compress_1/run5 2>$TEST_DIR/compress_1/run5_std_err.txt >$TEST_DIR/compress_1/run5_std_out.txt
) &

(
[ -d "$TEST_DIR/compress_01" ] || mkdir $TEST_DIR/compress_01
cargo run --release -- --no-rayon -d 3 -n 500 -s 25 --cell-list --pressure inf --pcell 0.002 --axial 0.0 --shear 0.0 --moves 50000 --sweeps 5000 --adjust --adjust-gap 0.1 --uniform-moves --logfiles $TEST_DIR/compress_01/run1 -o $TEST_DIR/compress_01/run1 2>$TEST_DIR/compress_01/run1_std_err.txt >$TEST_DIR/compress_01/run1_std_out.txt

cargo run --release -- --no-rayon -d 3 -n 500 -s 25 --cell-list --pressure inf --pcell 0.002 --axial 0.0 --shear 0.0 --moves 50000 --sweeps 5000 --adjust --adjust-gap 0.1 --uniform-moves --logfiles $TEST_DIR/compress_01/run2 -o $TEST_DIR/compress_01/run2 2>$TEST_DIR/compress_01/run2_std_err.txt >$TEST_DIR/compress_01/run2_std_out.txt

cargo run --release -- --no-rayon -d 3 -n 500 -s 25 --cell-list --pressure inf --pcell 0.002 --axial 0.0 --shear 0.0 --moves 50000 --sweeps 5000 --adjust --adjust-gap 0.1 --uniform-moves --logfiles $TEST_DIR/compress_01/run3 -o $TEST_DIR/compress_01/run3 2>$TEST_DIR/compress_01/run3_std_err.txt >$TEST_DIR/compress_01/run3_std_out.txt

cargo run --release -- --no-rayon -d 3 -n 500 -s 25 --cell-list --pressure inf --pcell 0.002 --axial 0.0 --shear 0.0 --moves 50000 --sweeps 5000 --adjust --adjust-gap 0.1 --uniform-moves --logfiles $TEST_DIR/compress_01/run4 -o $TEST_DIR/compress_01/run4 2>$TEST_DIR/compress_01/run4_std_err.txt >$TEST_DIR/compress_01/run4_std_out.txt

cargo run --release -- --no-rayon -d 3 -n 500 -s 25 --cell-list --pressure inf --pcell 0.002 --axial 0.0 --shear 0.0 --moves 50000 --sweeps 5000 --adjust --adjust-gap 0.1 --uniform-moves --logfiles $TEST_DIR/compress_01/run5 -o $TEST_DIR/compress_01/run5 2>$TEST_DIR/compress_01/run5_std_err.txt >$TEST_DIR/compress_01/run5_std_out.txt
) &

(
[ -d "$TEST_DIR/compress_001" ] || mkdir $TEST_DIR/compress_001
cargo run --release -- --no-rayon -d 3 -n 500 -s 25 --cell-list --pressure inf --pcell 0.002 --axial 0.0 --shear 0.0 --moves 50000 --sweeps 5000 --adjust --adjust-gap 0.01 --uniform-moves --logfiles $TEST_DIR/compress_001/run1 -o $TEST_DIR/compress_001/run1 2>$TEST_DIR/compress_001/run1_std_err.txt >$TEST_DIR/compress_001/run1_std_out.txt

cargo run --release -- --no-rayon -d 3 -n 500 -s 25 --cell-list --pressure inf --pcell 0.002 --axial 0.0 --shear 0.0 --moves 50000 --sweeps 5000 --adjust --adjust-gap 0.01 --uniform-moves --logfiles $TEST_DIR/compress_001/run2 -o $TEST_DIR/compress_001/run2 2>$TEST_DIR/compress_001/run2_std_err.txt >$TEST_DIR/compress_001/run2_std_out.txt

cargo run --release -- --no-rayon -d 3 -n 500 -s 25 --cell-list --pressure inf --pcell 0.002 --axial 0.0 --shear 0.0 --moves 50000 --sweeps 5000 --adjust --adjust-gap 0.01 --uniform-moves --logfiles $TEST_DIR/compress_001/run3 -o $TEST_DIR/compress_001/run3 2>$TEST_DIR/compress_001/run3_std_err.txt >$TEST_DIR/compress_001/run3_std_out.txt

cargo run --release -- --no-rayon -d 3 -n 500 -s 25 --cell-list --pressure inf --pcell 0.002 --axial 0.0 --shear 0.0 --moves 50000 --sweeps 5000 --adjust --adjust-gap 0.01 --uniform-moves --logfiles $TEST_DIR/compress_001/run4 -o $TEST_DIR/compress_001/run4 2>$TEST_DIR/compress_001/run4_std_err.txt >$TEST_DIR/compress_001/run4_std_out.txt

cargo run --release -- --no-rayon -d 3 -n 500 -s 25 --cell-list --pressure inf --pcell 0.002 --axial 0.0 --shear 0.0 --moves 50000 --sweeps 5000 --adjust --adjust-gap 0.01 --uniform-moves --logfiles $TEST_DIR/compress_001/run5 -o $TEST_DIR/compress_001/run5 2>$TEST_DIR/compress_001/run5_std_err.txt >$TEST_DIR/compress_001/run5_std_out.txt
) &
