#!/bin/bash

cargo build --release

export RAYON_NUM_THREADS=3
export RUST_BACKTRACE=full

# https://www.cyberciti.biz/faq/howto-check-if-a-directory-exists-in-a-bash-shellscript/
# https://stackoverflow.com/questions/13553173/whats-the-meaning-of-the-operator-in-linux-shell
# Also ubuntu default crontab
[ -d "test_eq_sphere" ] || mkdir test_eq_sphere

(
[ -d "test_eq_sphere/true_0_5" ] || mkdir test_eq_sphere/true_0_5
# Set up true equilibrium system, using auto adjust to get in the vicinity
# Then do a measurement run, with fixed step sizes determined by last sweep
cargo run --release -- --log-volume-step --no-clamp -d 3 -n 500 --cell-list -s 20.0 --pressure 0.5 --pcell 0.002 --moves 50000 --sweeps 5000 --adjust --logfiles test_eq_sphere/true_0_5/equilibration -o test_eq_sphere/true_0_5/equilibration 2>test_eq_sphere/true_0_5/eq_std_err.txt >test_eq_sphere/true_0_5/eq_std_out.txt
cargo run --release -- --log-volume-step --no-clamp --cell-list --initfile test_eq_sphere/true_0_5/equilibration_final.dat --sweeps 1000 --schedulefile test_eq_sphere/true_0_5/equilibration_schedule.json --logfiles test_eq_sphere/true_0_5/measure -o test_eq_sphere/true_0_5/measure 2>test_eq_sphere/true_0_5/me_std_err.txt >test_eq_sphere/true_0_5/me_std_out.txt

# Now do the same, but use the Torquato ASC method, which involves a linearization of deformations

[ -d "test_eq_sphere/linear_0_5" ] || mkdir test_eq_sphere/linear_0_5
# Set up linear equilibrium system, using auto adjust to get in the vicinity
# Then do a measurement run, with fixed step sizes determined by last sweep
cargo run --release -- -d 3 -n 500 --cell-list -s 20.0 --pressure 0.5 --pcell 0.002 --moves 50000 --sweeps 5000 --adjust --logfiles test_eq_sphere/linear_0_5/equilibration -o test_eq_sphere/linear_0_5/equilibration 2>test_eq_sphere/linear_0_5/eq_std_err.txt >test_eq_sphere/linear_0_5/eq_std_out.txt
cargo run --release -- --cell-list --initfile test_eq_sphere/linear_0_5/equilibration_final.dat --sweeps 1000 --schedulefile test_eq_sphere/linear_0_5/equilibration_schedule.json --logfiles test_eq_sphere/linear_0_5/measure -o test_eq_sphere/linear_0_5/measure 2>test_eq_sphere/linear_0_5/me_std_err.txt >test_eq_sphere/linear_0_5/me_std_out.txt

# Do one final set, which uses the linear steps, but turns off shear and axial deformations
[ -d "test_eq_sphere/isolinear_0_5" ] || mkdir test_eq_sphere/isolinear_0_5
# Set up isolinear equilibrium system, using auto adjust to get in the vicinity
# Then do a measurement run, with fixed step sizes determined by last sweep
cargo run --release -- --axial 0.0 --shear 0.0 -d 3 -n 500 --cell-list -s 20.0 --pressure 0.5 --pcell 0.002 --moves 50000 --sweeps 5000 --adjust --logfiles test_eq_sphere/isolinear_0_5/equilibration -o test_eq_sphere/isolinear_0_5/equilibration 2>test_eq_sphere/isolinear_0_5/eq_std_err.txt >test_eq_sphere/isolinear_0_5/eq_std_out.txt
cargo run --release -- --cell-list --initfile test_eq_sphere/isolinear_0_5/equilibration_final.dat --sweeps 1000 --schedulefile test_eq_sphere/isolinear_0_5/equilibration_schedule.json --logfiles test_eq_sphere/isolinear_0_5/measure -o test_eq_sphere/isolinear_0_5/measure 2>test_eq_sphere/isolinear_0_5/me_std_err.txt >test_eq_sphere/isolinear_0_5/me_std_out.txt
) &

# Do the same, for 3 more pressures. Chosen with help of Erpenbeck and Wood (1984)

(
# p = 1.0
[ -d "test_eq_sphere/true_1" ] || mkdir test_eq_sphere/true_1
cargo run --release -- --log-volume-step --no-clamp -d 3 -n 500 --cell-list -s 20.0 --pressure 1 --pcell 0.002 --moves 50000 --sweeps 5000 --adjust --logfiles test_eq_sphere/true_1/equilibration -o test_eq_sphere/true_1/equilibration 2>test_eq_sphere/true_1/eq_std_err.txt >test_eq_sphere/true_1/eq_std_out.txt
cargo run --release -- --log-volume-step --no-clamp --cell-list --initfile test_eq_sphere/true_1/equilibration_final.dat --sweeps 1000 --schedulefile test_eq_sphere/true_1/equilibration_schedule.json --logfiles test_eq_sphere/true_1/measure -o test_eq_sphere/true_1/measure 2>test_eq_sphere/true_1/me_std_err.txt >test_eq_sphere/true_1/me_std_out.txt

[ -d "test_eq_sphere/linear_1" ] || mkdir test_eq_sphere/linear_1
cargo run --release -- -d 3 -n 500 --cell-list -s 20.0 --pressure 1 --pcell 0.002 --moves 50000 --sweeps 5000 --adjust --logfiles test_eq_sphere/linear_1/equilibration -o test_eq_sphere/linear_1/equilibration 2>test_eq_sphere/linear_1/eq_std_err.txt >test_eq_sphere/linear_1/eq_std_out.txt
cargo run --release -- --cell-list --initfile test_eq_sphere/linear_1/equilibration_final.dat --sweeps 1000 --schedulefile test_eq_sphere/linear_1/equilibration_schedule.json --logfiles test_eq_sphere/linear_1/measure -o test_eq_sphere/linear_1/measure 2>test_eq_sphere/linear_1/me_std_err.txt >test_eq_sphere/linear_1/me_std_out.txt

[ -d "test_eq_sphere/isolinear_1" ] || mkdir test_eq_sphere/isolinear_1
cargo run --release -- --axial 0.0 --shear 0.0 -d 3 -n 500 --cell-list -s 20.0 --pressure 1 --pcell 0.002 --moves 50000 --sweeps 5000 --adjust --logfiles test_eq_sphere/isolinear_1/equilibration -o test_eq_sphere/isolinear_1/equilibration 2>test_eq_sphere/isolinear_1/eq_std_err.txt >test_eq_sphere/isolinear_1/eq_std_out.txt
cargo run --release -- --cell-list --initfile test_eq_sphere/isolinear_1/equilibration_final.dat --sweeps 1000 --schedulefile test_eq_sphere/isolinear_1/equilibration_schedule.json --logfiles test_eq_sphere/isolinear_1/measure -o test_eq_sphere/isolinear_1/measure 2>test_eq_sphere/isolinear_1/me_std_err.txt >test_eq_sphere/isolinear_1/me_std_out.txt
) &

(
# p = 1.5
[ -d "test_eq_sphere/true_1_5" ] || mkdir test_eq_sphere/true_1_5
cargo run --release -- --log-volume-step --no-clamp -d 3 -n 500 --cell-list -s 20.0 --pressure 1.5 --pcell 0.002 --moves 50000 --sweeps 5000 --adjust --logfiles test_eq_sphere/true_1_5/equilibration -o test_eq_sphere/true_1_5/equilibration 2>test_eq_sphere/true_1_5/eq_std_err.txt >test_eq_sphere/true_1_5/eq_std_out.txt
cargo run --release -- --log-volume-step --no-clamp --cell-list --initfile test_eq_sphere/true_1_5/equilibration_final.dat --sweeps 1000 --schedulefile test_eq_sphere/true_1_5/equilibration_schedule.json --logfiles test_eq_sphere/true_1_5/measure -o test_eq_sphere/true_1_5/measure 2>test_eq_sphere/true_1_5/me_std_err.txt >test_eq_sphere/true_1_5/me_std_out.txt

[ -d "test_eq_sphere/linear_1_5" ] || mkdir test_eq_sphere/linear_1_5
cargo run --release -- -d 3 -n 500 --cell-list -s 20.0 --pressure 1.5 --pcell 0.002 --moves 50000 --sweeps 5000 --adjust --logfiles test_eq_sphere/linear_1_5/equilibration -o test_eq_sphere/linear_1_5/equilibration 2>test_eq_sphere/linear_1_5/eq_std_err.txt >test_eq_sphere/linear_1_5/eq_std_out.txt
cargo run --release -- --cell-list --initfile test_eq_sphere/linear_1_5/equilibration_final.dat --sweeps 1000 --schedulefile test_eq_sphere/linear_1_5/equilibration_schedule.json --logfiles test_eq_sphere/linear_1_5/measure -o test_eq_sphere/linear_1_5/measure 2>test_eq_sphere/linear_1_5/me_std_err.txt >test_eq_sphere/linear_1_5/me_std_out.txt

[ -d "test_eq_sphere/isolinear_1_5" ] || mkdir test_eq_sphere/isolinear_1_5
cargo run --release -- --axial 0.0 --shear 0.0 -d 3 -n 500 --cell-list -s 20.0 --pressure 1.5 --pcell 0.002 --moves 50000 --sweeps 5000 --adjust --logfiles test_eq_sphere/isolinear_1_5/equilibration -o test_eq_sphere/isolinear_1_5/equilibration 2>test_eq_sphere/isolinear_1_5/eq_std_err.txt >test_eq_sphere/isolinear_1_5/eq_std_out.txt
cargo run --release -- --cell-list --initfile test_eq_sphere/isolinear_1_5/equilibration_final.dat --sweeps 1000 --schedulefile test_eq_sphere/isolinear_1_5/equilibration_schedule.json --logfiles test_eq_sphere/isolinear_1_5/measure -o test_eq_sphere/isolinear_1_5/measure 2>test_eq_sphere/isolinear_1_5/me_std_err.txt >test_eq_sphere/isolinear_1_5/me_std_out.txt
) &

(
# p = 2.0
[ -d "test_eq_sphere/true_2" ] || mkdir test_eq_sphere/true_2
cargo run --release -- --log-volume-step --no-clamp -d 3 -n 500 --cell-list -s 20.0 --pressure 2 --pcell 0.002 --moves 50000 --sweeps 5000 --adjust --logfiles test_eq_sphere/true_2/equilibration -o test_eq_sphere/true_2/equilibration 2>test_eq_sphere/true_2/eq_std_err.txt >test_eq_sphere/true_2/eq_std_out.txt
cargo run --release -- --log-volume-step --no-clamp --cell-list --initfile test_eq_sphere/true_2/equilibration_final.dat --sweeps 1000 --schedulefile test_eq_sphere/true_2/equilibration_schedule.json --logfiles test_eq_sphere/true_2/measure -o test_eq_sphere/true_2/measure 2>test_eq_sphere/true_2/me_std_err.txt >test_eq_sphere/true_2/me_std_out.txt

[ -d "test_eq_sphere/linear_2" ] || mkdir test_eq_sphere/linear_2
cargo run --release -- -d 3 -n 500 --cell-list -s 20.0 --pressure 2 --pcell 0.002 --moves 50000 --sweeps 5000 --adjust --logfiles test_eq_sphere/linear_2/equilibration -o test_eq_sphere/linear_2/equilibration 2>test_eq_sphere/linear_2/eq_std_err.txt >test_eq_sphere/linear_2/eq_std_out.txt
cargo run --release -- --cell-list --initfile test_eq_sphere/linear_2/equilibration_final.dat --sweeps 1000 --schedulefile test_eq_sphere/linear_2/equilibration_schedule.json --logfiles test_eq_sphere/linear_2/measure -o test_eq_sphere/linear_2/measure 2>test_eq_sphere/linear_2/me_std_err.txt >test_eq_sphere/linear_2/me_std_out.txt

[ -d "test_eq_sphere/isolinear_2" ] || mkdir test_eq_sphere/isolinear_2
cargo run --release -- --axial 0.0 --shear 0.0 -d 3 -n 500 --cell-list -s 20.0 --pressure 2 --pcell 0.002 --moves 50000 --sweeps 5000 --adjust --logfiles test_eq_sphere/isolinear_2/equilibration -o test_eq_sphere/isolinear_2/equilibration 2>test_eq_sphere/isolinear_2/eq_std_err.txt >test_eq_sphere/isolinear_2/eq_std_out.txt
cargo run --release -- --cell-list --initfile test_eq_sphere/isolinear_2/equilibration_final.dat --sweeps 1000 --schedulefile test_eq_sphere/isolinear_2/equilibration_schedule.json --logfiles test_eq_sphere/isolinear_2/measure -o test_eq_sphere/isolinear_2/measure 2>test_eq_sphere/isolinear_2/me_std_err.txt >test_eq_sphere/isolinear_2/me_std_out.txt
) &
