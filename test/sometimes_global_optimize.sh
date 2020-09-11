# This script sometimes produces something near the global optimum, at least on
# a program close but not exactly commit 589f305bfdfbd594d41cfc351e73850b9d68862f
# I also used it to choose the symmetric_invalidation() test case.
# Not sure this was the exact command, but it is close
cargo run --release -- --adjust --cell-list --pcell 0.1 -n 10 -s 10 --moves 4000 --sweeps 1000 --pressure inf
