# Because of
# https://github.com/rust-lang/cargo/issues/7169
# this is a hack to get a proxy for whether
# --locked was run on a cargo install
# Also ensures we run cargo clean first, to
# make sure all files are fresh and environment
# constants get built in properly
build :
	cargo clean; cargo build --features using_make --release

install :
	cargo clean; cargo install --features using_make --locked --force --path .
