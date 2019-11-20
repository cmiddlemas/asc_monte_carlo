build :
	USING_MAKE=true cargo build --release

install :
	USING_MAKE=true cargo install --locked --force --path . 
