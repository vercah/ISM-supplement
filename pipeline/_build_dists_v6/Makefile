.PHONY: all clean

all: _build_dists_v6

_build_dists_v6:
	RUSTFLAGS="-C target-cpu=native" cargo build --release
	ln -sf target/release/build_dists_v6 _build_dists_v6

clean:
	cargo clean
	rm -f _build_dists_v6
