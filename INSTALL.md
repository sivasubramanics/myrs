# Installation & Building

## Prerequisites (For Cross-Compiling to Linux)

If you are on macOS and want to compile static binaries for Linux, install `zig` and `cargo-zigbuild`:

```bash
brew install zig
cargo install cargo-zigbuild
rustup target add x86_64-unknown-linux-musl
```

---

## Building

### Native Build (macOS Apple Silicon / Intel)

```bash
cargo build --release
```
The compiled binary will be located at:
`./target/release/myrs`

---

### Cross-Compiling for Linux (`x86_64` Static Binary)

Produces a fully static binary that runs on any Linux distribution without `glibc` dependency issues:
```bash
cargo zigbuild --release --target x86_64-unknown-linux-musl
```
The compiled binary will be located at:
`./target/x86_64-unknown-linux-musl/release/myrs`
