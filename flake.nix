{
  description = "Rust project with BLAS/LAPACK using Nix flakes";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
    rust-overlay.url = "github:oxalica/rust-overlay";
  };

  outputs = { self, nixpkgs, flake-utils, rust-overlay }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs {
          inherit system;
          overlays = [ rust-overlay.overlays.default ];
        };

        rust = pkgs.rust-bin.stable.latest.default;

      in {
        # ⛏ Build Rust project
        packages.default = pkgs.rustPlatform.buildRustPackage {
          pname = "hf_rust";
          version = "0.1.0";

          src = ./.;

          cargoLock = {
            lockFile = ./Cargo.lock;
          };

          nativeBuildInputs = [
            pkgs.pkg-config
            pkgs.gcc
          ];

          buildInputs = [
            pkgs.openblas
            pkgs.lapack
          ];

          # Optional for linking
          LD_LIBRARY_PATH = "${pkgs.openblas}/lib:${pkgs.lapack}/lib";

          # Optional: helps linking with BLAS in crates like ndarray-linalg
          RUSTFLAGS = "-L${pkgs.openblas}/lib -L${pkgs.lapack}/lib";
        };

        # 🧪 Development shell
        devShell = pkgs.mkShell {
          buildInputs = [
            rust
            pkgs.pkg-config
            pkgs.openblas
            pkgs.lapack
            pkgs.gcc
          ];

          shellHook = ''
            export RUST_BACKTRACE=1
            export LAPACK_LIB=openblas
            export BLAS_LIB=openblas
          '';
        };
      });
}
