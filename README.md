# HartreeFockRustProgram

This is a reimplementation of the HartreeFockPythonProgram in Rust from
[NickelAndCopper's](https://youtube.com/playlist?list=PL-hN8vfRaU7jSFHN1ZSAMNe_2nXhwAmzM&si=ANjI8kPn-5v_3Kvs)
YouTube Playlist (also, here is the
[GitHub](https://github.com/nickelandcopper/HartreeFockPythonProgram)
for the HatreeFockPythonProgram).

Currently the code is set up with the STO-3G basis set for Hydrogen,
based on the values available from the [Basis Set
Exhange](https://www.basissetexchange.org/basis/6-31g/format/json/?version=1&elements=1)

STO-3G Reference: [Ref](https://www.basissetexchange.org/references/sto-6g/format/txt/?version=1&elements=1)  

## Build Instructions

Dependencies:
- Rust toolchain:
  - cargo
  - rustc
- gsl - GNU Scientific Library
- OpenBLAS


Clone the repo:

``` bash
git clone https://github.com/jlheflin/HartreeFockRustProgram.git
```

### Rust Build

Run cargo build:

``` bash
cd ./HatreeFockRustProgram
cargo build
```

Run the program:

``` bash
cargo run
# OR
./target/debug/hf_rust
```

## Build Result

The output should be the following:

``` bash
Total Energy: -1.0659994615565611
```

# References
Pritchard, Benjamin P., Doaa Altarawy, Brett Didier, Tara D. Gibson, and
Theresa L. Windus. 2019. <span>“New Basis Set Exchange: An Open,
up-to-Date Resource for the Molecular Sciences Community.”</span>
<em>Journal of Chemical Information and Modeling</em> 59 (11): 4814–20.
<a
href="https://doi.org/10.1021/acs.jcim.9b00725">https://doi.org/10.1021/acs.jcim.9b00725</a>.
