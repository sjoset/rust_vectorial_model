# Development:
With nix flakes enabled, use ```nix develop``` to enter development shell

## Testing
Use ```cargo test``` to run available code tests

## Building
cargo build --release

## Running
cargo run --release --bin [rust_vect_oneshot, total_fragments_vs_rh] [program arguments]

or run the binaries available after building in target/release

# Usage
## rust_vect_oneshot
The main vectorial model binary.

### Command line
Takes two command line arguments:
#### Input config file
a yaml file containing the model configuration
An example configuration yaml file is provided in `example_input/parameters_example.yaml`.
#### Output file
the filename to use for model output.

## total_fragments_vs_rh
Binary for running a given model configuration over a range of heliocentric distances.

### Command line
#### Input config file
a yaml file containing the model configuration
An example configuration yaml file is provided in `example_input/parameters_example.yaml`.
#### Output file
the filename to use for model output.
#### rh_start
Starting point of heliocentric grid, in AU
#### rh_end
Ending point of heliocentric grid, in AU
#### rh_gridsize
How many heliocentric distances to run along this range

### Output
The output is a file with format of

r_helio  total_fragments
-------  ---------------
1.0      1.0e30
1.1      1.1e30
...
etc.
