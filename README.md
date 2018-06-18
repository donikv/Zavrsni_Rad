# Landau-Vishkin-Nussinov Algorithm for Pair-wise Sequence Alignment

This is a fast implementation of the Landau-Vishkin-Nussinov algorithm for pair-wise alignment of sequences.

## Building
This library uses CMAKE to build libraries (static and shared) and binaries (apps and tests).
Execute following commands to build it using CMAKE:

1. `cd build`
2. `cmake -D CMAKE_BUILD_TYPE=Release ..`
3. `make`

This will create binaries in `build/bin/` directory and libraries (static and shared) in `build/lib/` directory.
You can run `./test.sh` to confirm that the installation was successful.

Also there is a simple build script provided in the repository, which is executed by running `./build.sh` from the root directory.

## Usage
This implementation is available as a standalone library, but it also offers a simple executable named main that can be found in the `build/bin/` directory. This executable is run from the command line, where the first argument should be `SHW` for prefix alignment, `NW` for global, and `HW` for infix alignment. The second and third arguments represent the path to the files containing FASTA representations of the pattern and text for which the alignment should be calculated. Next argument denotes the number of repetitions, and the final argument represents a boolean value which tells the program whether or not it should also reconstruct the alignment path. This is returned in the CIGAR format.

## Publication
Detailed explanation of this implementation is available in the `Paper` directory of this repository, but unfortunately the paper is only available in Croatian.

## Useful links
### Original paper
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC339353/pdf/nar00270-0052.pdf

### Edlib
https://github.com/Martinsos/edlib
