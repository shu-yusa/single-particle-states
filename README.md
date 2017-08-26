# Single Particle State
Code for calculating nuclear single particle states with the Woods-Saxon potential.

## Prerequisites
`make` and `gfortran` need to be installed.

## Compile & Execute
Compile the source codes by the following.
```bash
make
```
This will produce an executable file `a.out`. Execute the calculation by
```bash
./a.out
```

## Input file
`input_scat` is the input file for the calculation. You can configure reaction details and output directory.

## Details of source codes and modules
* global_constant.f90 ... Definition of constants.
* input_data.f90 ... Read input parameters from the input file.
* potentials.f90 ... Class for calculating potential.
* sp_basis.f90 ... Class for calculating single particle states.
* main.f90 ... Main program.
* print_array.f90 ... Utility functions for showing array on the terminal.

## Output
Binding energies are printed on the screen respectively for neutron and proton.
