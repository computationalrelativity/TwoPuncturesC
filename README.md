# TwoPunctures

Create initial data for two puncture black holes using a single domain
spectral method following

```
Marcus Ansorg, Bernd Brügmann, Wolfgang Tichy,
"A single-domain spectral method for black hole puncture data",
PRD 70, 064011 (2004), arXiv:gr-qc/0404056.
```

This is a stand-alone version of `Cactus`/`EinsteinToolkit` C code:

```
https://bitbucket.org/einsteintoolkit/einsteininitialdata/src/master/TwoPunctures/
```

Main modifs from `Cactus`/`EinsteinToolkit`

 * Switched off code for matter terms (`sources` array is always to 0)
 * Simplified some prototypes in `TP_Newton.c` 
 * Changed filenames and reduced to one header file
 * Parameters managed in a minimal/dummy way, see `TP_Utilities.c`
 * Commented out the line calling `Set_Initial_Guess_for_u()`

Please cite the [original paper](https://arxiv.org/abs/gr-qc/0404056) and [this paper](https://arxiv.org/abs/2101.08289) if you use this code. 

## HowTo

 * Compile (using `gcc`) with:
   ```make clean && make```
   
 * Compile (using `icc`) with:
   ```make clean && make CC=icc```
 
 * Additional flags may be similarily passed (see `Makefile`).

Shared library will automatically be generated.
 
## Interface with GR-Athena++

Use the `Makefile` here- as is.
Linking without omp suppression in the header and compile Athena++.

See [this issue](https://bitbucket.org/bernuzzi/twopuncturesc/issues/3/athena-interface-and-omp)
 
## History

 * SB git issues activated from [v1.0.0](https://bitbucket.org/bernuzzi/twopuncturesc/commits/tag/v1.0.0)
 * FZ 09/2020 checked mem leaks (`fb7d1c0`)
 * SB 12/2019 changed the way `derivs` data are allocated, deallocate, and passed around. Extended parameter functionality, added output routines for derivs and bam's puncture_ps.
 * BD 09/2019 created lib
 * SB 01/2019 started
