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
 
TODO

 * Correct internal memory leaks (Valgrind)
 * Verify stand alone version against `Cactus`/`EinsteinToolkit`
 * Verify Cartesian grid interpolator against `Cactus`/`EinsteinToolkit` and/or the stand alone code

SB 01/2019,
BD 09/2019
