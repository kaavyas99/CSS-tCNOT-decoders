# Transversal CNOT decoders for a pair of surface codes

This is C++ code to implement 

1) a 2SCQM with MWPM decoding

2) a tCNOT with single-update and ordered MWPM decoding

## Features

- This code is a Pauli frame simulator designed to quickly Monte-Carlo sample simulations of a pair of $d \times d$ rotated surface codes measured for $2d$ rounds, either with a transversal CNOT applied after round $d$ or not.

- It currently supports two-qubit Pauli errors, randomly chosen at uniform from $\{ I,X,Y,Z\} ^{\otimes 2} \ \{I \otimes I\}$, along with erasures after two-qubit gates.

- The noise model, number of measurement rounds, and gate schedules are adjustable within the code. 

- The base decoder is minimum-weight perfect matching (MWPM), with Djikstra inserted in order to account for edge weights adjusted by erasures. 

## Compilation

The code can be compiled with the provided makefile using the standard `make` command. It uses gcc 13.2.0 

The underlying MWPM package being used is the Kolmogorov implementation of Blossom V. The relevant files have been removed to comply with the package's distribution guidelines. To use this code, download the Blossom V files and copy them into the root directory; copy `GeomPerfectMatching.h` from `GEOM` into the root directory.

## Usage

./run takes in 11 input parameters:

./run 

1: `operation`

2: `d_min`   3: `d_max` 4: `d_step` 

5: `p_min` 6: `p_max` 7: `p_step` 

8: `R_e` 9:`erasure_form` 

10:`no_of_samples` 11:`rand_seed`

Where 

- the `operation` can be `d` (decoupled surface codes measured for $2d$ rounds), `s` (tCNOT with single-update decoding), or `o` (tCNOT with ordered decoding)

- $d$ represents a range of system sizes (integers) that are run from ($d_{min}$, $d_{max}$ , $d_{step}$ )
 
- $p$ represents physical error rates (floats less than 1) that are run from ($p_{min}$, $p_{max}$ , $p_{step}$ )

- $R_e$ is the erasure fraction, i.e. The Pauli error rate is $(1-p) * R_e$ and the mutually exclusive erasure error rate is $p * R_e$

- The `erasure_form` can be `u` (unbiased), `n` (biased with native gates), or `b` (biased with bias-preserving gates).

- `no_of_samples` determines the number of Monte Carlo samples collected for each ($d$, $p$ ) combination

- `rand_seed` determines the starting random seed for the sampling


For example,

```
./run d 3 4 2 0.001 0.002 0.001 0 u 1000 0
```

Will return 

```
simulation paramters: (transversal_CNOT=0), (decoder_variant=d)
3 0.001 0 u 1 0 2 0 3 1000
```

where the numbers after `u` are  
```
ctrlXFailCount  ctrlZFailCount trgtXFailCount trgtZFailCount FailCount NoOfSamples
```

## Citation

If you want to cite this code in an academic work, please cite our upcoming arXiv preprint.

