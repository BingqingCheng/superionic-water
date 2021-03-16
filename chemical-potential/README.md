mu-mlp-cl-*.dat: chemical potentials of the phase at the classical, MLP level

mu-mlp2pbe-*.dat: the correction term mu^{PBE} - mu^{MLP} computed from free energy perturbation.

mu-nqe-*.dat: the difference between the mu with nuclear quantum effects and the classical mu.

To get the PBE chemical potential with nuclear quantum effects, one can do:
mu-mlp-cl-*.dat + mu-mlp2pbe-*.dat + mu-nqe-*.dat
