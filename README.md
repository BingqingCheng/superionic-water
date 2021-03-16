# Supplementary Data set for "Predicting the phase behaviors of superionic water atplanetary conditions".

## Data

EOS: tabulated equation of states
			
diffusion-coefficient: diffusion coefficients for oxygen and hydrogen for superionic and liquid water	

phase-boundary: melting lines (Tm) and ice-superionic transition lines (Ts).

low-T-ice: structures of the low-temperature ice phases used for benchmarking DFT.

quantum-kinetic-energy: kinetic energies of quantum-mechanical nuclei (Ek).

chemical-potential: tabulated chemical potentials for the superionic phases (bcc, fcc, hcp) and liquid water.	

## Input files

nnp-md-example: LAMMPS input for molecular dynamics simulations.

nnp-pimd-example: i-PI + LAMMPS input for path-integral molecular dynamics simulations.

nnp-sl-umbrella-example: LAMMPS + PLUMED input for running umbrella sampling simulations on coexistence systems.

## Machine-learning potential
nnp: specifications of the neural-network potential for high-pressure water based on PBE DFT. Competitable with N2P2 or RuNNer.

training-set: training set of the MLP, N2P2/RuNNer format. In atomic units (Hartree/bohr).
