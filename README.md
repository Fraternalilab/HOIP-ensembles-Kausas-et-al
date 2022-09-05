# Characterisation of HOIP RBR E3 ligase conformational dynamics using integrative modelling 

## Description

The code repository contains experimental data, computational observables, analysis code and conformational ensemble generation pipelines used in the generation of figures of paper.

## List of contents:

- **experimental_data**: contains experimental SAXS scattering curve and ATSAS analysis output.
- **MD_ensemble**: contains representative clusters of the MD trajectory, collective variable (radius of gyration and interdomain distances) and theoretical scattering profiles.
- **MaxPars_ensemble**: contains a script pipeline used for generating a MaxPars ensemble and analysis output files.
- **MaxEnt_ensemble**: contains scripts to setup, run and analyse BioEn optimisation to generate MaxEnt ensemble. In addition, contains output files containing optimized weights and scattering profiles for each theta and extracted significant states.
- **figures.ipynb**: A Jupyter notebook containing the code for reproducing main text and supplementary figures.

## Zenodo repository:

Zenodo [repository](https://zenodo.org/record/7041795) contains the rest of raw simulation and ensemble data. The repository follows same folder structure as the GitHub repository.

- **MD_ensemble/MD_traj**: contains the combined 1 us simulation .xtc trajectory of the extended HOIP RBR domain. 
- **MaxEnt_ensemble/input.files**: contains _BioEn_ simulation and output files, such as formatted experimental scattering .dat file, BioEn pickle input file with set simulation parameters and formatted simulated scattering data for each trajectory frame of extended HOIP RBR domain. The simulation and output files were generated using _MaxEnt_ensemble/0X_ scripts.
- **MaxEnt_ensemble/output.files/significant_states**: contains a trajectory .xtc of MaxEnt significant conformers.
- **MaxEnt_ensemble/combinations**: contains _BioEn_ input and output files for different ensemble combinations of individual simulations.

## Required packages

- For construction of MaxPars ensemble, please install [_bayesaxs_](https://github.com/mariuskausas/bayesaxs) package.
- For construction of MaxEnt ensemble, please install [BioEn](https://github.com/bio-phys/BioEn) package.
