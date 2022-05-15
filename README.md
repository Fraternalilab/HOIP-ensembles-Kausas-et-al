# Characterisation of HOIP RBR E3 ligase conformational dynamics using integrative modelling 

## Description

The code repository contains experimental data, computational observables, analysis code and conformational ensemble generation pipelines used in the generation of figures of paper.

## List of contents:

- **experimental_data**: contains experimental SAXS scattering curve and ATSAS analysis output.
- **MD_ensemble**: contains representative clusters of the MD trajectory, collective variable (radius of gyration and interdomain distances) and theoretical scattering profiles.
- **MaxPars_ensemble**: contains a script pipeline used for generating a MaxPars ensemble and analysis output files.
- **MaxEnt_ensemble**: contains scripts to setup, run and analyse BioEn optimisation to generate MaxEnt ensemble. In addition, contains output files containing optimized weights and scattering profiles for each $\theta$ and extracted significant states.
- **figures.ipynb**: A Jupyter notebook containing the code for reproducing main text and supplementary figures.

## Required packages

- For construction of MaxPars ensemble, please install [_bayesaxs_ package](https://github.com/mariuskausas/bayesaxs).
- For construction of MaxEnt ensemble, please install [BioEn](https://github.com/bio-phys/BioEn)
