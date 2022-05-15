# Load libraries
import numpy as np
import matplotlib.pyplot as plt
import bayesaxs as bs

# Initialize HDBSCAN clusterer
clustering = bs.HDBSCAN()
clustering.set_title("HOIP_extended")

# # Load the trajectory
clustering.load_traj(top_path="MD_ensemble/MD_traj/HOIP_extended.pdb", traj_path="MD_ensemble/MD_traj/HOIP_extended_combined_MD_simulations_1us.xtc")

# Select CA atoms of RBR subdomains for HBSCAN clustering
ca_atoms_R1 = clustering.get_traj().topology.select("name CA and residue {} to {}".format(699, 751))
ca_atoms_IBR = clustering.get_traj().topology.select("name CA and residue {} to {}".format(797, 841))
ca_atoms_R2 = clustering.get_traj().topology.select("name CA and residue {} to {}".format(869, 935))

# Combine CA atom indices
ca_atom_selection = np.hstack((ca_atoms_R1, ca_atoms_IBR, ca_atoms_R2))

# Perform clustering on pairwise interdomain distances
clustering.fit_predict(metric="distances", atom_selection=ca_atom_selection)

# Save cluster labels
np.save("HOIP_extended_cluster_labels", clustering.get_cluster_labels())

# Load cluster labels
clustering.load_cluster_labels("HOIP_extended_cluster_labels.npy")

# Save traj clusters and extract cluster leaders
clustering.save_traj_clusters()
clustering.save_cluster_leaders(metric="distances", atom_selection=ca_atom_selection)

# Initialize experimental curve
curve = bs.Curve()
curve.load_txt("experimental_data/HOIP_experimental_SAXS.dat")

# Scattering clustering
analysis = bs.Scatter()
analysis.set_title("HOIP_extended")
analysis.load_cluster_leaders("./HOIP_extended_cluster_leaders/")

analysis.calc_scattering(curve)
analysis.load_fits("HOIP_extended_fits/")

analysis.calc_pairwise_chi_matrix()
analysis.cluster_fits(method="ward", metric="euclidean", cutoff_value=0.15)
analysis.calc_representative_fits()

# Bayesian Markov chain Monte Carlo sampling
repfits = bs.Scatter()
repfits.load_representative_fits("HOIP_extended_repfits/")

sampler = bs.Sampler()
sampler.load_curves(repfits.get_representative_fits())
n_curves = len(sampler.get_curves())
states = sampler.inference_single_basis(n_states=n_curves,
                                        step="metropolis",
                                        num_samples=50000,
                                        chains=1,
                                        burn=1000)
bs.save_pickle("HOIP_extended_basis_set", states)
