#!/usr/bin/python3

import glob
import operator
import numpy as np
import mdtraj as mdt
import matplotlib.pyplot as plt
# import bayesaxs as bs

###########################

#    Global variables     #

###########################

title = "HOIP"
theta = 100
std = 1
path_to_top = "HOIP_extended.pdb"
path_to_traj = "HOIP_extended_combined_MD_simulations_1us.xtc"


###########################

#   Significant states    #

###########################


def load_txt(path_to_txt):
    """ Load .txt file."""
    return np.loadtxt(path_to_txt)


def get_significant_indices(path_to_weights, std):
    """ Extract indices of significant weights."""
    wopt = load_txt(path_to_weights)
    std_cutoff = wopt.mean() + std * wopt.std()
    significant_idx = np.where(wopt >= std_cutoff)[0]
    return wopt, significant_idx


def visualise_significant_weights(wopt, std):
    """ Visualise significant states using standard deviation cutoff."""
    tick_params = dict(labelsize=22, length=10, width=1)

    fs = 20
    fig = plt.figure(figsize=[8, 5])
    ax = fig.add_subplot(111)

    plt.plot(wopt, label="Optimized weights", color="k", zorder=1)
    plt.hlines(wopt.mean(), xmin=0, xmax=wopt.shape[0], label="Mean", color="tab:blue", zorder=2, linewidth=2)
    plt.hlines(wopt.mean() + std * wopt.std(),
            xmin=0,
            xmax=wopt.shape[0],
            label= "{} $\\times$ $\\sigma$".format(std),
            color="tab:orange",
            zorder=2,
            linewidth=2)

    ax.tick_params(**tick_params)

    ax.legend(fontsize=tick_params["labelsize"] - 10, frameon=False)

    ax.set_xlabel("Indices", fontsize=fs)
    ax.set_ylabel("Optimised weights", fontsize=fs)

    plt.tight_layout()
    plt.savefig("HOIP_extended_significant_states_std_{}.png".format(std), dpi=300)


def save_significant_traj(path_to_traj, path_to_top, significant_idx):
    """ Save significant trajectory states."""
    traj = mdt.load(path_to_traj, top=path_to_top)
    significant_traj = traj[significant_idx]
    significant_traj.save_xtc("MaxEnt_significant_traj_frames.xtc")
    return
    
    
wopt, significant_idx = get_significant_indices(path_to_weights="bioen_wopt_theta_{}_{}.txt".format(theta, title), std=std)
visualise_significant_weights(wopt=wopt, std=std)
save_significant_traj(path_to_traj=path_to_traj, path_to_top=path_to_top, significant_idxsignificant_idx)


###################################

#   Cluster significant states    #

###################################

# Initialize HDBSCAN clusterer
clustering = bs.HDBSCAN()
clustering.set_title(title)

# Load the trajectory
clustering.load_traj(top_path=path_to_top, traj_path="MaxEnt_significant_traj_frames.xtc")

# Select CA atoms of RBR subdomains for HBSCAN clustering
# According to indices in .pdb topology file, "name CB" selects CA atoms, and not "name CA" 
ca_atoms_R1 = clustering.get_traj().topology.select("name CB and residue {} to {}".format(699, 751))
ca_atoms_IBR = clustering.get_traj().topology.select("name CB and residue {} to {}".format(797, 841))
ca_atoms_R2 = clustering.get_traj().topology.select("name CB and residue {} to {}".format(869, 935))

# Combine CA atom indices
ca_atom_selection = np.hstack((ca_atoms_R1, ca_atoms_IBR, ca_atoms_R2))

# Perform clustering on pairwise interdomain distances
clustering.fit_predict(metric="distances", atom_selection=ca_atom_selection)

# Save cluster labels
np.save("MaxEnt_cluster_labels", clustering.get_cluster_labels())

# Save traj clusters and extract cluster leaders
clustering.save_traj_clusters()
clustering.save_cluster_leaders(metric="distances", atom_selection=ca_atom_selection)


def calc_rgyr(pdb):
    """ Calculate a radius of gyration for a protein molecule."""
    traj = mdt.load(pdb, top=pdb)
    rgyr = mdt.compute_rg(traj)[0] * 10 # Convert to Angstroms
    return rgyr


def get_top_labels(labels):
    """ Get percentages of top clusters."""
    
    # Remove noise
    labels = labels[np.invert(labels == -1)]

    # Count unique instances and convert to percentage
    unique, counts = np.unique(labels, return_counts=True)
    percentage = np.divide(counts, labels.shape[0])

    # Dictionary of labels and percentages
    label_percentage = dict(zip(unique, percentage))

    # Get a sorted labels and percentages dictionary from higest to lowest
    sorted_label_percentage = list(reversed(sorted(label_percentage, key=label_percentage.get)))

    # Get the top
    top_labels = {}
    
    with open('MaxEnt_ensemble_weights.txt', 'w') as f:
    
        for label in sorted_label_percentage:
            percentage = np.round(label_percentage[label], 5)
            path_to_cluster_leader = "{}_cluster_leaders/cluster_leader_{}.pdb".format(title, label)
            f.write("Label: {}, Percentage: {}, Rg: {} \n".format(label, percentage, calc_rgyr(path_to_cluster_leader)))
            top_labels[label] = percentage    
    
    return top_labels

# Extract top labels
top_labels = get_top_labels(np.load("MaxEnt_ensemble_cluster_labels.npy"))


###################################

#  Visualize significant states   #

###################################

# Generate pymol script
with open('visualize_significant_clusters.pml', 'w') as f:

    # Misc options
    f.write("# Setting misc options \n")
    f.write("hide all \n")
    f.write("bg_color white \n")
    f.write("unset opaque_background \n")
    f.write("unset depth_cue \n")
    f.write("\n")
    
    # Load cluster leaders
    f.write("# Loading cluster leaders \n")
    for idx, label in enumerate(top_labels.keys()):
        f.write("load cluster_leader_{}.pdb \n".format(label))
    f.write("\n")
    
    # Additional misc options
    f.write("# More misc options \n")
    f.write("as cartoon \n")
    f.write("set ray_trace_gain, 0 \n")
    f.write("set ray_trace_mode, 1 \n")
    # f.write("set ray_trace_color, black\n")
    f.write("\n")
    
    # Select domains
    f.write("# Domain selection \n")
    f.write("select RING1, resi 699-751 \n")
    f.write("select L1, resi 752-795 \n")
    f.write("select IBR, resi 796-841 \n")
    f.write("select L2, resi 842-867 \n")
    f.write("select RING2, resi 868-1071 \n")
    f.write("select ZINCS, name ZN \n")
    f.write("select C885, resi 885 \n")
    f.write("\n")

    # Domain colors
    f.write("# Domain colors \n")
    f.write("color skyblue, RING1 \n")
    f.write("color gray80, L1 \n")
    f.write("color tv_orange, IBR \n")
    f.write("color gray80, L2 \n")
    f.write("color forest, RING2 \n")
    f.write("color tv_yellow, C885 \n")
    f.write("show spheres, C885  \n")
    f.write("show spheres, ZINCS \n")
    f.write("\n")
    
    # Define cluster colors
    f.write("# Cluster colors \n")
    max_cluster = max(top_labels.items(), key=operator.itemgetter(1))[0]
    
    # Align each cluster leader
    f.write("# Align each cluster leader on IBR \n")
    f.write("select alignment, model cluster_leader_{} and resi 796-841 \n".format(max_cluster))
    for idx, label in enumerate(top_labels.keys()):
        f.write("align cluster_leader_{}, alignment \n".format(label))
    f.write("\n")
