#!/usr/bin/python3

import numpy as np
from scipy import stats
import MDAnalysis as mda
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def get_top_clusters(labels, n_top):
    """ Extract top cluster percentage."""
    
    # Remove HBDSCAN noise
    labels = labels[np.invert(labels == -1)]

    # Count unique instances and convert to percentage
    unique, counts = np.unique(labels, return_counts=True)
    percentage = np.divide(counts, labels.shape[0])

    # Dictionary of labels and percentages
    label_percentage = dict(zip(unique, percentage))

    # Get a sorted labels and percentages dictionary from higest to lowest
    sorted_label_percentage = list(reversed(sorted(label_percentage, key=label_percentage.get)))

    # Get the top
    top_labels = []
    for label in sorted_label_percentage[:n_top]:
        print("Label: {}, Percentage: {}".format(label, np.round(label_percentage[label], 5)))
        top_labels.append(label)    
    
    return top_labels


def domain_dist(domains, domain1, domain2):
    """ Calculate distance between COM of domains."""

    cm = dict((name, dom.center_of_mass()) for name, dom in domains.items())

    return np.linalg.norm(cm[domain1] - cm[domain2])


def calculate_RBR_dist(top):
    """ Calculate distance between three RBR subdomains using their COM for a trajectory."""

    # Initialize the system
    u = mda.Universe(top)

    # Defined domains for three RBR subdomains as a dict
    domains = {
        'RING1': u.select_atoms("resid 699:751 and (backbone or name CB)"),
        'IBR': u.select_atoms("resid 797:841 and (backbone or name CB)"),
        'RING2': u.select_atoms("resid 869:935 and (backbone or name CB)")}

    # Distances (d1, d2, d3) to be stored in an numpy array
    d1d2_distances = np.array([(domain_dist(domains, "RING1", "IBR"), domain_dist(domains, "IBR", "RING2")) for ts in u.trajectory])

    return d1d2_distances[0]


def get_histogram(xval, yval, nbins=100, weights=None):
    """ Compute 2D histogram of two given variables."""
    z, xedge, yedge = np.histogram2d(xval, yval, bins=nbins, weights=weights)
    x = 0.5 * (xedge[:-1] + xedge[1:])
    y = 0.5 * (yedge[:-1] + yedge[1:])
    return z.T, x, y


def density(z):
    """ Compute a probability density function."""
    return z / float(z.sum())


def free_energy(z):
    """ Log-transform of the probability density function."""
    prob = density(z)
    free_energy = np.inf * np.ones(shape=z.shape)
    nonzero = prob.nonzero()
    free_energy[nonzero] = -np.log(prob[nonzero])
    return free_energy


def convert_fes(cv1, cv2, bins, weights=None):
    """ 2D histogram two collective variables and prepare input for FES plot.

    The log values are taken of two collective variables and 2D histogram is arranged.
    Weights can be provided to weight the resulted 2D histogram. x and y should be provided
    as (n,) numpy arrays.

    """
    ## Generate a 2D histogram of first two collective variables
    # Define range of collective variables
    xmin = np.floor(cv1.min())
    xmax = np.ceil(cv1.max())
    ymin = np.floor(cv2.min())
    ymax = np.ceil(cv2.max())

    # Calculate bin size for each collective variable
    xbin_size = np.abs(xmax - xmin) / float(bins)
    ybin_size = np.abs(ymax - ymin) / float(bins)
    
    xedges = np.arange(xmin, xmax, xbin_size)
    yedges = np.arange(ymin, ymax, ybin_size)

    # Histogram values into 2D array
    H, xedges, yedges = get_histogram(cv1, cv2, bins, weights=weights)
    f = free_energy(H)

    return f, xedges, yedges


def load_txt(path_to_txt):
    """ Load .txt file."""
    return np.loadtxt(path_to_txt)


def difference_fes(path_to_cv1, path_to_cv2, path_to_weights, bins, labels_dict_me, output_name):
    """ Plot a difference FES."""

    # Load collective variables and weights
    cv1 = load_txt(path_to_cv1)
    cv2 = load_txt(path_to_cv2)
    weights = load_txt(path_to_weights)

    # Convert CVs to a matrix
    H, xedges, yedges = convert_fes(cv1, cv2, bins)
    Hw, xedgesw, yedgesw = convert_fes(cv1, cv2, bins, weights)
    
    # Normalize to zero
    H -= H.min()
    Hw -= Hw.min()
    
    # Remove infs
    Hinf = H[~np.isinf(H)]
    Hwinf = Hw[~np.isinf(Hw)]
    
    # Set contour level base
    levels_base = np.round(Hinf.max())
    
    # Set range for non-weighted and weighted matrices
    if Hinf.min() < Hwinf.min():
        vmin = np.floor(Hinf.min())
    else:
        vmin = np.floor(Hwinf.min())
        
    if Hinf.max() > Hwinf.max():
        vmax = np.ceil(Hinf.max())
    else:
        vmax = np.ceil(Hwinf.max())
        
    # Set your own vmax
    vmax = 5
    
    # Set font and figure dimensions
    fs = 20
    fig = plt.figure()
    fig, (ax0, ax1, ax2) = plt.subplots(figsize=[22, 6], ncols=3)
    
    # Non-weighted
    im0 = ax0.contourf(H, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
                        cmap="Spectral", # GnBu_r
                        levels=np.linspace(vmin, int(vmax), int(levels_base * 2 + 1)),
                        corner_mask=True,
                        alpha=0.75)
    im0c = ax0.contour(H, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
                        cmap="gist_gray",
                        levels=np.linspace(vmin, int(vmax), int(levels_base * 2 + 1)),
                        linewidths=0.5,
                        corner_mask=True,
                        alpha=1)
    
    ax0.set_title("$-\ln p_{MD}$", fontsize=fs+4)

    cbar0 = fig.colorbar(im0, ax=ax0)
    cbar0.ax.tick_params(labelsize=20)
    cbar0.set_ticks(np.arange(vmin, vmax + 1, 1))

    ax0.tick_params(labelsize=20)
    ax0.set_ylabel("$D_{IBR-RING2}$ $(\AA)$", fontsize=fs)
    ax0.set_xlabel("$D_{RING1-IBR}$ $(\AA)$", fontsize=fs)

    # Weighted
    im1 = ax1.contourf(Hw, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
                        cmap="Spectral", # GnBu_r
                        levels=np.linspace(vmin, int(vmax), int(levels_base * 2 + 1)),
                        corner_mask=True,
                        alpha=0.75)
    im1c = ax1.contour(Hw, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
                        cmap="gist_gray",
                        levels=np.linspace(vmin, int(vmax), int(levels_base * 2 + 1)),
                        linewidths=0.5,
                        corner_mask=True,
                        alpha=1)

    ax1.set_title("$-\ln p_{MaxEnt}$", fontsize=fs+4)
        
    custom_lines = [Line2D([0], [0], color="deepskyblue", lw=4)]

    cbar1 = fig.colorbar(im1, ax=ax1)
    cbar1.ax.tick_params(labelsize=20)
    cbar1.set_ticks(np.arange(vmin, vmax + 1, 1))

    ax1.tick_params(labelsize=20)    
    ax1.set_ylabel("$D_{IBR-RING2}$ $(\AA)$", fontsize=fs)
    ax1.set_xlabel("$D_{RING1-IBR}$ $(\AA)$", fontsize=fs)
    
    # Difference matrix
    diff = Hw - H
    
    # Remove nans
    diffnan = diff[~np.isnan(diff)]
    
    # Set range for difference matrix
    if np.abs(diffnan.min()) > np.abs(diffnan.max()):
        vmin_max_diff = np.ceil(np.abs(diffnan.min()))
    else:
        vmin_max_diff = np.ceil(np.abs(diffnan.max()))

    
    im2 = ax2.contourf(diff, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
                        cmap="seismic_r", # GnBu_r
                        levels=np.linspace(-vmin_max_diff, int(vmin_max_diff), int(levels_base * 4 + 1)),
                        corner_mask=True,
                        alpha=0.75)
    im2c = ax2.contour(diff, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
                        cmap="gist_gray",
                        levels=np.linspace(-vmin_max_diff, int(vmin_max_diff), int(levels_base * 4 + 1)),
                        linewidths=0.75,
                        corner_mask=True,
                        alpha=1)
    
    # Plot cluster labels
    for key in labels_dict_me.keys():
        
        x = labels_dict_me[key][0]
        y = labels_dict_me[key][1]

        ax1.plot(x, y, 'o', markersize=15, color="deepskyblue", markeredgewidth=2, markeredgecolor='k')
        ax1.annotate(key,
                    (x, y),
                    xytext=(4, 6),
                    textcoords='offset points',
                    fontsize=20,
                    weight="bold")
        
        ax2.plot(x, y, 'o', markersize=15, color="deepskyblue", markeredgewidth=2, markeredgecolor='k')
        ax2.annotate(key,
                    (x, y),
                    xytext=(4, 6),
                    textcoords='offset points',
                    fontsize=20,
                    weight="bold")

    ax2.set_title("$\Delta G$ MaxEnt - MD", fontsize=fs+4)
    ax2.set_ylabel("$D_{IBR-RING2}$ $(\AA)$", fontsize=fs)
    ax2.set_xlabel("$D_{RING1-IBR}$ $(\AA)$", fontsize=fs)

    cbar2 = fig.colorbar(im2, ax=ax2)
    cbar2.ax.tick_params(labelsize=20)
    cbar2.set_ticks(np.arange(-vmin_max_diff, vmin_max_diff + 1, vmin_max_diff / 2))
    ax2.tick_params(labelsize=20)
    
    plt.tight_layout()
    
    plt.savefig(output_name + ".png", dpi=300)
    plt.show()


# Calculate bioen labels
bioen_labels = [0, 15, 12, 17]
bioen_updated_labels = [1, 2, 3, 4]

def top_label_parameter_bioen(top_labels):
    parameters = []
    for label in top_labels:
        parameters.append(calculate_RBR_dist("./output.files/significant_states/HOIP_cluster_leaders/cluster_leader_{}.pdb".format(label)))
    return np.array(parameters)

top_param_bioen = top_label_parameter_bioen(bioen_labels)
top_dict_bioen = dict(zip(bioen_updated_labels, top_param_bioen))

# Plot the difference FES
difference_fes(path_to_cv1="d1.txt",
            path_to_cv2="d2.txt",
            bins=100,
            path_to_weights="bioen_wopt_theta_100_HOIP_extended.txt",
            labels_dict_me=top_dict_bioen,
            output_name ="HOIP_extended_fes_difference")
            