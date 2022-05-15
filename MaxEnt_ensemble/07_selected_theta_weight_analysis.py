#!/usr/bin/python3

import glob 
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns
import scipy.stats as stats


def load_txt(path_to_txt):
    """ Load .txt."""
    return np.loadtxt(path_to_txt)


# Load observables
rgyr = load_txt("rgyr.txt")
d1 = load_txt("d1.txt")
d2 = load_txt("d2.txt")
d3 = load_txt("d3.txt")
chi2_skl = load_txt("chi2_skl.txt")
wopt_cumsum = load_txt("bioen_wopt_cumsum_theta_100.txt")


# Plot chi2 vs Skl
# ================

fs = 20
fig = plt.figure(figsize=[4, 3.5])
ax = fig.add_subplot(111)

theta = 100

scatter = ax.scatter(chi2_skl[:,1], chi2_skl[:,0],
            marker='o',
            s=20,
            c=np.arange(8),
            facecolors="none",
            linewidths=2,
            cmap="Spectral")

labels = ["$\\theta = 10^6$",
          "$\\theta = 10^5$",
          "$\\theta = 10^4$",
          "$\\theta = 10^3$",
          "$\\theta = 10^2$", 
          "$\\theta = 10$",
          "$\\theta = 1$",
          "$\\theta = 0$"][::-1]

for idx, point in enumerate(chi2_skl):
        
        x = point[1]
        y = point[0]
        
        if labels[idx] == "$\\theta = 10^5$":
            ax.annotate(labels[idx],
                     (x, y),
                     xytext=(4, -4),
                     textcoords='offset points',
                     fontsize=14)
        
        elif labels[idx] == "$\\theta = 10^2$":
            ax.annotate(labels[idx],
                     (x, y),
                     xytext=(4, 10),
                     textcoords='offset points',
                     fontsize=14)
        
        elif labels[idx] == "$\\theta = 1$":
            ax.annotate(labels[idx],
                     (x, y),
                     xytext=(8, 6),
                     textcoords='offset points',
                     fontsize=14)
        
        else:
            ax.annotate(labels[idx],
                     (x, y),
                     xytext=(4, 6),
                     textcoords='offset points',
                     fontsize=14)
        
theta_idx = 3
        
x_side = 0.3
y_side = 0.08
rect = matplotlib.patches.Rectangle((chi2_skl[theta_idx:theta_idx + 1, 1] - x_side / 2,
                                    chi2_skl[theta_idx:theta_idx + 1, 0] - y_side / 2),
                                    x_side, 
                                    y_side,
                                    linewidth=1.5,
                                    edgecolor='k',
                                    facecolor='none')
ax.add_patch(rect)

plt.hlines(1, -1, 10, linewidths=2, linestyle="--", color="k")

ax.set_xlabel(r'$S_{\mathrm{KL}}$', fontsize=fs)
ax.set_ylabel(r'$\chi^{2}$', fontsize=fs)

ax.set_xlim(-0.2, 5.4)
ax.set_ylim(0.95, 2.25)

ax.tick_params(labelsize=16)

ax.set_xticks([0, 1, 2, 3, 4, 5])
ax.set_yticks([1, 1.5, 2])

custom_lines = [Line2D([0], [0], color="k", lw=2)]
ax.legend(custom_lines, ["$S_{KL} = $" + "{}".format(round(chi2_skl[theta_idx:theta_idx + 1, 1][0], 2))],
          loc="upper right",
          fontsize=14,
          frameon=False)

plt.tight_layout()
plt.savefig("HOIP_extended_selected_chi2_skl.png", dpi=300)


# Plot cumulative sum for the optimised weights at a given theta
# ==============================================================

fs = 20
fig = plt.figure(figsize=[4, 3.5])
ax = fig.add_subplot(111)

shape = wopt_cumsum.shape[0]
x = np.arange(1, shape + 1)
ref = np.cumsum(np.ones(shape) / shape)

plt.plot(x, ref, linewidth=3, color="k", linestyle="--", label="Reference")
plt.plot(x, wopt_cumsum, linewidth=3, color="tab:cyan", label="$\\theta=100$")

ax.set_xlabel("Weight index", fontsize=fs)
ax.set_ylabel("Cumulative weights", fontsize=fs)

ax.set_ylim(0, 1)

ax.tick_params(labelsize=16)

ax.set_xticks([0, 25000, 50000])
ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])

ax.legend(loc="lower right", fontsize=14, frameon=False)

plt.tight_layout()
plt.savefig("HOIP_extended_theta_100_wopt_cumsum.png", dpi=300)


# Plot experimental scattering vs MD and MaxEnt optimised ensembles
# =================================================================


def get_exp(path_to_dat, i0, skiprows):
    """ Load I(0) normalized experimental scattering data."""
    data = np.loadtxt(path_to_dat, skiprows=skiprows)
    q = data[:,0]
    exp = data[:,1] / i0
    sigma = data[:,2] / i0
    return q, exp, sigma


# Define i0 for normalization
i0_n = 0.015

# Load experimental data
q, exp, sigma = get_exp(path_to_dat="HOIP_experimental_SAXS.dat",
                        i0=i0_n,
                        skiprows=0)


def get_fit(path_to_fit, skiprows, i0):
    """ Get theoretical fit values."""
    data = np.loadtxt(path_to_fit, skiprows=skiprows)
    fit = data[:,3] / i0
    return fit


# Load fits
fit_md = get_fit("HOIP_extended_combined_MD_simulations_1us_averaged_crysol_fit.fit", skiprows=0, i0=i0_n)
fit_bioen = get_fit("bioen_scattering_wopt_theta_100.dat", skiprows=0, i0=i0_n)
fits = [fit_md, fit_bioen]

# Define residuals
residuals_md = (exp - fit_md) / sigma
residuals_bioen = (exp - fit_bioen) / sigma
residuals = [residuals_md, residuals_bioen]

# Plot a fit between experimental and theoretical scattering data
tick_params = dict(labelsize=22, length=10, width=1)
marker_style = dict(c='tab:blue', marker='o', s=25, alpha=0.4, edgecolors="k")

fig = plt.figure(figsize=[8, 5])

titles = ["MD", "MaxEnt, $\\theta = 100$"]
colors=["tab:red", "deepskyblue"]

ax1 = plt.subplot2grid((4, 3), (0, 0), colspan=3, rowspan=3)

ax1.plot(q, exp, color="k", linewidth=3, zorder=2, label="Experimental SAXS")
ax1.fill_between(q, exp-sigma, exp+sigma, color='tab:grey', alpha=1, zorder=1)

for idx, fit in enumerate(fits):
    ax1.plot(q, fit, color=colors[idx], alpha=1, linewidth=3, label=titles[idx])

ax1.semilogy()
ax1.tick_params(labelsize=22, axis="y", length=10)
ax1.set_xticklabels([])
ax1.set_xlim(0, 0.28)
ax1.set_ylabel('$I(q)/I(0)$', fontsize=22)
ax1.legend(fontsize=14, frameon=False)


ax2 = plt.subplot2grid((4, 3), (3, 0), colspan=3)

for idx, residual in enumerate(residuals):
    ax2.scatter(q, residual, s=3, color=colors[idx], zorder=1, marker='o', alpha=1)
    
ax2.axhline(y=0, xmin=0, xmax=1, ls='--', color="k", zorder=2, linewidth=2)
ax2.tick_params(**tick_params)
ax2.set_xlim(0, 0.28)
ax2.set_xlabel('$q$ $(\AA^{-1})$', fontsize=22)
ax2.set_yticks([-4, 4])
ax2.set_ylabel('$(I\Delta)/\sigma_{exp}$', fontsize=22)

plt.tight_layout()
plt.savefig("HOIP_extended_experimental_scattering_vs_MD_and_MaxEnt.png", dpi=300, bbox_inches='tight')


# Plot MaxEnt optimised distribution of radius of gyration values
# ===============================================================

fs = 22
fig = plt.figure(figsize=[8, 5])
ax = fig.add_subplot(111)

density = stats.kde.gaussian_kde(rgyr)
density_weighted = stats.kde.gaussian_kde(rgyr, weights=weight)
x = np.arange(20, 50, .1)

ax.vlines(30.33, 0, 1, linestyle="--", linewidth=3, color="k", label=r"Experimental $\langle R_g \rangle$")
ax.vlines(rgyr.mean(), 0, 1, linestyle="--", linewidth=3, color="tab:red", label=r"MD $\langle R_g \rangle$")
ax.vlines((rgyr * weight).sum(), 0, 1, linestyle="--", linewidth=3, color="deepskyblue", label=r"MaxEnt $\langle R_g \rangle$, $\theta$ = 100")

ax.plot(x, density(x),
          color="tab:red",
          linewidth=3)
ax.fill_between(x, density(x), color="tab:red", alpha=0.1)

ax.plot(x, density_weighted(x),
          color="deepskyblue",
          linewidth=3)
ax.fill_between(x, density_weighted(x), color="deepskyblue", alpha=0.1)

ax.set_xlim(20, 50)
ax.set_ylim(0, 0.2)

ax.legend(fontsize=14, frameon=False)

custom_lines = [Line2D([0], [0], color="k", lw=4),
                Line2D([0], [0], color="tab:red", lw=4),
                Line2D([0], [0], color="deepskyblue", lw=4)]

ax.legend(custom_lines, [r"Experimental $\langle R_g \rangle$", r"MD $\langle R_g \rangle$", r"MaxEnt $\langle R_g \rangle$, $\theta$ = 100"], fontsize=14, frameon=False)

ax.set_xticks([20, 25, 30, 35, 40, 45, 50])

ax.tick_params(**tick_params)

ax.set_xlabel("$R_g$ ($\AA$)", fontsize=fs)
ax.set_ylabel('Probability', fontsize=fs)

plt.tight_layout()
plt.savefig("HOIP_extended_maxent_optimized_rgyr.png", dpi=300, bbox_inches = "tight")


# Plot MaxEnt optimised distribution of RING1-IBR interdomain distances
# =====================================================================

tick_params = dict(labelsize=22, length=10, width=1)

fs = 22
fig = plt.figure(figsize=[8, 5])
ax = fig.add_subplot(111)

theta = 100
weight_file_name = glob.glob("bioen_wopt_theta_{}_*.txt".format(theta))[0]
weight = load_txt(weight_file_name)                                                                      
avg_opt = np.sum(d1 * weight)

density = stats.kde.gaussian_kde(d1)
density_weighted = stats.kde.gaussian_kde(d1, weights=weight)
x = np.arange(20, 70, .1)

ax.plot(x, density(x),
          color="tab:red",
          linewidth=3)
ax.fill_between(x, density(x), color="tab:red", alpha=0.1)
ax.vlines(np.mean(d1), 0, 1,
          linestyle="--",
          color="tab:red",
          linewidth=3) 

ax.plot(x, density_weighted(x),
          color="tab:cyan",
          linewidth=3)
ax.fill_between(x, density_weighted(x), color="tab:cyan", alpha=0.1)
ax.vlines(avg_opt, 0, 1,
          linestyle="--",
          color="tab:cyan",
          linewidth=3) 

ax.set_xlim(20, 70)
ax.set_ylim(0, 0.2)
ax.tick_params(**tick_params)

custom_lines = [Line2D([0], [0], color="tab:red", lw=3),
               Line2D([0], [0], color="tab:cyan", lw=3)]
labels = [r"Simulated $\langle D_{R1-IBR} \rangle$", r"Optimized $\langle D_{{R1-IBR}} \rangle, \theta$ = {}".format(theta)]
ax.legend(custom_lines, labels, fontsize=14, frameon=False, loc="upper right")

ax.set_xlabel(r'D$_{R1-IBR}$ ($\AA$)', fontsize=fs)
ax.set_ylabel(r'Probability', fontsize=fs)

plt.tight_layout()
plt.savefig("HOIP_extended_maxent_optimized_d1.png", dpi=300)


# Plot MaxEnt optimised distribution of IBR-RING2 interdomain distances
# =====================================================================

tick_params = dict(labelsize=22, length=10, width=1)

fs = 22
fig = plt.figure(figsize=[8, 5])
ax = fig.add_subplot(111)

theta = 100
weight_file_name = glob.glob("bioen_wopt_theta_{}_*.txt".format(theta))[0]
weight = load_txt(weight_file_name)                                                                      
avg_opt = np.sum(d2 * weight)

density = stats.kde.gaussian_kde(d2)
density_weighted = stats.kde.gaussian_kde(d2, weights=weight)
x = np.arange(20, 85, .1)

ax.plot(x, density(x),
          color="tab:red",
          linewidth=3)
ax.fill_between(x, density(x), color="tab:red", alpha=0.1)
ax.vlines(np.mean(d2), 0, 1,
          linestyle="--",
          color="tab:red",
          linewidth=3) 

ax.plot(x, density_weighted(x),
          color="tab:cyan",
          linewidth=3)
ax.fill_between(x, density_weighted(x), color="tab:cyan", alpha=0.1)
ax.vlines(avg_opt, 0, 1,
          linestyle="--",
          color="tab:cyan",
          linewidth=3) 

ax.set_xlim(20, 85)
ax.set_ylim(0, 0.15)
ax.tick_params(**tick_params)

ax.set_yticks([0, 0.05, 0.1, 0.15])

custom_lines = [Line2D([0], [0], color="tab:red", lw=3),
               Line2D([0], [0], color="tab:cyan", lw=3)]
labels = [r"Simulated $\langle D_{IBR-R2} \rangle$", r"Optimized $\langle D_{{IBR-R2}} \rangle, \theta$ = {}".format(theta)]
ax.legend(custom_lines, labels, fontsize=14, frameon=False, loc="upper right")

ax.set_xlabel(r'D$_{IBR-R2}$ ($\AA$)', fontsize=fs)
ax.set_ylabel(r'Probability', fontsize=fs)

plt.tight_layout()
plt.savefig("HOIP_extended_maxent_optimized_d2.png", dpi=300)


# Plot MaxEnt optimised distribution of RING1-RING2 interdomain distances
# =======================================================================

tick_params = dict(labelsize=22, length=10, width=1)

fs = 22
fig = plt.figure(figsize=[8, 5])
ax = fig.add_subplot(111)

theta = 100
weight_file_name = glob.glob("bioen_wopt_theta_{}_*.txt".format(theta))[0]
weight = load_txt(weight_file_name)                                                                      
avg_opt = np.sum(d3 * weight)

density = stats.kde.gaussian_kde(d3)
density_weighted = stats.kde.gaussian_kde(d3, weights=weight)
x = np.arange(10, 120, .1)

ax.plot(x, density(x),
          color="tab:red",
          linewidth=3)
ax.fill_between(x, density(x), color="tab:red", alpha=0.1)
ax.vlines(np.mean(d3), 0, 1,
          linestyle="--",
          color="tab:red",
          linewidth=3) 

ax.plot(x, density_weighted(x),
          color="tab:cyan",
          linewidth=3)
ax.fill_between(x, density_weighted(x), color="tab:cyan", alpha=0.1)
ax.vlines(avg_opt, 0, 1,
          linestyle="--",
          color="tab:cyan",
          linewidth=3) 

ax.set_xlim(10, 120)
ax.set_ylim(0, 0.1)
ax.tick_params(**tick_params)

ax.set_yticks([0, 0.05, 0.1])

custom_lines = [Line2D([0], [0], color="tab:red", lw=3),
               Line2D([0], [0], color="tab:cyan", lw=3)]
labels = [r"Simulated $\langle D_{R1-R2} \rangle$", r"Optimized $\langle D_{{R1-R2}} \rangle, \theta$ = {}".format(theta)]
ax.legend(custom_lines, labels, fontsize=14, frameon=False, loc="upper right")

ax.set_xlabel(r'D$_{R1-R2}$ ($\AA$)', fontsize=fs)
ax.set_ylabel(r'Probability', fontsize=fs)

plt.tight_layout()
plt.savefig("HOIP_extended_maxent_optimized_d3.png", dpi=300)
