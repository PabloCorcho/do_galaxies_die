#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot Fig. 7 from Corcho-Caballero et al. (2021).

Created on Thu Feb 25 11:38:17 2021
@author: pablo
"""

from read_suite import DataSet
import numpy as np
from matplotlib import pyplot as plt
# from astropy.io import fits

from glob import glob

# from scipy.stats import binned_statistic
from scipy.interpolate import interp1d

# from matplotlib import cm
# from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

data_sets_names = {
    'MAGNETICUM.fits': 'MAG500',
    'tng300-1.fits': 'IllustrisTNG300',
    'tng100-1.fits': 'IllustrisTNG100',
    'SIMBA_m100n1024_151.fits': 'SIMBA150',
    'SIMBA_m50n512_151.fits': 'SIMBA75',
    'SIMBA_m25n512_151.fits': 'SIMBA35',
    'RefL0100N1504.fits': 'EAGLE100',
    'RefL0050N0752.fits': 'EAGLE50',
    'RefL0025N0752.fits': 'EAGLE25'
    }


simcolors = ['darkgreen', 'darkgreen', 'darkgreen',
             'darkorange', 'darkorange', 'darkorange',
             'crimson', 'crimson', 'crimson']

# =============================================================================
# SIMULATIONS
# =============================================================================

data_sets = glob('simu_data/*.fits')

# =============================================================================
# SDSS & GAMA data
# =============================================================================

extra_path = 'obs_data/derived_data/'
gama_all_lgm_ssfr_pdf_vmax = np.loadtxt(extra_path+'gama_fundamentalplane_vmax.txt')
all_lgm_ssfr_pdf_vmax = np.loadtxt(extra_path+'sdss_fundamentalplane_vmax.txt')


ssfr, dssfr = np.loadtxt(extra_path+'gama_prior_ssfr.txt', usecols=(0, 1),
                         unpack=True)
lgm, dlgm = np.loadtxt(extra_path+'gama_prior_lgm.txt', usecols=(0, 1),
                       unpack=True)
print('2D Histrograms with:\n - {} Mass bins\n - {} sSFR bins'.format(
    lgm.size, ssfr.size))

lin_ssfr = 10**ssfr

print('Computing conditional probability dp(ssfr|M)dssfr')
norm_lgm_ssfr_vmax = np.sum(all_lgm_ssfr_pdf_vmax*dssfr[:, np.newaxis], axis=0)
norm_gama_lgm_ssfr_vmax = np.sum(
    gama_all_lgm_ssfr_pdf_vmax*dssfr[:, np.newaxis], axis=0)
all_lgm_ssfr_pdf_vmax /= norm_lgm_ssfr_vmax[np.newaxis, :]
gama_all_lgm_ssfr_pdf_vmax /= norm_gama_lgm_ssfr_vmax[np.newaxis, :]

# Cumulative PDF
F_lgm_ssfr_vmax = np.cumsum(all_lgm_ssfr_pdf_vmax*dssfr[:, np.newaxis], axis=0)
gama_F_lgm_ssfr_vmax = np.cumsum(gama_all_lgm_ssfr_pdf_vmax*dssfr[:, np.newaxis], axis=0)


sdss_qg_fraction = []
gama_qg_fraction = []

for i in range(len(lgm)):
    sdss_cum_pdf = F_lgm_ssfr_vmax[:, i]
    obs_quench = np.interp(-11, ssfr, sdss_cum_pdf)
    sdss_qg_fraction.append(obs_quench)

    gama_cum_pdf = gama_F_lgm_ssfr_vmax[:, i]
    obs_quench = np.interp(-11, ssfr, gama_cum_pdf)
    gama_qg_fraction.append(obs_quench)


# =============================================================================
# %% Quenched/passive fraction (FIG. 7)
# =============================================================================

fig, axs = plt.subplots(ncols=3, nrows=3, figsize=(6, 6),
                        sharex=True, sharey=True,
                        gridspec_kw={'wspace': 0.1, 'hspace': 0.1})

all_axs = axs.flatten()

for ith, data_i in enumerate(data_sets_names.keys()):
    print('Data set ', ith+1)

    # ax = axs[0]
    simu_i = DataSet('simu_data/'+data_i)

    simu_i.mass_filter(min_mass=10**8.5, max_mass=10**11.5)
    simu_i.replace_ssfr0(5e-15)

    M = simu_i.mass
    sSFR = simu_i.ssfr

    qg = np.where(sSFR == 5e-15)[0]
    M_q = M[qg]

    mass_bin_edges = np.linspace(8.4, 11.6, min(50, sSFR.size//50))
    mass_bins = (mass_bin_edges[:-1]+mass_bin_edges[1:])/2

    h_mass, _ = np.histogram(np.log10(M), bins=mass_bin_edges)
    h_q, _ = np.histogram(np.log10(M_q), bins=mass_bin_edges)

    H, xbins, ybins = np.histogram2d(
        np.log10(M), np.log10(sSFR),
        range=[[8.4, 11.6], [-16, -6]], bins=15)

    H = H/np.sum(H*np.diff(ybins)[np.newaxis, :], axis=1)[:, np.newaxis]
    H_cond = np.cumsum(H*np.diff(ybins)[np.newaxis, :], axis=1)

    xbins = (xbins[1:]+xbins[:-1])/2
    ybins = (ybins[1:]+ybins[:-1])/2

    obs_qg_frac = interp1d(ybins, H_cond, axis=1)(-11)

    # PLOT

    ax = all_axs[ith]

    ax.fill_between(xbins, obs_qg_frac, alpha=.5,
                    label=data_sets_names[data_i],
                    color=simcolors[ith],
                    hatch="\\", edgecolor='k')

    ax.fill_between(mass_bins, h_q/h_mass, color=simcolors[ith], alpha=1,
                    hatch="/", edgecolor='k')
    ax.plot(lgm, gama_qg_fraction, 'k--', lw=3)
    ax.plot(lgm, sdss_qg_fraction, 'k', lw=3)
    ax.set_xlim(8.5, 11.6)
    ax.tick_params(direction='in', top=True, right=True)
    ax.annotate(data_sets_names[data_i], xy=(.05, .85), xycoords='axes fraction')
    if ith >= 6:
        ax.set_xlabel(r'$\log(M/M_\odot)$')


ax = ax = all_axs[1]

legend_elements = [Line2D([0], [0], color='k', lw=4, label='SDSS (passive)'),
                   Line2D([0], [0], ls='--', color='k', label='GAMA (passive)',
                          markerfacecolor='g', markersize=15),
                   Patch(facecolor='white', edgecolor='k', hatch="\\\\\\\\",
                         label=r'$sSFR\leq 10^{-11}$ yr$^{-1}$ (Passive)'),
                   Patch(facecolor='white', edgecolor='k', hatch="////",
                         label=r'$sSFR=0$ (Quenched)')]
ax.legend(handles=legend_elements, loc=(-0.95, 1.1), ncol=2,
          framealpha=1, edgecolor='k')
fig.savefig('figures/quenched_fraction_obs_approach.png', bbox_inches='tight')

# ...
