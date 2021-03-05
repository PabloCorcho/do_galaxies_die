#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Plot Fig. 4b from Corcho-Caballero et al. (2021)."""

import numpy as np
from matplotlib import pyplot as plt

import pandas as pd

data_path = 'obs_data/gama_sample.csv'
data = pd.read_csv(data_path)

lgm_p50 = data['lgm_tot_p50']
mass_mask = np.where(lgm_p50 >= 8.5)[0]


w = data['weights'][mass_mask].values

ssfr_p2p5 = data['specsfr_tot_p2p5'][mass_mask].values
ssfr_p16 = data['specsfr_tot_p16'][mass_mask].values
ssfr_p50 = data['specsfr_tot_p50'][mass_mask].values
ssfr_p84 = data['specsfr_tot_p84'][mass_mask].values
ssfr_p97p5 = data['specsfr_tot_p97p5'][mass_mask].values

lgm_p2p5 = data['lgm_tot_p2p5'][mass_mask].values
lgm_p16 = data['lgm_tot_p16'][mass_mask].values
lgm_p50 = data['lgm_tot_p50'][mass_mask].values
lgm_p84 = data['lgm_tot_p84'][mass_mask].values
lgm_p97p5 = data['lgm_tot_p97p5'][mass_mask].values

sfr_p2p5 = ssfr_p2p5 + lgm_p2p5
sfr_p16 = ssfr_p16 + lgm_p16
sfr_p50 = ssfr_p50 + lgm_p50
sfr_p84 = ssfr_p84 + lgm_p84
sfr_p97p5 = ssfr_p97p5 + lgm_p97p5

sfr_pcnt = np.array([sfr_p2p5, sfr_p16, sfr_p50, sfr_p84, sfr_p97p5])
ssfr_pcnt = np.array([ssfr_p2p5, ssfr_p16, ssfr_p50, ssfr_p84, ssfr_p97p5])
lgm_pcnt = np.array([lgm_p2p5, lgm_p16, lgm_p50, lgm_p84, lgm_p97p5])

percent = np.array([2.5, 16, 50, 84, 97.5])/100
percentiles = np.linspace(0, 1, 20)

sfr_pcnt = np.array([sfr_p2p5, sfr_p16, sfr_p50, sfr_p84, sfr_p97p5])
ssfr_pcnt = np.array([ssfr_p2p5, ssfr_p16, ssfr_p50, ssfr_p84, ssfr_p97p5])
lgm_pcnt = np.array([lgm_p2p5, lgm_p16, lgm_p50, lgm_p84, lgm_p97p5])


# %% GAMA UNCERTAINTIES (Fig. 4-b)

left_sigma = ssfr_p50-ssfr_p16
right_sigma = ssfr_p84-ssfr_p50

mean_sigma = (left_sigma + right_sigma)/2

h, sigma_edges = np.histogram(mean_sigma, bins=30)

cumulative_frac = np.cumsum(h[::-1])

cumulative_frac= cumulative_frac/cumulative_frac[-1]
sigma_bins = (sigma_edges[:-1]+sigma_edges[1:])/2

quenched_frac = np.interp(0.3, sigma_bins, cumulative_frac[::-1])

H, xedges, yedges = np.histogram2d(left_sigma, right_sigma, bins=20, 
                                   range=[[-1,1.5],[-1,1.5]],
                                   density=True)
xbins = (xedges[:-1]+xedges[1:])/2
ybins = (yedges[:-1]+yedges[1:])/2
H = (H-H.min())/(H.max()-H.min())

fig = plt.figure(figsize=(5,4))
ax = fig.add_subplot(111)
# ax.contourf(xbins, ybins, H.T, levels=[0.05, 1], colors='k', alpha=0.1)
mappable = ax.scatter(left_sigma, right_sigma, s=3, c=ssfr_p50, vmin=-13,
           vmax=-8, cmap='jet_r', alpha=0.3)
ax.contourf(xbins, ybins, H.T, levels=[0.05, 1], colors='k', alpha=0.2,
            # hatches=['\\']
            )
cbar = plt.colorbar(mappable=mappable, ax=ax, label=r'$\log(sSFR_{50})$')
cbar.set_alpha(1)
cbar.draw_all()

ax.annotate('b',xy=(-.15, 1.05), xycoords='axes fraction',
            ha='left',
            color='k',
            fontweight='bold',
            fontsize=15)

ax.plot([0, 1.5], [0, 1.5], 'k--')
epsilon = np.linspace(0,1)
ax.plot(-np.log10(1-epsilon), np.log10(1+epsilon), 'k:')
# ax.plot([0, 1.5], [0.2, 1.7], 'k:')
# ax.plot([0, 1.5], [-.2, 1.3], 'k:')
ax.set_xlabel(r'$\log(sSFR_{50})-\log(sSFR_{16})$')
ax.set_ylabel(r'$\log(sSFR_{84})-\log(sSFR_{50})$')
ax.grid(b=True)
ax.set_xlim(0,1.5)
ax.set_ylim(0,1.5)
plt.savefig('figures/gama_uncertainties.png', bbox_inches='tight')

# ...