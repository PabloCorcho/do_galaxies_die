#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Plot Fig. 4a from Corcho-Caballero et al. (2021)."""

import numpy as np
from matplotlib import pyplot as plt

from astropy.io import fits
import pandas as pd

# =============================================================================
# LOADING DATA
# =============================================================================

data_path = 'obs_data/SDSS_z01_PabloCorcho.fit'
abs_photometry_path = 'obs_data/SDSS_absvalues.csv'

data = fits.open(data_path)
photo = pd.read_csv(abs_photometry_path)

lgm_p50 = data[1].data['lgm_tot_p50']
mass_mask = np.where(lgm_p50 >= 8.5)[0]

photo = pd.read_csv(abs_photometry_path)

vmax = photo['Vmax'].values
w = 1/vmax
w = w[mass_mask]

ssfr_p2p5 = data[1].data['specsfr_tot_p2p5'][mass_mask]
ssfr_p16 = data[1].data['specsfr_tot_p16'][mass_mask]
ssfr_p50 = data[1].data['specsfr_tot_p50'][mass_mask]
ssfr_p84 = data[1].data['specsfr_tot_p84'][mass_mask]
ssfr_p97p5 = data[1].data['specsfr_tot_p97p5'][mass_mask]

sfr_p2p5 = data[1].data['sfr_tot_p2p5'][mass_mask]
sfr_p16 = data[1].data['sfr_tot_p16'][mass_mask]
sfr_p50 = data[1].data['sfr_tot_p50'][mass_mask]
sfr_p84 = data[1].data['sfr_tot_p84'][mass_mask]
sfr_p97p5 = data[1].data['sfr_tot_p97p5'][mass_mask]

lgm_p2p5 = data[1].data['lgm_tot_p2p5'][mass_mask]
lgm_p16 = data[1].data['lgm_tot_p16'][mass_mask]
lgm_p50 = data[1].data['lgm_tot_p50'][mass_mask]
lgm_p84 = data[1].data['lgm_tot_p84'][mass_mask]
lgm_p97p5 = data[1].data['lgm_tot_p97p5'][mass_mask]

# %% SDSS UNCERTAINTIES (Fig. 4-a)

left_sigma = ssfr_p50-ssfr_p16
right_sigma = ssfr_p84-ssfr_p50

mean_sigma = (left_sigma + right_sigma)/2

good_ssfr = (ssfr_p50>-15)&(ssfr_p50<-6)
mask1 = (left_sigma<0.3)&(right_sigma<0.3)&(np.abs(left_sigma-right_sigma)<0.2)&good_ssfr
mask2 = (left_sigma<0.5)&(right_sigma<0.5)&(np.abs(left_sigma-right_sigma)<0.2)&good_ssfr

H, xedges, yedges = np.histogram2d(left_sigma, right_sigma, bins=20, 
                                   range=[[-1,1.5],[-1,1.5]],
                                   density=True)
xbins = (xedges[:-1]+xedges[1:])/2
ybins = (yedges[:-1]+yedges[1:])/2
H = (H-H.min())/(H.max()-H.min())

h, sigma_edges = np.histogram(mean_sigma, bins=30)

cumulative_frac = np.cumsum(h[::-1])

cumulative_frac= cumulative_frac/cumulative_frac[-1]
sigma_bins = (sigma_edges[:-1]+sigma_edges[1:])/2

quenched_frac = np.interp(0.3, sigma_bins, cumulative_frac[::-1])

fig = plt.figure(figsize=(5,4))
ax = fig.add_subplot(111)

mappable = ax.scatter(left_sigma, right_sigma, s=3, c=ssfr_p50, vmin=-13,
           vmax=-8, cmap='jet_r', alpha=0.1)
ax.contourf(xbins, ybins, H.T, levels=[0.05, 1], colors='k', alpha=0.2,
            # hatches=['\\']
            )
cbar = plt.colorbar(mappable=mappable, ax=ax, label=r'$\log(sSFR_{50})$')
cbar.set_alpha(1)
cbar.draw_all()

ax.annotate('a',xy=(-.15, 1.05), xycoords='axes fraction',
            ha='left',
            color='k',
            fontweight='bold',
            fontsize=15)

ax.plot([0, 1.5], [0, 1.5], 'k--')
epsilon = np.linspace(0,1)
ax.plot(-np.log10(1-epsilon), np.log10(1+epsilon), 'k:')
ax.set_xlabel(r'$\log(sSFR_{50})-\log(sSFR_{16})$')
ax.set_ylabel(r'$\log(sSFR_{84})-\log(sSFR_{50})$')
ax.grid(b=True)
ax.set_xlim(0,1.5)
ax.set_ylim(0,1.5)
plt.savefig('figures/sdss_uncertainties.png', bbox_inches='tight')

# ...