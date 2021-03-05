#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Plot Figs. 3 and 8 from Corcho-Caballero et al. (2021)."""

from __future__ import print_function, division
from read_suite import DataSet
import numpy as np
from matplotlib import pyplot as plt


# %% ==========================================================================
#  Survey detection (FIG. 8)
# =============================================================================

ssfr_distri_gama = np.loadtxt('obs_data/derived_data/gama_ssfr_pdf_vmax_redshift.txt')
ssfr_distri_gama_active = np.loadtxt('obs_data/derived_data/gama_ssfr_pdf_vmax_redshift_active.txt')
ssfr_distri_gama_passive = np.loadtxt('obs_data/derived_data/gama_ssfr_pdf_vmax_redshift_passive.txt')

ssfr_gama, dssfr = np.loadtxt('obs_data/derived_data/gama_prior_ssfr.txt',
                              usecols=(0, 1), unpack=True)

data_sets_names = {
                    'tng300-1.fits': 'IllustrisTNG300',
                    'tng100-1.fits': 'IllustrisTNG100',
                    'RefL0100N1504.fits': 'EAGLE100',
                    'RefL0050N0752.fits': 'EAGLE50',
                    'RefL0025N0752.fits': 'EAGLE25'
                  }

simcolors = ['g', 'lime', 'purple', 'crimson', 'fuchsia']


fig = plt.figure(figsize=(6, 6))
ax = fig.add_subplot(111)

for ith, data_i in enumerate(data_sets_names.keys()):
    print('Data set ', ith+1)

    simu_i = DataSet('simu_data/'+data_i)

    simu_i.mass_filter(min_mass=10**8.5, max_mass=10**11.5)
    simu_i.replace_ssfr0(1e-30)

    simu_i.get_photometry()
    u = simu_i.u[simu_i.mass_mask]
    g = simu_i.r[simu_i.mass_mask]
    r = simu_i.r[simu_i.mass_mask]

    d_u_max_sdss = 10**((19-u)/5+1)/1e6
    d_g_max_sdss = 10**((19-g)/5+1)/1e6
    d_r_max_sdss = 10**((17.7-r)/5+1)/1e6

    d_max_sdss = np.min([d_u_max_sdss, d_g_max_sdss, d_r_max_sdss], axis=0)

    z_max_sdss = 70*d_max_sdss/3e5

    d_u_max_gama = 10**((23-u)/5+1)/1e6
    d_g_max_gama = 10**((23-g)/5+1)/1e6
    d_r_max_gama = 10**((19.8-r)/5+1)/1e6

    d_max_gama = np.min([d_u_max_gama, d_g_max_gama, d_r_max_gama], axis=0)
    z_max_gama = 70*d_max_gama/3e5

    M = simu_i.mass
    sSFR = simu_i.ssfr

    ssfr_mask = sSFR == 1e-30

    H_sdss, z_bins = np.histogram(z_max_sdss[ssfr_mask])
    z_bins_sdss = (z_bins[:-1] + z_bins[1:]) / 2
    cum_H_sdss = np.cumsum(H_sdss[::-1])
    cum_H_sdss = cum_H_sdss/cum_H_sdss[-1]

    H_gama, z_bins = np.histogram(z_max_gama[ssfr_mask])
    z_bins_gama = (z_bins[:-1] + z_bins[1:]) / 2
    cum_H_gama = np.cumsum(H_gama[::-1])
    cum_H_gama = cum_H_gama/cum_H_gama[-1]

    ax.plot(z_bins_sdss[::-1], cum_H_sdss,
            label=data_sets_names[data_i]+' (SDSS)', ls='-',
            color=simcolors[ith], lw=3)

    ax.plot(z_bins_gama[::-1], cum_H_gama,
            label=data_sets_names[data_i]+' (GAMA)', ls='--',
            color=simcolors[ith], lw=3)
    ax.set_yscale('log')
    ax.set_xlim(0.1, 0)

ax.tick_params(which='both', direction='in', right=True, top=True)
ax.legend(framealpha=1, edgecolor='k', fontsize=9)

ax.grid(b=True)
ax.set_ylim(1e-3, 2)
ax.set_ylabel(r'Cumulative fraction of quenched galaxies')
ax.set_xlabel(r'z$_{max}$ for solid detection')
fig.savefig('figures/detected_quench_galaxies.png', bbox_inches='tight')

# =============================================================================
# %% sSFR distribution (FIG 3.)
# =============================================================================
simcolors = [
                'darkgreen',
                'darkgreen',
                'crimson',
                'crimson',
                'crimson']

fig = plt.figure(figsize=(6, 6))
ax = fig.add_subplot(212)

ax.fill_between(ssfr_gama, ssfr_distri_gama_active, color='k', alpha=0.3,
                label='GAMA (active)')
ax.fill_between(ssfr_gama, ssfr_distri_gama_passive, color='k', alpha=0.3,
                lw=4, label='GAMA (passive)')

ax.annotate('b', xy=(-.15, 0.9), xycoords='axes fraction',
            ha='left',
            color='k',
            fontweight='bold',
            fontsize=15)

fig2 = plt.figure(figsize=(6, 6))
ax2 = fig.add_subplot(211)
ax2.fill_between(ssfr_gama, np.cumsum(dssfr*ssfr_distri_gama),
                 color='k', alpha=0.3,
                 label='GAMA')
ax2.annotate('a', xy=(-.15, 0.9), xycoords='axes fraction',
             ha='left',
             color='k',
             fontweight='bold',
             fontsize=15)

data_sets_names = {
                    'tng300-1.fits': 'IllustrisTNG300',
                    'tng100-1.fits': 'IllustrisTNG100',
                    'RefL0100N1504.fits': 'EAGLE100',
                    'RefL0050N0752.fits': 'EAGLE50',
                    'RefL0025N0752.fits': 'EAGLE25'
                  }
simcolors = [
                'darkgreen',
                'darkgreen',
                'crimson',
                'crimson',
                'crimson'
             ]
sims_ls = ['--', ':', '-', '--', ':']


ax2.minorticks_on()
ax.minorticks_on()

for ith, data_i in enumerate(data_sets_names.keys()):
    print('Data set ', ith+1)

    # ax = axs[0]
    simu_i = DataSet('simu_data/'+data_i)

    simu_i.mass_filter(min_mass=10**8.5, max_mass=10**11.5)
    simu_i.replace_ssfr0(1e-16)

    M = simu_i.mass
    sSFR = simu_i.ssfr

    H, xbins = np.histogram(np.log10(sSFR),
                            range=[-16, -7], bins=50, density=True)
    xbin = (xbins[1:]+xbins[:-1])/2

    ax.plot(xbin, H, color=simcolors[ith], lw=3, label=data_sets_names[data_i],
            ls=sims_ls[ith]
            )

    ax2.plot(xbin, np.cumsum(H*np.diff(xbins)), color=simcolors[ith],
             lw=3, ls=sims_ls[ith],
             label=data_sets_names[data_i])

ax.tick_params(which='both', direction='in', right=True, top=True)
ax.set_ylim(1e-3, 2)
ax.set_yscale('log')
ax.set_xlabel(r'$\log_{10}(sSFR)$ [1/yr]')
ax.set_ylabel(r'$\frac{dp(sSFR)}{d\log_{10}(sSFR)}$', fontsize=12)
ax.legend(loc=(.1, .5), fontsize=7, framealpha=0)
ax.set_xlim(-16.1, -7.5)
ax2.set_xlim(-16.1, -7.5)
ax2.tick_params(labelbottom=False)

ax2.set_ylabel(r'Cumulaive probability')
ax2.set_ylim(0, 1.05)
ax2.tick_params(which='both', direction='in', right=True, top=True)
ax2.legend(loc='lower right', fontsize=7, framealpha=0)
fig.subplots_adjust(hspace=0)
fig.savefig('figures/ssfr_distrib_z005.png', bbox_inches='tight')

# ...
