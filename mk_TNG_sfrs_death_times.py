#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Plot Figs. 5 and 6 from Corcho-Caballero et al. (2021)."""

from astropy.io import fits
import numpy as np
from matplotlib import pyplot as plt

from matplotlib.lines import Line2D
from matplotlib.patches import Patch


simus = ['TNG100-1', 'TNG300-1']

for simu in simus:

    data_table = 'simu_data/'+simu+'_sfrs.fits'

    hdul = fits.open(data_table)[1]

    ids = hdul.data['ID']
    instantaneous_sfr = hdul.data['inst_sfr']
    instantaneous_sfr_2r = hdul.data['inst_sfr_tworad']
    sfr_100 = hdul.data['100myr_sfr']
    sfr_300 = hdul.data['300myr_sfr']
    sfr_1000 = hdul.data['1000myr_sfr']
    sfr_100_2r = hdul.data['100myr_sfr_tworad']
    sfr_300_2r = hdul.data['300myr_sfr_tworad']
    sfr_1000_2r = hdul.data['1000myr_sfr_tworad']
    dead_time = hdul.data['dead_time']

    total_stellar_mass = hdul.data['total_stellar_mass']
    lgm = np.log10(total_stellar_mass)
    total_stellar_mass_2r = hdul.data['total_stellar_mass_tworad']
    lgm_2r = np.log10(total_stellar_mass_2r)

    qg_inst = np.where(instantaneous_sfr == 0)[0]
    qg_inst_2r = np.where(instantaneous_sfr_2r == 0)[0]
    qg_100 = np.where(sfr_100 == 0)[0]
    qg_300 = np.where(sfr_300 == 0)[0]
    qg_1000 = np.where(sfr_1000 == 0)[0]
    qg_100_2r = np.where(sfr_100_2r == 0)[0]
    qg_300_2r = np.where(sfr_300_2r == 0)[0]
    qg_1000_2r = np.where(sfr_1000_2r == 0)[0]


# =============================================================================
# DEATH TIMES (FIG. 6)
# =============================================================================

    plt.figure(figsize=(5, 4))
    plt.hist(dead_time[qg_inst], bins=np.linspace(0, 13.7),
             log=False,
             cumulative=True,
             histtype='step', color='k', density=True,
             lw=2,
             label=r'$SFR_{inst}$')
    plt.hist(dead_time[qg_inst_2r], bins=np.linspace(0, 13.7),
             log=False,
             cumulative=True,
             histtype='step', ls='--', color='k', density=True,
             lw=2,
             label=r'$SFR_{inst}$')
    plt.hist(dead_time[qg_100], bins=np.linspace(0, 13.7),
             log=False,
             cumulative=True,
             histtype='step', color='blue', density=True,
             lw=2,
             label=r'$SFR_{100}$')
    plt.hist(dead_time[qg_100_2r], bins=np.linspace(0, 13.7),
             log=False,
             cumulative=True,
             histtype='step', ls='--', color='blue', density=True,
             lw=2,
             label=r'$SFR_{100}$')
    plt.hist(dead_time[qg_300], bins=np.linspace(0, 13.7),
             log=False,
             cumulative=True,
             histtype='step', color='green', density=True,
             lw=2,
             label=r'$SFR_{100}$')
    plt.hist(dead_time[qg_300_2r], bins=np.linspace(0, 13.7),
             log=False,
             cumulative=True,
             histtype='step', ls='--', color='green', density=True,
             lw=2,
             label=r'$SFR_{100}$')
    plt.hist(dead_time[qg_1000], bins=np.linspace(0, 13.7),
             log=False,
             cumulative=True,
             histtype='step', color='red', density=True,
             lw=2,
             label=r'$SFR_{1000}$')
    plt.hist(dead_time[qg_1000_2r], bins=np.linspace(0, 13.7),
             log=False,
             cumulative=True,
             histtype='step', ls='--', color='red', density=True,
             lw=2,
             label=r'$SFR_{1000}$')

    # panel = {'TNG100-1': 'a', 'TNG300-1': 'b'}
    # plt.annotate(panel[simu], xy=(-.2, .95), xycoords='axes fraction',
    #              weight="bold", fontsize=15)

    plt.minorticks_on()
    plt.tick_params(which='both', direction='in', top=True, right=True)
    plt.ylim(1e-3, 1.1)
    plt.xlim(1e-1, 13)
    plt.xlabel(r'$\tau_{death}$ [Gyr]', fontsize=12)
    plt.ylabel('Cumulative fraction', fontsize=12)

    custom_lines = [Patch(facecolor='k', edgecolor='k',
                          label=r'$SFR_{inst}$'),
                    Patch(facecolor='b', edgecolor='b',
                          label=r'$SFR_{100}$'),
                    Patch(facecolor='g', edgecolor='g',
                          label=r'$SFR_{300}$'),
                    Patch(facecolor='r', edgecolor='r',
                          label=r'$SFR_{1000}$'),
                    Line2D([0], [0], color='k', lw=2, label='All'),
                    Line2D([0], [0], color='k', ls='--', lw=2, label=r'$2R_e$')
                    ]

    plt.legend(handles=custom_lines, framealpha=0, loc='lower right')
    plt.savefig('figures/dead_times_'+simu+'.png', bbox_inches='tight')


# =============================================================================
# sSFR distrib (Fig. 5)
# =============================================================================

    ssfr_bins = np.hstack([0, np.logspace(-14., -7, 50)])
    plt.figure(figsize=(6, 6))

    plt.hist(instantaneous_sfr/total_stellar_mass, log=True, bins=ssfr_bins,
             color='k',
             histtype='step', lw=2, label=r's$SFR_{inst}$', alpha=1)
    plt.hist(instantaneous_sfr_2r/total_stellar_mass_2r, log=True, bins=ssfr_bins,
             color='k',
             histtype='step', lw=2, ls='--', label=r'$sSFR_{inst-2R_e}$', alpha=1)

    plt.hist(sfr_100/total_stellar_mass, log=True, bins=ssfr_bins,
             color='blue',
             histtype='step', lw=2, label=r'$sSFR_{100}$', alpha=1)
    plt.hist(sfr_100_2r/total_stellar_mass_2r, log=True, bins=ssfr_bins,
             color='blue',
             histtype='step', lw=2, ls='--', label=r'$sSFR_{100-2R_e}$', alpha=1)

    plt.hist(sfr_300/total_stellar_mass, log=True, bins=ssfr_bins,
             color='green',
             histtype='step', lw=2, label=r'$sSFR_{300}$', alpha=1)
    plt.hist(sfr_300_2r/total_stellar_mass_2r, log=True, bins=ssfr_bins,
             color='green',
             histtype='step', lw=2, ls='--', label=r'$sSFR_{300}$', alpha=1)

    plt.hist(sfr_1000/total_stellar_mass, log=True, bins=ssfr_bins,
             color='red',
             histtype='step', lw=2, label=r'$sSFR_{1000}$', alpha=1)
    plt.hist(sfr_1000_2r/total_stellar_mass_2r, log=True, bins=ssfr_bins,
             color='red',
             histtype='step', lw=2, ls='--', label=r'$sSFR_{300}$', alpha=1)

    # panel = {'TNG100-1': 'a', 'TNG300-1': 'b'}
    # plt.annotate(panel[simu], xy=(-.15, .95), xycoords='axes fraction',
    #              weight="bold", fontsize=15)

    plt.xscale('symlog', linthresh=1e-14)
    plt.xlabel(r'$sSFR~[yr^{-1}]$')
    plt.ylabel(r'Galaxies per bin')
    plt.tick_params(which='both', direction='in', top=True, right=True)

    custom_lines = [Patch(facecolor='k', edgecolor='k',
                          label=r'$sSFR_{inst}$'),
                    Patch(facecolor='b', edgecolor='b',
                          label=r'$sSFR_{100}$'),
                    Patch(facecolor='g', edgecolor='g',
                          label=r'$sSFR_{300}$'),
                    Patch(facecolor='r', edgecolor='r',
                          label=r'$sSFR_{1000}$'),
                    Line2D([0], [0], color='k', lw=2, label='All'),
                    Line2D([0], [0], color='k', ls='--', lw=2, label=r'$2R_e$')
                    ]
    plt.xlim(0, 5e-7)
    plt.legend(handles=custom_lines, framealpha=0)
    plt.savefig('figures/ssfrs_distrib_illustris'+simu+'.png', bbox_inches='tight')


# ...
