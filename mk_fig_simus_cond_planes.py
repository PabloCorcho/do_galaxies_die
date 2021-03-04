#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 11:19:31 2021

@author: pablo
"""

import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits


"""
This script computes the plane dlog(ssfr)(ssfr|M)/dlog(ssfr) both density and
cumulative. 

"""

# =============================================================================
# SDSS & GAMA probability planes
# =============================================================================

extra_path = 'obs_data/derived_data/'
gama_all_lgm_ssfr_pdf_vmax = np.loadtxt(extra_path+'gama_fundamentalplane_vmax.txt') 
all_lgm_ssfr_pdf_vmax = np.loadtxt(extra_path+'sdss_fundamentalplane_vmax.txt')

ssfr, dssfr = np.loadtxt(extra_path+'gama_prior_ssfr.txt', usecols=(0,1),
                         unpack=True) 
lgm, dlgm = np.loadtxt(extra_path+'gama_prior_lgm.txt', usecols=(0,1),
                         unpack=True) 
print('2D Histrograms with:\n - {} Mass bins\n - {} sSFR bins'.format(lgm.size, ssfr.size))


print('Computing conditional probability dp(ssfr|M)dssfr') 
norm_lgm_ssfr_vmax = np.sum(all_lgm_ssfr_pdf_vmax*dssfr[:, np.newaxis], axis=0)    
norm_gama_lgm_ssfr_vmax = np.sum(gama_all_lgm_ssfr_pdf_vmax*dssfr[:, np.newaxis], axis=0)    
all_lgm_ssfr_pdf_vmax /= norm_lgm_ssfr_vmax[np.newaxis, :]
gama_all_lgm_ssfr_pdf_vmax /= norm_gama_lgm_ssfr_vmax[np.newaxis, :]

# Cumulative PDF
F_lgm_ssfr_vmax = np.cumsum(all_lgm_ssfr_pdf_vmax*dssfr[:, np.newaxis], axis=0)
gama_F_lgm_ssfr_vmax = np.cumsum(gama_all_lgm_ssfr_pdf_vmax*dssfr[:, np.newaxis], axis=0)


# %%
# =============================================================================
# LOAD SIMULATIONS
# =============================================================================

data_sets_names = {
                     'MAGNETICUM.fits':'MAG500',
                     'tng300-1.fits':'IllustrisTNG300',
                    'tng100-1.fits':'IllustrisTNG100',                    
                    'SIMBA_m100n1024_151.fits':'SIMBA150',                    
                    'SIMBA_m50n512_151.fits':'SIMBA75',
                    'SIMBA_m25n512_151.fits':'SIMBA35',
                    'RefL0100N1504.fits':'EAGLE100',
                    'RefL0050N0752.fits':'EAGLE50',                    
                    'RefL0025N0752.fits':'EAGLE25'
                  }


simcolors = ['darkgreen', 'darkgreen', 'darkgreen',
             'darkorange', 'darkorange', 'darkorange',
             'crimson', 'crimson',
              'crimson']

sim_ls = ['-', '--', ':', '-', '--', ':', '-', '--', ':']

simcolors = simcolors[:len(data_sets_names)]
# =============================================================================
# %% CONDITIONAL AND DENSITY PLANES (FIG. 1 & FIG. 2)
# =============================================================================

# conditional distribution in mass bins figure
fig1, axs1 = plt.subplots(nrows=3, ncols=3, figsize=(6,6), 
                          gridspec_kw={'hspace':0.1, 'wspace':0.1})
fig2, axs2 = plt.subplots(nrows=3, ncols=3, figsize=(6,6), 
                          gridspec_kw={'hspace':0.1, 'wspace':0.1})



ssfr_bin_edges = np.arange(-16, -7, .2)
ssfr_bins = (ssfr_bin_edges[:-1] + ssfr_bin_edges[1:])/2
lgm_bin_edges = np.linspace(8.5, 11.5, 10)
lgm_bins = (lgm_bin_edges[1:]+lgm_bin_edges[:-1])/2


conturf_colors = ['indianred', 'darkorange', 'blue', 'mediumpurple']

fig3, axs3 = plt.subplots(nrows=3, ncols=3, sharex=True, sharey=True,
                          figsize=(8, 6))
axs3 = axs3.flatten()


# positions =[0,1,2,]

for ith, data_i in enumerate(data_sets_names.keys()):
    
    hdul = fits.open('simu_data/simu_derived_data/'+data_sets_names[data_i]+
                     '_planes.fits')
    
    density_plane = hdul[2].data
    cumulative_plane = hdul[3].data
    cumulative_plane_filtered = hdul[4].data
    
    M = hdul[1].data['mass'] # masses of all galaxies
    # masses of galaxies with sfr>0
    M_filtered = hdul[1].data['mass_filtered'][:cumulative_plane_filtered.shape[0]]
    
    ssfr_bins = hdul[1].data['ssfr_bins'][:density_plane.shape[1]]
    
    print('Data set ', ith+1)
    print('--> '+data_sets_names[data_i], M.size, M_filtered.size)

    
    # PLOT                                                                                                                                                                                                      
    
    ax = axs1.flatten()[ith]
    mappable = ax.contourf(M, ssfr_bins, 
                           cumulative_plane.T,
                           levels=[0.05, 0.16, 0.5, 0.84, 0.95], 
                           colors=conturf_colors)
    
    
    CS = ax.contour(lgm, ssfr, F_lgm_ssfr_vmax, levels=[0.05, 0.16, 0.5, 0.84, 0.95], 
            colors='k', linewidths=[1, 1.5, 3.5, 1.5, 1], linestyles=['solid'])
    
    # ax.clabel(CS, inline=True, fontsize=7)
    
    CS = ax.contour(lgm, ssfr, gama_F_lgm_ssfr_vmax, levels=[0.05, 0.16, 0.5, 0.84, 0.95], 
            colors='k', linewidths=[1, 1.5,3.5,1.5, 1], linestyles=['dashed'])
    
    
    ax.annotate(data_sets_names[data_i], xy=(.05, .05), 
                xycoords='axes fraction', va='bottom',
                ha='left', fontsize=12)
    
    # ax.clabel(CS, inline=True, fontsize=7)
    ax.set_ylim(-16., -8.5)
    ax.set_xlim(8.5, 11.5)
    ax.set_xticks([9,10,11])
    
    # plt.colorbar(mappable, ax=ax)
    
    ax.tick_params(direction='in', top=True, right=True, labelleft=False,
                   labelbottom=False, width=1, length=6)
    if (ith==0)|(ith==3)|(ith==6):
        ax.tick_params(labelleft=True)
        ax.set_ylabel(r'$\log(sSFR/yr)$')
    if (ith==6)|(ith==7)|(ith==8):
        ax.tick_params(labelbottom=True)
        ax.set_xlabel(r'$\log(M/M_\odot)$')
        
    # ax.annotate(data_sets_names[data_i], xy=(.96, .94), xycoords='axes fraction', va='top',
    #             ha='right', fontsize=12)
   
    
    # conditional figure -----------------------------------------------------
    
    for jth in range(lgm_bins.size):
        mask = (M>lgm_bin_edges[jth])&(M<lgm_bin_edges[jth+1])
        
        ax = axs3[jth]
        
        mask_obs = (lgm>lgm_bin_edges[jth])&(lgm<lgm_bin_edges[jth+1])

        ax.fill_between(x=ssfr, y1=np.nanmean(all_lgm_ssfr_pdf_vmax[:, mask_obs], axis=1),
                        y2=0, hatch='/', alpha=0.1/len(data_sets_names), color='k')
        
        ax.fill_between(x=ssfr, y1=np.nanmean(gama_all_lgm_ssfr_pdf_vmax[:, mask_obs], axis=1),
                        y2=0, hatch='o', alpha=0.1/len(data_sets_names), color='k')
        
        ax.plot(ssfr_bins, np.nanmean(density_plane[mask, :], axis=0),
                color=simcolors[ith], label=data_sets_names[data_i], 
                ls=sim_ls[ith], alpha=0.95, lw=2)
        
        
        ax.set_yscale('symlog', linthresh=0.001)
        ax.set_ylim(0, np.nanmax(density_plane))
        ax.tick_params(direction='in', top=True, right=True)
        
        if data_sets_names[data_i] == 'EAGLE25':
            ax.annotate( '{:04.2f}'.format(lgm_bins[jth]), xy=(.96,.90),
                        xycoords='axes fraction', ha='right', va='top')                        
            if jth ==0:
                ax.set_ylabel(r'$\frac{dp(log(sSFR)|M)}{dlog(sSFR)}$', fontsize=14)
                
            if (jth==6)|(jth==7)|(jth==8):
                ax.set_xlabel(r'$\log(sSFR/yr)$', fontsize=14)
                
        ax.set_yticks([1e-3, 1e-2, 1e-1, 1])
        ax.set_xticks([-16, -14, -12, -10, -8])
        # ax.set_yticklabels(['0.001', '0.01', '0.1', '1'])
    if data_sets_names[data_i] == 'EAGLE25':
        axs3[1].legend(loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.8),
                        fontsize=12, framealpha=1, edgecolor='k')
    #------------------------------------------------------------------------
    
    # Filtered (sfr=0) plane
    ax = axs2.flatten()[ith]
   
    mappable = ax.contourf(M_filtered, ssfr_bins, 
                           cumulative_plane_filtered.T,
                           levels=[0.05, 0.16, 0.5, 0.84, 0.95], 
                           colors=conturf_colors)
    
    
    
    CS = ax.contour(lgm, ssfr, F_lgm_ssfr_vmax, levels=[0.05, 0.16, 0.5, 0.84, 0.95], 
            colors='k', linewidths=[1, 1.5,3.5,1.5, 1], linestyles=['solid'])
    # ax.clabel(CS, inline=True, fontsize=7)
    
    CS = ax.contour(lgm, ssfr, gama_F_lgm_ssfr_vmax, levels=[0.05, 0.16, 0.5, 0.84, 0.95], 
            colors='k', linewidths=[1, 1.5,3.5,1.5, 1], linestyles=['dashed'])
    
    
    # ax.clabel(CS, inline=True, fontsize=7)
    
    ax.set_xlim(8.5, 11.5)
    ax.set_ylim(-16., -8.5)
    ax.set_xticks([9,10,11])
    
    
    ax.tick_params(direction='in', top=True, right=True, labelleft=False,
                   labelbottom=False, width=1, length=6)
    # if positions[ith][1][0]==0:
    #     ax.tick_params(labelleft=True)
    #     ax.set_ylabel(r'$\log(sSFR/yr)$')
        
    ax.annotate(data_sets_names[data_i], xy=(.05, .05), 
                xycoords='axes fraction', va='bottom',
                ha='left', fontsize=12)
    
   
    if (ith==0)|(ith==3)|(ith==6):
        ax.tick_params(labelleft=True)
        ax.set_ylabel(r'$\log(sSFR/yr)$')
    if (ith==6)|(ith==7)|(ith==8):
        ax.tick_params(labelbottom=True)
        ax.set_xlabel(r'$\log(M/M_\odot)$')
   
ax = fig1.add_subplot(111)
cbar_ax = ax.inset_axes([0.15, 1.10, 0.7, 0.02])
ax.axis('off')
cb = plt.colorbar(mappable, cax=cbar_ax, orientation='horizontal')
cb.ax.tick_params(labelsize=14)
cb.set_label(label='Cumulative probability', labelpad=-50)
ax.annotate('a', xy=(-0.1, 1.05), xycoords='axes fraction', weight='bold',
            fontsize=16)
ax = fig2.add_subplot(111)
cbar_ax = ax.inset_axes([0.15, 1.10, 0.7, 0.02])
ax.axis('off')
cb = plt.colorbar(mappable, cax=cbar_ax, orientation='horizontal')
cb.ax.tick_params(labelsize=14)
cb.set_label(label='Cumulative probability', labelpad=-50)
ax.annotate('b', xy=(-0.1, 1.05), xycoords='axes fraction', weight='bold',
            fontsize=16)

    
print('Process completed')    

fig3.subplots_adjust(wspace=0.05, hspace=0.05)

fig1.savefig('figures/cprob_all_simus.png', bbox_inches='tight')
fig2.savefig('figures/cprob_all_simus_flitered.png', bbox_inches='tight')

fig3.savefig('figures/cprob_mass_bins.png', bbox_inches='tight')

