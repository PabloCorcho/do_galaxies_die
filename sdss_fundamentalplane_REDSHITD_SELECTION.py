#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 25 17:50:35 2020

@author: pablo
"""


import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
import pandas as pd


# =============================================================================
# This script computes the probability density distribution of log(M)-log(sSFR) 
# plane for SDSS data restricted to z<0.05
# =============================================================================

data_path = 'obs_data/SDSS_z01_PabloCorcho.fit'
abs_photometry_path = 'obs_data/SDSS_absvalues.csv'

data = fits.open(data_path)
photo = pd.read_csv(abs_photometry_path)

lgm_p50 = data[1].data['lgm_tot_p50']
mass_mask = np.where((lgm_p50>=8.5)&(lgm_p50<=11.5))[0]

photo = pd.read_csv(abs_photometry_path)

vmax = photo['Vmax'].values
w = 1/vmax
w = w[mass_mask]
phot_objid = photo['objid'][mass_mask]

objid = data[1].data['objid'][mass_mask]


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

Z = data[1].data['Z'][mass_mask]

sfr_pcnt = np.array([sfr_p2p5, sfr_p16, sfr_p50, sfr_p84, sfr_p97p5])
ssfr_pcnt = np.array([ssfr_p2p5, ssfr_p16, ssfr_p50, ssfr_p84, ssfr_p97p5])
lgm_pcnt = np.array([lgm_p2p5, lgm_p16, lgm_p50, lgm_p84, lgm_p97p5])

percent = np.array([2.5, 16, 50, 84, 97.5])/100


ssfrs = np.linspace(-15, -7, 81)
dssfr = np.diff(ssfrs)
ssfr = (ssfrs[:-1]+ssfrs[1:])/2    

sfrs = np.linspace(-6, 2, 81)
dsfr = np.diff(sfrs)
sfr = (sfrs[:-1]+sfrs[1:])/2    

lgms = np.linspace(7.5, 12.5, 51)
dlgm = np.diff(lgms)
lgm = (lgms[:-1]+lgms[1:])/2


all_sfr_pdf_vmax = np.zeros_like(sfr)
all_sfr_pdf = np.zeros_like(sfr)

all_ssfr_pdf_vmax = np.zeros_like(ssfr)
all_ssfr_pdf = np.zeros_like(ssfr)
all_ssfr_pdf_vmax_redshift = np.zeros_like(ssfr)

all_lgm_pdf_vmax = np.zeros_like(lgm)
all_lgm_pdf = np.zeros_like(lgm)

all_lgm_ssfr_pdf_vmax = np.zeros((ssfr.size, lgm.size))
all_lgm_ssfr_pdf = np.zeros((ssfr.size, lgm.size))

all_lgm_sfr_pdf_vmax = np.zeros((sfr.size, lgm.size))
all_lgm_sfr_pdf = np.zeros((sfr.size, lgm.size))

total_weight = 0
bad = 0

for i in range(len(ssfr_p50)):
# for i in range(10):
    print(i)
    
    cumpdf = np.interp(ssfrs, ssfr_pcnt[:, i], percent)    
    pdf = np.diff(cumpdf)
    dssfr = np.diff(ssfrs)
    ssfr_pdf = pdf/dssfr    
    ssfr_pdf /= np.sum(ssfr_pdf*dssfr) + 1e-30
    
    cumpdf = np.interp(lgms, lgm_pcnt[:, i], percent)    
    pdf = np.diff(cumpdf)
    dlgm = np.diff(lgms)
    lgm_pdf = pdf/dlgm    
    lgm_pdf /= np.sum(lgm_pdf*dlgm) + 1e-30
    
    cumpdf = np.interp(sfrs, sfr_pcnt[:, i], percent)    
    pdf = np.diff(cumpdf)
    dsfr = np.diff(sfrs)
    sfr_pdf = pdf/dsfr        
    all_sfr_pdf += sfr_pdf
    all_sfr_pdf_vmax += sfr_pdf*w[i]
    
    all_lgm_sfr_pdf += sfr_pdf[:, np.newaxis]*lgm_pdf[np.newaxis, :]
    all_lgm_sfr_pdf_vmax += sfr_pdf[:, np.newaxis]*lgm_pdf[np.newaxis, :]*w[i]
        
    if (np.sum(ssfr_pdf*dssfr)<0.98)|(np.sum(lgm_pdf*dlgm)<0.98):
        bad += 1
        continue
    
    
    if Z[i]<=0.05:
        all_ssfr_pdf_vmax_redshift += ssfr_pdf*w[i]
        total_weight += w[i]
    
print('Bad galaxies in the sample', bad)

plt.figure()
plt.semilogy(ssfr, all_ssfr_pdf_vmax_redshift/total_weight, c='k')
plt.hist(ssfr_p50, weights=w, density=True, range=[-15, -6], histtype='step', bins=30, color='r')
plt.hist(ssfr_p50, density=True, range=[-15, -6], histtype='step', bins=30, color='b')
plt.ylim(1e-3, 1)

np.savetxt('obs_data/derived_data/sdss_ssfr_pdf_vmax_redshift.txt', 
           all_ssfr_pdf_vmax_redshift/total_weight) 

# ...
