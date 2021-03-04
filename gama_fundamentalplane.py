import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
import pandas as pd


# =============================================================================
# This script computes the probability density distribution of log(M)-log(sSFR) 
# plane for GAMA data
# =============================================================================

data_path = 'obs_data/gama_sample.csv'
data = pd.read_csv(data_path)

lgm_p50 = data['lgm_tot_p50']
mass_mask = np.where((lgm_p50>=8.5)&(lgm_p50<=11.5))[0]

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

all_lgm_pdf_vmax = np.zeros_like(lgm)
all_lgm_pdf = np.zeros_like(lgm)

all_lgm_sfr_pdf_vmax = np.zeros((sfr.size, lgm.size))
all_lgm_sfr_pdf = np.zeros((sfr.size, lgm.size))

all_lgm_ssfr_pdf_vmax = np.zeros((ssfr.size, lgm.size))
all_lgm_ssfr_pdf = np.zeros((ssfr.size, lgm.size))

total_weight = 0
bad = 0
for i in range(len(ssfr_p50)):
    print(i)
    
    cumpdf = np.interp(ssfrs, ssfr_pcnt[:, i], percent)    
    pdf = np.diff(cumpdf)    
    ssfr_pdf = pdf/dssfr    
    ssfr_pdf /= np.sum(ssfr_pdf*dssfr) + 1e-30
    
    
            
    cumpdf = np.interp(lgms, lgm_pcnt[:, i], percent)    
    pdf = np.diff(cumpdf)
    lgm_pdf = pdf/dlgm    
    lgm_pdf /= np.sum(lgm_pdf*dlgm) + 1e-30
    
    cumpdf = np.interp(sfrs, sfr_pcnt[:, i], percent)    
    pdf = np.diff(cumpdf)    
    sfr_pdf = pdf/dsfr        
    all_sfr_pdf += sfr_pdf
    all_sfr_pdf_vmax += sfr_pdf*w[i]
    
    all_lgm_sfr_pdf += sfr_pdf[:, np.newaxis]*lgm_pdf[np.newaxis, :]
    all_lgm_sfr_pdf_vmax += sfr_pdf[:, np.newaxis]*lgm_pdf[np.newaxis, :]*w[i]
        
    if (np.sum(ssfr_pdf*dssfr)<0.98)|(np.sum(lgm_pdf*dlgm)<0.98):
        bad += 1
        continue
    
    total_weight += w[i]
    
    if lgm_pdf.all() is not np.nan:
        all_lgm_pdf += lgm_pdf
        all_lgm_pdf_vmax += lgm_pdf*w[i]            
    
        if ssfr_pdf.all() is not np.nan:        
            all_ssfr_pdf += ssfr_pdf        
            all_ssfr_pdf_vmax += ssfr_pdf*w[i]

            all_lgm_ssfr_pdf += ssfr_pdf[:, np.newaxis]*lgm_pdf[np.newaxis, :]
            all_lgm_ssfr_pdf_vmax += ssfr_pdf[:, np.newaxis]*lgm_pdf[np.newaxis, :]*w[i]
            
        else: 
            print('Bad sSFR')
    else: 
        print('Bad mass')

print('Bad galaxies in the sample', bad)


N_massbin, _ = np.histogram(lgm_p50, bins=lgms)


np.savetxt('obs_data/derived_data/gama_galaxies_permassbin.txt', 
           N_massbin)

np.savetxt('obs_data/derived_data/gama_fundamentalplane.txt', 
           all_lgm_ssfr_pdf) 
np.savetxt('obs_data/derived_data/gama_fundamentalplane_vmax.txt', 
           all_lgm_ssfr_pdf_vmax/total_weight) 
np.savetxt('obs_data/derived_data/gama_mass_pdf.txt', 
           all_lgm_pdf) 
np.savetxt('obs_data/derived_data/gama_mass_pdf_vmax.txt', 
           all_lgm_pdf_vmax/total_weight) 
np.savetxt('obs_data/derived_data/gama_ssfr_pdf.txt', 
           all_ssfr_pdf) 
np.savetxt('obs_data/derived_data/gama_ssfr_pdf_vmax.txt', 
           all_ssfr_pdf_vmax/total_weight) 
np.savetxt('obs_data/derived_data/gama_fundamentalplaneSFR.txt', 
           all_lgm_sfr_pdf) 
np.savetxt('obs_data/derived_data/gama_fundamentalplaneSFR_vmax.txt', 
           all_lgm_sfr_pdf_vmax/total_weight) 

np.savetxt('obs_data/derived_data/gama_prior_lgm.txt', 
           np.array([lgm, dlgm]).T) 
np.savetxt('obs_data/derived_data/gama_prior_ssfr.txt', 
           np.array([ssfr, dssfr]).T) 


