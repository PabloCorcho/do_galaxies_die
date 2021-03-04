from read_suite import DataSet, Compute_plane
import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits

from glob import glob

from scipy.stats import binned_statistic

from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

from os.path import isdir
from os import mkdir

"""
This script computes the plane dlog(ssfr)(ssfr|M)/dlog(ssfr) for each 
simulation suite. 
"""

# %%
# =============================================================================
# SIMULATIONS
# =============================================================================
# data_sets = glob('simu_data/*.fits')

data_sets = [
                'simu_data/MAGNETICUM.fits',
                'simu_data/tng300-1.fits',
               'simu_data/tng100-1.fits',                                             
               'simu_data/SIMBA_m100n1024_151.fits',
              'simu_data/SIMBA_m50n512_151.fits',
              'simu_data/SIMBA_m25n512_151.fits',
              'simu_data/RefL0100N1504.fits',
              'simu_data/RefL0050N0752.fits',
              'simu_data/RefL0025N0752.fits'
             ]
data_sets_names = [
                       'MAG500',
                      'IllustrisTNG300',
                    'IllustrisTNG100',                                     
                     'SIMBA150',                    
                    'SIMBA75',
                    'SIMBA35',
                    'EAGLE100',
                    'EAGLE50',
                    'EAGLE25'                    
                    ]


if not isdir('simu_data/simu_derived_data'):
    mkdir('simu_data/derived_data')
    print('Folder created: simu_data/simu_derived_data')
    

# %%

ssfr_bin_edges = np.arange(-16, -7, 0.2)

ssfr_bins = (ssfr_bin_edges[:-1] + ssfr_bin_edges[1:])/2
    

lgm_bin_edges = np.arange(7.5, 12.5, 0.2)
lgm_bins = (lgm_bin_edges[1:]+lgm_bin_edges[:-1])/2


for ith, data_i in enumerate(data_sets):
    print('Data set ', ith+1)


    simu_i = DataSet(data_i)    
    simu_i.mass_filter(min_mass=10**8.5, max_mass=10**11.5)
    simu_i.replace_ssfr0(1e-16)
    
    Computer = Compute_plane(simu_i)
    
    Computer.define_ssfr_bins(ssfr_bin_edges)
    Computer.running_plane(semi_width=min(300, len(simu_i.mass)//50))

    M = np.log10(np.sort(simu_i.mass))
    density_plane = Computer.density_plane
    cumulative_plane = Computer.plane
            
    c1 = fits.Column(name='mass', array=M, format='E')    
    c2 = fits.Column(name='ssfr_bins', array=ssfr_bins, format='E')
        
    hdr = fits.Header()
    hdr['COMMENT1'] = "Derived data from simulation: "+data_sets_names[ith]     
    hdr['COMMENT2'] = "Masses are sorted"
    hdr['COMMENT3'] = "Mass units Msun"    
    hdr['COMMENT4'] = "sSFR units 1/yr"  
    hdr['COMMENT5'] = "Image 1: density plane dlog(ssfr)(ssfr|M)dlog(ssfr)"  
    hdr['COMMENT6'] = "Image 2: cumulative plane dlog(ssfr)(ssfr|M)dlog(ssfr)"  

    empty_primary = fits.PrimaryHDU(header=hdr)

    image1 = fits.ImageHDU(density_plane)
    image2 = fits.ImageHDU(cumulative_plane)

    ##########################################################################
    
    simu_i = DataSet(data_i)
    
    simu_i.mass_filter(min_mass=10**8.5, max_mass=10**11.5)
    simu_i.sfr_filter(min_sfr=1e-5)
        
    Computer = Compute_plane(simu_i)
    
    Computer.define_ssfr_bins(ssfr_bin_edges)
    Computer.running_plane(semi_width=min(300, len(simu_i.mass)//50))
    
    M = np.log10(np.sort(simu_i.mass))
        
    c3 = fits.Column(name='mass_filtered', array=M, format='E')    
    

    density_plane = Computer.density_plane
    cumulative_plane = Computer.plane
    
    table_hdu = fits.BinTableHDU.from_columns([c1, c2, c3])
    
    hdr['COMMENT6'] = "Image 3: Same as Image 2 without sfr=0 galaxies"  
    image3 = fits.ImageHDU(cumulative_plane)


    hdu = fits.HDUList([empty_primary, table_hdu, image1, image2, image3])

    hdu.writeto('simu_data/simu_derived_data/'+data_sets_names[ith]+'_planes.fits', 
            overwrite=True)
    
print('Process completed')    

# ...
