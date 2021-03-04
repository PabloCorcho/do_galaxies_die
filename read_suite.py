import numpy as np
from astropy.io import fits


class DataSet(object):
    def __init__(self, path, mode='tot'):
        
        self.data_path = path
        
        hdul = fits.open(path)
        self.table = hdul[1]
        # self.mode = mode
        
        self.get_sfr(mode=mode)
        self.get_mass(mode=mode)
        self.get_ssfr()
        
    def get_sfr(self, mode):
        self.sfr = self.table.data[mode+'_sfr']
        
    def get_mass(self, mode):
        self.mass = self.table.data[mode+'_mass']
    
    def get_ssfr(self):
        self.ssfr = self.sfr/self.mass
    
    def get_photometry(self):
        self.u = self.table.data['u']
        self.g = self.table.data['g']
        self.r = self.table.data['r']
        
    def mass_filter(self, min_mass=0, max_mass=1e14):
        mass_filter = (self.mass>=min_mass)&(self.mass<=max_mass)
        self.mass = self.mass[mass_filter]
        self.sfr = self.sfr[mass_filter]
        self.ssfr = self.ssfr[mass_filter]
        self.mass_mask = mass_filter
        
    def sfr_filter(self, min_sfr=0, max_sfr=1e4):
        sfr_filter = (self.sfr>=min_sfr)&(self.sfr<=max_sfr)
        self.mass = self.mass[sfr_filter]
        self.sfr = self.sfr[sfr_filter]
        self.ssfr = self.ssfr[sfr_filter]
    
    def replace_ssfr0(self, min_ssfr):
        self.ssfr[self.ssfr==0] = min_ssfr
        

class Compute_plane(object):
    
    def __init__(self, DataSet):
        self.DS = DataSet
        
    def define_ssfr_bins(self, bin_edges, log=True):
        self.ssfr_bin_edges = bin_edges
        self.log = log
        
    def running_plane(self, semi_width=50):
        M = self.DS.mass
        sSFR = self.DS.ssfr
        mass_sort = np.argsort(M)
        M = M[mass_sort]
        sSFR = sSFR[mass_sort]
        
        self.plane = np.zeros((M.size, len(self.ssfr_bin_edges)-1))
        self.density_plane = np.zeros((M.size, len(self.ssfr_bin_edges)-1))
        
        if self.log:
            sSFR = np.log10(sSFR)
        
        for ith in range(M.size):
            print('Running % completion -->', ith/M.size*100)
            if ith > semi_width:                
                ssfr_i = sSFR[ith-semi_width:ith+semi_width]            
            else:
                ssfr_i = sSFR[0:2*ith]
                
            h, _ = np.histogram(ssfr_i, bins=self.ssfr_bin_edges, density=True)
            
            H = np.cumsum(h)
            
            H = H/H[-1]
            
            self.density_plane[ith, :] = h
            self.plane[ith, :] = H
        
    def running_quenched_fraction(self, semi_width=50):
        M = self.DS.mass
        sSFR = self.DS.ssfr
        mass_sort = np.argsort(M)
        M = M[mass_sort]
        sSFR = sSFR[mass_sort]
        
        self.quench_frac = np.zeros((M.size))
        self.quench_value = np.nanmin(sSFR)
                
        if self.log:
            sSFR = np.log10(sSFR)
            self.quench_value = np.nanmin(sSFR)
            
        for ith in range(M.size):
            print('Running % completion -->', ith/M.size*100)
            if ith > semi_width:                
                ssfr_i = sSFR[ith-semi_width:ith+semi_width]            
            else:
                ssfr_i = sSFR[0:2*ith]
            if ith ==0:
                self.quench_frac[0] = np.nan
                continue
                            
            
            frac = len(ssfr_i[ssfr_i==self.quench_value])/len(ssfr_i)
            
            self.quench_frac[ith] = frac
            

if __name__ == '__main__':
    
    path = 'simu_data/tng300-1.fits'
    
    my_set = DataSet(path)    
    my_set.mass_filter(min_mass=10**8.5, max_mass=10**11.5)
    
    
    