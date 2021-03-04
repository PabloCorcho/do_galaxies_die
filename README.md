# do_galaxies_die
This repository comprises all the relevant scripts necessary to reproduce the results presented in 'Do galaxies die? Different views from simulations and observations in the Local Universe."

## Brief instructions. 

### Generate Figs. 1 & 2. 
  - Run **compute_simus_planes.py**. This code will compute the probability maps corresponding to each simulation included in **simu_data**.
  - Run **mk_fig_simus_cond_planes.py**.
  
### Generate Figs. 3 & 8. 
  - Run **mk_fig_detection.py**. 

### Generate Fig. 4.
 - Run **mk_fig_sdss_uncertainties.py**.
 - Run **mk_fig_gama_uncertainties.py**.
 
### Generate Figs. 5 & 6.
 - Run **mk_TNG_sfrs_death_times.py**.
 
### Generate Fig. 7.
 - Run **mk_fig_quenched_fractions.py**.
 
### Generate Fig. 8.
 - Run **mk_fig_detection.py**.
 
## All figures.
If the probability maps for each simulations are already storaged, it is possible to compute all the figures by running **make_all_figures.py**
