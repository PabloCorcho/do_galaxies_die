# do_galaxies_die

## Overview
This repository comprises all the relevant scripts necessary to reproduce the results presented in "Do galaxies die? Different views from simulations and observations in the Local Universe".

### Installation

```
git clone https://github.com/PabloCorcho/do_galaxies_die
cd do_galaxies_die
```
Run
```
python compute_simus_planes.py
```
This code will compute the probability maps corresponding to each simulation included in **simu_data**. 1.7 Gb of memory will be necessary to store the output results.

Then run
```
python make_all_figures.py
```
This code will generate all the figures presented in Corcho et al. 2021

### Individual figures

If you want to replot a single figure, you can directly run the following scripts:

- Figs. 1 & 2: 
```
python mk_fig_simus_cond_planes.py
```
- Figs. 3 & 8: 
```
python mk_fig_detection.py
```

- Fig. 4: 
```
python mk_fig_sdss_uncertainties.py
python mk_fig_gama_uncertainties.py
```

- Figs. 5 & 6: 
```
python mk_TNG_sfrs_death_times.py
```

- Fig. 7: 
```
python mk_fig_quenched_fractions.py
```

# System Requirements
## Hardware requirements

These scripts only require a standard computer with enough RAM to support the in-memory operations.

## Software requirements 

All scripts have been tested on the following systems:
+ Linux: Ubuntu 20.04

### Python dependencies 

```
os
numpy
scipy
astropy
pandas
matplotlib
```	

# License

This project is covered under the **GNU General Public License v3.0**



