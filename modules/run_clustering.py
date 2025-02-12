# Libraries

# Standard
import numpy as np
import matplotlib.pyplot as plt
import os

# Matplotlib figures to annotate
from matplotlib.patches import Rectangle
from matplotlib.patches import Ellipse
from matplotlib.patches import Circle

# Astropy
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u

# Spectral cubes
from spectral_cube import SpectralCube

import math
import scipy

# External modules
from module_data_path import cube_data_path, plot_data_path, fits_data_path, mask_data_path, catalog_data_path
from module_utils import rms, smooth, cube_mom0, cube_mom8, cube_smoothing, plot_mom8, plot_mom8_comparison
from module_clustering import make_clustering, make_catalog, make_plot_clusters, make_mask, catalog_mask_drop
from module_data_spectra import make_spectra

stages = [4]

# cube smoothing
def stage1():

    # data, plots and fits files directory path
    data_path = cube_data_path()
    fits_path = fits_data_path()
    plots_path = plot_data_path()

    # original data cubes
    data_path_12co = os.path.join(data_path, 'dr21_12co_cube.fits')
    data_path_13co = os.path.join(data_path, 'dr21_13co_cube.fits')
    data_path_c18o = os.path.join(data_path, 'dr21_c18o_cube.fits')

    # Source prefix
    prefix_source = 'dr21'
    
    # Smoothing (for all emissions)
    cube_smoothing(data_path=data_path_12co, output_path=fits_path, prefix_source=prefix_source, prefix_emission='12co', efficiency=1.0, kernel_px=1, write_fits=True)
    cube_smoothing(data_path=data_path_13co, output_path=fits_path, prefix_source=prefix_source, prefix_emission='13co', efficiency=1.0, kernel_px=1, write_fits=True)
    cube_smoothing(data_path=data_path_c18o, output_path=fits_path, prefix_source=prefix_source, prefix_emission='c18o', efficiency=1.0, kernel_px=1, write_fits=True)

    # Smoothed cubes path
    data_sm_path_12co = os.path.join(fits_path, f'{prefix_source}_12co_smoothed.fits')
    data_sm_path_13co = os.path.join(fits_path, f'{prefix_source}_13co_smoothed.fits')
    data_sm_path_c18o = os.path.join(fits_path, f'{prefix_source}_c18o_smoothed.fits')

    # Make mom8 fits files (using smoothed cubes)
    cube_mom8(cube_path=data_sm_path_12co, velmin=-50.0, velmax=50.0, output_path=os.path.join(fits_path, f'{prefix_source}_12co_mom8.fits'), write_fits=True)
    cube_mom8(cube_path=data_sm_path_13co, velmin=-50.0, velmax=50.0, output_path=os.path.join(fits_path, f'{prefix_source}_13co_mom8.fits'), write_fits=True)
    cube_mom8(cube_path=data_sm_path_c18o, velmin=-50.0, velmax=50.0, output_path=os.path.join(fits_path, f'{prefix_source}_c18o_mom8.fits'), write_fits=True)

    # Mom8 fits path
    data_mom8_path_12co = os.path.join(fits_path, f'{prefix_source}_12co_mom8.fits')
    data_mom8_path_13co = os.path.join(fits_path, f'{prefix_source}_13co_mom8.fits')
    data_mom8_path_c18o = os.path.join(fits_path, f'{prefix_source}_c18o_mom8.fits')

    # Plot mom8 (for all emissions)
    plot_mom8(path=data_mom8_path_12co, output_path=plots_path, prefix_source=prefix_source, prefix_emission='12co', gamma=1.0, vmin=0.0, vmax=45.0)
    plot_mom8(path=data_mom8_path_13co, output_path=plots_path, prefix_source=prefix_source, prefix_emission='13co', gamma=1.0, vmin=0.0, vmax=25.0)
    plot_mom8(path=data_mom8_path_c18o, output_path=plots_path, prefix_source=prefix_source, prefix_emission='c18o', gamma=1.0, vmin=0.0, vmax=5.0)
    
# clustering with astrodendro, extract first cluster catalog and masks    
def stage2():

    # data, plots and fits files directory path
    data_path = cube_data_path()
    fits_path = fits_data_path()
    plots_path = plot_data_path()
    mask_path = mask_data_path()
    catalog_path = catalog_data_path()


    # original data cubes
    data_path_12co = os.path.join(data_path, 'dr21_12co_cube.fits')
    data_path_13co = os.path.join(data_path, 'dr21_13co_cube.fits')
    data_path_c18o = os.path.join(data_path, 'dr21_c18o_cube.fits')

    # Source prefix
    prefix_source = 'dr21'

    # Smoothed cubes path
    data_sm_path_12co = os.path.join(fits_path, f'{prefix_source}_12co_smoothed.fits')
    data_sm_path_13co = os.path.join(fits_path, f'{prefix_source}_13co_smoothed.fits')
    data_sm_path_c18o = os.path.join(fits_path, f'{prefix_source}_c18o_smoothed.fits')

    # Mom8 fits path
    data_mom8_path_12co = os.path.join(fits_path, f'{prefix_source}_12co_mom8.fits')
    data_mom8_path_13co = os.path.join(fits_path, f'{prefix_source}_13co_mom8.fits')
    data_mom8_path_c18o = os.path.join(fits_path, f'{prefix_source}_c18o_mom8.fits')

    # Astrodendro hyperparameters
    T_rms = 0.35 # in K (corrected for main beam efficiency)
    T_min = 3.0 * T_rms
    T_delta = 2.0 * T_rms
    n_vox = 16

    make_clustering(cube_path=data_sm_path_c18o, catalog_path=catalog_path, mask_path=mask_path, T_min=T_min, T_delta=T_delta, n_vox=n_vox, prefix_source=prefix_source, prefix_emission='c18o', catalog=True, mask=True)
    
    make_plot_clusters(mom_path=data_mom8_path_c18o, catalog_path=catalog_path, mask_path=mask_path, plots_path=plots_path, prefix_source=prefix_source, prefix_emission='c18o', gamma=1.0, vmin=0.0, vmax=5.0)

# drop bad indexes (by eye) before and after
def stage3():
    # fits, plots, catalog and mask files directory path
    fits_path = fits_data_path()
    plots_path = plot_data_path()
    catalog_path = catalog_data_path()
    mask_path = mask_data_path()

    # Source prefix
    prefix_source = 'dr21'

    # Mom8 fits path
    data_mom8_path_12co = os.path.join(fits_path, f'{prefix_source}_12co_mom8.fits')
    data_mom8_path_13co = os.path.join(fits_path, f'{prefix_source}_13co_mom8.fits')
    data_mom8_path_c18o = os.path.join(fits_path, f'{prefix_source}_c18o_mom8.fits')

    # Drop list
    drop_list = [13]

    plot_mom8_comparison(mom_path=data_mom8_path_c18o, plots_path=plots_path, catalog_path=catalog_path, prefix_source=prefix_source, prefix_emission='c18o', dropped=False, gamma=1.0, vmin=0.0, vmax=5.0)
    catalog_mask_drop(catalog_path=catalog_path, mask_path=mask_path, drop_list=drop_list, prefix_source=prefix_source, prefix_emission='c18o')
    plot_mom8_comparison(mom_path=data_mom8_path_c18o, plots_path=plots_path, catalog_path=catalog_path, prefix_source=prefix_source, prefix_emission='c18o', dropped=True, gamma=1.0, vmin=0.0, vmax=5.0)

# spectra extraction
def stage4():
    
    # data, plots, fits, catalog and mask files directory path
    data_path = cube_data_path()
    fits_path = fits_data_path()
    plots_path = plot_data_path()
    mask_path = mask_data_path()
    catalog_path = catalog_data_path()


    # original data cubes
    data_path_12co = os.path.join(data_path, 'dr21_12co_cube.fits')
    data_path_13co = os.path.join(data_path, 'dr21_13co_cube.fits')
    data_path_c18o = os.path.join(data_path, 'dr21_c18o_cube.fits')

    # Source prefix
    prefix_source = 'dr21'

    # Smoothed cubes path
    data_sm_path_12co = os.path.join(fits_path, f'{prefix_source}_12co_smoothed.fits')
    data_sm_path_13co = os.path.join(fits_path, f'{prefix_source}_13co_smoothed.fits')
    data_sm_path_c18o = os.path.join(fits_path, f'{prefix_source}_c18o_smoothed.fits')

    # Mom8 fits path
    data_mom8_path_12co = os.path.join(fits_path, f'{prefix_source}_12co_mom8.fits')
    data_mom8_path_13co = os.path.join(fits_path, f'{prefix_source}_13co_mom8.fits')
    data_mom8_path_c18o = os.path.join(fits_path, f'{prefix_source}_c18o_mom8.fits')

    make_spectra(cube_path=data_path, catalog_path=catalog_path, mask_path=mask_path, plots_path=plots_path,
                 prefix_source=prefix_source,prefix_emission='c18o', prefix_cube='12co', efficiency=1.0, height=4.0, distance=15)
    make_spectra(cube_path=data_path, catalog_path=catalog_path, mask_path=mask_path, plots_path=plots_path,
                 prefix_source=prefix_source, prefix_emission='c18o', prefix_cube='13co', efficiency=1.0, height=2.0, distance=15)
    make_spectra(cube_path=data_path, catalog_path=catalog_path, mask_path=mask_path, plots_path=plots_path,
                 prefix_source=prefix_source, prefix_emission='c18o', prefix_cube='c18o', efficiency=1.0, height=0.75, distance=5)



if __name__ == '__main__': 
    
    if 1 in stages:
        stage1()
    elif 2 in stages:
        stage2()
    elif 3 in stages:
        stage3()
    elif 4 in stages:
        stage4()

