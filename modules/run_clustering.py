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
from module_data_path import cube_data_path, plot_data_path, fits_data_path
from module_utils import rms, smooth, cube_mom0, cube_mom8, cube_smoothing, plot_mom8

def main():

    # data, plots and fits files mainpath
    data_path = cube_data_path()
    fits_path = fits_data_path()
    plots_path = plot_data_path()

    #filenames
    data_path_12co = os.path.join(data_path, 'dr21_12co_cube.fits')
    data_path_13co = os.path.join(data_path, 'dr21_13co_cube.fits')
    data_path_c18o = os.path.join(data_path, 'dr21_c18o_cube.fits')

    # Prefix
    prefix_source = 'dr21'
    
    # Smoothing (for all emissions)
    cube_smoothing(data_path=data_path_12co, output_path=fits_path, prefix_source=prefix_source, prefix_emission='12co', kernel_px=1, write_fits=True)
    cube_smoothing(data_path=data_path_13co, output_path=fits_path, prefix_source=prefix_source, prefix_emission='13co', kernel_px=1, write_fits=True)
    cube_smoothing(data_path=data_path_c18o, output_path=fits_path, prefix_source=prefix_source, prefix_emission='c18o', kernel_px=1, write_fits=True)

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
    

if __name__ == '__main__':
    main()

