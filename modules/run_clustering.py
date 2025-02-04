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
from module_utils import rms, smooth, cube_mom0, cube_mom8, cube_smoothing

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
    source_prefix = 'dr21'
    
    # Smoothing (for all emissions)
    cube_smoothing(data_path=data_path_12co, output_path=fits_path, source_prefix=source_prefix, emission_prefix='12co', kernel_px=1, write_fits=True)
    cube_smoothing(data_path=data_path_13co, output_path=fits_path, source_prefix=source_prefix, emission_prefix='13co', kernel_px=1, write_fits=True)
    cube_smoothing(data_path=data_path_c18o, output_path=fits_path, source_prefix=source_prefix, emission_prefix='c18o', kernel_px=1, write_fits=True)

if __name__ == '__main__':
    main()

