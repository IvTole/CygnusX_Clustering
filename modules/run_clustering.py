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
from module_utils import rms, smooth, cube_mom0, cube_mom8

# Data import


# plots and fits files mainpath
data_path = cube_data_path()
fits_path = fits_data_path()
plots_path = plot_data_path()