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
from module_utils import rms, smooth, cube_mom0, cube_mom8



