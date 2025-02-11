# Libraries

# Standard
import numpy as np
import pandas as pd
import os

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.colors import PowerNorm
from matplotlib.ticker import LogFormatter
from matplotlib.patches import Rectangle
from matplotlib.patches import Ellipse
from matplotlib.patches import Circle

## Astropy
from astropy.io import fits
from astropy import wcs
from astropy import units as u
from astropy.wcs import WCS

# Lmfit
from lmfit.models import GaussianModel

# Data cubes
from spectral_cube import SpectralCube

# uncertainties
from uncertainties import ufloat
from uncertainties import unumpy as unp

def make_spectra(cube_path, catalog_path, mask_path, plots_path, prefix_source, prefix_emission, prefix_cube):
    
    catalog = pd.read_csv(os.path.join(catalog_path, f"{prefix_source}_catalog_{prefix_emission}.csv"))
    print(f'Open catalog for {prefix_source}, line = {prefix_emission}')
    mask = np.load(os.path.join(mask_path, f'{prefix_source}_{prefix_emission}_masks.npy'))
    print(f'Open mask for {prefix_source}, line = {prefix_emission}')
    cube = SpectralCube.read(os.path.join(cube_path, f'{prefix_source}_{prefix_cube}_cube.fits'))
    print(f'Extracting spectra for {prefix_source}, line = {prefix_cube}')


