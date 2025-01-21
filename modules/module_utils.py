# Libraries

# Standard
import numpy as np
import scipy
import math

# Astropy
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u

# Data Cube
from spectral_cube import SpectralCube

def rms(image):
    """""
    Returns root mean square error (rms) of an image (2d array).

    Parameters:
        image(np.darray):The 2d array used to calculate the rms.

    Returns:
        rms(float):The rms of the image.   
    """
    
    rms = np.sqrt((np.mean(image**2.0)))
    return rms

# Image smoothing using a Gaussian Kernel
def smooth(image, kern_px=1):
    """
    Returns a smoothed image in the first HDU of the input file.

    Parameters:
        image(2d np.darray): image to be smoothed.
        kern_px: FWHM of kernel in pixels.

    Returns:
        f1(2d np.darray): Smoothed 2d array.
    """

    f1=scipy.ndimage.gaussian_filter(image, kern_px/(2*math.sqrt(2*math.log(2))))
    return f1

def cube_mom8(cube,velmin,velmax):
    """
    Returns the moment 8 (max intensity) image of a data cube.

    Parameters:
        cube(SpectralCube): data cube from which the moment is computed.
        velmin(float): min value of the spectral range (in km/s)
        velmax(float): max value of the spectral range (in km/s)

    Returns:
        moment(SpectralCube 2d): Moment 8 image (2d) of the data cube.
    """

    cube_slab = cube.spectral_slab(velmin *u.km / u.s, velmax *u.km / u.s)
    moment = cube_slab.max(axis = 0)
    return moment

def cube_mom0(cube,velmin,velmax):
    """
    Returns the moment 0 (mean intensity) image of a data cube.

    Parameters:
        cube(SpectralCube): data cube from which the moment is computed.
        velmin(float): min value of the spectral range (in km/s)
        velmax(float): max value of the spectral range (in km/s)

    Returns:
        moment(SpectralCube 2d): Moment 0 image (2d) of the data cube.
    """

    cube_slab = cube.spectral_slab(velmin*u.km/u.s, velmax*u.km/u.s)
    moment = cube_slab.with_spectral_unit(u.km/u.s).moment(order=0)
    return moment