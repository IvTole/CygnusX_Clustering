# Libraries

# Standard
import numpy as np
import scipy
import math
import os

import matplotlib.pyplot as plt 
from matplotlib.colors import LogNorm
from matplotlib.colors import PowerNorm
from matplotlib.ticker import LogFormatter

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

def cube_mom8(cube_path,velmin,velmax,output_path,write_fits=False):
    """
    Returns the moment 8 (max intensity) image of a data cube.

    Parameters:
        cube(SpectralCube): data cube from which the moment is computed.
        velmin(float): min value of the spectral range (in km/s)
        velmax(float): max value of the spectral range (in km/s)

    Returns:
        moment(SpectralCube 2d): Moment 8 image (2d) of the data cube.
    """

    cube = SpectralCube.read(cube_path)
    cube.allow_huge_operations = True
    cube_slab = cube.spectral_slab(velmin *u.km / u.s, velmax *u.km / u.s)
    moment = cube_slab.max(axis = 0)

    if write_fits:
        moment.write(output_path, overwrite=True)

    return moment

def cube_mom0(cube_path,velmin,velmax,output_path,write_fits=False):
    """
    Returns the moment 0 (mean intensity) image of a data cube.

    Parameters:
        cube(SpectralCube): data cube from which the moment is computed.
        velmin(float): min value of the spectral range (in km/s)
        velmax(float): max value of the spectral range (in km/s)

    Returns:
        moment(SpectralCube 2d): Moment 0 image (2d) of the data cube.
    """

    cube = SpectralCube.read(cube_path)
    cube_slab = cube.spectral_slab(velmin*u.km/u.s, velmax*u.km/u.s)
    moment = cube_slab.with_spectral_unit(u.km/u.s).moment(order=0)

    if write_fits:
        moment.write(output_path, overwrite=True)
        
    return moment

def cube_smoothing(data_path, output_path, prefix_source, prefix_emission, kernel_px=1, write_fits=False):
    hdu = fits.open(data_path)[0]
    for v in range(0,hdu.data.shape[0]):
        hdu.data[v,:,:] = smooth(hdu.data[v,:,:],kern_px=kernel_px)
    print('Smoothing done for:', data_path)
    
    if write_fits:
        print('Writing new cube in following path:', output_path)
        hdu.writeto(os.path.join(output_path,prefix_source+'_'+prefix_emission+'_smoothed.fits'),
                overwrite = True)

def plot_mom8(path, output_path, prefix_source, prefix_emission, gamma=1.0, vmin=0.0, vmax=25.0):
    hdu = fits.open(path)[0]

    fig = plt.figure()

    ax = fig.add_subplot(111, projection = WCS(hdu.header))

    im = ax.imshow(hdu.data, cmap = 'RdBu_r',
                   norm = PowerNorm(gamma=gamma, vmin=vmin, vmax=vmax))

    ### Axis parameters ###
    lat = ax.coords['glat']
    lat.set_axislabel('Galactic Latitude', size = 12, alpha = 1.0)
    lat.set_ticks(width = 1, spacing = 0.1 * u.deg)
    lat.set_ticklabel(size = 12, exclude_overlapping=True)
    lat.display_minor_ticks(True)

    lon = ax.coords['glon']
    lon.set_axislabel('Galactic Longitude', size = 12, alpha = 1.0)
    lon.set_ticks(width = 1, spacing = 0.1 * u.deg)
    lon.set_ticklabel(size = 12, exclude_overlapping=True)
    lon.display_minor_ticks(True)

    ### Annotations ###
    ax.annotate(prefix_source + ', ' + prefix_emission + ' Peak Temperature', xy = (5,5), xytext = (5, 5), color='black',
            fontsize = 8, bbox = dict(boxstyle = "round", fc = "w", alpha = 0.0))

    ### Colorbar ###
    cbar = plt.colorbar(im, pad=.01)
    cbar.set_label(r'$T_{\rm MB}^{\rm \ peak}$', labelpad = 4, y = 0.5, rotation=90, size = 14)
    #cbar.ax.tick_params(labelsize=14)
    #cbar.ax.locator_params(nbins=6)

    plt.savefig(os.path.join(output_path, prefix_source + '_' + prefix_emission + '_mom8.pdf'),
                bbox_inches = 'tight')
    plt.close()
    