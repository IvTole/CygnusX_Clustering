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

def make_spectra(cube_path, catalog_path, mask_path, plots_path, prefix_source, prefix_emission, prefix_cube, efficiency=1.0):
    
    catalog = pd.read_csv(os.path.join(catalog_path, f"{prefix_source}_catalog_{prefix_emission}_dropped.csv"))
    print(f'Open catalog for {prefix_source}, line = {prefix_emission}')
    mask = np.load(os.path.join(mask_path, f'{prefix_source}_{prefix_emission}_masks_dropped.npy'))
    print(f'Open mask for {prefix_source}, line = {prefix_emission}')
    cube = SpectralCube.read(os.path.join(cube_path, f'{prefix_source}_{prefix_cube}_cube.fits'))
    print(f'Extracting spectra for {prefix_source}, line = {prefix_cube}')

    # Extraction (by index)
    clump_index = 0
    fig_all = plt.figure(figsize=(15,15))     # initialize fig for plot with all spectra

    for index in range(0,len(mask)):
    
        print(f'clump_idx = {clump_index}')

        for j in range(0,mask[index].shape[0]):
            for k in range(0,mask[index].shape[1]):
                if mask[index][j,k] == True:
                    x_init = j
                    y_init = k
                    #print('True')
                    break
            else:
                continue
            break

        spec = cube[:,x_init,y_init]
        x = spec.spectral_axis
        x = x.to_value()
        y = spec.to_value()

        x = np.array(x)
        y = np.array(y)

        # I ignore some indexes, in case there are nan values
        min_ch = 10
        max_ch = -10
        x = x[min_ch:max_ch]
        y = y[min_ch:max_ch]

        count_npix = 1

        for j in range(0,mask[index].shape[0]):
            for k in range(0,mask[index].shape[1]):
        
                if (j==x_init and k==y_init):
                    continue
                else:
            
                    if mask[index][j,k] == True:
                        spec = cube[:,j,k]
                        xt = spec.spectral_axis
                        xt = xt.to_value()
                        yt = spec.to_value()
        
                        xt = np.array(xt)
                        yt = np.array(yt)
        
                        xt = xt[min_ch:max_ch]
                        yt = yt[min_ch:max_ch]
        
                        y = y + yt
    
                        count_npix += 1
    
        print(f'n_vox = {count_npix}')

        y_p = [((1/count_npix)*ys/efficiency) for ys in y]
        x_p = [value / 1000 for value in x] # in km/s

        # Gaussian Fit
        max_index = y_p.index(max(y_p)) # index of the maximum value, used as velocity parameter

        gmodel = GaussianModel(prefix='p1_')
        params = gmodel.make_params(p1_amplitude=10, p1_center=x_p[max_index], p1_sigma=0.5)
        result = gmodel.fit(y_p, params, x=x_p)
        #print(result.fit_report())

        # Plot (individual)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlim(-25,25)
        ax.set_xlabel(r'V$_{\rm LSR}$ [km/s]')
        ax.set_ylabel(r'T$_{\rm MB}$')
        ax.set_title(f'Clump {clump_index}, {prefix_cube}')
        ax.plot(x_p, y_p, drawstyle='steps-mid', c='black')
        ax.plot(x_p, result.best_fit, 'r')
        fig.savefig(os.path.join(plots_path, f'{prefix_cube}_spectra', f'Clump_{clump_index}_{prefix_source}.pdf'), bbox_inches='tight')
        print(f'Saving figure for clump {clump_index}')
        plt.close(fig)

        # Plot (all)
        ax = fig_all.add_subplot(int(len(catalog)/3)+1,3,clump_index+1)
        ax.set_xlim(-25,25)
        ax.set_ylim(-1.5,max(y_p) + 0.5)
        ax.set_ylabel(' ')
        ax.set_xlabel(' ')
        ax.plot(x_p,y_p,drawstyle='steps-mid', c='black', label = f'Cl_{clump_index}_{prefix_cube}')
        ax.plot(x_p, result.best_fit, 'r')
        ax.legend(loc='upper left', fontsize=6)

        clump_index += 1

    fig_all.text(0.5, 0.08, r'V$_{\rm LSR}$', ha='center', size = 20)
    fig_all.text(0.08, 0.5, r'T$_{\rm MB}$', va='center', rotation='vertical', size=20)
    
    fig_all.savefig(os.path.join(plots_path,f'{prefix_cube}_spectra', f'Clump_{prefix_source}_all.pdf'),
                bbox_inches = 'tight')
    print(f'Saving figure for all spectra of {prefix_cube} for {prefix_source}')
    plt.close(fig_all)
