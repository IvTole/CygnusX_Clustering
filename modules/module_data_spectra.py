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
from astropy.io.fits.verify import VerifyWarning
import warnings
warnings.simplefilter('ignore', category=VerifyWarning)


# Lmfit
from lmfit.models import GaussianModel

# scipy
from scipy.signal import find_peaks

# Data cubes
from spectral_cube import SpectralCube

# uncertainties
from uncertainties import ufloat
from uncertainties import unumpy as unp

from module_utils import Gauss_area

def make_spectra(cube_path, catalog_path, mask_path, plots_path, prefix_source, prefix_emission, prefix_cube, height, distance, efficiency=1.0):
    
    catalog = pd.read_csv(os.path.join(catalog_path, f"{prefix_source}_catalog_{prefix_emission}_dropped.csv"))
    print(f'Open catalog for {prefix_source}, line = {prefix_emission}')
    mask = np.load(os.path.join(mask_path, f'{prefix_source}_{prefix_emission}_masks_dropped.npy'))
    print(f'Open mask for {prefix_source}, line = {prefix_emission}')
    cube = SpectralCube.read(os.path.join(cube_path, f'{prefix_source}_{prefix_cube}_cube.fits'))
    print(f'Extracting spectra for {prefix_source}, line = {prefix_cube}')

    # Extraction (by index)
    clump_index = 0
    fig_all = plt.figure(figsize=(15,15))     # initialize fig for plot with all spectra
    dfs = [] # list of dataframes for fit parameters

    for index in range(0,len(mask)):
    
        print(f'clump_idx = {clump_index}')

        x_p, y_p, n_vox = spectra_extraction(cube=cube, mask=mask[index], efficiency=efficiency)

        # Gaussian Fit
        # find peaks in emission
        peaks, _ = find_peaks(y_p, height=height, distance=distance)
        ngaus = len(peaks)
        print(f'ngaus={ngaus}')

        # initial parameters using peak and index
        gaussian_params = []
        for i, peak in enumerate(peaks):
            amp = y_p[peak]  # Peak height as amplitude
            cen = x_p[peak]  # Peak position as center
            wid = 1  # Initial guess for width (adjustable)
            gaussian_params.extend([amp, cen, wid])
        gaussian_params = tuple(gaussian_params)

        gaussian_result = gaussian_fit(*gaussian_params,ngaus=ngaus,x=x_p,y=y_p)

        df_params = export_lmfit_results(result=gaussian_result, id=clump_index, n_vox=n_vox)
        dfs.append(df_params)

        # Gaussian Fit
        

        #gmodel = GaussianModel(prefix='p1_')
        #result = gmodel.fit(y_p, params, x=x_p)
        ##params = gmodel.make_params(p1_amplitude=10, p1_center=x_p[max_index], p1_sigma=0.5)
        #print(result.fit_report())

        # Plot (individual)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlim(-25,25)
        ax.set_xlabel(r'V$_{\rm LSR}$ [km/s]')
        ax.set_ylabel(r'T$_{\rm MB}$')
        ax.set_title(f'Clump {clump_index}, {prefix_cube}')
        ax.plot(x_p, y_p, drawstyle='steps-mid', c='black')
        ax.plot(x_p, gaussian_result.best_fit, 'r')
        fig.savefig(os.path.join(plots_path, f'{prefix_cube}_spectra', f'Clump_{clump_index}_{prefix_source}.pdf'), bbox_inches='tight')
        print(f'Saving figure for clump {clump_index}')
        plt.close(fig)

        # Plot (all)
        ax = fig_all.add_subplot(int(len(catalog)/3)+1,3,clump_index+1)
        ax.set_xlim(-45,45)
        ax.set_ylim(-1.5,max(y_p) + 0.5)
        ax.set_ylabel(' ')
        ax.set_xlabel(' ')
        ax.plot(x_p,y_p,drawstyle='steps-mid', c='black', label = f'Cl_{clump_index}_{prefix_cube}')
        ax.plot(x_p, gaussian_result.best_fit, 'r')
        ax.legend(loc='upper left', fontsize=6)

        clump_index += 1

    fig_all.text(0.5, 0.08, r'V$_{\rm LSR}$', ha='center', size = 20)
    fig_all.text(0.08, 0.5, r'T$_{\rm MB}$', va='center', rotation='vertical', size=20)
    
    fig_all.savefig(os.path.join(plots_path,f'{prefix_cube}_spectra', f'{prefix_source}_{prefix_cube}_spectra_all.pdf'),
                bbox_inches = 'tight')
    print(f'Saving figure for all spectra of {prefix_cube} for {prefix_source}')
    plt.close(fig_all)

    # export csv for fit parameters
    df_params_concat = pd.concat(dfs, axis=0)
    df_params_concat.to_csv(os.path.join(catalog_path,f'{prefix_source}_{prefix_cube}_fit_params.csv'),
                            header = True,
                            index = False)


def spectra_extraction(cube, mask, efficiency):

    for j in range(0,mask.shape[0]):
        for k in range(0,mask.shape[1]):
            if mask[j,k] == True:
                x_init = j
                y_init = k
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

    # Ignore some indexes, in case there are nan values
    min_ch = 10
    max_ch = -10
    x = x[min_ch:max_ch]
    y = y[min_ch:max_ch]

    count_npix = 1

    for j in range(0,mask.shape[0]):
        for k in range(0,mask.shape[1]):
    
            if (j==x_init and k==y_init):
                continue
            else:
        
                if mask[j,k] == True:
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

    y_p = [((1/count_npix)*ys/efficiency) for ys in y] # corrected for mb
    x_p = [value / 1000 for value in x] # in km/s

    return x_p, y_p, count_npix

def gaussian_fit(*args,ngaus,x,y):

    gmodel = None
    for i in range(1, ngaus + 1):
        prefix = f'p{i}_'
        if gmodel is None:
            gmodel = GaussianModel(prefix=prefix)
        else:
            gmodel = gmodel + GaussianModel(prefix=prefix)

    params = gmodel.make_params()
    for i in range(1, ngaus + 1):
        prefix = f'p{i}_'
        amplitude, center, sigma = args[(i - 1) * 3:(i - 1) * 3 + 3]
        params[f'{prefix}amplitude'].set(value=amplitude)
        params[f'{prefix}center'].set(value=center)
        params[f'{prefix}sigma'].set(value=sigma)

    result = gmodel.fit(y, params, x=x)

    return result

def export_lmfit_results(result, id, n_vox):
    """
    Extracts Gaussian fit parameters from an lmfit result object and exports them to a CSV file.

    Parameters:
    - result : lmfit ModelResult -> The fitted lmfit result object

    Returns:
    - df_params : pandas DataFrame -> Fit parameters
    """

    fit_data = []

    # Extract the number of Gaussian components from the result
    gaussians = [key.split('_')[0] for key in result.params if "amplitude" in key]
    unique_gaussians = sorted(set(gaussians), key=lambda x: int(x[1:]))  # Sort by index

    for i, g_prefix in enumerate(unique_gaussians):
        fit_data.append({
            "Clump_id":id,
            "Gaussian": i + 1,
            "N_vox": n_vox,
            "Amplitude": result.params[f"{g_prefix}_amplitude"].value,
            "Amplitude_err": result.params[f"{g_prefix}_amplitude"].stderr,
            "Center": result.params[f"{g_prefix}_center"].value,
            "Center_err": result.params[f"{g_prefix}_center"].stderr,
            "Sigma": result.params[f"{g_prefix}_sigma"].value,
            "Sigma_err": result.params[f"{g_prefix}_sigma"].stderr,
            "FWHM": result.params[f"{g_prefix}_fwhm"].value,
            "FWHM_err": result.params[f"{g_prefix}_fwhm"].stderr,
            "Height": result.params[f"{g_prefix}_height"].value,
            "Height_err": result.params[f"{g_prefix}_height"].stderr,
            "Area": Gauss_area(ufloat(result.params[f"{g_prefix}_height"].value, result.params[f"{g_prefix}_height"].stderr), ufloat(result.params[f"{g_prefix}_fwhm"].value, result.params[f"{g_prefix}_fwhm"].stderr)).n,
            "Area_err": Gauss_area(ufloat(result.params[f"{g_prefix}_height"].value, result.params[f"{g_prefix}_height"].stderr), ufloat(result.params[f"{g_prefix}_fwhm"].value, result.params[f"{g_prefix}_fwhm"].stderr)).s
        })
    df_params = pd.DataFrame(fit_data)

    return df_params