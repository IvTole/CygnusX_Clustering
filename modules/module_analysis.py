import numpy as np
import pandas as pd
import os

## Astropy
from astropy.io import fits
from astropy import wcs
from astropy import units as u
from astropy.wcs import WCS

# uncertainties
from uncertainties import ufloat
from uncertainties import unumpy as unp

def radius(distance,theta_pix,n_sky):
    '''
    Obtains a radius with distance, pixel scale (in radians), and number of projected pixels on sky
    '''       
    res = distance * theta_pix * np.sqrt(n_sky / np.pi)
    return res

def T_ex(T12_mb):
    res = (5.53) / (unp.log(1 + ((5.53) / (T12_mb + 0.83)) ))
    return res

def M_LTE(Tex,I):
    X18 = 5.9E06
    res = 13.2 * (X18/5.9E06) * Tex * unp.exp(5.27/Tex)* I
    return res

def number_density(M,R): # expressed in cm-3
    # 1pc = 3.086e+18
    # 1 solar mass = 1.989e+33 g
    mu = 1.36
    mp = 1.67E-24
    res = (3*M) / (4 * np.pi * 2 * mu * mp * (R**3.0))
    res = res * 1.989e+33 # solar mass to grams
    res = res / ((3.086e+18)**3.0) # pc to cm
    return res

def M_vir(R,sigma):
    k=0
    G = 4.301E-03 # Gravitational constant in km2 Mpc M_sun-1 s-2
    res = ((3*(5 - 2*k)) / (G*(3-k))) * R * (sigma**2.0)
    return res

def cube_info(cube_path, prefix_source, prefix_emission):
    hdu = fits.open(os.path.join(cube_path,f'{prefix_source}_{prefix_emission}_cube.fits'))[0]

    ### Cube scale information
    # reference pixel
    pix_ref_lon = hdu.header['CRPIX1']
    pix_ref_lat = hdu.header['CRPIX2']
    pix_ref_vel = hdu.header['CRPIX3']

    # center
    center_lon = hdu.header['CRVAL1']
    center_lat = hdu.header['CRVAL2']
    center_vel = hdu.header['CRVAL3']

    # coord increase
    pix_incr_lon = hdu.header['CDELT1']
    pix_incr_lat = hdu.header['CDELT2']
    pix_incr_vel = hdu.header['CDELT3']

    # pix max range
    pix_max_lon = hdu.header['NAXIS1']
    pix_max_lat = hdu.header['NAXIS2']
    pix_max_vel = hdu.header['NAXIS3']
    
    ## min max values
    min_lon = center_lon - (pix_incr_lon*(pix_ref_lon))
    #print('min longitude =', min_lon, ' deg')
    max_lon = center_lon + (pix_incr_lon*(pix_max_lon - pix_ref_lon))
    #print('max longitude =', max_lon, ' deg')
    min_lat = center_lat - (pix_incr_lat*(pix_ref_lat))
    #print('min latitude =', min_lat, ' deg')
    max_lat = center_lat + (pix_incr_lat*(pix_max_lat - pix_ref_lat))
    #print('max latitude =', max_lat, ' deg')
    min_vel = center_vel - (pix_incr_vel*(pix_ref_vel))
    #print('min velocity =', min_vel, ' m/s')
    max_vel = center_vel + (pix_incr_vel*(pix_max_vel - pix_ref_vel))
    #print('max velocity =', max_vel, ' m/s')

    return min_lon, max_lon, pix_incr_lon, min_lat, max_lat, pix_incr_lat, min_vel, max_vel, pix_incr_vel



def make_analysis(cube_path, catalog_path, mask_path, prefix_source, prefix_emission, prefix_cube, source_distance):

    min_lon, max_lon, pix_incr_lon, min_lat, max_lat, pix_incr_lat, min_vel, max_vel, pix_incr_vel = cube_info(cube_path=cube_path, prefix_source=prefix_source, prefix_emission=prefix_cube)

    cat_fit_12 = pd.read_csv(os.path.join(catalog_path, f"{prefix_source}_12co_fit_params.csv"))
    print(f'Open catalog(fit parameters) for {prefix_source}, line = 12co')
    cat_fit_18 = pd.read_csv(os.path.join(catalog_path, f"{prefix_source}_c18o_fit_params.csv"))
    print(f'Open catalog(fit parameters) for {prefix_source}, line = c18o')
    cat_cl_18 = pd.read_csv(os.path.join(catalog_path, f"{prefix_source}_catalog_{prefix_emission}_dropped.csv"))
    print(f'Open catalog(astrodendro clusters, dropped) for {prefix_source}, line = {prefix_emission}')

    mask = np.load(os.path.join(mask_path, f'{prefix_source}_{prefix_emission}_masks_dropped.npy'))
    print(f'Open mask for {prefix_source}, line = {prefix_emission}')

    cols = ['Id','l','b','Cen_c18o','FWHM_c18o','Height_c18o','Height_12co','I_c18o','R_cl','T_ex','M_LTE','n_H2','M_vir']
    phys_params_dict = {col: [] for col in cols} # initialize dictionary

    clump_indexes = sorted(cat_fit_18['Clump_id'].unique())

    for index in clump_indexes:
        cat_fit_sub_12 = cat_fit_12[cat_fit_12['Clump_id']==index]
        cat_fit_sub_18 = cat_fit_18[cat_fit_18['Clump_id']==index]
        cat_cl_sub_18 = cat_cl_18[cat_cl_18['index']==index]

        # Find the index of the closest Gaussian in 12co in velocity
        reference_velocity = cat_fit_sub_18.iloc[0]["Center"]
        closest_idx = (cat_fit_sub_12["Center"] - reference_velocity).abs().idxmin()
        cat_fit_sub_12 = cat_fit_sub_12.loc[[closest_idx]]

        # Filling the last catalog

        theta_pix = abs(pix_incr_lat) * 0.0174533 # pixel physical size in (deg -> rad)

        # Id
        phys_params_dict['Id'].append(index)

        # Position
        value = round(cat_cl_sub_18['x_cen'].iloc[0] * pix_incr_lon + min_lon, 4)
        phys_params_dict['l'].append(value)
        value = round(cat_cl_sub_18['y_cen'].iloc[0] * pix_incr_lat + min_lat, 4)
        phys_params_dict['b'].append(value)

        # Velocity
        value = str(round(cat_fit_sub_18['Center'].iloc[0],2))+'('+str(round(cat_fit_sub_18['Center_err'].iloc[0],2))+')'
        phys_params_dict['Cen_c18o'].append(value)

        # FWHM
        value = str(round(cat_fit_sub_18['FWHM'].iloc[0],2))+'('+str(round(cat_fit_sub_18['FWHM_err'].iloc[0],2))+')'
        phys_params_dict['FWHM_c18o'].append(value)

        # Height
        value = str(round(cat_fit_sub_18['Height'].iloc[0],2))+'('+str(round(cat_fit_sub_18['Height_err'].iloc[0],2))+')'
        phys_params_dict['Height_c18o'].append(value)
        value = str(round(cat_fit_sub_12['Height'].iloc[0],2))+'('+str(round(cat_fit_sub_12['Height_err'].iloc[0],2))+')'
        phys_params_dict['Height_12co'].append(value)

        # Integrated intensity
        int_intensity = ufloat(cat_fit_sub_18['Area'].iloc[0],cat_fit_sub_18['Area_err'].iloc[0]) * (theta_pix * source_distance)**2.0 * cat_fit_sub_18['N_sky'].iloc[0]
        value = str(round(int_intensity.n,2))+'('+str(round(int_intensity.s,2))+')'
        phys_params_dict['I_c18o'].append(value)

        # Radius
        rad = radius(distance=source_distance,theta_pix=theta_pix,n_sky=cat_fit_sub_18['N_sky'].iloc[0])
        value = round(rad.n,2)
        phys_params_dict['R_cl'].append(value)

        # Excitation temperature
        temp = T_ex(ufloat(cat_fit_sub_12['Height'].iloc[0],cat_fit_sub_12['Height_err'].iloc[0]))
        value = str(round(temp.n,2))+'('+str(round(temp.s,2))+')'
        phys_params_dict['T_ex'].append(value)

        # LTE Mass
        mass = M_LTE(temp, int_intensity)/100 # in units of 10^{2} solar mass
        value = str(round(mass.n,2))+'('+str(round(mass.s,2))+')'
        phys_params_dict['M_LTE'].append(value)

        # Number density H2
        density = number_density(mass, rad)
        value = str(round(density.n,2))+'('+str(round(density.s,2))+')'
        phys_params_dict['n_H2'].append(value)

        # Virial mass
        dv = ufloat(cat_fit_sub_18['FWHM'].iloc[0],cat_fit_sub_18['FWHM_err'].iloc[0])
        vir = M_vir(rad, dv/(2.0 * np.sqrt(2.0 * np.log(2.0)))) / 100 # in units of 10^{2} solar mass
        value = str(round(vir.n,2))+'('+str(round(vir.s,2))+')'
        phys_params_dict['M_vir'].append(value)

    phys_params_df = pd.DataFrame(phys_params_dict)
    phys_params_df.to_csv(os.path.join(catalog_path,f"{prefix_source}_clumps_catalog.csv"),
                          header = True,
                          index = False)
    print(f'Saving final {prefix_emission} clump catalog for {prefix_source}')

    latex_table = phys_params_df.to_latex(index=False)
    with open(os.path.join(catalog_path,f"{prefix_source}_clumps_catalog_latex.txt"),"w") as f:
        f.write(latex_table)
    print(f'LaTeX catalog table saved for {prefix_source}')