import numpy as np
from astropy.io import fits
import os
from matplotlib.patches import Ellipse
from scipy.optimize import curve_fit
from scipy.integrate import quad

import pickle

###############################################################################
#SBP fitting imports
import sep
import h5py
import numpy as np

import os

import matplotlib as mpl
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)

from astropy.modeling import models, fitting

#---------------------------------------------------------------------------#
#User imports
import sys
sys.path.append('/Users/fardila/Documents/Github/kungpao')
from kungpao.galsbp import galSBP
from kungpao.display import display_single, random_cmap

###############################################################################
#SETUP IRAF STUFF
# For Kungpao
x_images = '/Users/fardila/anaconda/envs/hsc_hydro/iraf/bin.macosx/x_images.e'
x_ttools = '/Users/fardila/anaconda/envs/hsc_hydro/iraf_extern/tables/bin.macosx/x_ttools.e'
x_isophote = '/Users/fardila/anaconda/envs/hsc_hydro/iraf_extern/stsdas/bin.macosx/x_isophote.e'
#---------------------------------------------------------------------------#

#---------------------------------------------------------------------------#
# About the Colormaps
IMG_CMAP = plt.get_cmap('viridis')
IMG_CMAP.set_bad(color='black')

SEG_CMAP = random_cmap(ncolors=512, background_color=u'white')
SEG_CMAP.set_bad(color='white')
SEG_CMAP.set_under(color='white')

from pyraf import iraf

iraf.tables()
iraf.stsdas()
iraf.analysis()
iraf.isophote()

iraf.unlearn('ellipse')
iraf.unlearn('bmodel')
###############################################################################

def get_pixel_scale(file):
    f = h5py.File(file, 'r')
    map_size = f['config'].attrs['map_range_min']
    n_pixels = f['config'].attrs['map_npixel']
    pixel_scale=2 * (map_size/n_pixels)

    return pixel_scale

def mu_iso(iso, pixel_scale):
    return 10**(np.log10(iso['intens'] / (pixel_scale**2.0))+ np.log10(0.7 ** 2.0))

def mu_extrap(iso, ini_r=50, final_r=100):
    power_law_iso = np.log10(powerlaw(iso['sma_kpc'],*fit_power_law_to_iso(iso, ini_r, final_r)))
    return 10**(power_law_iso+ np.log10(0.7 ** 2.0))

def get_median_profile(isos, pixel_scale, quantity = 'intens', rmin=0.05, rmax=4.7, nbin=150):
    """Get the median profiles."""
    sma_common = np.linspace(rmin, rmax, nbin)

    if quantity == 'intens':
        mu = np.nanmedian(np.stack([interp1d((gal['sma'] * pixel_scale) ** 0.25,
                                               np.log10(gal[quantity] / (pixel_scale ** 2)),
                                               bounds_error=False,
                                               fill_value=np.nan,
                                               kind='slinear')(sma_common)
                               for gal in isos]), axis=0)
    elif quantity == 'growth_ori':
        mu = np.nanmedian(np.stack([interp1d((gal['sma'] * pixel_scale) ** 0.25,
                                               np.log10(gal[quantity]),
                                               bounds_error=False,
                                               fill_value=np.nan,
                                               kind='slinear')(sma_common)
                               for gal in isos]), axis=0)
    elif quantity == 'ratio':

        mu = np.nanmedian(np.stack([interp1d((gal['sma'] * pixel_scale) ** 0.25,
                                               mu_extrap(gal)/mu_iso(gal, pixel_scale),
                                               bounds_error=False,
                                               fill_value=np.nan,
                                               kind='slinear')(sma_common)
                               for gal in isos]), axis=0)

    return sma_common, mu

def load_pkl(filename):
    try:
        import cPickle as pickle
    except:
        warnings.warn("## cPickle is not available!!")
        import pickle

    if os.path.isfile(filename):
        pklFile = open(filename, 'rb')
        data = pickle.load(pklFile)
        pklFile.close()

        return data
    else:
        warnings.warn("## Can not find %s, return None" % filename)
        return None

def open_pkl(file_name):
    pkl = open(file_name,'rb')
    array = pickle.load(pkl)
    pkl.close()
    return array


def save_to_fits(image, name):
    """
    Save a 2-D array as fits image.
    """
    hdu = fits.PrimaryHDU(image)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(name, overwrite=True)

    return


def show_illustris(img_cen, img_sat, img_icl, img_all, pixel_scale):
    """
    Show the images of Illustris simulated galaxy.
    """
    fig = plt.figure(figsize=(12, 12))
    fig.subplots_adjust(hspace=0.0, wspace=0.0,
                        top=1.0, right=1.0,
                        left=0.0, bottom=0.0)

    ax1 = plt.subplot(2, 2, 1)
    ax1 = display_single(img_all, ax=ax1,
                         contrast=0.15,
                         scale_bar_length=20.0,
                         scale_bar_loc='right',
                         stretch='log10',
                         pixel_scale=1.0,
                         physical_scale=pixel_scale,
                         color_bar=True)

    ax2 = plt.subplot(2, 2, 2)
    ax2 = display_single(img_cen, ax=ax2,
                         contrast=0.10,
                         scale_bar_length=20.0,
                         scale_bar_loc='right',
                         stretch='log10',
                         pixel_scale=1.0,
                         physical_scale=pixel_scale,
                         color_bar=True)

    ax3 = plt.subplot(2, 2, 3)
    ax3 = display_single(img_sat, ax=ax3,
                         contrast=0.12,
                         scale_bar_length=20.0,
                         scale_bar_loc='right',
                         stretch='log10',
                         pixel_scale=1.0,
                         physical_scale=pixel_scale,
                         color_bar=True)

    ax4 = plt.subplot(2, 2, 4)
    ax4 = display_single(img_icl, ax=ax4,
                         contrast=0.01,
                         scale_bar_length=20.0,
                         scale_bar_loc='right',
                         stretch='log10',
                         pixel_scale=1.0,
                         physical_scale=pixel_scale,
                         color_bar=True)

    return fig

def fit_isophotes(image, iso, iso_bin, prefix, suffix, stage, pixel_scale):
    # Here, we use the stsdas.isophote.analysis.bmodel function to reconstruct
    # a 2-D model using the isophote information
    # So that we can subtract it from the original image, and see how well it does.
    try:
        os.remove(prefix + suffix + '_ellip_' + str(stage) + '.fits')
    except Exception:
        pass

    iraf.bmodel(parent=prefix + suffix + '.fits',
                table=iso_bin,
                output=prefix + suffix + '_ellip_' + str(stage) + '.fits',
                minsma=0.0,
                highar='no')

    img_ellip = fits.open(prefix + suffix + '_ellip_'+str(stage)+'.fits')[0].data

    fig = plt.figure(figsize=(8, 8))
    fig.subplots_adjust(hspace=0.0, wspace=0.0,
                        left=0.0, bottom=0.0,
                        top=1.0, right=1.0)
    ax1 = fig.add_subplot(1, 1, 1)

    ax1 = display_single(image - img_ellip,
                         ax=ax1,
                         contrast=0.15,
                         scale_bar_length=20.0,
                         scale_bar_loc='right',
                         stretch='linear',
                         pixel_scale=1.0,
                         physical_scale=pixel_scale,
                         color_bar=True)

    # Overplot a subsample of isophotes on the image
    iso_ellip = galSBP.convIso2Ell(iso)
    for ii, e in enumerate(iso_ellip):
        if (ii % 3 == 0):
            ax1.add_artist(e)
            e.set_clip_box(ax1.bbox)
            e.set_alpha(0.6)
            e.set_edgecolor('r')
            e.set_facecolor('none')
            e.set_linewidth(2.0)

    return iso_ellip


def check_profile(iso):
    # Normally it is better to also check the profiles for geometry
    fig = plt.figure(figsize=(8, 7))
    fig.subplots_adjust(hspace=0.0, wspace=0.0,
                        left=0.13, bottom=0.10,
                        top=0.97, right=0.97)

    ax1 = plt.subplot(4, 1, 1)
    ax1.grid(linewidth=2.0, linestyle='--', alpha=0.5)
    ax1.plot(iso['sma'], iso['x0'], linewidth=3.0)
    ax1.set_ylabel(r'$\mathrm{X}_0$', fontsize=20)

    ax2 = plt.subplot(4, 1, 2)
    ax2.grid(linewidth=2.0, linestyle='--', alpha=0.5)
    ax2.plot(iso['sma'], iso['y0'], linewidth=3.0)
    ax2.set_ylabel(r'$\mathrm{Y}_0$', fontsize=20)

    ax3 = plt.subplot(4, 1, 3)
    ax3.grid(linewidth=2.0, linestyle='--', alpha=0.5)
    ax3.plot(iso['sma'], iso['ell'], linewidth=3.0)
    ax3.set_ylabel(r'$e$', fontsize=20)

    ax4 = plt.subplot(4, 1, 4)
    ax4.grid(linewidth=2.0, linestyle='--', alpha=0.5)
    ax4.plot(iso['sma'], iso['pa'], linewidth=3.0)
    ax4.set_ylabel(r'$\mathrm{PA}$', fontsize=20)
    ax4.set_xlabel(r'$\mathrm{SMA}$', fontsize=25)

###############################################################################
#get_masses imports
def oneD_mass(galaxy_iso, radius_px):
    mass=np.interp(radius_px,galaxy_iso['sma'], galaxy_iso['growth_ori'])
    return np.log10(mass)

def find_closest(objects, x0=100., y0=100.):
    xs=objects['x']
    ys=objects['y']
    distances=np.sqrt(((xs-x0)**2) + ((ys-y0)**2))

    closest_index=np.argmin(distances)

    return objects[closest_index]

# def powerlaw(x, m, c, c0):
#     return c0 + x**m * c
def powerlaw(x, m, c):
    return x**m * c

def ellipse_perimeter(a,b):
    '''perimeter of an ellipse from Ramanujan's approximation'''
    perimeter = np.pi * ( 3*(a+b) - np.sqrt( (3*a + b) * (a + 3*b) ) )
    return perimeter

def mass_per_r(r, power_law_params, e):
    q=1-e
    b=r*q
    return powerlaw(r,*power_law_params)*ellipse_perimeter(r,b)

def fit_power_law_to_iso(iso, ini_r, final_r):

    x=iso['sma_kpc'][(iso['sma_kpc']>ini_r) & (iso['sma_kpc']<final_r)]
    y=iso['intens_kpc'][(iso['sma_kpc']>ini_r) & (iso['sma_kpc']<final_r)]

    p_fit_power, _ = curve_fit(powerlaw, x, y, p0=[-2,10**10])

    return p_fit_power

def extrapolated_1D_mass(iso, max_r):

    #fit power law between 50 and 100kpc
    p_fit_power = fit_power_law_to_iso(iso, 50, 100)
    #print p_fit_power

    #integrate between 100 and max_r to get mass
    e = iso['ell'][1]
    mass_beyond100, mass_err =quad(mass_per_r, 100, max_r, args=tuple([p_fit_power, e]))

    #add to galSBP mass to get total mass
    mass_100 = 10**(oneD_mass(iso, 100))
    power_law_fit_mass = np.log10(mass_100+mass_beyond100)

    return power_law_fit_mass

def get_mass_maps(sim_file, gal_n=0):

    # Load general simulation and galaxy properties
    f = h5py.File(sim_file, 'r')
    cat_sh_mstar = np.array(f['catsh_SubhaloMassType'][:, 4])

    cen_insitu = np.array(f['map_star_rho_insitu_xy'])
    cen_exsitu = np.array(f['map_star_rho_exsitu_xy'])
    map_stars_cen = cen_exsitu + cen_insitu

    fuzz_insitu = np.array(f['map_star_rho_fuzz_insitu_xy'])
    fuzz_exsitu = np.array(f['map_star_rho_fuzz_exsitu_xy'])
    map_stars_fuzz = fuzz_exsitu + fuzz_insitu

    sats_insitu = np.array(f['map_star_rho_oshs_insitu_xy'])
    sats_exsitu = np.array(f['map_star_rho_oshs_exsitu_xy'])
    map_stars_sats = sats_exsitu + sats_insitu

    map_size = f['config'].attrs['map_range_min']
    n_pixels = f['config'].attrs['map_npixel']
    pixel_scale=2 * (map_size/n_pixels)

    f.close()

    #make maps
    img_cen = map_stars_cen[gal_n] * (pixel_scale ** 2) # Central
    img_sat = map_stars_sats[gal_n] * (pixel_scale ** 2) # Satellites
    img_icl = map_stars_fuzz[gal_n] * (pixel_scale ** 2) # Diffuse
    img_cen_sat = (img_cen + img_sat)           # Central + Satellites
    img_cen_icl = (img_cen + img_icl)           # Central + Satellites
    img_all = (img_cen + img_sat + img_icl)           # Central + Satellites + Diffuse

    #catalog mass
    m_cat = np.log10(cat_sh_mstar[gal_n])

    return img_cen, img_cen_icl, pixel_scale, m_cat

def get_masses_iso(sim_file, sim_name, resolution, intMode='mean', components='cen', gal_n=0):

    # Load maps
    mass_map_cen, mass_map_cen_icl, pixel_scale, m_cat = get_mass_maps(sim_file, gal_n=gal_n)

    #postage mass
    m_post = np.log10(np.sum(mass_map_cen))
    m_post_icl = np.log10(np.sum(mass_map_cen_icl))


    #ouput maps
    maps_location='/Users/fardila/Documents/GitHub/HSC_vs_hydro/Figures/fits_files/{0}/'.format(resolution)

    file_name=sim_name+'_'+str(gal_n)+'_xy'
    fits_prefix = maps_location + file_name

    #components
    if components == 'cen':
        save_to_fits(mass_map_cen, fits_prefix + '_cen.fits')
        data=mass_map_cen
    elif components == 'cen+icl':
        save_to_fits(mass_map_cen_icl, fits_prefix + '_cen+icl.fits')
        data=mass_map_cen_icl
    else:
        raise ValueError('only cen or cen+icl allowed for now')

    # save_to_fits(mass_map_cen, fits_prefix + '_cen.fits')
    # save_to_fits(img_cen_sat, fits_prefix + '_cen_sat.fits')
    # save_to_fits(img_cen_icl, fits_prefix + '_cen_icl.fits')
    # save_to_fits(img_all, fits_prefix + '_all.fits')

    suffix='_'+components

    #central pixels
    x0=len(data)/2.
    y0=len(data)/2.

    ###########################################################################
    #ellipse information
    #########################
    if resolution == 'highres':
        #get ellipse information from iso file of quick
        iso_quick = load_pkl('/Users/fardila/Documents/GitHub/HSC_vs_hydro/Figures/fits_files/quick/{0}_{1}_xy_cen_ellip_3.pkl'.format(sim_name, gal_n))
        q = 1- iso_quick['ell'][-1]
        theta = iso_quick['pa'][-1]* np.pi /180.
    elif resolution == 'quick':
        #get background
        bkg = sep.Background(data, bw=10, bh=10, fw=5, fh=5)
        bkg_subtraced_data = data - bkg

        thresh = 50 * bkg.globalrms
        objects = sep.extract(bkg_subtraced_data, thresh, minarea = 100,
                              deblend_nthresh=24, deblend_cont=0.1)

        #find object closest to image center
        obj = find_closest(objects, x0=x0, y0=y0)

        #ellipse parameters
        theta = obj['theta']
        q = obj['b']/ obj['a']
    #########################


    a_10, a_30, a_100 = (10. / pixel_scale), (30. / pixel_scale), (100. / pixel_scale)
    b_10, b_30, b_100 =  a_10 * q, a_30 * q, a_100 * q



    # plot background-subtracted image
    m, s = np.mean(data), np.std(data)
    fig, ax = plt.subplots()
    im = ax.imshow(data, interpolation='nearest', cmap=plt.get_cmap('viridis'),
                   vmin=m-s, vmax=m+s, origin='lower')

    # plot an ellipse for each object
    e_30 = Ellipse(xy=(x0, y0),
                 width=a_30,
                 height=b_30,
                 angle=theta * 180. / np.pi)
    e_30.set_facecolor('none')
    e_30.set_edgecolor('red')
    ax.add_artist(e_30)

    e_100 = Ellipse(xy=(x0, y0),
                 width=a_100,
                 height=b_100,
                 angle=theta * 180. / np.pi)
    e_100.set_facecolor('none')
    e_100.set_edgecolor('red')
    ax.add_artist(e_100)


    plt.savefig('/Users/fardila/Documents/GitHub/HSC_vs_hydro/Figures/ellipses/{0}/{1}'.format(resolution, file_name))
    plt.clf()

    ###########################################################################
    #2D masses
    flux_10, fluxerr_10, flag_10 = sep.sum_ellipse(data, x0, y0,
                                                   a_10, b_10, theta)
    flux_30, fluxerr_30, flag_30 = sep.sum_ellipse(data, x0, y0,
                                                   a_30, b_30, theta)
    flux_100, fluxerr_100, flag_100 = sep.sum_ellipse(data, x0, y0,
                                                      a_100, b_100, theta)

    ###########################################################################
    #1D masses from galSBP
#     try:
    iso, iso_bin = galSBP.galSBP(maps_location+file_name+suffix+'.fits',
                                     galX=x0,
                                     galY=y0,
                                     galQ=q,
                                     galPA=theta* 180. / np.pi,
                                     maxSma=250,
                                     iniSma=50.0,
                                     stage=3,
                                     intMode=intMode,
                                     ellipStep=0.05,
                                     pix=pixel_scale,
                                     zpPhoto=0.0,
                                     isophote=x_isophote,
                                     xttools=x_ttools,
                                     recenter=True,
                                     savePng=False,
                                     verbose=True,
                                     uppClip=3.0,
                                     lowClip=3.0,
                                     nClip=2)


    ###########################################################################
    iso['sma_kpc'] = iso['sma'] * pixel_scale
    iso['intens_kpc']=iso['intens'] / (pixel_scale**2)

    m_1d_10, m_1d_30, m_1d_100 = oneD_mass(iso, 10.), \
                                oneD_mass(iso, 30.), \
                                oneD_mass(iso, 100.)

    #integrated mass from extrapolation
    extrap_mass = extrapolated_1D_mass(iso, 500)

#     except:
#         iso,m_1d_10, m_1d_30, m_1d_100, extrap_mass  = -99.99, -99.99, -99.99, -99.99, -99.99


    m_2d_10, m_2d_30, m_2d_100 = np.log10(flux_10), \
                                np.log10(flux_30), \
                                np.log10(flux_100)


    masses = [m_cat, m_post, m_post_icl, m_1d_10, m_1d_30, m_1d_100, m_2d_10,
            m_2d_30, m_2d_100, extrap_mass]

    return iso, masses


###############################################################################


#OLD FUNCTIONS
# def get_masses(sim_file, sim_name, gal_n=0):
#
#     # Load general simulation and galaxy properties
#     f = h5py.File(sim_file, 'r')
#     cat_sh_mstar = np.array(f['cat_sh_mstar'])
#
#     try:
#         map_stars = np.array(f['map_stars'])
#     except:
#         map_stars_insitu = np.array(f['map_stars_insitu'])
#         map_stars_exsitu = np.array(f['map_stars_exsitu'])
#         map_stars = map_stars_exsitu + map_stars_insitu
#
#     map_size = f.attrs['stellar_map_size']
#     n_pixels = f.attrs['stellar_map_np']
#
#     pixel_scale=2 * (map_size/n_pixels)
#     f.close()
#
#     #make maps
#     img_cen = map_stars[gal_n, 0, 1] * (pixel_scale ** 2) # Central
#     img_sat = map_stars[gal_n, 1, 1] * (pixel_scale ** 2) # Satellites
#     img_icl = map_stars[gal_n, 2, 1] * (pixel_scale ** 2) # Diffuse
#     img_cen_sat = (img_cen + img_sat)           # Central + Satellites
#     img_cen_icl = (img_cen + img_icl)           # Central + Satellites
#     img_all = (img_cen + img_sat + img_icl)           # Central + Satellites + Diffuse
#
#     #convert the image into unit of stellar mass instead of mass density
#     m_cat = np.log10(cat_sh_mstar[gal_n])
#     m_post = np.log10(np.sum(img_cen))
#     m_post_icl = np.log10(np.sum(img_cen_icl))
#
#     #ouput maps
#     maps_location='/Users/fardila/Documents/GitHub/HSC_vs_hydro/Figures/fits_files/'
#     file_name=sim_name+'_'+str(gal_n)+'_xy'
#     fits_prefix = maps_location + file_name
#     save_to_fits(img_cen, fits_prefix + '_cen.fits')
#     save_to_fits(img_cen_sat, fits_prefix + '_cen_sat.fits')
#     save_to_fits(img_cen_icl, fits_prefix + '_cen_icl.fits')
#     save_to_fits(img_all, fits_prefix + '_all.fits')
#
#     data=img_cen
#     suffix='_cen'
#
#     ###########################################################################
#     #get background
#     bkg = sep.Background(data, bw=10, bh=10, fw=5, fh=5)
#     bkg_subtraced_data = data - bkg
#
#     thresh = 500 * bkg.globalrms
#     objects = sep.extract(bkg_subtraced_data, thresh, minarea = 100,
#                           deblend_nthresh=24, deblend_cont=0.1)
#
#     #find object closest to image center
#     obj = find_closest(objects)
#
#     #ellipse parameters
#     theta = obj['theta']
#     q = obj['b']/ obj['a']
#
#     a_10, a_30, a_100 = (10. / pixel_scale), (30. / pixel_scale), (100. / pixel_scale)
#     b_10, b_30, b_100 =  a_10 * q, a_30 * q, a_100 * q
#
#
#
#     # plot background-subtracted image
#     m, s = np.mean(data), np.std(data)
#     fig, ax = plt.subplots()
#     im = ax.imshow(data, interpolation='nearest', cmap=plt.get_cmap('viridis'),
#                    vmin=m-s, vmax=m+s, origin='lower')
#
#     # plot an ellipse for each object
#     e_30 = Ellipse(xy=(obj['x'], obj['y']),
#                  width=a_30,
#                  height=b_30,
#                  angle=theta * 180. / np.pi)
#     e_30.set_facecolor('none')
#     e_30.set_edgecolor('red')
#     ax.add_artist(e_30)
#
#     e_100 = Ellipse(xy=(obj['x'], obj['y']),
#                  width=a_100,
#                  height=b_100,
#                  angle=theta * 180. / np.pi)
#     e_100.set_facecolor('none')
#     e_100.set_edgecolor('red')
#     ax.add_artist(e_100)
#
#
#     plt.savefig('/Users/fardila/Documents/GitHub/HSC_vs_hydro/Figures/ellipses/'+file_name)
#     plt.clf()
#
#     ###########################################################################
#     #2D masses
#     flux_10, fluxerr_10, flag_10 = sep.sum_ellipse(data, 100., 100.,
#                                                    a_10, b_10, theta)
#     flux_30, fluxerr_30, flag_30 = sep.sum_ellipse(data, 100., 100.,
#                                                    a_30, b_30, theta)
#     flux_100, fluxerr_100, flag_100 = sep.sum_ellipse(data, 100., 100.,
#                                                       a_100, b_100, theta)
#
#     ###########################################################################
#     #1D masses from galSBP
#     try:
#         iso, iso_bin = galSBP.galSBP(maps_location+file_name+suffix+'.fits',
#                                          galX=100.,
#                                          galY=100.,
#                                          galQ=q,
#                                          galPA=theta* 180. / np.pi,
#                                          maxSma=250,
#                                          iniSma=50.0,
#                                          stage=3,
#                                          intMode='median',
#                                          ellipStep=0.05,
#                                          pix=pixel_scale,
#                                          zpPhoto=0.0,
#                                          isophote=x_isophote,
#                                          xttools=x_ttools,
#                                          recenter=True,
#                                          savePng=False,
#                                          verbose=True)
#
#
#         ###########################################################################
#         iso['sma_kpc'] = iso['sma'] * pixel_scale
#         iso['intens_kpc']=iso['intens'] / (pixel_scale**2)
#
#         m_1d_10, m_1d_30, m_1d_100 = oneD_mass(iso, 10.), \
#                                     oneD_mass(iso, 30.), \
#                                     oneD_mass(iso, 100.)
#
#         #integrated mass from extrapolation
#         extrap_mass = extrapolated_1D_mass(iso, 500)
#
#     except:
#         iso,m_1d_10, m_1d_30, m_1d_100, extrap_mass  = -99.99, -99.99, -99.99, -99.99, -99.99
#
#
#     m_2d_10, m_2d_30, m_2d_100 = np.log10(flux_10), \
#                                 np.log10(flux_30), \
#                                 np.log10(flux_100)
#
#
#     masses = [m_cat, m_post, m_post_icl, m_1d_10, m_1d_30, m_1d_100, m_2d_10,
#             m_2d_30, m_2d_100, extrap_mass]
#
#     return iso, masses


# def get_iso(sim_file, sim_name, components='cen', gal_n=0):
#     #central pixels
#     x0=150.
#     y0=150.
#
#     # Load maps
#     mass_map_cen, mass_map_cen_icl, pixel_scale, m_cat = get_mass_maps(sim_file, gal_n=gal_n)
#
#     #postage mass
#     m_post = np.log10(np.sum(mass_map_cen))
#     m_post_icl = np.log10(np.sum(mass_map_cen_icl))
#
#
#     #ouput maps
#     maps_location='/Users/fardila/Documents/GitHub/HSC_vs_hydro/Figures/fits_files/quick_800/'
#
#     file_name=sim_name+'_'+str(gal_n)+'_xy'
#     fits_prefix = maps_location + file_name
#
#     if components == 'cen':
#         save_to_fits(mass_map_cen, fits_prefix + '_cen.fits')
#         data=mass_map_cen
#     elif components == 'cen+icl':
#         save_to_fits(mass_map_cen_icl, fits_prefix + '_cen+icl.fits')
#         data=mass_map_cen_icl
#     else:
#         raise ValueError('only cen or cen+icl allowed for now')
#     # save_to_fits(img_cen_sat, fits_prefix + '_cen+sat.fits')
#     # save_to_fits(img_all, fits_prefix + '_all.fits')
#
#
#     suffix='_'+components
#
#     ###########################################################################
#     #get background
#     bkg = sep.Background(data, bw=10, bh=10, fw=5, fh=5)
#     bkg_subtraced_data = data - bkg
#
#     thresh = 50 * bkg.globalrms
#     objects = sep.extract(bkg_subtraced_data, thresh, minarea = 100,
#                           deblend_nthresh=24, deblend_cont=0.1)
#
#     #find object closest to image center
#     obj = find_closest(objects, x0=x0, y0=y0)
#
#     #ellipse parameters
#     theta = obj['theta']
#     q = obj['b']/ obj['a']
#
#     a_10, a_30, a_100 = (10. / pixel_scale), (30. / pixel_scale), (100. / pixel_scale)
#     b_10, b_30, b_100 =  a_10 * q, a_30 * q, a_100 * q
#
#     ###########################################################################
#     #1D masses from galSBP
#     try:
#         iso, iso_bin = galSBP.galSBP(maps_location+file_name+suffix+'.fits',
#                                          galX=x0,
#                                          galY=y0,
#                                          galQ=q,
#                                          galPA=theta* 180. / np.pi,
#                                          maxSma=250,
#                                          iniSma=50.0,
#                                          stage=3,
#                                          intMode='median',
#                                          ellipStep=0.05,
#                                          pix=pixel_scale,
#                                          zpPhoto=0.0,
#                                          isophote=x_isophote,
#                                          xttools=x_ttools,
#                                          recenter=True,
#                                          savePng=False,
#                                          verbose=True)
#
#
#         ###########################################################################
#         iso['sma_kpc'] = iso['sma'] * pixel_scale
#         iso['intens_kpc']=iso['intens'] / (pixel_scale**2)
#
#     except ValueError:
#         iso = -99.99
#
#
#     return iso

# def get_masses_iso(sim_file, sim_name, intMode='median', components='cen', gal_n=0):
#     #central pixels
#     x0=150.
#     y0=150.
#
#     # Load maps
#     mass_map_cen, mass_map_cen_icl, pixel_scale, m_cat = get_mass_maps(sim_file,
#      gal_n=gal_n)
#
#     #postage mass
#     m_post = np.log10(np.sum(mass_map_cen))
#     m_post_icl = np.log10(np.sum(mass_map_cen_icl))
#
#
#     #ouput maps
#     maps_location='/Users/fardila/Documents/GitHub/HSC_vs_hydro/Figures/fits_files/quick_800/'
#
#     file_name=sim_name+'_'+str(gal_n)+'_xy'
#     fits_prefix = maps_location + file_name
#
#     #components
#     if components == 'cen':
#         save_to_fits(mass_map_cen, fits_prefix + '_cen.fits')
#         data=mass_map_cen
#     elif components == 'cen+icl':
#         save_to_fits(mass_map_cen_icl, fits_prefix + '_cen+icl.fits')
#         data=mass_map_cen_icl
#     else:
#         raise ValueError('only cen or cen+icl allowed for now')
#
#     # save_to_fits(mass_map_cen, fits_prefix + '_cen.fits')
#     # save_to_fits(img_cen_sat, fits_prefix + '_cen_sat.fits')
#     # save_to_fits(img_cen_icl, fits_prefix + '_cen_icl.fits')
#     # save_to_fits(img_all, fits_prefix + '_all.fits')
#
#     suffix='_'+components
#
#     ###########################################################################
#     #get background
#     bkg = sep.Background(data, bw=10, bh=10, fw=5, fh=5)
#     bkg_subtraced_data = data - bkg
#
#     thresh = 50 * bkg.globalrms
#     objects = sep.extract(bkg_subtraced_data, thresh, minarea = 100,
#                           deblend_nthresh=24, deblend_cont=0.1)
#
#     #find object closest to image center
#     obj = find_closest(objects, x0=x0, y0=y0)
#
#     #ellipse parameters
#     theta = obj['theta']
#     q = obj['b']/ obj['a']
#
#     a_10, a_30, a_100 = (10. / pixel_scale), (30. / pixel_scale), (100. / pixel_scale)
#     b_10, b_30, b_100 =  a_10 * q, a_30 * q, a_100 * q
#
#
#
#     # plot background-subtracted image
#     m, s = np.mean(data), np.std(data)
#     fig, ax = plt.subplots()
#     im = ax.imshow(data, interpolation='nearest', cmap=plt.get_cmap('viridis'),
#                    vmin=m-s, vmax=m+s, origin='lower')
#
#     # plot an ellipse for each object
#     e_30 = Ellipse(xy=(obj['x'], obj['y']),
#                  width=a_30*2,
#                  height=b_30*2,
#                  angle=theta * 180. / np.pi)
#     e_30.set_facecolor('none')
#     e_30.set_edgecolor('red')
#     ax.add_artist(e_30)
#
#     e_100 = Ellipse(xy=(obj['x'], obj['y']),
#                  width=a_100*2,
#                  height=b_100*2,
#                  angle=theta * 180. / np.pi)
#     e_100.set_facecolor('none')
#     e_100.set_edgecolor('red')
#     ax.add_artist(e_100)
#
#
#     plt.savefig('/Users/fardila/Documents/GitHub/HSC_vs_hydro/Figures/ellipses/quick_800/'+file_name)
#     plt.clf()
#
#     ###########################################################################
#     #2D masses
#     flux_10, fluxerr_10, flag_10 = sep.sum_ellipse(data, x0, y0,
#                                                    a_10, b_10, theta)
#     flux_30, fluxerr_30, flag_30 = sep.sum_ellipse(data, x0, y0,
#                                                    a_30, b_30, theta)
#     flux_100, fluxerr_100, flag_100 = sep.sum_ellipse(data, x0, y0,
#                                                       a_100, b_100, theta)
#
#     ###########################################################################
#     #1D masses from galSBP
#     try:
#         iso, iso_bin = galSBP.galSBP(maps_location+file_name+suffix+'.fits',
#                                          galX=x0,
#                                          galY=y0,
#                                          galQ=q,
#                                          galPA=theta* 180. / np.pi,
#                                          maxSma=250,
#                                          iniSma=50.0,
#                                          stage=3,
#                                          intMode=intMode,
#                                          ellipStep=0.05,
#                                          pix=pixel_scale,
#                                          zpPhoto=0.0,
#                                          isophote=x_isophote,
#                                          xttools=x_ttools,
#                                          recenter=True,
#                                          savePng=False,
#                                          verbose=True)
#
#
#         ###########################################################################
#         iso['sma_kpc'] = iso['sma'] * pixel_scale
#         iso['intens_kpc']=iso['intens'] / (pixel_scale**2)
#
#         m_1d_10, m_1d_30, m_1d_100 = oneD_mass(iso, 10.), \
#                                     oneD_mass(iso, 30.), \
#                                     oneD_mass(iso, 100.)
#
#         #integrated mass from extrapolation
#         extrap_mass = extrapolated_1D_mass(iso, 800)
#
#     except:
#         iso,m_1d_10, m_1d_30, m_1d_100, extrap_mass  = -99.99, -99.99, -99.99, -99.99, -99.99
#
#
#     m_2d_10, m_2d_30, m_2d_100 = np.log10(flux_10), \
#                                 np.log10(flux_30), \
#                                 np.log10(flux_100)
#
#
#     masses = [m_cat, m_post, m_post_icl, m_1d_10, m_1d_30, m_1d_100, m_2d_10,
#             m_2d_30, m_2d_100, extrap_mass]
#
#     return iso, masses

# def get_median_profile(isos, pixel_scale, quantity = 'intens', rmin=0.05, rmax=4.7, nbin=150):
#     """Get the median profiles."""
#     sma_common = np.linspace(rmin, rmax, nbin)
#
#     if quantity == 'intens':
#         mu = np.nanmedian(np.stack([interp1d((gal['sma'] * pixel_scale) ** 0.25,
#                                                np.log10(gal[quantity] / (pixel_scale ** 2)),
#                                                bounds_error=False,
#                                                fill_value=np.nan,
#                                                kind='slinear')(sma_common)
#                                for gal in isos]), axis=0)
#     elif quantity == 'growth_ori':
#         mu = np.nanmedian(np.stack([interp1d((gal['sma'] * pixel_scale) ** 0.25,
#                                                np.log10(gal[quantity]),
#                                                bounds_error=False,
#                                                fill_value=np.nan,
#                                                kind='slinear')(sma_common)
#                                for gal in isos]), axis=0)
#
#
#     return sma_common, mu
