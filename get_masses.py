from __future__ import print_function, division, absolute_import

import pickle
from matplotlib.patches import Ellipse


###############################################################################
#SBP fitting imports

import sep
import h5py
import numpy as np
from functions import *

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

import time
time0=time.time()

###############################################################################
def oneD_mass(galaxy_iso, radius_px):
    mass=np.interp(radius_px,galaxy_iso['sma'], galaxy_iso['growth_ori'])
    return np.log10(mass)

def find_closest(objects, x0=100., y0=100.):
    xs=objects['x']
    ys=objects['y']
    distances=np.sqrt(((xs-x0)**2) + ((ys-y0)**2))

    closest_index=np.argmin(distances)
    
    return objects[closest_index]

def get_masses(sim_file, sim_name, gal_n=0):

    # Load general simulation and galaxy properties
    f = h5py.File(sim_file, 'r')
    cat_sh_mstar = np.array(f['cat_sh_mstar'])

    try:
        map_stars = np.array(f['map_stars'])
    except:
        map_stars_insitu = np.array(f['map_stars_insitu'])
        map_stars_exsitu = np.array(f['map_stars_exsitu'])
        map_stars = map_stars_exsitu + map_stars_insitu

    map_size = f.attrs['stellar_map_size']
    n_pixels = f.attrs['stellar_map_np']

    pixel_scale=2 * (map_size/n_pixels)
    f.close()

    #make maps
    img_cen = map_stars[gal_n, 0, 1] * (pixel_scale ** 2) # Central
    img_sat = map_stars[gal_n, 1, 1] * (pixel_scale ** 2) # Satellites
    img_icl = map_stars[gal_n, 2, 1] * (pixel_scale ** 2) # Diffuse
    img_cen_sat = (img_cen + img_sat)           # Central + Satellites
    img_cen_icl = (img_cen + img_icl)           # Central + Satellites
    img_all = (img_cen + img_sat + img_icl)           # Central + Satellites + Diffuse

    #convert the image into unit of stellar mass instead of mass density
    m_cat = np.log10(cat_sh_mstar[gal_n])
    m_post = np.log10(np.sum(img_cen))
    m_post_icl = np.log10(np.sum(img_cen_icl))

    #ouput maps
    maps_location='/Users/fardila/Documents/GitHub/HSC_vs_hydro/Figures/fits_files/'
    file_name=sim_name+'_'+str(gal_n)+'_xy'
    fits_prefix = maps_location + file_name
    save_to_fits(img_cen, fits_prefix + '_cen.fits')
    save_to_fits(img_cen_sat, fits_prefix + '_cen_sat.fits')
    save_to_fits(img_cen_icl, fits_prefix + '_cen_icl.fits')
    save_to_fits(img_all, fits_prefix + '_all.fits')

    data=img_cen
    suffix='_cen'

    ###########################################################################
    #get background
    bkg = sep.Background(data, bw=10, bh=10, fw=5, fh=5)
    bkg_subtraced_data = data - bkg

    thresh = 500 * bkg.globalrms
    objects = sep.extract(bkg_subtraced_data, thresh, minarea = 100,
                          deblend_nthresh=24, deblend_cont=0.1)

    #find object closest to image center
    obj = find_closest(objects)

    #ellipse parameters
    theta = obj['theta'][0]
    q = obj['b'][0] / obj['a'][0]

    a_10, a_30, a_100 = (10. / pixel_scale), (30. / pixel_scale), (100. / pixel_scale)
    b_10, b_30, b_100 =  a_10 * q, a_30 * q, a_100 * q



    # plot background-subtracted image
    m, s = np.mean(data), np.std(data)
    fig, ax = plt.subplots()
    im = ax.imshow(data, interpolation='nearest', cmap=plt.get_cmap('viridis'),
                   vmin=m-s, vmax=m+s, origin='lower')

    # plot an ellipse for each object
    e_30 = Ellipse(xy=(obj['x'][0], obj['y'][0]),
                 width=a_30,
                 height=b_30,
                 angle=theta * 180. / np.pi)
    e_30.set_facecolor('none')
    e_30.set_edgecolor('red')
    ax.add_artist(e_30)

    e_100 = Ellipse(xy=(obj['x'][0], obj['y'][0]),
                 width=a_100,
                 height=b_100,
                 angle=theta * 180. / np.pi)
    e_100.set_facecolor('none')
    e_100.set_edgecolor('red')
    ax.add_artist(e_100)


    plt.savefig('/Users/fardila/Documents/GitHub/HSC_vs_hydro/Figures/ellipses/'+file_name)
    plt.clf()

    ###########################################################################
    #2D masses
    flux_10, fluxerr_10, flag_10 = sep.sum_ellipse(data, 100., 100.,
                                                   a_10, b_10, theta)
    flux_30, fluxerr_30, flag_30 = sep.sum_ellipse(data, 100., 100.,
                                                   a_30, b_30, theta)
    flux_100, fluxerr_100, flag_100 = sep.sum_ellipse(data, 100., 100.,
                                                      a_100, b_100, theta)

    ###########################################################################
    #1D masses from galSBP
    iso, iso_bin = galSBP.galSBP(maps_location+file_name+suffix+'.fits',
                                     galX=100.,
                                     galY=100.,
                                     galQ=q,
                                     galPA=theta* 180. / np.pi,
                                     maxSma=250,
                                     iniSma=50.0,
                                     stage=3,
                                     intMode='median',
                                     ellipStep=0.05,
                                     pix=pixel_scale,
                                     zpPhoto=0.0,
                                     isophote=x_isophote,
                                     xttools=x_ttools,
                                     recenter=True,
                                     savePng=False,
                                     verbose=True)


    ###########################################################################
    iso['sma_kpc'] = iso['sma'] * pixel_scale
    iso['intens_kpc']=iso['intens'] / (pixel_scale**2)

    m_1d_10, m_1d_30, m_1d_100 = oneD_mass(iso, 10./pixel_scale), \
                                oneD_mass(iso, 30./pixel_scale), \
                                oneD_mass(iso, 100./pixel_scale)

    m_2d_10, m_2d_30, m_2d_100 = np.log10(flux_10[0]), \
                                np.log10(flux_30[0]), \
                                np.log10(flux_100[0])



    masses = [m_cat, m_post, m_post_icl, m_1d_10, m_1d_30, m_1d_100, m_2d_10,
            m_2d_30, m_2d_100]

    return iso, masses



###############################################################################
#data
Illustris_file = '/Users/fardila/Documents/GitHub/HSC_vs_hydro/Data/Illustris/galaxies_orig_11.2.hdf5'
TNG_file = '/Users/fardila/Documents/GitHub/HSC_vs_hydro/Data/TNG/galaxies_tng75_11.2.hdf5'

###############################################################################
#run on Illustris
isos_illustris=[]
masses_illustris=[]

for i in range(339):
    print('^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^')
    print('^^^^^^^^GALAXY '+str(i)+'^^^^^^^^^^^^^^')
    print('^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^')

    try:
        iso, masses = get_masses(Illustris_file,'Illustris',gal_n=i)

    except:
        iso=-99.99
        masses=-99.99


    isos_illustris.append(iso)
    masses_illustris.append(masses)

#save as pickles
pkl_isos = open('/Users/fardila/Documents/GitHub/HSC_vs_hydro/Data/Illustris/Illustris_isos.pkl','wb')
pickle.dump(isos_illustris,pkl_isos)
pkl_isos.close()

pkl_masses = open('/Users/fardila/Documents/GitHub/HSC_vs_hydro/Data/Illustris/Illustris_masses.pkl','wb')
pickle.dump(masses_illustris,pkl_masses)
pkl_masses.close()

###############################################################################
#run on TNG
isos_tng=[]
masses_tng=[]
for i in range(235):

    print('^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^')
    print('^^^^^^^^GALAXY '+str(i)+'^^^^^^^^^^^^^^')
    print('^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^')

    try:
        iso, masses = get_masses(TNG_file,'TNG',gal_n=i)

    except:
        iso=-99.99
        masses=-99.99



    isos_tng.append(iso)
    masses_tng.append(masses)

#save as pickles
pkl_isos = open('/Users/fardila/Documents/GitHub/HSC_vs_hydro/Data/TNG/TNG_isos.pkl','wb')
pickle.dump(isos_tng,pkl_isos)
pkl_isos.close()

pkl_masses = open('/Users/fardila/Documents/GitHub/HSC_vs_hydro/Data/TNG/TNG_masses.pkl','wb')
pickle.dump(masses_tng,pkl_masses)
pkl_masses.close()

print(time.time()-time0, ' seconds')
