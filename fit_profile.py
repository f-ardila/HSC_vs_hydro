def fit_profile(sim_file, pixel_scale, gal_n=0, cen=True, sat=False, icl=False, plots=True):

    # Load general simulation and galaxy properties
    f = h5py.File(sim_file, 'r')
    snap_a = f.attrs['snap_a'] #scale factor!
    snap_z = f.attrs['snap_z'] #redshift!
    n_pixels = f.attrs['stellar_map_np']
    center_pixel=n_pixels/2.


    n_galaxies = len(f['cat_sh_id'])
    sh_idx = np.array(f['cat_sh_id'])
    grp_idx = np.array(f['cat_grp_id'])
    cat_is_primary = np.array(f['cat_grp_is_primary'])
    cat_sh_mstar = np.array(f['cat_sh_mstar'])
    cat_sh_pos_bound = np.array(f['cat_sh_pos_bound'])
    cat_sh_halfmrad_stars = np.array(f['cat_sh_halfmrad_stars'])

    try:
        map_stars = np.array(f['map_stars'])
    except:
        map_stars_insitu = np.array(f['map_stars_insitu'])
        map_stars_exsitu = np.array(f['map_stars_exsitu'])
        map_stars = map_stars_exsitu + map_stars_insitu
    map_size = f.attrs['stellar_map_size']
    n_pixels = f.attrs['stellar_map_np']
    f.close()

    #make maps
    img_cen = map_stars[gal_n, 0, 1] * (pixel_scale ** 2) # Central
    img_sat = map_stars[gal_n, 1, 1] * (pixel_scale ** 2) # Satellites
    img_icl = map_stars[gal_n, 2, 1] * (pixel_scale ** 2) # Diffuse
    img_cen_sat = (img_cen + img_sat)           # Central + Satellites
    img_cen_icl = (img_cen + img_icl)           # Central + Satellites
    img_all = (img_cen + img_sat + img_icl)           # Central + Satellites + Diffuse

    #convert the image into unit of stellar mass instead of mass density
    log_mstar = np.log10(cat_sh_mstar[0])
    log_mcen = np.log10(np.sum(img_cen))

    #ouput maps
    maps_location='/Users/fardila/Documents/GitHub/HSC_vs_hydro/notebooks/felipe_test/maps/'
    fits_prefix = maps_location+'illustris_1_xy'
    save_to_fits(img_cen, fits_prefix + '_cen.fits')
    save_to_fits(img_cen_sat, fits_prefix + '_cen_sat.fits')
    save_to_fits(img_cen_icl, fits_prefix + '_cen_icl.fits')
    save_to_fits(img_all, fits_prefix + '_all.fits')

    #determine which map to use
    if cen:
        if sat:
            if icl:
                image=img_all
                print('Central+Satellites+ICL')
                suffix='_all'
            else:
                image=img_cen_sat
                print('Central+Satellites')
                suffix='_cen_sat'
        elif icl:
            image=img_cen_icl
            print('Central+ICL')
            suffix='_cen_icl'
        else:
            image=img_cen
            print('Central only')
            suffix='_cen'
    else:
        raise ValueError('Are you sure you don\'t want to use the central???')

    # Measure the background,
    # Here on the image has no noise and has a lot of diffuse features
    # so measure sky using a very small box (bw, bh values) helps us remove
    # a lot of the them, make the detection easier
    # See: http://sep.readthedocs.io/en/v1.0.x/api/sep.Background.html#sep.Background
    # For more details
    bkg = sep.Background(image, bw=10, bh=10, fw=5, fh=5)
    print("# Mean Sky / RMS Sky = %10.5f / %10.5f" % (bkg.globalback, bkg.globalrms))


    # Object detection after subtracting the background
    # Since there is no noise on the image, we will just use the global RMS
    # from the background estimate as error
    # Here, we use very high threshold, less aggressive deblending method to
    # make the detection focuses on the big object
    objs, seg = sep.extract(image - bkg.back(),
                                    20.0,
                                    err=bkg.globalrms,
                                    minarea=1000,
                                    deblend_nthresh=24,
                                    deblend_cont=0.1,
                                    segmentation_map=True)

    # And this is how you turn the segmentation image into a mask
    # You can remove certain object from the segmentation map first
    seg_mask = (seg > 0)
    print("# Detect %d objects" % len(objs))

#     #################################################################################################
#     #stage 1: Free center and geometry
#     #################################################################################################
#     print('****************STAGE 1****************')
#     stage = 1
#     iso_1, iso_1_bin = galSBP.galSBP(maps_location+'illustris_1_xy'+suffix+'.fits',
#                                              galX=objs[0]['x'],
#                                              galY=objs[0]['y'],
#                                              galQ=(objs[0]['b'] / objs[0]['a']),
#                                              galPA=(objs[0]['theta'] * 180.0 / np.pi),
#                                              maxSma=150,
#                                              iniSma=10.0,
#                                              stage=1,
#                                              intMode='median',
#                                              ellipStep=0.1,
#                                              pix=pixel_scale,
#                                              zpPhoto=0.0,
#                                              isophote=x_isophote,
#                                              xttools=x_ttools,
#                                              savePng=False, verbose=True)
#     print('# Output file : %s' % iso_1_bin)
#     print('# Total stellar mass from the profile: logM = %7.4f' % (iso_1['mag_tot'][0] / -2.5))

#     if plots:
#         iso_ellip = fit_isophotes(image, iso_1, iso_1_bin, fits_prefix, suffix, stage, pixel_scale)
#         check_profile(iso_1)


    #################################################################################################
    #stage 2: fixed center, and let the geometry to be free
    #################################################################################################
    print('****************STAGE 2****************')
    stage=2

    iso_2, iso_2_bin = galSBP.galSBP(maps_location+'illustris_1_xy'+suffix+'.fits',
                                             galX=center_pixel,
                                             galY=center_pixel,
                                             galQ=0.6,
                                             galPA=-50.0,
                                             maxSma=220,
                                             iniSma=10.0,
                                             stage=2,
                                             intMode='median',
                                             ellipStep=0.05,
                                             pix=pixel_scale,
                                             zpPhoto=0.0,
                                             harmonics='1 2 3 4',
                                             isophote=x_isophote,
                                             xttools=x_ttools,
                                             recenter=False,
                                             savePng=False,
                                             verbose=True)
    print('# Output file : %s' % iso_2_bin)
    print('# Total stellar mass from the profile: logM = %7.4f' % (iso_2['mag_tot'][0] / -2.5))

    if plots:
        iso_ellip = fit_isophotes(image, iso_2, iso_2_bin, fits_prefix, suffix, stage, pixel_scale)
        check_profile(iso_2)


    ###########################################################################
    #stage 3: fix everything
    ###########################################################################
    print('****************STAGE 3****************')
    stage=3

    iso_3, iso_3_bin = galSBP.galSBP(maps_location+'illustris_1_xy'+suffix+'.fits',
                                             galX=center_pixel,
                                             galY=center_pixel,
                                             galQ=iso_2['avg_q'][0],
                                             galPA=iso_2['avg_pa'][0],
                                             maxSma=250,
                                             iniSma=50.0,
                                             stage=3,
                                             intMode='median',
                                             ellipStep=0.05,
                                             pix=pixel_scale,
                                             zpPhoto=0.0,
                                             harmonics='1 2 3 4',
                                             isophote=x_isophote,
                                             xttools=x_ttools,
                                             recenter=True,
                                             savePng=False,
                                             verbose=True)
    print('# Output file : %s' % iso_3_bin)
    print('# Total stellar mass from the profile: logM = %7.4f' % (iso_3['mag_tot'][0] / -2.5))

    if plots:
        iso_ellip = fit_isophotes(image, iso_3, iso_3_bin, fits_prefix, suffix, stage, pixel_scale)
        check_profile(iso_3)

    ###########################################################################

    return iso_3, iso_3_bin

def oneD_profile(iso):

    #fit Sersic profiles
    fitter = fitting.LevMarLSQFitter()
    r_mask = (iso['sma'] * 3.0 > 6.0)
    xx = iso['sma'][r_mask]
    yy = iso['intens'][r_mask] / 9.0

    # Let try to fit a Sersic function
    ser1_init = models.Sersic1D(amplitude=np.nanmedian(iso['intens'] / 9.0),
                                r_eff=6.5, n=6.0)
    ser1_fit = fitter(ser1_init, xx, yy)
    print("# 1-Sersic model:")
    print(ser1_fit.r_eff * 3.0, ser1_fit.n)

    # Let try to fit 2-Sersic functions then
    ser2_init = models.Sersic1D(amplitude=np.nanmedian(iso['intens'] / 8.0),
                                r_eff=2.5, n=3.0) + \
                models.Sersic1D(amplitude=np.nanmedian(iso['intens'] / 12.0),
                                r_eff=20.5, n=2.0)
    ser2_fit = fitter(ser2_init, xx, yy)
    print("\n# 2-Sersic model:")
    print(ser2_fit.r_eff_0 * 3.0, ser2_fit.n_0)
    print(ser2_fit.r_eff_1 * 3.0, ser2_fit.n_1)

    # Let try to fit 3-Sersic functions then
    ser3_init = models.Sersic1D(amplitude=np.nanmedian(iso['intens'] / 6.0),
                                r_eff=2.5, n=2.0) + \
                models.Sersic1D(amplitude=np.nanmedian(iso['intens'] / 10.0),
                                r_eff=12.5, n=1.5) + \
                models.Sersic1D(amplitude=np.nanmedian(iso['intens'] / 15.0),
                                r_eff=25.5, n=1.0)
    ser3_fit = fitter(ser3_init, xx, yy)
    print("\n# 3-Sersic model:")
    print(ser3_fit.r_eff_0 * 3.0, ser3_fit.n_0)
    print(ser3_fit.r_eff_1 * 3.0, ser3_fit.n_1)
    print(ser3_fit.r_eff_2 * 3.0, ser3_fit.n_2)


    #Plot the 1-D mass density profile, the curve-of-growth of mass, and the residuals

    # I label the 6.0 kpc (twice the pixel size, not resolved within) and 100 kpc radius
    fig = plt.figure(figsize=(8, 7))
    fig.subplots_adjust(hspace=0.0, wspace=0.0,
                        left=0.13, bottom=0.10,
                        top=0.97, right=0.97)

    ax1 = plt.subplot(3, 1, 1)
    ax1.grid(linewidth=2.0, linestyle='--', alpha=0.5)
    ax1.plot((iso['sma'] * 3.0) ** 0.25,
             np.log10(iso['intens'] / 9.0), linewidth=3.0,
             label=r'$\mathrm{Central}$')
    ax1.plot((iso['sma'] * 3.0) ** 0.25,
             np.log10(ser1_fit(iso['sma'])), linewidth=3.0, alpha=0.5,
             linestyle='--', label=r'$\mathrm{1\ Sersic}$')
    ax1.plot((iso['sma'] * 3.0) ** 0.25,
             np.log10(ser2_fit(iso['sma'])), linewidth=3.0, alpha=0.5,
             linestyle='--', label=r'$\mathrm{2\ Sersic}$')
    ax1.plot((iso['sma'] * 3.0) ** 0.25,
             np.log10(ser3_fit(iso['sma'])), linewidth=3.0, alpha=0.5,
             linestyle='--', label=r'$\mathrm{3\ Sersic}$')
    ax1.set_ylabel(r'$\mu_{\star}\ (M_{\star}/\mathrm{kpc}^2)$', fontsize=20)
    ax1.axvline(100.0 ** 0.25, linestyle='--', linewidth=3.0, alpha=0.6)
    ax1.axvline(6.0 ** 0.25, linestyle='-', linewidth=3.0, alpha=0.6, c='r')
    ax1.set_xlim(1.0, 4.5)
    ax1.set_ylim(4.0, 11.0)
    ax1.legend(fontsize=14)

    ax2 = plt.subplot(3, 1, 2)
    ax2.grid(linewidth=2.0, linestyle='--', alpha=0.5)
    ax2.plot([0.0, 4.5], [0.0, 0.0])
    ax2.plot((iso['sma'] * 3.0) ** 0.25,
             (np.log10(iso['intens'] / 9.0) -
              np.log10(ser1_fit(iso['sma']))),
             linewidth=3.0, alpha=0.5,
             linestyle='--', label=r'$\mathrm{1\ Sersic}$')
    ax2.plot((iso['sma'] * 3.0) ** 0.25,
             (np.log10(iso['intens'] / 9.0) -
              np.log10(ser2_fit(iso['sma']))),
             linewidth=3.0, alpha=0.5,
             linestyle='--', label=r'$\mathrm{2\ Sersic}$')
    ax2.plot((iso['sma'] * 3.0) ** 0.25,
             (np.log10(iso['intens'] / 9.0) -
              np.log10(ser3_fit(iso['sma']))),
             linewidth=3.0, alpha=0.5,
             linestyle='--', label=r'$\mathrm{3\ Sersic}$')
    ax2.set_ylabel(r'$\mathrm{Residuals}$', fontsize=20)
    ax2.set_xlabel(r'$\mathrm{SMA/kpc}^{1/4}$', fontsize=25)
    ax2.axvline(100.0 ** 0.25, linestyle='--', linewidth=3.0, alpha=0.6)
    ax2.axvline(6.0 ** 0.25, linestyle='-', linewidth=3.0, alpha=0.6, c='r')
    ax2.set_xlim(1.0, 4.5)
    #ax2.set_ylim(-2.0, 2.0)


    ax3 = plt.subplot(3, 1, 3)
    ax3.grid(linewidth=2.0, linestyle='--', alpha=0.5)
    ax3.plot((iso['sma'] * 3.0) ** 0.25,
             np.log10(iso['growth_ori']), linewidth=3.0)
    ax3.set_ylabel(r'$\mathrm{Curve\ of\ growth}\ (M_{\star})$', fontsize=20)
    ax3.set_xlabel(r'$\mathrm{SMA/kpc}^{1/4}$', fontsize=25)
    ax3.axvline(100.0 ** 0.25, linestyle='--', linewidth=3.0, alpha=0.6)
    ax3.axvline(6.0 ** 0.25, linestyle='-', linewidth=3.0, alpha=0.6, c='r')
    ax3.set_xlim(1.0, 4.5)
    
