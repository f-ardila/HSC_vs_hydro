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
