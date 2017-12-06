import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
plt.rc('text', usetex=True)

from functions import load_pkl
import glob

################################################################################
#HSC
hsc_data='/Users/fardila/Documents/Github/HSC_vs_hydro/Data/HSC/'
hscAvgProf0 = load_pkl(hsc_data+"hscAvgProf0.pkl")
hscAvgProf1 = load_pkl(hsc_data+"hscAvgProf1.pkl")
hscAvgProf2 = load_pkl(hsc_data+"hscAvgProf2.pkl")

rm0_sl, rm0_ml, rm0_aml = hscAvgProf0['all'], hscAvgProf0['med'], hscAvgProf0['avg']
rm1_sl, rm1_ml, rm1_aml = hscAvgProf1['all'], hscAvgProf1['med'], hscAvgProf1['avg']
rm2_sl, rm2_ml, rm2_aml = hscAvgProf2['all'], hscAvgProf2['med'], hscAvgProf2['avg']

# Universal RSMA array
RSMA_COMMON = np.arange(0.4, 4.2, 0.01)

# These are the median stellar mass density profiles for HSC galaxies at 0.3 < z < 0.5
# in three mass bins
# rm0 : 11.4 < logM_100kpc < 11.6
# rm1 : 11.6 < logM_100kpc < 11.8
# rm2 : 11.8 < logM_100kpc < 12.0
# They are on a common radius array, and we use (r ** 0.25) as radius

#####################################################
#import files
data_loc='/Users/fardila/Documents/GitHub/HSC_vs_hydro/Figures/fits_files/'
illustris_pickels=glob.glob('/Users/fardila/Documents/GitHub/HSC_vs_hydro/Figures/fits_files/Illustris*3.pkl')
tng_pickels=glob.glob('/Users/fardila/Documents/GitHub/HSC_vs_hydro/Figures/fits_files/TNG*3.pkl')




###############################################################################
#plot
fig, ax1 = plt.subplots(figsize=(15, 10))
ax2 = ax1.twiny()
new_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b',
                '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

for illustris_pkl, tng_pkl in zip(illustris_pickels,tng_pickels):

    iso_illustris = load_pkl(illustris_pkl)
    iso_tng = load_pkl(tng_pkl)

    #m_star= isos[i][1]

    ax1.plot((iso_illustris['sma']) ** 0.25,
         (iso_illustris['sbp_cor'] / -2.5) + np.log10(0.7 ** 2.0), linewidth=2.0, c='b', alpha=0.5)
             #label='Illustris Galaxy '+str(i)+': $M_{\star} = $'+str(round(m_star,1)))
    ax1.plot((iso_tng['sma']) ** 0.25,
         (iso_tng['sbp_cor'] / -2.5) + np.log10(0.7 ** 2.0), linewidth=2.0, c='r', alpha=0.5)
             #label='Illustris Galaxy '+str(i)+': $M_{\star} = $'+str(round(m_star,1)))


## Median profiles from HSC
ax1.plot(RSMA_COMMON, rm0_aml[2], linestyle='--', linewidth=4.0,
         alpha=0.5, zorder=8, label='HSC: $11.4 < \log{M_{100kpc}} < 11.6$')
ax1.plot(RSMA_COMMON, rm1_aml[2], linestyle='--', linewidth=4.0,
         alpha=0.5, zorder=8, label='HSC: $11.6 < \log{M_{100kpc}} < 11.8$')
ax1.plot(RSMA_COMMON, rm2_aml[2], linestyle='--', linewidth=4.0,
         alpha=0.5, zorder=8, label='HSC: $11.8 < \log{M_{100kpc}} < 12.0$')


ax1.set_xlim(0.9, 4.5)
ax1.set_ylim(4, 10)

#add twin x axis in kpc
x1, x2 = ax1.get_xlim()
ax2.set_xlim(x1, x2)
ax2.figure.canvas.draw()
ax2.xaxis.set_ticks([2**0.25, 5**0.25, 10**0.25, 50**0.25, 100**0.25, 200**0.25, 300**0.25])
ax2.xaxis.set_ticklabels([2, 5, 10, 50, 100, 200, 300])

ax1.tick_params(axis='both', which='major', labelsize=15)
ax2.tick_params(axis='both', which='major', labelsize=15)
ax1.set_xlabel(r'$\mathrm{SMA/kpc}^{1/4}$', fontsize=25)
ax2.set_xlabel(r'$\mathrm{kpc}$', fontsize=25)
ax1.set_ylabel(r'$\mu_{\star}\ (M_{\star}/\mathrm{kpc}^2)$', fontsize=20)
ax1.axvline(100.0 ** 0.25, linestyle='--', linewidth=3.0, alpha=0.6)
ax1.axvline(6.0 ** 0.25, linestyle='-', linewidth=3.0, alpha=0.6, c='r')

illustris_line = mlines.Line2D([], [], color='b', markersize=15, alpha=0.5, label='Illustris' )
tng_line = mlines.Line2D([], [], color='r', markersize=15, alpha=0.5, label='TNG' )
hsc_line1 = mlines.Line2D([], [], markersize=15, alpha=0.5, linestyle='--', linewidth=4.0,
                            c=new_colors[0], label='HSC: $11.4 < \log{M_{100kpc}} < 11.6$' )
hsc_line2 = mlines.Line2D([], [], markersize=15, alpha=0.5, linestyle='--', linewidth=4.0,
                            c=new_colors[1], label='HSC: $11.6 < \log{M_{100kpc}} < 11.8$' )
hsc_line3 = mlines.Line2D([], [], markersize=15, alpha=0.5, linestyle='--', linewidth=4.0,
                            c=new_colors[2], label='HSC: $11.8 < \log{M_{100kpc}} < 12.0$')


ax1.legend(handles=[illustris_line, tng_line, hsc_line1, hsc_line2, hsc_line3], fontsize=25)

plt.show()

print(time.time()-time0, ' seconds')
