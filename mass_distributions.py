from __future__ import print_function, division, absolute_import

import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
plt.rc('text', usetex=True)




###############################################################################
#data
Illustris_file = '/Users/fardila/Documents/GitHub/HSC_vs_hydro/Data/Illustris/galaxies_orig_11.2.hdf5'
TNG_file = '/Users/fardila/Documents/GitHub/HSC_vs_hydro/Data/TNG/galaxies_tng75_11.2.hdf5'


###############################################################################
# Illustris
f = h5py.File(Illustris_file, 'r')
cat_sh_mstars = np.array(f['cat_sh_mstar'])
f.close()
m_stars_illustris=np.log10(cat_sh_mstars)
mean_illustris=round(np.mean(m_stars_illustris),2)
median_illustris=round(np.median(m_stars_illustris),2)
std_illustris=round(np.std(m_stars_illustris),2)

#tng
f = h5py.File(TNG_file, 'r')
cat_sh_mstars = np.array(f['cat_sh_mstar'])
f.close()
m_stars_tng=np.log10(cat_sh_mstars)
mean_tng=round(np.mean(m_stars_tng),2)
median_tng=round(np.median(m_stars_tng),2)
std_tng=round(np.std(m_stars_tng),2)

###############################################################################
#plot

plt.hist([m_stars_illustris, m_stars_tng], 500, label=['Illustris: Mean: '+ str(mean_illustris)+', Median: '+ str(median_illustris) + ', StDev: ' + str(std_illustris), 'TNG: Mean: '+ str(mean_tng)+', Median: '+ str(median_tng) + ', StDev: ' + str(std_tng)],
        cumulative=True, histtype='step', normed=True)

plt.legend(loc='lower center')
plt.xlabel('$\log{M_{\star}}$')
plt.show()
