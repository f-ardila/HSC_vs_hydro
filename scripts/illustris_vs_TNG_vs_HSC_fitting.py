from __future__ import print_function, division, absolute_import

from fit_profile import *

import numpy as np
import os
import pickle



import time
time0=time.time()

###############################################################################
#data
Illustris_file = '/Users/fardila/Documents/GitHub/HSC_vs_hydro/Data/Illustris/galaxies_orig_11.2.hdf5'
TNG_file = '/Users/fardila/Documents/GitHub/HSC_vs_hydro/Data/TNG/galaxies_tng75_11.2.hdf5'

###############################################################################
#run on Illustris
isos_illustris=[]
for i in range(339):
    print('^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^')
    print('^^^^^^^^GALAXY '+str(i)+'^^^^^^^^^^^^^^')
    print('^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^')

    try:
        iso_cen, iso_cen_bin, log_mstar, log_mcen = fit_profile(Illustris_file,
                                            'Illustris',gal_n=i, plots=False)

    except:
        continue
        

    isos_illustris.append([iso_cen,log_mstar, log_mcen])

#save as pickle
pkl = open('/Users/fardila/Documents/GitHub/HSC_vs_hydro/Data/Illustris/Illustris.pkl','wb') 
pickle.dump(isos_illustris,pkl)   
pkl.close()

###############################################################################
#run on TNG
isos_tng=[]
for i in range(235):

    print('^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^')
    print('^^^^^^^^GALAXY '+str(i)+'^^^^^^^^^^^^^^')
    print('^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^')

    try:
        iso_cen, iso_cen_bin, log_mstar, log_mcen= fit_profile(TNG_file,'TNG',\
                                                            gal_n=i, plots=False)

    except:
        continue



    isos_tng.append([iso_cen,log_mstar, log_mcen])

#save as pickle
pkl = open('/Users/fardila/Documents/GitHub/HSC_vs_hydro/Data/TNG/TNG.pkl','wb') 
pickle.dump(isos_tng,pkl)   
pkl.close()

print(time.time()-time0, ' seconds')


