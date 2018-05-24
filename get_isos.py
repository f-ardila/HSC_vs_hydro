from __future__ import print_function, division, absolute_import

import pickle
from functions import *

import time
time0=time.time()

#data
Illustris_file_quick = '/Users/fardila/Documents/GitHub/HSC_vs_hydro/Data/Illustris/galaxies_stellarmaps_orig_11.2.hdf5'
TNG_file_quick = '/Users/fardila/Documents/GitHub/HSC_vs_hydro/Data/TNG/galaxies_stellarmaps_tng75_11.2.hdf5'

#which components?
comp = 'cen+icl'
###############################################################################
#run on Illustris
isos_illustris=[]

for i in range(339):


    print('^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^')
    print('^^^^^^^^GALAXY '+str(i)+'^^^^^^^^^^^^^^')
    print('^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^')

    try:
        iso= get_iso(Illustris_file_quick,'Illustris', components=comp, gal_n=i)

    except:
        iso=-99.99

    isos_illustris.append(iso)

#save as pickles
pkl_isos = open('/Users/fardila/Documents/GitHub/HSC_vs_hydro/Data/Illustris/Illustris_isos_quick_'+ comp +'.pkl','wb')
pickle.dump(isos_illustris,pkl_isos)
pkl_isos.close()

###############################################################################
#run on TNG
isos_tng=[]

for i in range(235):

    print('^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^')
    print('^^^^^^^^GALAXY '+str(i)+'^^^^^^^^^^^^^^')
    print('^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^')

    try:
        iso= get_iso(TNG_file_quick,'TNG', components=comp, gal_n=i)

    except:
        iso=-99.99

    isos_tng.append(iso)

#save as pickles
pkl_isos = open('/Users/fardila/Documents/GitHub/HSC_vs_hydro/Data/TNG/TNG_isos_quick_'+ comp +'.pkl','wb')
pickle.dump(isos_tng,pkl_isos)
pkl_isos.close()

print(time.time()-time0, ' seconds')
