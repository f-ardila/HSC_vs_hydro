from __future__ import print_function, division, absolute_import

from functions import *

import time
time0=time.time()

resolution=sys.argv[1]

#data
Illustris_file_orig = '/Users/fardila/Documents/GitHub/HSC_vs_hydro/Data/Illustris/galaxies_orig_11.2.hdf5'
TNG_file_orig = '/Users/fardila/Documents/GitHub/HSC_vs_hydro/Data/TNG/galaxies_tng75_11.2.hdf5'

Illustris_file_quick = '/Users/fardila/Documents/GitHub/HSC_vs_hydro/Data/Illustris/galaxies_stellarmaps_orig_11.2.hdf5'
TNG_file_quick = '/Users/fardila/Documents/GitHub/HSC_vs_hydro/Data/TNG/galaxies_stellarmaps_tng75_11.2.hdf5'

TNG_file_highres = '/Users/fardila/Documents/GitHub/HSC_vs_hydro/Data/TNG/galaxies_stellarmaps_tng75_11.2_highres.hdf5'
Illustris_file_highres = '/Users/fardila/Documents/GitHub/HSC_vs_hydro/Data/Illustris/galaxies_stellarmaps_orig_11.2_highres.hdf5'


if resolution=='quick':
    Illustris_file = Illustris_file_quick
    TNG_file = TNG_file_quick
elif resolution == 'highres':
    Illustris_file = Illustris_file_highres
    TNG_file = TNG_file_highres

###############################################################################
#run on Illustris
isos_illustris=[]
masses_illustris=[]

for i in range(339):

    print('^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^')
    print('^^^^^^^^GALAXY '+str(i)+'^^^^^^^^^^^^^^')
    print('^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^')

    try:
        iso, masses = get_masses_iso(Illustris_file,'Illustris', resolution, gal_n=i)

    except ValueError:
        iso=-99.99
        masses=-99.99

    isos_illustris.append(iso)
    masses_illustris.append(masses)

#save as pickles
pkl_isos = open('/Users/fardila/Documents/GitHub/HSC_vs_hydro/Data/Illustris/Illustris_isos_{0}.pkl'.format(resolution),'wb')
pickle.dump(isos_illustris,pkl_isos)
pkl_isos.close()

pkl_masses = open('/Users/fardila/Documents/GitHub/HSC_vs_hydro/Data/Illustris/Illustris_masses_{0}.pkl'.format(resolution),'wb')
pickle.dump(masses_illustris,pkl_masses)
pkl_masses.close()

###############################################################################
# run on TNG
isos_tng=[]
masses_tng=[]

for i in range(235):

    print('^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^')
    print('^^^^^^^^GALAXY '+str(i)+'^^^^^^^^^^^^^^')
    print('^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^')

    try:
        iso, masses = get_masses_iso(TNG_file,'TNG', resolution, gal_n=i, intMode='mean')

    except:
        iso=-99.99
        masses=-99.99

    isos_tng.append(iso)
    masses_tng.append(masses)

#save as pickles
pkl_isos = open('/Users/fardila/Documents/GitHub/HSC_vs_hydro/Data/TNG/TNG_isos_{0}.pkl'.format(resolution),'wb')
pickle.dump(isos_tng,pkl_isos)
pkl_isos.close()

pkl_masses = open('/Users/fardila/Documents/GitHub/HSC_vs_hydro/Data/TNG/TNG_masses_{0}.pkl'.format(resolution),'wb')
pickle.dump(masses_tng,pkl_masses)
pkl_masses.close()

print(time.time()-time0, ' seconds')
