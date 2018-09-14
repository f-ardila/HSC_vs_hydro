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
    print('^^^^^ILLUSTRIS GALAXY '+str(i)+'^^^^^^')
    print('^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^')

    try:
        iso = get_iso(Illustris_file,'Illustris', resolution, intMode='mean', components='cen', gal_n=i)
        masses = get_masses(iso, Illustris_file,'Illustris', resolution, rs=[300,500,800], gal_n=i)
    except ValueError:
        iso=-99.99
        masses=-99.99

    isos_illustris.append(iso)
    masses_illustris.append(masses)

#save as pickles
outfile_loc = '/Users/fardila/Documents/GitHub/HSC_vs_hydro/Data/Illustris/'

pkl_masses = open(outfile_loc+'Illustris_masses_{0}.pkl'.format(resolution),'wb')
pickle.dump(vstack(masses_illustris),pkl_masses)
pkl_masses.close()

pkl_isos = open(outfile_loc+'Illustris_isos_{0}.pkl'.format(resolution),'wb')
pickle.dump(isos_illustris,pkl_isos)
pkl_isos.close()

###############################################################################
# run on TNG
isos_tng=[]
masses_tng=[]

for i in range(235):

    print('^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^')
    print('^^^^^^^^^TNG GALAXY '+str(i)+'^^^^^^^^^')
    print('^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^')

    try:

        iso = get_iso(TNG_file,'TNG', resolution, intMode='mean', components='cen', gal_n=i)
        masses = get_masses(iso, TNG_file,'TNG', resolution, rs=[300,500,800], gal_n=i)

    except:
        iso=-99.99
        masses=-99.99

    isos_tng.append(iso)
    masses_tng.append(masses)

#save as pickles
outfile_loc = '/Users/fardila/Documents/GitHub/HSC_vs_hydro/Data/TNG/'

pkl_masses = open(outfile_loc+'TNG_masses_{0}.pkl'.format(resolution),'wb')
pickle.dump(vstack(masses_tng),pkl_masses)
pkl_masses.close()

pkl_isos = open(outfile_loc+'TNG_isos_{0}.pkl'.format(resolution),'wb')
pickle.dump(isos_tng,pkl_isos)
pkl_isos.close()

print(time.time()-time0, ' seconds')
