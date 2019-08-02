from functions import *
import argparse
import h5py
import numpy as np

import time
time0=time.time()

#arguments
# parser = argparse.ArgumentParser()
# # parser.add_argument('-masses', action='store_true', help="calculate only masses")
# # parser.add_argument('-isos', action='store_true', help="calculate only isos")
# parser.add_argument("resolution", help="resolution of maps at which to make measurements",
#                     choices=['quick', 'highres'])
# args = parser.parse_args()
#
# #resolution argument
# resolution=args.resolution


components = 'cen'

xray_file = '/Users/fardila/Documents/GitHub/HSC_vs_hydro/CC/xray_clusterfinder_tng300_100.hdf5'
maps_file_highres = '/Users/fardila/Documents/GitHub/HSC_vs_hydro/CC/galaxies_stellarmaps_hr_11.2_tng300_072.hdf5'
maps_file_quick= '/Users/fardila/Documents/GitHub/HSC_vs_hydro/CC/galaxies_stellarmaps_lr_11.2_tng300_072.hdf5'


f = h5py.File(xray_file, 'r')
xray_ids = np.array(f['catgrp_id'])
f.close()


f = h5py.File(maps_file_highres, 'r')
all_map_indexes = range(len(np.array(f['catgrp_is_primary'])))
mask = np.in1d(np.array(f['catgrp_id']), xray_ids) & np.array(f['catgrp_is_primary'])
maps_indexes_to_use = np.array(all_map_indexes)[mask]

print(maps_indexes_to_use)
f.close()



###############################################################################
resolution='quick'
isos_tng = [get_iso(maps_file_quick, 'TNG300', resolution, intMode='mean',
                    components=components, gal_n=i) for i in maps_indexes_to_use]
masses_tng = [get_masses(isos_tng[i], TNG_file, 'TNG', resolution,
                         rs=[300, 500, 800], gal_n=i) for i in maps_indexes_to_use]

# save as pickles
outfile_loc = '/Users/fardila/Documents/GitHub/HSC_vs_hydro/CC/'
save_pkl(outfile_loc+'CC_TNG300_isos_{0}_{1}.pkl'.format(components, resolution), isos_tng)
save_pkl(outfile_loc+'CC_TNG300_masses_{0}_{1}.pkl'.format(components,
                                                     resolution), vstack(masses_tng))
                                                     
###############################################################################
resolutions = 'highres'
isos_tng = [get_iso(maps_file_highres, 'TNG300', resolution, intMode='mean',
                    components=components, gal_n=i) for i in maps_indexes_to_use]
masses_tng = [get_masses(isos_tng[i], TNG_file, 'TNG', resolution,
                         rs=[300, 500, 800], gal_n=i) for i in maps_indexes_to_use]

# save as pickles
outfile_loc = '/Users/fardila/Documents/GitHub/HSC_vs_hydro/CC/'
save_pkl(outfile_loc+'CC_TNG300_isos_{0}_{1}.pkl'.format(components, resolution), isos_tng)
save_pkl(outfile_loc+'CC_TNG300_masses_{0}_{1}.pkl'.format(components,
                                                     resolution), vstack(masses_tng))
