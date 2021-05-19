from functions import *
import argparse

import time
time0=time.time()

#arguments
parser = argparse.ArgumentParser()
# parser.add_argument('-masses', action='store_true', help="calculate only masses")
# parser.add_argument('-isos', action='store_true', help="calculate only isos")
parser.add_argument("resolution", help="resolution of maps at which to make measurements",
                    choices=['quick', 'highres'])
args = parser.parse_args()

#resolution argument
resolution=args.resolution
if resolution=='quick':
    Illustris_file = Illustris_file_quick
    TNG_file = TNG_file_quick
elif resolution == 'highres':
    Illustris_file = Illustris_file_highres
    TNG_file = TNG_file_highres

###############################################################################
#run on Illustris
outfile_loc = '/Users/fardila/Documents/GitHub/HSC_vs_hydro/Data/Illustris/'
isos_illustris = load_pkl(outfile_loc+'Illustris_isos_{0}.pkl'.format(resolution))
masses_illustris = [get_masses(isos_illustris[i], Illustris_file,'Illustris', resolution, rs=[300,500,800], gal_n=i) for i in range(339)]

#save as pickles
save_pkl(outfile_loc+'Illustris_masses_{0}.pkl'.format(resolution), vstack(masses_illustris))

###############################################################################
# run on TNG
outfile_loc = '/Users/fardila/Documents/GitHub/HSC_vs_hydro/Data/TNG/'
isos_tng = load_pkl(outfile_loc+'TNG_isos_{0}.pkl'.format(resolution))
masses_tng = [get_masses(isos_tng[i], TNG_file,'TNG', resolution, rs=[300,500,800], gal_n=i) for i in range(235)]

#save as pickles
save_pkl(outfile_loc+'TNG_masses_{0}.pkl'.format(resolution), vstack(masses_tng))

print(time.time()-time0, ' seconds')
