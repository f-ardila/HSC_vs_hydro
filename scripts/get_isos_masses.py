from functions import *
import argparse

import time
time0 = time.time()

# arguments
parser = argparse.ArgumentParser()
# parser.add_argument('-masses', action='store_true', help="calculate only masses")
# parser.add_argument('-isos', action='store_true', help="calculate only isos")
parser.add_argument("resolution", help="resolution of maps at which to make measurements",
                    choices=['quick', 'highres'])

parser.add_argument("-sats", action='store_true', help="include satellites in addition to centrals")
parser.add_argument("-icl", action='store_true', help="include ICL (fuzz) in addition to centrals")

parser.add_argument("-TNG", action='store_true', help="include TNG only")
parser.add_argument("-TNG300", action='store_true', help="include TNG300 only")
parser.add_argument("-Illustris", action='store_true', help="include Illustris only")
args = parser.parse_args()

# resolution argument
resolution = args.resolution
if resolution == 'quick':
    Illustris_file = Illustris_file_quick
    TNG_file = TNG_file_quick
elif resolution == 'highres':
    Illustris_file = Illustris_file_highres
    TNG_file = TNG_file_highres

# components arguments
components = 'cen'
if args.sats:
    if args.icl:
        components = 'all'
    else:
        components += '+sats'
elif args.icl:
    components += '+icl'

if not args.TNG and not args.Illustris and not args.TNG300:
    args.TNG = True
    args.Illustris = True

###############################################################################
# run on Illustris
if args.Illustris:
    isos_illustris = [get_iso(Illustris_file, 'Illustris', resolution,
                              intMode='mean', components=components, gal_n=i) for i in range(339)]
    masses_illustris = [get_masses(isos_illustris[i], Illustris_file, 'Illustris', resolution, rs=[
                                   300, 500, 800], gal_n=i) for i in range(339)]

    # save as pickles
    outfile_loc = '/Users/fardila/Documents/GitHub/HSC_vs_hydro/Data/Illustris/'
    save_pkl(
        outfile_loc+'Illustris_isos_{0}_{1}.pkl'.format(components, resolution), isos_illustris)
    save_pkl(
        outfile_loc+'Illustris_masses_{0}_{1}.pkl'.format(components, resolution), vstack(masses_illustris))

###############################################################################
# run on TNG
if args.TNG:
    isos_tng = [get_iso(TNG_file, 'TNG', resolution, intMode='mean',
                        components=components, gal_n=i) for i in range(235)]
    masses_tng = [get_masses(isos_tng[i], TNG_file, 'TNG', resolution,
                             rs=[300, 500, 800], gal_n=i) for i in range(235)]

    # save as pickles
    outfile_loc = '/Users/fardila/Documents/GitHub/HSC_vs_hydro/Data/TNG/'
    save_pkl(outfile_loc+'TNG_isos_{0}_{1}.pkl'.format(components, resolution), isos_tng)
    save_pkl(outfile_loc+'TNG_masses_{0}_{1}.pkl'.format(components,
                                                         resolution), vstack(masses_tng))

###############################################################################
# run on TNG
if args.TNG300:
    components = 'cen'
    resolution = 'quick'
    sim_name = 'TNG300'
    TNG_file = TNG300_file_quick
    f = h5py.File(TNG_file, 'r')
    n = len(f['catgrp_Group_M_Crit200'])
    centrals = np.array(f['catgrp_is_primary'])
    f.close()
    # n = 1
    print(len(np.arange(n)[centrals]))
    # isos_tng = [get_iso(TNG_file, sim_name, resolution, intMode='mean',
    #                     components=components, gal_n=i) for i in [390,391]]
    isos_tng = [get_iso(TNG_file, sim_name, resolution, intMode='mean',
                        components=components, gal_n=i) for i in np.arange(n)[centrals]]
    # masses_tng = [get_masses(isos_tng[i], TNG_file, sim_name, resolution,
    #                          rs=[300, 500, 800], gal_n=i) for i in range(len())]

    # save as pickles
    outfile_loc = '/Users/fardila/Documents/GitHub/HSC_vs_hydro/Data/TNG/'
    save_pkl(outfile_loc+'{0}_isos_{1}_{2}.pkl'.format(sim_name, components, resolution), isos_tng)
    # save_pkl(outfile_loc+'TNG_masses_{0}_{1}.pkl'.format(components,
    #                                                      resolution), vstack(masses_tng))

print(time.time()-time0, ' seconds')
