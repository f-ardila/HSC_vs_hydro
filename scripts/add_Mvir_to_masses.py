import numpy as np
import h5py
import pickle

from colossus.halo import mass_adv
from colossus.cosmology import cosmology


def open_pkl(file_name):
    pkl = open(file_name,'rb')
    array = pickle.load(pkl)
    pkl.close()
    return array

def save_pkl(file, array):
    pkl = open(file,'wb')
    pickle.dump(array,pkl)
    pkl.close()

def get_Mvir(sim_map_file):
    f = h5py.File(sim_map_file, 'r')
    m_200c = np.log10(np.array(f['catgrp_Group_M_Crit200']))
    z = f['config'].attrs['snap_z']
    #convert to Mvir
    m_vir, r_vir, c_vir = mass_adv.changeMassDefinitionCModel(10**m_200c,
                                                                z, '200c',
                                                                'vir',
                                                                profile='nfw',
                                                                c_model='diemer19')
    log_m_vir = np.log10(m_vir)
    f.close()

    return log_m_vir, m_200c

#open pickels
tng_masses_file = '/Users/fardila/Documents/GitHub/HSC_vs_hydro/Data/TNG/TNG_masses_highres.pkl'
tng_map_highres_file = '/Users/fardila/Documents/GitHub/HSC_vs_hydro/Data/TNG/galaxies_stellarmaps_tng75_11.2_highres.hdf5'

illustris_masses_file = '/Users/fardila/Documents/GitHub/HSC_vs_hydro/Data/Illustris/Illustris_masses_highres.pkl'
illustris_map_highres_file = '/Users/fardila/Documents/GitHub/HSC_vs_hydro/Data/Illustris/galaxies_stellarmaps_orig_11.2_highres.hdf5'

#array of masses
tng_masses = open_pkl(tng_masses_file)
illustris_masses = open_pkl(illustris_masses_file)

#get halo mass
cosmology.setCosmology('planck15')
tng_m_vir, tng_m200c = get_Mvir(tng_map_highres_file)
cosmology.setCosmology('illustris')
illustris_m_vir, illustris_m200c = get_Mvir(illustris_map_highres_file)

#add mvir to array
tng_masses['m_200c'] = tng_m200c
illustris_masses['m_200c'] = illustris_m200c

#delete previous virial masses wwhich were converted from m220c and missing an h
try:
    tng_masses.remove_columns(['m_vir', 'm_vir(converted)'])
    illustris_masses.remove_columns(['m_vir', 'm_vir(converted)'])
except:
    pass

################################################################################
#get true virial masses from group catalog

h_assumed = 0.7
h_TNG = 0.6774
h_illustris = 0.704
#TNG
files_directory = '/Users/fardila/Documents/GitHub/HSC_vs_hydro/Data/TNG/group_catalogs/'
MTopHat200_file = 'fof_subhalo_tab_072.Group.Group_M_TopHat200.hdf5'
f = h5py.File(files_directory+MTopHat200_file, 'r')
tng_MTopHat200 = np.array(f['Group'][u'Group_M_TopHat200'])* 1e10 / h_assumed
f.close()

#get ids
f = h5py.File(tng_map_highres_file, 'r')
tng_cat_group_id = np.array(f['catgrp_id'])
f.close()

#Illustris
files_directory = '/Users/fardila/Documents/GitHub/HSC_vs_hydro/Data/Illustris/group_catalogs/'
MTopHat200_file = 'groups_108.Group.Group_M_TopHat200.hdf5'
f = h5py.File(files_directory+MTopHat200_file, 'r')
Illustris_MTopHat200 = np.array(f['Group'][u'Group_M_TopHat200'])* 1e10 / h_assumed
f.close()

#get ids
f = h5py.File(illustris_map_highres_file, 'r')
illustris_cat_group_id = np.array(f['catgrp_id'])
f.close()

#add Mtophat to array
tng_masses['m_tophat200'] = [np.log10(tng_MTopHat200[id]) for id in tng_cat_group_id]
illustris_masses['m_tophat200'] = [np.log10(Illustris_MTopHat200[id]) for id in illustris_cat_group_id]

################################################################################
#save
save_pkl(tng_masses_file, tng_masses)
save_pkl(illustris_masses_file, illustris_masses)
