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

    return log_m_vir

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
tng_m_vir = get_Mvir(tng_map_highres_file)
cosmology.setCosmology('illustris')
illustris_m_vir = get_Mvir(illustris_map_highres_file)

#add mvir to array
tng_masses['m_vir'] = tng_m_vir
illustris_masses['m_vir'] = illustris_m_vir

#save
save_pkl(tng_masses_file+'new', tng_masses)
save_pkl(illustris_masses_file+'new', illustris_masses)
