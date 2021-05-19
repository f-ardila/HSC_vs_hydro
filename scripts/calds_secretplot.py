from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table, Column
import numpy as np
from halotools.mock_observables import return_xyz_formatted_array
from halotools.mock_observables import delta_sigma
from halotools.mock_observables import tpcf
from astropy.cosmology import WMAP7 as cosmo

# Massive black 2
# Lbox = 100 h^-1 Mpc
# WMAP 7 cosmology

# This is using the new DeltaSigma module
# March 2017
# New DS module computes DS using projected mass distirbution
# To see how to use this, look at the html page under this directory

def main():

        #--------------------------
        #--------------------------
        
        Lbox = 100.0   # 100/h
        period = np.array([Lbox,Lbox,Lbox])

        #numpy.logspace(start, stop, num=50, endpoint=True, base=10.0, dtype=None)[source]
        # This is in Mpc/h comoving
        rp_bins = np.logspace(-1.5,1.5,20)
        print rp_bins
        
        # ---- DARK MATTER
        dm_pos = np.fromfile('./position_data/dm/position',dtype=(np.float64,3))

        # These are in units of kpc, put into Mpc
        x=dm_pos[:,0]/1000.0
        y=dm_pos[:,1]/1000.0
        z=dm_pos[:,2]/1000.0
                
        # Format the array for halotools
        dmpos = return_xyz_formatted_array(x, y, z)
        
        ntotparticles = 5754585088.0
        dm_downsampling_factor = ntotparticles/10000000.0

        # Constant particle mass
        dm_particle_mass = 0.00110449e+10

        # ---- GAS
        gas_pos = np.fromfile('./position_data/gas/position',dtype=(np.float64,3))
        
        # particle mass array
        gas_pm_array = np.fromfile('./position_data/gas/mass',dtype=np.float64)

        # Proper units !!
        gas_pm_array = gas_pm_array*1e10
                
        # These are in units of kpc, put into Mpc
        x=gas_pos[:,0]/1000.0
        y=gas_pos[:,1]/1000.0
        z=gas_pos[:,2]/1000.0
                
        # Format the array for halotools
        gaspos = return_xyz_formatted_array(x, y, z)
        
        ntotparticles = 5620192089.0
        gas_downsampling_factor = ntotparticles/10000000.0

        # ---- STARS
        star_pos = np.fromfile('./position_data/star/position',dtype=(np.float64,3))
        
        # particle mass array
        star_pm_array = np.fromfile('./position_data/star/mass',dtype=np.float64)

        # Proper units !!
        star_pm_array = star_pm_array*1e10
                
        # These are in units of kpc, put into Mpc
        x=star_pos[:,0]/1000.0
        y=star_pos[:,1]/1000.0
        z=star_pos[:,2]/1000.0
                
        # Format the array for halotools
        starpos = return_xyz_formatted_array(x, y, z)
        
        ntotparticles = 569213225.0
        star_downsampling_factor = ntotparticles/10000000.0
        
        #--------------------------
        # Now do the galaxy samples
        #--------------------------
        
        # Read in the Galaxy File
        gal_list = np.loadtxt('/Users/alexie/Work/HSC/CompareHydroProfiles/DataFromMBII/mblack2_mhalo_matched_highM100.asc')
        x=gal_list[:,4]
        y=gal_list[:,5]
        z=gal_list[:,6]
        # These are in units of kpc, put into Mpc
        x=x/1000.0
        y=y/1000.0
        z=z/1000.0
        # Format the array for halotools
        galpos = return_xyz_formatted_array(x, y, z)

        r, ds = delta_sigma(galpos, dmpos, dm_particle_mass, dm_downsampling_factor, rp_bins, Lbox)
        ds = ds/1e12
        table = {'rmpc': r,'deltasigma': ds}   
        ascii.write(table, 'ds_dm_mblack2_mhalo_matched_highM100.dat', formats={'rmpc': '%.5f', 'deltasigma': '%.5f'},Writer=ascii.CommentedHeader)
        print 'ds', ds

        r, ds = delta_sigma(galpos, gaspos, gas_pm_array, gas_downsampling_factor, rp_bins, Lbox)
        ds = ds/1e12
        table = {'rmpc': r,'deltasigma': ds}   
        ascii.write(table, 'ds_gas_mblack2_mhalo_matched_highM100.dat', formats={'rmpc': '%.5f', 'deltasigma': '%.5f'},Writer=ascii.CommentedHeader)
        print 'ds', ds

        r, ds = delta_sigma(galpos, starpos, star_pm_array, star_downsampling_factor, rp_bins, Lbox)
        ds = ds/1e12
        table = {'rmpc': r,'deltasigma': ds}   
        ascii.write(table, 'ds_star_mblack2_mhalo_matched_highM100.dat', formats={'rmpc': '%.5f', 'deltasigma': '%.5f'},Writer=ascii.CommentedHeader)
        print 'ds', ds

        #--------------------------
        #--------------------------
        
        # Read in the Galaxy File
        gal_list = np.loadtxt('/Users/alexie/Work/HSC/CompareHydroProfiles/DataFromMBII/mblack2_mhalo_matched_lowM100.asc')
        x=gal_list[:,4]
        y=gal_list[:,5]
        z=gal_list[:,6]
        # These are in units of kpc, put into Mpc
        x=x/1000.0
        y=y/1000.0
        z=z/1000.0
        # Format the array for halotools
        galpos = return_xyz_formatted_array(x, y, z)

        r, ds = delta_sigma(galpos, dmpos, dm_particle_mass, dm_downsampling_factor, rp_bins, Lbox)
        ds = ds/1e12
        table = {'rmpc': r,'deltasigma': ds}   
        ascii.write(table, 'ds_dm_mblack2_mhalo_matched_lowM100.dat', formats={'rmpc': '%.5f', 'deltasigma': '%.5f'},Writer=ascii.CommentedHeader)
        print 'ds', ds

        r, ds = delta_sigma(galpos, gaspos, gas_pm_array, gas_downsampling_factor, rp_bins, Lbox)
        ds = ds/1e12
        table = {'rmpc': r,'deltasigma': ds}   
        ascii.write(table, 'ds_gas_mblack2_mhalo_matched_lowM100.dat', formats={'rmpc': '%.5f', 'deltasigma': '%.5f'},Writer=ascii.CommentedHeader)
        print 'ds', ds

        r, ds = delta_sigma(galpos, starpos, star_pm_array, star_downsampling_factor, rp_bins, Lbox)
        ds = ds/1e12
        table = {'rmpc': r,'deltasigma': ds}   
        ascii.write(table, 'ds_star_mblack2_mhalo_matched_lowM100.dat', formats={'rmpc': '%.5f', 'deltasigma': '%.5f'},Writer=ascii.CommentedHeader)
        print 'ds', ds

        #--------------------------
        # M10
        #--------------------------
        
        # Read in the Galaxy File
        gal_list = np.loadtxt('/Users/alexie/Work/HSC/CompareHydroProfiles/DataFromMBII/mblack2_mhalo_matched_highM10.asc')
        x=gal_list[:,4]
        y=gal_list[:,5]
        z=gal_list[:,6]
        # These are in units of kpc, put into Mpc
        x=x/1000.0
        y=y/1000.0
        z=z/1000.0
        # Format the array for halotools
        galpos = return_xyz_formatted_array(x, y, z)

        r, ds = delta_sigma(galpos, dmpos, dm_particle_mass, dm_downsampling_factor, rp_bins, Lbox)
        ds = ds/1e12
        table = {'rmpc': r,'deltasigma': ds}   
        ascii.write(table, 'ds_dm_mblack2_mhalo_matched_highM10.dat', formats={'rmpc': '%.5f', 'deltasigma': '%.5f'},Writer=ascii.CommentedHeader)
        print 'ds', ds

        r, ds = delta_sigma(galpos, gaspos, gas_pm_array, gas_downsampling_factor, rp_bins, Lbox)
        ds = ds/1e12
        table = {'rmpc': r,'deltasigma': ds}   
        ascii.write(table, 'ds_gas_mblack2_mhalo_matched_highM10.dat', formats={'rmpc': '%.5f', 'deltasigma': '%.5f'},Writer=ascii.CommentedHeader)
        print 'ds', ds

        r, ds = delta_sigma(galpos, starpos, star_pm_array, star_downsampling_factor, rp_bins, Lbox)
        ds = ds/1e12
        table = {'rmpc': r,'deltasigma': ds}   
        ascii.write(table, 'ds_star_mblack2_mhalo_matched_highM10.dat', formats={'rmpc': '%.5f', 'deltasigma': '%.5f'},Writer=ascii.CommentedHeader)
        print 'ds', ds

        #--------------------------
        #--------------------------
        
        # Read in the Galaxy File
        gal_list = np.loadtxt('/Users/alexie/Work/HSC/CompareHydroProfiles/DataFromMBII/mblack2_mhalo_matched_lowM10.asc')
        x=gal_list[:,4]
        y=gal_list[:,5]
        z=gal_list[:,6]
        # These are in units of kpc, put into Mpc
        x=x/1000.0
        y=y/1000.0
        z=z/1000.0
        # Format the array for halotools
        galpos = return_xyz_formatted_array(x, y, z)

        r, ds = delta_sigma(galpos, dmpos, dm_particle_mass, dm_downsampling_factor, rp_bins, Lbox)
        ds = ds/1e12
        table = {'rmpc': r,'deltasigma': ds}   
        ascii.write(table, 'ds_dm_mblack2_mhalo_matched_lowM10.dat', formats={'rmpc': '%.5f', 'deltasigma': '%.5f'},Writer=ascii.CommentedHeader)
        print 'ds', ds

        r, ds = delta_sigma(galpos, gaspos, gas_pm_array, gas_downsampling_factor, rp_bins, Lbox)
        ds = ds/1e12
        table = {'rmpc': r,'deltasigma': ds}   
        ascii.write(table, 'ds_gas_mblack2_mhalo_matched_lowM10.dat', formats={'rmpc': '%.5f', 'deltasigma': '%.5f'},Writer=ascii.CommentedHeader)
        print 'ds', ds

        r, ds = delta_sigma(galpos, starpos, star_pm_array, star_downsampling_factor, rp_bins, Lbox)
        ds = ds/1e12
        table = {'rmpc': r,'deltasigma': ds}   
        ascii.write(table, 'ds_star_mblack2_mhalo_matched_lowM10.dat', formats={'rmpc': '%.5f', 'deltasigma': '%.5f'},Writer=ascii.CommentedHeader)
        print 'ds', ds








        
if __name__=='__main__':
    main()
