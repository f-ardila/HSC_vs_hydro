import numpy as np
import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

###################################################################################################

def main():

	plotStellarMaps('Data/galaxies_star_maps.hdf5')

	return

###################################################################################################

def plotStellarMaps(fn, hosts = True, sats = True, verbose = True):

	map_min = 1E4
	map_max = 1E10
	cmap = plt.cm.afmhot

	# Load general simulation and galaxy properties
	f = h5py.File(fn, 'r')
	snap_a = f.attrs['snap_a']
	n_galaxies = len(f['cat_sh_id'])
	sh_idx = np.array(f['cat_sh_id'])
	grp_idx = np.array(f['cat_grp_id'])
	cat_is_primary = np.array(f['cat_grp_is_primary'])
	cat_sh_mstar = np.array(f['cat_sh_mstar'])
	cat_sh_pos_bound = np.array(f['cat_sh_pos_bound'])
	cat_sh_halfmrad_stars = np.array(f['cat_sh_halfmrad_stars'])
	map_stars = np.array(f['map_stars'])
	map_size = f['map_stars_size'].value
	n_pixels = f['map_stars_npixel'].value
	f.close()

	print('Found %d galaxies, scale factor %.2f. Maps are %.1f kpc in size, %d pixels.' \
		% (n_galaxies, snap_a, map_size, n_pixels))

	# Plot galaxies
	for i in range(n_galaxies):

		if (not hosts and cat_is_primary[i]) or (not sats and not cat_is_primary[i]):
			continue

		if verbose:
			print('Galaxy %d, Subhalo ID %8d, Group ID %8d, host %d, M* %.1e, xyz [%6.2f %6.2f %6.2f]' \
				% (i, sh_idx[i], grp_idx[i], cat_is_primary[i], cat_sh_mstar[i],
				cat_sh_pos_bound[i, 0] * 0.001, cat_sh_pos_bound[i, 1] * 0.001, cat_sh_pos_bound[i, 2] * 0.001))

		# Plot
		n_dir = 3
		n_samples = 3
		sample_labels = [r'$\mathrm{Galaxy}$', r'$\mathrm{Group\ members}$', r'$\mathrm{Unbound}$']
		panel_size = 3.5
		cbar_size = 0.2
		cbar_label_size = 0.7
		space_left = 1.0
		space = 0.2
		label_fontsize = 14
		label_top = 0.93
		label_x = 0.02
		fwidth = space_left + n_dir * (panel_size + space)
		fheight = space_left + n_samples * (panel_size + space) + cbar_size + cbar_label_size + space
		fig = plt.figure(figsize = (fwidth, fheight))
		gs = gridspec.GridSpec(1 + n_samples, n_dir, height_ratios = [cbar_size / panel_size, 1.0, 1.0, 1.0])
		plt.subplots_adjust(left = space_left / fwidth, right = 1.0 - space / fwidth,
						bottom = space_left / fheight, top = 1.0 - space / fheight - cbar_label_size / fheight,
						hspace = space / panel_size, wspace = space / panel_size)
		panels = []
		for j in range(n_dir):

			dir_x, dir_y = direction(j)
			label_dirx = 'xyz'[dir_x]
			label_diry = 'xyz'[dir_y]

			panels.append([])
			panels[j].append(fig.add_subplot(gs[0, j]))
			plt.gca().set_xticklabels([])
			plt.gca().set_yticklabels([])
			if j != 1:
				plt.axis('off')

			for k in range(n_samples):
				panels[j].append(fig.add_subplot(gs[1 + k, j]))

				if j == 0:
					plt.ylabel(r'${\rm kpc}$')
				else:
					plt.gca().set_yticklabels([])
				if k == 2:
					plt.xlabel(r'${\rm kpc}$')
				else:
					plt.gca().set_xticklabels([])

				q_label = r'$\log_{10}\ \rho_*\ (M_{\odot}/{\rm kpc}^2)$'
				plt.xlim(-map_size, map_size)
				plt.ylim(-map_size, map_size)
				map = map_stars[i, k, j, :, :]
				map[map < map_min] = map_min
				map[map > map_max] = map_max
				map = np.log10(map)
				vmin = np.log10(map_min)
				vmax = np.log10(map_max)
				norm = mpl.colors.Normalize(vmin = vmin, vmax = vmax)
				plt.imshow(map.T, extent = [-map_size, map_size, -map_size, map_size],
						interpolation = 'nearest', cmap = cmap, vmin = vmin, vmax = vmax)

				label_color = 'w'
				if j == 0:
					plt.text(label_x, label_top, sample_labels[k],
						transform = plt.gca().transAxes, color = label_color, fontsize = label_fontsize)
				if k == 0:
					plt.text(0.87, label_top, r'$%c%c$' % (label_dirx, label_diry),
						transform = plt.gca().transAxes, color = label_color, fontsize = label_fontsize)

				# Add color bar
				if j == 1 and k == 0:
					plt.sca(panels[1][0])
					ax = plt.gca()
					ticks = np.arange(vmin, vmax + 1.0, 1.0)
					cb = mpl.colorbar.ColorbarBase(ax, orientation = 'horizontal', cmap = cmap,
								norm = norm, ticks = ticks, format = r'$%.0f$')
					cb.set_label(q_label, rotation = 0, labelpad = 8)
					cb.ax.xaxis.set_ticks_position('top')
					cb.ax.xaxis.set_label_position('top')
					cb.ax.xaxis.set_tick_params(pad = 5)

		# Add labels
		plt.sca(panels[0][0])
		label_color = 'k'
		plt.text(0.0, 3.0, r'$\mathrm{Subhalo\ idx:}\ %d$' % (sh_idx[i]),
				transform = plt.gca().transAxes, color = label_color, fontsize = label_fontsize)
		plt.text(0.0, 1.5, r'$\mathrm{Group\ idx:}\ %d$' % (grp_idx[i]),
				transform = plt.gca().transAxes, color = label_color, fontsize = label_fontsize)
		if cat_is_primary[i]:
			label = r'$\mathrm{Central}$'
		else:
			label = r'$\mathrm{Satellite}$'
		plt.text(0.0, 0.0, label,
				transform = plt.gca().transAxes, color = label_color, fontsize = label_fontsize)

		plt.text(0.5, 3.0, r'$M^* = %.1f$' % np.log10(cat_sh_mstar[i]),
				transform = plt.gca().transAxes, color = label_color, fontsize = label_fontsize)
		plt.text(0.5, 1.5, r'$R^*_{1/2} = %.1f$' % cat_sh_halfmrad_stars[i],
				transform = plt.gca().transAxes, color = label_color, fontsize = label_fontsize)

		plt.savefig('Figures/galaxy_%04d.pdf' % i)
		plt.clf()

	return

###################################################################################################

def direction(idx):

	if idx == 0:
		dir_x = 0
		dir_y = 1
	elif idx == 1:
		dir_x = 0
		dir_y = 2
	else:
		dir_x = 1
		dir_y = 2

	return dir_x, dir_y

###################################################################################################
# Trigger
###################################################################################################

if __name__ == "__main__":
	main()
