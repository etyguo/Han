#!/usr/bin/python

"""Plots all the cluster properties of interest for the singular, polytropic sphere model"""

import sys

from Nbody6IO import *

endtime = 20.
perc = np.arange( 0.1, 1., 0.1 )
colors = None
#The y-axis limits for the Lagrangian radii plot
yl = [0.02, 20]
yl2 = [0.3, 1.2]
xl = [1.e-2, 20]
mlow = 0.001
nlow = 0.001
plot_suff = '.eps'

try:
    fname = sys.argv[1]
except IndexError:
    fname = 'OUT3'

try:
    snaps = plot_surface_density_profiles(fname, save='mass_surface_density_profiles.png', legend=True, lowcut=mlow)
except IOError:
    sys.exit('No NBODY6 output file found or specified!')
    
plot_surface_density_profiles(fname, ylab=r'N $pc^{-2}$', save='number_surface_density_profiles' + plot_suff, legend=True, m_or_n='n', lowcut=nlow, snaps=snaps)

plot_mf(fname, save='imf' + plot_suff, logn=False)
plot_mf(fname, save='imf_logn' + plot_suff, logn=True)
plot_rf(fname, save='irf' + plot_suff, logn=False)
plot_vf(fname, save='ivf' + plot_suff, logn=False)

plot_lagrangian_radii_and_meanmass_in_time(fname, endtime, perc, colors, legend=True, save='lagrangian_radii_evolution' + plot_suff, ylim=yl, core=True, const_mp=True, xlim=xl, snapend=snaps[-1], ylim2=yl2)

plot_npairs_in_time(fname, endtime, save='npair_evolution' + plot_suff, snapend=snaps[-1])

plot_nstars_in_time(fname, endtime, save='ntot_evo' + plot_suff, snapend=snaps[-1])

plot_mtot_in_time(fname, endtime, save='mtot_evo' + plot_suff, snapend=snaps[-1])

plot_mean_masses_in_radii_in_time(fname, save='mean_mass_in_logr' + plot_suff, lgnd=True, snaps=snaps)

plot_most_massive(fname, endtime, save='most_massive_r' + plot_suff, snapend=snaps[-1])

track_most_massive(fname, endtime, save='most_massive_r' + plot_suff, snapend=snaps[-1])

try:
    plot_semi_dist(0, save='semi_dist_0' + plot_suff)
except IOError:
    print('No binaries here.')
