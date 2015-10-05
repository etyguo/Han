import struct
import numpy as np
import matplotlib.pyplot as plt

from matplotlib import rc_file
rc_file('/etc/matplotlibrc')

#from phoenix_functions import plot_bar_no_inlines

def plot_bar_no_inlines(hist, bins, fig=None, ax=None, style='k-'):
    """
    Plots a bar plot with no inner lines between adjacent bars.
    """
    if fig is None and ax is None:
        fig = plt.figure()

    npoints = 2 * np.size(bins) + 1
    x = np.zeros(npoints)
    y = np.zeros(npoints)

    j = 0
    for i in np.arange(npoints):
        if i == npoints - 1:
            x[i] = bins[j-1]
        elif i == 0:
            x[i] = bins[0]
        elif i % 2 == 1:
            x[i] = bins[j]
            j += 1
        elif i % 2 == 0:
            x[i] = bins[j]

    j = 0
    for i in np.arange(npoints):
        if i == 0:
            y[i] = 0
        elif i == npoints - 1:
            continue
        elif i == npoints - 2:
            y[i] = hist[j-1]
        elif i % 2 == 1:
            y[i] = hist[j]
        elif i % 2 == 0:
            y[i] = hist[j]
            j += 1

    if ax is None:
        plt.plot(x, y, style)
    else:
        ax.plot(x, y, style)

    return fig

myr_to_s = 3.15569e13
pc_to_km = 3.08567758e13
pc_to_au = 206265.
#G in SI units
G = 6.673e-11
days_to_s = 86400.
days_to_yrs = 0.00274
msun_to_kg = 1.988e30

#class HolgersNbody4Snapshot(object):

  #def add_fixed_record(self, names, type_lens, type_chars, idx= None):
    #if not idx is None:
      #self.current_idx= idx

    #self.current_idx+= self.rec_header
    #rec= {}
    #for name, tlen, tchar in zip(names, type_lens, type_chars):
      #rec[name]= (self.current_idx, tlen, tchar)
      #self.current_idx+= tlen
    #self.current_idx+= self.rec_header

    #return rec

  #def __init__(self, file, int_len=4, float_len=4, rec_header=8):
    #self.file= open(filename)
    #self.int_len= int_len
    #self.float_len= float_len
    #self.rec_header= 8
    #self.current_idx= idx
    #self.record_header= {}
    #for name in ['ntot', 'model', 'nrun']:

def nbody6_snapshot_len(n_stars):
  # header for first record (WRITE statement)
  snap_len= 4
  # 4 integers; first one is star number
  snap_len+= 4*4
  # "footer" for first record
  snap_len+= 4
  # header for main record
  snap_len+= 4
  # float array that contains some global information about the snapshot
  snap_len+= 4*20
  # the star masses 
  snap_len+= 4*n_stars
  # the positions of the stars (x1, y1, z1, x2, y2, z2, ...)
  snap_len+= 4*3*n_stars
  # their velocities
  snap_len+= 4*3*n_stars
  # densities 
  #snap_len+= 4*n_stars
  # potential
  #snap_len+= 4*n_stars
  # the star id's
  snap_len+= 4*n_stars
  # the star types
  #snap_len+= 4*n_stars
  # the footer
  snap_len+= 4

  return snap_len

def read_nbody6_header(f):
  pos= f.tell()
  chunk= f.read(32)
  f.seek(pos)
  ntot, model, nrun, nk= struct.unpack('=llll', chunk[4:4+16])

  return (ntot, model, nrun, nk)

def get_nbody6_header(raw_snap):
  ntot, model, nrun, nk= struct.unpack('=llll', raw_snap[4:4+16])

  return (ntot, model, nrun, nk)

def get_nbody6_ntot(raw_snap):
  """
  Returns the total number of stars (whether in binaries or not).
  """
  ntot= get_nbody6_header(raw_snap)[0]
  # Subtract the number of center of mass particles.
  nb= get_as_dict(raw_snap)['NPAIRS']

  return ntot-nb
  

def get_as_dict(raw_snap):
  as_keys= [ 'TTOT', 'NPAIRS', 'RBAR', 'ZMBAR', 'RTIDE', 'TIDAL4',
       'RDENS', 'TTOT/TCR', 'TSCALE', 'VSTAR', 'RC', 'NC',
       'VC', 'RHOM', 'CMAX', 'RSCALE', 'RSMIN', 'DMIN1' ]
  a_s= get_as_array(raw_snap)
  #populate the dict 
  as_dict= dict(zip(as_keys[:6], a_s[:6]))
  as_dict['RDENS']= a_s[6:9]
  as_dict.update(zip(as_keys[7:17], a_s[9:19]))
  as_dict['DMIN1']= a_s[-1]
  as_dict['NPAIRS']= int(as_dict['NPAIRS'])

  return as_dict

def count_nbody6_snapshots(filename):
  f= open(filename)
  f.seek(0, 2)
  end_pos= f.tell()
  f.seek(0)

  cur= 0; n_snap=0
  while cur< end_pos:
    ntot= read_nbody6_header(f)[0]
    snaplen= nbody6_snapshot_len(ntot)
    cur+= snaplen
    n_snap+= 1
    f.seek(snaplen, 1)

  if not cur==end_pos:
    print "cur", cur, 'end_pos', end_pos  

  f.close()
  return n_snap



def read_nbody6_snapshot_raw(filename, snap_no):
  """
  Reads a single snapshot at the snap_no'th position.
  """
  f=open(filename)
  pos=0;
  ntot= read_nbody6_header(f)[0]
  snaplen= nbody6_snapshot_len(ntot)
  for i in np.arange(snap_no):
    pos+= snaplen
    f.seek(pos)
    ntot= read_nbody6_header(f)[0]
    snaplen= nbody6_snapshot_len(ntot)

  snap= f.read(snaplen)
  f.close()
  return snap

def nbody6_snapshots_raw(filename):
  f= open(filename)
  nsnaps= count_nbody6_snapshots(filename)
  for i in np.arange(nsnaps):
    ntot= read_nbody6_header(f)[0]
    snaplen= nbody6_snapshot_len(ntot)
    yield f.read(snaplen)

  f.close()
  return

def as_offset():
  offset = 4 + 4*4 + 4 + 4

  return offset

def mass_offset(raw_snap):
  nk= get_nbody6_header(raw_snap)[-1]

  return as_offset()+4*nk

def xs_offset(raw_snap):
  ntot= get_nbody6_header(raw_snap)[0]

  return mass_offset(raw_snap)+ntot*4

def vs_offset(raw_snap):
  ntot= get_nbody6_header(raw_snap)[0]

  return xs_offset(raw_snap)+ntot*4*3

def rho1_offset(raw_snap):
  ntot= get_nbody6_header(raw_snap)[0]

  return vs_offset(raw_snap)+ ntot*4*3

def phi1_offset(raw_snap):
  ntot= get_nbody6_header(raw_snap)[0]

  return rho1_offset(raw_snap)+ ntot*4

def name_offset(raw_snap):
  ntot= get_nbody6_header(raw_snap)[0]

  return vs_offset(raw_snap)+ ntot*4*3

def kstar_offset(raw_snap):
  ntot= get_nbody6_header(raw_snap)[0]

  return name_offset(raw_snap)+ ntot*4


def get_mass_array(raw_snap):
  ntot= get_nbody6_header(raw_snap)[0]
  mass= np.array(struct.unpack(ntot*'f', raw_snap[mass_offset(raw_snap):xs_offset(raw_snap)]))

  return mass


def get_xs_array(raw_snap):
  ntot= get_nbody6_header(raw_snap)[0]
  xs= np.array(struct.unpack(3*ntot*'f', raw_snap[xs_offset(raw_snap):vs_offset(raw_snap)]))

  return xs.reshape(ntot, 3)

def get_vs_array(raw_snap):
  ntot = get_nbody6_header(raw_snap)[0]
  vs = np.array(struct.unpack(3*ntot*'f', raw_snap[vs_offset(raw_snap):name_offset(raw_snap)]))

  return vs.reshape(ntot, 3)


def get_star_masses(raw_snap):
  ntot= get_nbody6_ntot(raw_snap)
  mass_start= mass_offset(raw_snap)
  mass_end= mass_offset(raw_snap)+4*ntot
  mass= np.array(struct.unpack(ntot*'f', raw_snap[mass_start:mass_end]))

  return mass

def get_star_xs(raw_snap):
  ntot= get_nbody6_ntot(raw_snap)
  xs_start= xs_offset(raw_snap)
  xs_end= xs_offset(raw_snap)+4*ntot*3
  xs= np.array(struct.unpack(3*ntot*'f', raw_snap[xs_start:xs_end]))

  return xs.reshape(ntot, 3)

def get_as_array(raw_snap):
  nk=  get_nbody6_header(raw_snap)[-1]
  a_s= np.array(struct.unpack(nk*'f', raw_snap[as_offset():mass_offset(raw_snap)]))

  return a_s

def get_name_array(raw_snap):
  ntot=  get_nbody6_ntot(raw_snap)
  name= np.array(struct.unpack('='+ntot*'l', raw_snap[name_offset(raw_snap):kstar_offset(raw_snap)]))

  return name

def get_star_names(raw_snap):
  ntot=get_nbody6_header(raw_snap)[0]
  name_start= name_offset(raw_snap)
  name_end= name_offset(raw_snap)+ 4*ntot
  name= np.array(struct.unpack('='+ntot*'l', raw_snap[name_start:name_end]))

  return name

def get_kstar_array(raw_snap):
  nk=  get_nbody6_header(raw_snap)[-1]
  kstar_start= kstar_offset(raw_snap)
  ntot = get_nbody6_ntot(raw_snap)
  kstar_end= kstar_offset(raw_snap)+ 4*ntot
  kstar= np.array(struct.unpack('='+ntot*'l', raw_snap[kstar_start:kstar_end]))

  return kstar

def get_star_kstars(raw_snap):
  ntot= get_nbody6_ntot(raw_snap)
  kstar_start= kstar_offset(raw_snap)
  kstar_end= kstar_offset(raw_snap)+ 4*ntot
  kstar= np.array(struct.unpack('='+ntot*'l', raw_snap[kstar_start:kstar_end]))

  return kstar

def get_time(raw_snap):
  """
  Returns the time of the snapshot in physical units.
  """
  a= get_as_dict(raw_snap)
  
  return a['TTOT']*a['TSCALE']

def find_times(times, snapshotfile):
    snaps = np.arange(np.size(times))
    snapno = count_nbody6_snapshots(snapshotfile)

    found = 0
    diff = 9999
    for snap in np.arange(snapno):
        rawsnap = read_nbody6_snapshot_raw(snapshotfile, snap)
        t = get_time(rawsnap)

        if np.abs(t - times[found]) <= diff:
            diff = np.abs(t - times[found])
        elif np.abs(t - times[found]) > diff:
            snaps[found] = snap - 1
            found += 1
            diff = 9999

        if found >= len(times):
            break

    return snaps
            

def get_surface_radii(xs):
  """
  Returns the projected radial distance from the center in the xy-plane.
  """
  return np.sqrt(xs[ : , 0]**2 + xs[ : , 1]**2)

def get_radii(xs):
    return np.sqrt(xs[ : , 0]**2 + xs[ : , 1]**2 + xs[ : , 2]**2)

def get_bin_mids(bin_edges):
	"""Finds the midpoints of the bin edges given by numpy.histogram module.  This is essential 
	to plot the histogram data at the correct x coordinate.  
    """
	bin_mids = np.zeros(len(bin_edges)-1)
	for i in np.arange(len(bin_edges)-1):
		bin_mids[i] = (bin_edges[i+1]+bin_edges[i])/2.
	return bin_mids

def get_surface_density_from_mass_in_bins(binmass, binno, bin_edges, islog):
	"""
        Given the mass in each bin and total number in each bin it calculates the surface mass 
	densities, and surface number densities.  The errors are Poisson 1sigma.
    """
	dens = np.zeros(len(binmass))
	err_dens = np.zeros(len(binmass))
    
	for i in np.arange(len(binmass)):
		if islog==1:
			r1, r2 = 10.**bin_edges[i], 10.**bin_edges[i+1]
			dens[i] = binmass[i] / np.pi / (r2**2. - r1**2.)
		else:
			dens[i] = binmass[i] / np.pi / (bin_edges[i+1]**2. - bin_edges[i]**2.)
            
		err_dens[i] = dens[i] / np.sqrt(binno[i])

	return dens, err_dens

def get_surface_density_profile(snapshotfile, snapshotno, nobins=20, histrange=[-2, 1]):
    """
    snapshotfile: OUT3- the file which has all the snapshots.  
    snap_no: gives the snapshot index that needs to be loaded.  
    calculates the surface mass and number density profiles.  Can be easily 
    extended to luminosity density profile once the stellar properties are included.  
    Also returns the physical time for the snapshot.  
    """
    #reading a particular snapshot
    raw_snap = read_nbody6_snapshot_raw(snapshotfile, snapshotno)
    
    #Get the AS dict containing parameters of the snapshot
    as_dict = get_as_dict(raw_snap)
    
    #Get the positions, masses, and time from the snapshot. Note conversions using the AS dict.
    xs = get_star_xs(raw_snap) * as_dict['RBAR']
    mass = get_star_masses(raw_snap) * as_dict['ZMBAR']
    #get_time() already converts to Myr
    t = get_time(raw_snap)
    
    #Calculate the 2D radii (x and y)
    xs2d = get_surface_radii(xs)
    #Calculate the log of these radii
    xs2dlog = np.log10(xs2d)
    
    #if np.max(xs2dlog) > 1.:
    #    xs2dlogrange = [1., np.max(xs2dlog)]
    #else:
    #    xs2dlogrange = [np.min(xs2dlog), 1.]
    
    hist_mass, bin_edges = np.histogram(xs2dlog, bins=nobins, range=histrange, weights=mass)
    hist_no, bin_edges = np.histogram(xs2dlog, bins=nobins, range=histrange)
    
    bin_mids = get_bin_mids(bin_edges)
    
    dens, err_dens = get_surface_density_from_mass_in_bins(hist_mass, hist_no, bin_edges, 1)
    ndens, err_ndens = get_surface_density_from_mass_in_bins(hist_no, hist_no, bin_edges, 1)
	
    return t, bin_mids, dens, err_dens, ndens, err_ndens
	
def plot_surface_density_profile(snapshotfile, snapshotno, nobins=20, histrange=[-2, 1], color='k', fig=None, save=None, ylab=r'M$_{\odot}$ pc$^{-2}$', xlab=r'log(r/pc)', legend=None, m_or_n='m', lowcut=1e-3):
    t, mids, dens, err_dens, ndens, err_ndens = get_surface_density_profile(snapshotfile, snapshotno, nobins, histrange)

    if fig is None:
        plt.figure()
        plt.yscale('log')
        #plt.xscale('log')

    if m_or_n == 'm':
        for i in np.arange(np.size(dens)):
            if dens[i] - err_dens[i] <= 0.:
                err_dens[i] = dens[i] - lowcut
                
        plt.errorbar(mids, dens, yerr=err_dens, color=color)
    elif m_or_n == 'n':
        for i in np.arange(np.size(ndens)):
            if ndens[i] - err_ndens[i] <= 0.:
                err_ndens[i] = ndens[i] - lowcut
                
        plt.errorbar(mids, ndens, yerr=err_ndens, color=color)
        

    plt.ylabel(ylab)
    plt.xlabel(xlab)

    if legend is not None:
        plt.legend(legend)

    if save is not None:
        plt.savefig(save, dpi=200, clobber=True)

    return

def plot_surface_density_profiles(snapshotfile, times=[0., 1., 3., 10., 20.], colors=['k', 'm', 'b', 'g', 'r'], nobins=20, histrange=[-2, 1], fig=None, save=None, ylab=r'M$_{\odot}$ pc$^{-2}$', xlab=r'log(r/pc)', legend=None, m_or_n='m', lowcut=1e-3):
    snaps = find_times(times, snapshotfile)

    if legend is not None:
        leg = []
        for time in times:
            leg.append(str(int(time)) + ' Myr')
    else:
        leg = None

    fig = plt.figure()
    plt.yscale('log')
    #plt.xscale('log')

    for i, snap in enumerate(snaps):
        plot_surface_density_profile(snapshotfile, snap, nobins, histrange, colors[i], fig, save, ylab, xlab, leg, m_or_n, lowcut)

    return

def plot_mf(snapshotfile, snap=0, nobins=20, save=None, logn=True):
    rawsnap = read_nbody6_snapshot_raw(snapshotfile, snap)

    as_dict = get_as_dict(rawsnap)

    masses = get_star_masses(rawsnap) * as_dict['ZMBAR']

    logmass = np.log10(masses)

    hist, bins = np.histogram(logmass, nobins)

    plt.figure()

    if logn is True:
        plt.bar(bins[:-1], np.log10(hist), width=np.abs(bins[1]-bins[0]), facecolor='none', edgecolor='k')
        plt.ylabel(r'log(N)')
    else:
        plt.bar(bins[:-1], hist, width=np.abs(bins[1]-bins[0]), facecolor='none', edgecolor='k')
        plt.ylabel(r'N')
        
    plt.xlabel(r'log(M/M$_{\odot}$)')

    plt.xlim(bins[0], bins[-1])

    if save is not None:
        plt.savefig(save, dpi=200, clobber=True)

    return masses

def plot_rf(snapshotfile, snap=0, nobins=20, save=None, logn=True):
    rawsnap = read_nbody6_snapshot_raw(snapshotfile, snap)

    as_dict = get_as_dict(rawsnap)

    xs = get_star_xs(rawsnap) * as_dict['RBAR']
    rs = get_radii(xs)

    hist, bins = np.histogram(rs, nobins)

    plt.figure()

    if logn is True:
        plt.bar(bins[:-1], np.log10(hist), width=np.abs(bins[1]-bins[0]), facecolor='none', edgecolor='k')
        plt.ylabel(r'log(N)')
    else:
        plt.bar(bins[:-1], hist, width=np.abs(bins[1]-bins[0]), facecolor='none', edgecolor='k')
        plt.ylabel(r'N')
        
    plt.xlabel(r'r (pc)')

    plt.xlim(bins[0], bins[-1])

    if save is not None:
        plt.savefig(save, dpi=200, clobber=True)

    return rs

def plot_vf(snapshotfile, snap=0, nobins=20, save=None, logn=True):
    rawsnap = read_nbody6_snapshot_raw(snapshotfile, snap)

    as_dict = get_as_dict(rawsnap)

    vs = get_vs_array(rawsnap) * as_dict['RBAR'] / as_dict['TSCALE'] * pc_to_km / myr_to_s

    hist, bins = np.histogram(get_radii(vs), nobins)

    plt.figure()

    if logn is True:
        plt.bar(bins[:-1], np.log10(hist), width=np.abs(bins[1]-bins[0]), facecolor='none', edgecolor='k')
        plt.ylabel(r'log(N)')
    else:
        plt.bar(bins[:-1], hist, width=np.abs(bins[1]-bins[0]), facecolor='none', edgecolor='k')
        plt.ylabel(r'N')
        
    plt.xlabel(r'v (km/s)')

    plt.xlim(bins[0], bins[-1])

    if save is not None:
        plt.savefig(save, dpi=200, clobber=True)

    return vs


def lagrangian_radii(snapshotifle, snapno, percentiles=[0.05, 0.1, 0.25, 0.5, 0.75], mp=False):
    r_lagr = np.zeros(np.size(percentiles))
    
    rawsnap = read_nbody6_snapshot_raw(snapshotifle, snapno)

    as_dict = get_as_dict(rawsnap)

    xs = get_star_xs(rawsnap) * as_dict['RBAR']
    mass = get_star_masses(rawsnap) * as_dict['ZMBAR']
    rs = get_radii(xs)
    mtot = np.sum(mass)            

    sort = np.argsort(rs)

    per_cur = 0
    m_cur = 0.
    
    for i in sort:
        m_cur += mass[i]

        if mp is False and m_cur / mtot >= percentiles[per_cur]:
            r_lagr[per_cur] = rs[i]
            per_cur += 1
        elif mp is not False and m_cur >= mp[per_cur]:
            r_lagr[per_cur] = rs[i]
            per_cur += 1

        if per_cur >= np.size(percentiles):
            break

    return r_lagr, percentiles

def find_mp(snapshotfile, snapno=0,  percentiles=[0.05, 0.1, 0.25, 0.5, 0.75]):
    rawsnap = read_nbody6_snapshot_raw(snapshotfile, snapno)

    as_dict = get_as_dict(rawsnap)
    mass = get_star_masses(rawsnap) * as_dict['ZMBAR']
    mtot = np.sum(mass)
    
    mp = np.zeros(np.size(percentiles))

    for i, p in enumerate(percentiles):
        mp[i] = p * mtot

    return mp

def lagrangian_radii_in_time(snapshotfile, tend, percentiles=[0.05, 0.1, 0.25, 0.5, 0.75], const_mp=False):
    snapend = find_times([tend], snapshotfile)

    r_lagr = np.zeros((snapend + 1, np.size(percentiles)))
    t = np.zeros(snapend + 1)

    if const_mp is not False:
        mp = find_mp(snapshotfile, 0, percentiles)
    else:
        mp = False

    snap = 0
    while snap <= snapend:
        r_lagr[snap], perc = lagrangian_radii(snapshotfile, snap, percentiles, mp)
        rawsnap = read_nbody6_snapshot_raw(snapshotfile, snap)
        t[snap] = get_time(rawsnap)

        snap += 1

    return r_lagr, t, perc

def core_radius(snapshotfile, tend):
    snapend = find_times([tend], snapshotfile)

    rc = np.zeros(snapend + 1)
    t = np.zeros(snapend + 1)

    snap = 0
    while snap <= snapend:
        rawsnap = read_nbody6_snapshot_raw(snapshotfile, snap)
        t[snap] = get_time(rawsnap)
        asdict = get_as_dict(rawsnap)
        rc[snap] = asdict['RC'] * asdict['RBAR']

        snap += 1

    return rc, t

def plot_lagrangian_radii_in_time(snapshotfile, tend, percentiles=[0.05, 0.1, 0.25, 0.5, 0.75], plot_percentiles=[2, 3, 4], colors=['k', 'b', 'r'], legend=None, save=None, ylim=None, core=None, const_mp=False, xlim=None):
    r_lagr, t, perc = lagrangian_radii_in_time(snapshotfile, tend, percentiles, const_mp)

    plt.figure()
    plt.xscale('log')
    plt.yscale('log')

    if colors is not None:
        for j, i in enumerate(plot_percentiles):
            plt.plot(t, r_lagr[ : , i], color=colors[j])
    else:
        for j, i in enumerate(plot_percentiles):
            plt.plot(t, r_lagr[ : , i])        

    plt.xlabel(r't (Myr)')
    plt.ylabel(r'r (pc)')
    plt.xlim(t[0], t[-1])
    if ylim is not None:
        plt.ylim(ylim[0], ylim[1])
    if xlim is not None:
        plt.xlim(xlim[0], xlim[1])


    if core is not None:
        rc, t = core_radius(snapshotfile, tend)

        plt.plot(t, rc, 'k--')
        
    if legend is not None:
        lgnd = []
        for i in plot_percentiles:
            lgnd.append(str(int(percentiles[i] * 100)) + '\%')
            
        if core is not None:
            lgnd.append('Core')
            
        plt.legend(lgnd, loc='upper left')

    if save is not None:
        plt.savefig(save, dpi=200, clobber=True)

    return

def plot_npairs_in_time(snapshotfile, tend, save=None):
    snapend = find_times([tend], snapshotfile)

    npairs = np.zeros(snapend + 1)
    t = np.zeros(snapend + 1)

    snap = 0
    while snap <= snapend:
        rawsnap = read_nbody6_snapshot_raw(snapshotfile, snap)

        as_dict = get_as_dict(rawsnap)
        t[snap] = get_time(rawsnap)
        npairs[snap] = as_dict['NPAIRS']

        snap += 1

    plt.figure()
    plt.plot(t, npairs, color='k')
    plt.xlabel(r't (Myr)')
    plt.ylabel(r'N$_{bin}$')

    plt.ylim(np.min(npairs) - 1, np.max(npairs) + 1)

    if save is not None:
        plt.savefig(save, dpi=200, clobber=True)

    return

def plot_nstars_in_time(snapshotfile, tend, cutoff=0., save=None):
    snapend = find_times([tend], snapshotfile)

    nstars = np.zeros(snapend + 1)
    t = np.zeros(snapend + 1)

    snap = 0
    while snap <= snapend:
        rawsnap = read_nbody6_snapshot_raw(snapshotfile, snap)

        as_dict = get_as_dict(rawsnap)
        t[snap] = get_time(rawsnap)
        masses = get_star_masses(rawsnap) * as_dict['ZMBAR']
        nstars[snap] = np.size(np.where(masses > cutoff)[0])

        snap += 1

    plt.figure()
    plt.plot(t, nstars, color='k')
    plt.xlabel(r't (Myr)')
    if cutoff == 0.:
        plt.ylabel(r'N$_{tot}$')
    else:
        plt.ylabel(r'N(> 10 M$_{\odot}$)')
    plt.ylim(np.min(nstars) - 1, np.max(nstars) + 1)

    if save is not None:
        plt.savefig(save, dpi=200, clobber=True)

    return nstars

def plot_ntot_intime(snapshotfile, tend, save=None):
    snapend = find_times([tend], snapshotfile)

    n = np.zeros(snapend + 1)
    t = np.zeros(snapend + 1)

    snap = 0
    while snap <= snapend:
        rawsnap = read_nbody6_snapshot_raw(snapshotfile, snap)

        as_dict = get_as_dict(rawsnap)
        t[snap] = get_time(rawsnap)

        stars = get_star_masses(rawsnap)
        n[snap] = np.size(stars)

        snap += 1

    plt.figure()
    plt.plot(t, n, color='k')
    plt.xlabel(r't (Myr)')
    plt.ylabel(r'N$_{tot}$')

    plt.ylim(np.min(npairs) - 1, np.max(npairs) + 1)

    if save is not None:
        plt.savefig(save, dpi=200, clobber=True)

    return

def mean_mass_in_radii(snapshotfile, snap, rbinedges):
    mean_masses = np.zeros(np.size(rbinedges))
    mean_mass_errs = np.zeros(np.size(rbinedges))
    
    rawsnap = read_nbody6_snapshot_raw(snapshotfile, snap)

    as_dict = get_as_dict(rawsnap)

    xs = get_star_xs(rawsnap) * as_dict['RBAR']
    mass = get_star_masses(rawsnap) * as_dict['ZMBAR']
    rs = get_radii(xs)
    logrs = np.log10(rs)

    sort = np.argsort(logrs)

    binmasses = []
    bin_cur = 0
    for i in sort:
        if logrs[i] <= rbinedges[bin_cur]:
            binmasses.append(mass[i])
        else:
            mean_masses[bin_cur] = np.mean(binmasses)
            mean_mass_errs[bin_cur] = np.std(binmasses) / np.sqrt(np.size(binmasses))
            binmasses = []
            bin_cur += 1

        if bin_cur >= np.size(rbinedges):
            break

    return mean_masses, mean_mass_errs

def plot_mean_mass_in_radii(snapshotfile, snap, rbins, c='k', fig=None, save=None):
    mm, mmerr = mean_mass_in_radii(snapshotfile, snap, rbins[1:])

    rbinmids = get_bin_mids(rbins)

    if fig is None:
        plt.figure()

    plt.errorbar(rbinmids, mm, yerr=mmerr, color=c)

    plt.xlabel(r'log(r/pc)')
    plt.ylabel(r'$\langle$M$\rangle$ (M$_{\odot}$)')

    if save is not None:
        plt.savefig(save, dpi=200, clobber=True)

def plot_mean_masses_in_radii_in_time(snapshotfile, times=[0., 1., 3., 10., 20.], rbins=np.arange(-6, 2, 0.5), colors=['k', 'm', 'b', 'g', 'r'], save=None, lgnd=None):
    snaps = find_times(times, snapshotfile)

    fig = plt.figure()

    for i, snap in enumerate(snaps):
        plot_mean_mass_in_radii(snapshotfile, snap, rbins, colors[i], fig, save)

    if lgnd is not None:
        leg = []
        for time in times:
            leg.append(str(int(time)) + ' Myr')

        plt.legend(leg)

    if save is not None:
        plt.savefig(save, dpi=200, clobber=True)

    return

def find_most_massive(rawsnap, n=3, pairs=False):
    if pairs is not False:
        masses = get_mass_array(rawsnap)
    else: 
        masses = get_star_masses(rawsnap)


    sortargs = np.argsort(masses)

    return sortargs[ -n : ]

def plot_most_massive(snapshotfile, tend, n=3, pairs=False, colors=['r', 'b', 'g'], labels=['Third', 'Second', 'First'], save=None):
    snapend = find_times([tend], snapshotfile)

    most_massive_r = np.zeros((snapend + 1, n))
    hmr = np.zeros(snapend + 1)
    t = np.zeros(snapend + 1)
    most_massive_masses = np.zeros((snapend + 1, n))

    snap = 0
    while snap <= snapend:
        rawsnap = read_nbody6_snapshot_raw(snapshotfile, snap)

        as_dict = get_as_dict(rawsnap)
        t[snap] = get_time(rawsnap)

        most_massive = find_most_massive(rawsnap, n, pairs)

        if pairs is not False:
            xs = get_xs_array(rawsnap) * as_dict['RBAR']
            masses = get_mass_array(rawsnap) * as_dict['ZMBAR']
        else:
            xs = get_star_xs(rawsnap) * as_dict['RBAR']
            masses = get_star_masses(rawsnap) * as_dict['ZMBAR']
            
        rs = get_radii(xs)
        most_massive_r[snap] = rs[most_massive]

        most_massive_masses[snap] = masses[most_massive]
        
        hmr[snap] = lagrangian_radii(snapshotfile, snap, [0.5])[0]

        snap += 1

    plt.figure()
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r't (Myr)')
    plt.ylabel(r'r (pc)')
    plt.plot(t, hmr, 'k--', label=r'r$_h$')
    for i in np.arange(n):
        plt.plot(t, most_massive_r[ : , i ], colors[i], label=labels[i])

    if save is not None:
        plt.savefig(save, dpi=200, clobber=True)

    plt.legend(loc='upper left')

    plt.figure()

    for i in np.arange(n):
        plt.plot(t, most_massive_masses[ : , i ], colors[i], label=labels[i])

    plt.xlabel(r't (Myr)')
    plt.ylabel(r'Mass (M$_{\odot}$)')
    plt.legend(loc='upper left')
    plt.savefig(str(n) + '_most_massive.eps', dpi=200, clobber=True)    

    return most_massive_r, hmr, t, most_massive_masses

def track_most_massive(snapshotfile, tend, n=3, pairs=False, colors=['b', 'g', 'r'], labels=[r'$3^{rd}$', r'$2^{nd}$', r'$1^{st}$'], save=None):
    snapend = find_times([tend], snapshotfile)

    if n > 1:
        most_massive_r = np.zeros((snapend + 1, n))
        most_massive_masses = np.zeros((snapend + 1, n))
    else:
        most_massive_r = np.zeros((snapend + 1))
        most_massive_masses = np.zeros((snapend + 1))

    hmr = np.zeros(snapend + 1)
    t = np.zeros(snapend + 1)

    snap = 0
    while snap <= snapend:
        rawsnap = read_nbody6_snapshot_raw(snapshotfile, snap)

        most_massive = find_most_massive(rawsnap, n, pairs)
                        
        as_dict = get_as_dict(rawsnap)
        t[snap] = get_time(rawsnap)

        if pairs is not False:
            xs = get_xs_array(rawsnap) * as_dict['RBAR']
            masses = get_mass_array(rawsnap) * as_dict['ZMBAR']
        else:
            xs = get_star_xs(rawsnap) * as_dict['RBAR']
            masses = get_star_masses(rawsnap) * as_dict['ZMBAR']            
        
        rs = get_radii(xs)
        
        most_massive_r[snap] = rs[most_massive]

        most_massive_masses[snap] = masses[most_massive]

        hmr[snap] = lagrangian_radii(snapshotfile, snap, [0.5])[0]

        snap += 1

    plt.figure()
    plt.xscale('log')
    plt.xlim(0.1, 10)
    plt.ylim(0.01, 10)
    plt.yscale('log')
    plt.xlabel(r't (Myr)')
    plt.ylabel(r'r (pc)')
    plt.plot(t, hmr, 'k--', label=r'r$_h$')
    if n > 1:
        for i in np.arange(n):
            plt.plot(t, most_massive_r[ : , i ], colors[i], label=labels[i])
    else:
        plt.plot(t, most_massive_r[ : ], colors, label=labels)

    plt.legend(loc='lower right')

    if save is not None:
        plt.savefig(save, dpi=200, clobber=True)


    plt.figure()

    if n > 1:
        for i in np.arange(n):
            plt.plot(t, most_massive_masses[ : , i ], colors[i], label=labels[i])
    else:
        plt.plot(t, most_massive_masses[ : ], colors, label=labels)


    plt.xlabel(r't (Myr)')
    plt.ylabel(r'Mass (M$_{\odot}$)')
    plt.legend(loc='upper left')
    plt.savefig(str(n) + '_most_massive.png', dpi=200, clobber=True)    

    return most_massive_r, hmr, t, most_massive_masses, 
    
def plot_positions(snapshotfile, time, coord1=0, coord2=1, size=5, save=None):
    snap = find_times([time], snapshotfile)

    rawsnap = read_nbody6_snapshot_raw(snapshotfile, snap)

    as_dict = get_as_dict(rawsnap)
    t = get_time(rawsnap)

    xs = get_xs_array(rawsnap) * as_dict['RBAR']

    plt.figure()

    plt.scatter(xs[ : , coord1 ], xs[ : , coord2 ], color='r', s=size)
    if coord1 == 0:
        plt.xlabel('x (pc)')
    elif coord1 == 1:
        plt.xlabel('y (pc)')      
    elif coord1 == 2:
        plt.xlabel('z (pc)')
    if coord2 == 0:
        plt.ylabel('x (pc)')
    elif coord2 == 1:
        plt.ylabel('y (pc)')      
    elif coord2 == 2:
        plt.ylabel('z (pc)')

    if save is not None:
        plt.savefig(save, dpi=200, clobber=True)

    return

def plot_pos_most_massive(snapshotfile, time, n=3, coord1=0, coord2=1, size=5, scale=5, save=None):
    snap = find_times([time], snapshotfile)

    rawsnap = read_nbody6_snapshot_raw(snapshotfile, snap)

    as_dict = get_as_dict(rawsnap)
    t = get_time(rawsnap)

    xs = get_xs_array(rawsnap) * as_dict['RBAR']

    plt.figure()

    plt.scatter(xs[ : , coord1 ], xs[ : , coord2 ], color='r', s=size)
    if coord1 == 0:
        plt.xlabel('x (pc)')
    elif coord1 == 1:
        plt.xlabel('y (pc)')      
    elif coord1 == 2:
        plt.xlabel('z (pc)')
    if coord2 == 0:
        plt.ylabel('x (pc)')
    elif coord2 == 1:
        plt.ylabel('y (pc)')      
    elif coord2 == 2:
        plt.ylabel('z (pc)')

    most_massive = find_most_massive(rawsnap, n)

    plt.scatter(xs[most_massive, coord1], xs[most_massive, coord2], color='g', s=size*scale)

    if save is not None:
        plt.savefig(save, dpi=200, clobber=True)

    return

def veldisp_in_radii(snapshotfile, snap, rbinedges):
    veldisp = np.zeros(np.size(rbindedges))
    
    rawsnap = read_nbody6_snapshot_raw(snapshotfile, snap)

    as_dict = get_as_dict(rawsnap)

    xs = get_xs_array(rawsnap) * as_dict['RBAR']
    vels = get_vs_array()
    rs = get_radii(xs)
    logrs = np.log10(rs)

    sort = np.argsort(logrs)

    binvels = []
    bin_cur = 0
    for i in sort:
        if logrs[i] <= rbinedges[bin_cur]:
            binvels.append(vels[i])
        else:
            mean_masses[bin_cur] = np.mean(binmasses)
            mean_mass_errs[bin_cur] = np.std(binmasses) / np.sqrt(np.size(binmasses))
            binmasses = []
            bin_cur += 1

        if bin_cur >= np.size(rbinedges):
            break

    return mean_masses, mean_mass_errs

def print_mx(snapshotfile, snap, fname, objects=True):
    rawsnap = read_nbody6_snapshot_raw(snapshotfile, snap)

    as_dict = get_as_dict(rawsnap)

    starmasses = get_star_masses(rawsnap) * as_dict['ZMBAR']

    if objects is True:
        xs = get_xs_array(rawsnap) * as_dict['RBAR']
        mass = get_mass_array(rawsnap) * as_dict['ZMBAR']

        nbin = np.size(mass) - np.size(starmasses)

    f = open(fname, 'w')

    for i in np.arange(np.size(mass)):
        if objects is True and i == np.size(starmasses):
            f.write('\n' + str(mass[i]) + '\t' + str(xs[i, 0]) + '\t' + str(xs[i, 1]) + '\t' + str(xs[i, 2]) + '\n')
        else:
            f.write( str(mass[i]) + '\t' + str(xs[i, 0]) + '\t' + str(xs[i, 1]) + '\t' + str(xs[i, 2]) + '\n')

    f.close()

    return

def print_mxv(snapshotfile, snap, fname, objects=True):
    rawsnap = read_nbody6_snapshot_raw(snapshotfile, snap)

    as_dict = get_as_dict(rawsnap)

    starmasses = get_star_masses(rawsnap) * as_dict['ZMBAR']

    #Include binary centers of masses at the end of the file
    if objects is True:
        xs = get_xs_array(rawsnap) * as_dict['RBAR']
        mass = get_mass_array(rawsnap) * as_dict['ZMBAR']
        vs = get_vs_array(rawsnap) * as_dict['RBAR'] / as_dict['TSCALE'] * pc_to_km / myr_to_s

        nbin = np.size(mass) - np.size(starmasses)

    f = open(fname, 'w')

    for i in np.arange(np.size(mass)):
        if objects is True and i == np.size(starmasses):
            f.write('\n' + str(mass[i]) + '\t' + str(xs[i, 0]) + '\t' + str(xs[i, 1]) + '\t' + str(xs[i, 2]) + '\t' +  str(vs[i, 0]) + '\t' + str(vs[i, 1]) + '\t' + str(vs[i, 2]) + '\n')
        else:
            f.write( str(mass[i]) + '\t' + str(xs[i, 0]) + '\t' + str(xs[i, 1]) + '\t' + str(xs[i, 2]) + '\t' + str(vs[i, 0]) + '\t' + str(vs[i, 1]) + '\t' + str(vs[i, 2]) + '\n')

    f.close()

    return


def print_snap(snapshotfile, times=[0., 1., 3., 10., 20.], pre='mpos', suff='.dat'):
    snaps = find_times(times, snapshotfile)

    if pre == 'mpos':
        for i, snap in enumerate(snaps):
            print_mx(snapshotfile, snap, pre + str(int(times[i])) + suff)
    elif pre == 'mposvel':
        for i, snap in enumerate(snaps):
            print_mxv(snapshotfile, snap, pre + str(int(times[i])) + suff)

    return

def read_binout(snap, binfile='OUT9'):
    """
    Reads the binry output file for a specific snapshot.

    Parameters:
    -----------
    snap: int
        The snapshot to read
    binfile: string
        The file containing the binary properties for all snapshots. 
        Default is OUT9 since that is what NBODY6 automatically outputs thes to.

    Returns:
    --------
    p: float
        The period in days for each binary
    e: float
        The eccentricity for each binary.
    m1, m2: float
        The masses for each binary member in solar masses.
    """
    f = open(binfile, 'r')

    e = []
    m1 = []
    m2 = []
    p = []

    cursnap = -1
    i = 0
    while cursnap <= snap:
        line = f.readline()

        if line[0] == ' ' and i == 0:
            i += 1
        elif line[0] == ' ' and i == 1:
            i = 0
            cursnap += 1
        elif line[0] != ' ' and cursnap != snap:
            continue
        else:
            splt = line.split()
            e.append(float(splt[1]))
            m1.append(float(splt[4]))
            m2.append(float(splt[5]))
            p.append(float(splt[6]))

    return np.array(p), np.array(e), np.array(m1), np.array(m2)

def find_floats(x):
    """
    Finds only the floats in an array or list.

    Parameter:
    ----------
    x: array_like

    Returns:
    --------
    floats: indices of the float data types in x. E.g. the elements that aren't nan or inf.
    """
    floats = []

    for i in np.arange(np.size(x)):
        if -np.inf < x[i] < np.inf:
            floats.append(i)

    return floats

def semi_from_per2(per, m1, m2):
    """
    Calculates the semi-major axis in a binary using Kepler's Third Law.

    Paraneters:
    -----------
    per: float
        Period in days.
    m1, m2: float
        Masses in solar masses.

    Returns:
    --------
    The semi-major axis in AU.
    """
    per *= days_to_s
    m1 += msun_to_kg
    m2 += msun_to_kg
    
    return np.power(per**2 * G * (m1 + m2) / (4. * np.pi**2), 1./3.) / 1000. / pc_to_km * pc_to_au

def semi_from_per(per, m1, m2):
    """
    Calculates the semi-major axis in a binary using Kepler's Third Law.

    Paraneters:
    -----------
    per: float
        Period in days.
    m1, m2: float
        Masses in solar masses.

    Returns:
    --------
    The semi-major axis in AU.
    """
    per *= days_to_yrs
    
    return np.power(per**2 * (m1 + m2), 1./3.)

def plot_semi_dist(snap, snapshotfile='OUT3', binfile='OUT9', bins=np.arange(-3, 3.1, 0.1), save=None):
    rawsnap = read_nbody6_snapshot_raw(snapshotfile, snap)

    as_dict = get_as_dict(rawsnap)
    
    p, e, m1, m2 = read_binout(snap, binfile)

    #Calculate the semi-major axis in AU
    semi = semi_from_per(p, m1, m2)

    semi_hist, semi_bins = np.histogram(np.log10(semi), bins)

    fig = plot_bar_no_inlines(semi_hist, semi_bins)
    #plt.bar(semi_bins[:-1], semi_hist,width=np.abs(semi_bins[1]-semi_bins[0]), facecolor='none', edgecolor='k')
    plt.xlabel('log(a/AU)')
    plt.ylabel('N')
    plt.xlim(bins[0], bins[-1])

    if save is not None:
        plt.savefig(save, dpi=200, clobber=True)

    return
    

########################################################################################################################




def loop_over_lines(filename, bytesize=16*1024**2):
  f= open(filename)
  lines= f.readlines(bytesize)
  while lines:
    for l in lines:
      yield l
    del lines
    lines= f.readlines(bytesize)
  f.close()

  return

def goto_121_snapshot(filename, snap_no, bytesize=16*1024**2):
  # go to the start of the file
  lines= loop_over_lines(filename, bytesize)

  line=''
  # load a chunk into memory
  for cur_snap in range(snap_no):
    line= lines.next()
    while not line.strip().startswith('+++'):
      try:
        line= lines.next()
      except StopIteration:
        raise RuntimeError, 'End of file reached before reaching snapshot %i.'% snap_no
    if cur_snap>0: line= lines.next()

  return lines

def load_current_121_snapshot(lines_iter):
  try:
    time_line= lines_iter.next()
    snap_time= float(time_line)
    line= lines_iter.next()
    snap=[]
    while not line.strip().startswith('+++'):
      snap.append(map(float, line.split()))
      line=lines_iter.next()
    snap= array(snap)
    line= lines_iter.next()
  except StopIteration:
    raise RuntimeError, 'End of file reached before reaching end of snapshot.'

  return (snap_time, snap)

def load_121_snapshot(filename, snap_no, bytesize=16*1024**2):
  lines= goto_121_snapshot(filename, snap_no, bytesize)
  snap_time, snap= load_current_121_snapshot(lines)

  return (snap_time, snap)


def get_average_properties(no_of_snaps,start_snap,filename,M):
	"""input:
		no of snaps to take for averaging
		starting snap
		filename containing the snaps
		total mass at t=0
	return:
		overall cluster poperties at time t
		rc, M, N, rt, etc. 
	"""

	t_snaps = []
	rc_snaps = []
	rv_snaps = []
	M_snaps = []
	N_snaps = []
	rt_snaps = []
	Nc_snaps = []
	Mc_snaps = []

	for i in range(start_snap,start_snap+no_of_snaps,1):
		raw_snap = read_nbody6_snapshot_raw(filename,i)
		as_dict = get_as_dict(raw_snap)
		t_snaps.append(as_dict['TTOT']*as_dict['TSCALE'])
		rc_snaps.append(as_dict['RC']*as_dict['RSCALE'])
		rt_snaps.append(as_dict['RTIDE']*as_dict['RSCALE'])
		N_snaps.append(len(get_star_kstars(raw_snap)))
		M_snaps.append(sum(get_star_masses(raw_snap))*M)
		rv_snaps.append(as_dict['RBAR']*as_dict['RSCALE'])
		Nc_snaps.append(as_dict['NC'])
		m = get_star_masses(raw_snap)
		xs = get_xs_array(raw_snap)
		mc_tot = 0.
		dummy = 0
		try:
			for j in range(len(xs)):
    				r_2=0
    				for k in (0,1,2):
        				r_2 = r_2 + xs[j,k]**2
    				#print m[j]*as_dict['ZMBAR'], as_dict['ZMBAR'], rc_snaps[dummy]
    				r=r_2**0.5
    				if r<=rc_snaps[dummy]:
    					mc_tot = mc_tot + m[j]
		except IndexError:
			pass

		Mc_snaps.append(mc_tot*M)
		dummy = dummy+1


		
	t_myr, err_t_myr = mean(t_snaps), std(t_snaps)

	rc_pc, err_rc_pc = mean(rc_snaps), std(rc_snaps)
	
	rt_pc, err_rt_pc = mean(rt_snaps), std(rt_snaps)

	rv_pc, err_rv_pc = mean(rv_snaps), std(rv_snaps)
	
	M, err_M = mean(M_snaps), std(M_snaps)
	
	N, err_N = mean(N_snaps), std(N_snaps)

	Nc, err_Nc = mean(Nc_snaps), std(Nc_snaps)

	Mc, err_Mc = mean(Mc_snaps), std(Mc_snaps)

	print "t= %f +- %f rc= %f +- %f rt= %f +- %f rv= %f +- %f M= %f +- %f N= %f +- %f Nc= %f +- %f Mc= %f +- %f\n" %(t_myr, err_t_myr, rc_pc, err_rc_pc, rt_pc, err_rt_pc, rv_pc, err_rv_pc, M, err_M, N, err_N, Nc, err_Nc, Mc, err_Mc)

	return (t_myr, err_t_myr, rc_pc, err_rc_pc, rt_pc, err_rt_pc, rv_pc, err_rv_pc, M, err_M, N, err_N, Mc, err_Mc)


def sort_in_out(raw_snap):
	as_dict = get_as_dict(raw_snap)
	m = get_star_masses(raw_snap)
	xs = get_xs_array(raw_snap)
	rc = as_dict['RC']
	t = as_dict['TTOT']*as_dict['TSCALE']
	f=open('inside.dat','w')
	f1=open('outside.dat','w')
	f.write("#t = %f Myr %f NB rc = %f pc %f NB\n" %(t, as_dict['TTOT'], rc*as_dict['RSCALE'], rc))
	f.write("#1. r(NB) 2.m(NB) 3.id\n")
	f1.write("#t = %f Myr %f NB rc = %f pc %f NB\n" %(t, as_dict['TTOT'], rc*as_dict['RSCALE'], rc))
	f1.write("#1. r(NB) 2.m(NB) 3.id\n")
	inside_mass=[]
	outside_mass=[]
	try:
		for i in range(len(xs)):
    			r_2=0
    			for j in (0,1,2):
        			r_2 = r_2 + xs[i,j]**2
    			print m[i]*as_dict['ZMBAR'], as_dict['ZMBAR']
    			r=r_2**0.5
    			if r<2*rc:
        			f.write("%f %.12f %ld\n" %(r,m[i]*as_dict['ZMBAR'],i))
				inside_mass.append(log10(m[i]*as_dict['ZMBAR']))
    			else:
        			f1.write("%f %.12f %ld\n" %(r,m[i]*as_dict['ZMBAR'],i))
				outside_mass.append(log10(m[i]*as_dict['ZMBAR']))
	except IndexError:
		pass

	import gracePlot
	gpl=gracePlot.gracePlot()
	gpl.hold()
	inside_hd = histogram(inside_mass, 1000, range=(0.1,100.), normed=True)
	print inside_hd
	outside_hd = histogram(outside_mass, 1000, range=(0.1,100.), normed=True)

	gpl.plot(inside_hd[1][:-1], inside_hd[0])
	gpl.plot(outside_hd[1][:-1], outside_hd[0])

	return


def in_out_over_snaps(no_of_snaps,start_snap,filename):
	f=open('inside.dat','w')
	f1=open('outside.dat','w')
	t_snaps = []
	rc_snaps = []
	for i in range(start_snap,start_snap+no_of_snaps,1):
		raw_snap = read_nbody6_snapshot_raw(filename,i)
		as_dict = get_as_dict(raw_snap)
		t_snaps.append(as_dict['TTOT']*as_dict['TSCALE'])
		rc_snaps.append(as_dict['RC']*as_dict['RSCALE'])
	t_myr = mean(t_snaps)
	err_t_myr = std(t_snaps)
	rc_pc = mean(rc_snaps)
	err_rc_pc = std(rc_snaps)

	f.write("#t = %f +- %f Myr rc = %f +- %f pc\n" %(t_myr, err_t_myr, rc_pc, err_rc_pc))
	f.write("#1. r(pc) 2.m(Msun) 3.id\n")
	f1.write("#t = %f +- %f Myr rc = %f +- %f pc\n" %(t_myr, err_t_myr, rc_pc, err_rc_pc))
	f1.write("#1. r(pc) 2.m(Msun) 3.id\n")

	inside_mass_all=[]
	outside_mass_all=[]
	inside_m_snaps={}
	outside_m_snaps={}
	for k in range(start_snap,start_snap+no_of_snaps,1):
		raw_snap = read_nbody6_snapshot_raw(filename,k)
		as_dict = get_as_dict(raw_snap)
		m = get_star_masses(raw_snap)
		xs = get_xs_array(raw_snap)
		rc = as_dict['RC']
		t = as_dict['TTOT']*as_dict['TSCALE']
		print t
		inside_mass=[]
		outside_mass=[]
		try:
			for i in range(len(xs)):
    				r_2=0
    				for j in (0,1,2):
        				r_2 = r_2 + xs[i,j]**2
    				#print m[i]*as_dict['ZMBAR'], as_dict['ZMBAR']
    				r=r_2**0.5
    				if r<2*rc:
        				f.write("%f %.12f %ld\n" %(r*as_dict['RSCALE'],m[i]*as_dict['ZMBAR'],i))
					#inside_mass.append(log10(m[i]*as_dict['ZMBAR']))
					inside_mass.append(m[i]*as_dict['ZMBAR'])
					#inside_mass_all.append(log10(m[i]*as_dict['ZMBAR']))
					inside_mass_all.append(m[i]*as_dict['ZMBAR'])
					
    				else:
        				f1.write("%f %.12f %ld\n" %(r*as_dict['RSCALE'],m[i]*as_dict['ZMBAR'],i))
					#outside_mass.append(log10(m[i]*as_dict['ZMBAR']))
					#outside_mass_all.append(log10(m[i]*as_dict['ZMBAR']))
					outside_mass.append(m[i]*as_dict['ZMBAR'])
					outside_mass_all.append(m[i]*as_dict['ZMBAR'])

		except IndexError:
			pass

		inside_m_snaps[k] = inside_mass	
		outside_m_snaps[k] = outside_mass
		#print inside_m_snaps[k], outside_m_snaps[k]

	#print inside_m_snaps, outside_m_snaps
	#inside_hist_all, in_all_hist_edges = histogram(inside_mass_all, 100, range=(-1,2), normed=True)
	#outside_hist_all, out_all_hist_edges = histogram(outside_mass_all, 100, range=(-1,2), normed=True)
	inside_hist_all, in_all_hist_edges = histogram(inside_mass_all, 100, range=(0.1,100.), normed=True)
	outside_hist_all, out_all_hist_edges = histogram(outside_mass_all, 100, range=(0.1,100.), normed=True)

	#print inside_hist_all, outside_hist_all
	inside_hist = {}
	inside_edges = {}
	outside_hist = {}
	outside_edges = {}
	for i in inside_m_snaps.keys():
		#inside_hist[i], inside_edges[i] = histogram(inside_m_snaps[k], 100, range=(-1,2), normed=True)
		#outside_hist[i], outside_edges[i] = histogram(outside_m_snaps[k], 100, range=(-1,2), normed=True)
		inside_hist[i], inside_edges[i] = histogram(inside_m_snaps[k], 100, range=(0.1,100.), normed=True)
		outside_hist[i], outside_edges[i] = histogram(outside_m_snaps[k], 100, range=(0.1,100.), normed=True)
	#print inside_hist, outside_hist

	min = []
	pin = []
	dpin = []
	mout = []
	pout = []
	dpout = []
	for i in range(len(inside_hist_all)):
		binwidth = in_all_hist_edges[-1]/len(in_all_hist_edges)

		min.append(in_all_hist_edges[i]+0.5*binwidth)
		pin.append(inside_hist_all[i])
		errin = 0
		for j in inside_hist.keys():
			if inside_edges[j][i]==in_all_hist_edges[i]:
				errin = errin + (inside_hist_all[i] - inside_hist[j][i])**2
		dpin.append((errin**0.5)/len(inside_hist))

		mout.append(out_all_hist_edges[i]+0.5*binwidth)
		pout.append(outside_hist_all[i])
		errout = 0
		for j in outside_hist.keys():
			if outside_edges[j][i]==out_all_hist_edges[i]:
				errout = errout + (outside_hist_all[i] - outside_hist[j][i])**2
		dpout.append((errout**0.5)/len(outside_hist))

		#print "%.8f %.8f %.8f %.8f %.8f %.8f\n" %(min[i], pin[i], dpin[i], mout[i], pout[i], dpout[i])
	
	
	import gracePlot
	gpl=gracePlot.gracePlot()
	gpl.hold()
	gpl.plot(min, pin, dy=dpin, symbols=1, styles=1)
	gpl.plot(mout, pout, dy=dpout, symbols=1, styles=1)

	f.close()
	f1.close()

	return


def in_out_unequal_bins(no_of_snaps,start_snap,filename, bin_mem_in, bin_mem_out, insidefile, outsidefile):
	"""input:
		no of snaps to take for averaging
		starting snap
		filename containing the snaps
		no of members wanted in each bin: bins are divided with equal number of stars in each bin
	output:
		plots the mass function within 2rc of the cluster
		also prints out the masses and the time etc. in two different files
			r<2rc and r>2rc stellar masses
		insidefile: MF data inside r=2rc
		outsidefile: MF data outside r=2rc
	 """

	f=open('inside_3myr.dat','w')
	f1=open('outside_3myr.dat','w')
	t_snaps = []
	rc_snaps = []
	for i in range(start_snap,start_snap+no_of_snaps,1):
		raw_snap = read_nbody6_snapshot_raw(filename,i)
		as_dict = get_as_dict(raw_snap)
		t_snaps.append(as_dict['TTOT']*as_dict['TSCALE'])
		rc_snaps.append(as_dict['RC']*as_dict['RSCALE'])
	t_myr = mean(t_snaps)
	err_t_myr = std(t_snaps)
	rc_pc = mean(rc_snaps)
	err_rc_pc = std(rc_snaps)

	f.write("#t = %f +- %f Myr rc = %f +- %f pc\n" %(t_myr, err_t_myr, rc_pc, err_rc_pc))
	f.write("#1. r(pc) 2.m(Msun) 3.id\n")
	f1.write("#t = %f +- %f Myr rc = %f +- %f pc\n" %(t_myr, err_t_myr, rc_pc, err_rc_pc))
	f1.write("#1. r(pc) 2.m(Msun) 3.id\n")

	inside_mass_all=[]
	outside_mass_all=[]
	for k in range(start_snap,start_snap+no_of_snaps,1):
		raw_snap = read_nbody6_snapshot_raw(filename,k)
		as_dict = get_as_dict(raw_snap)
		m = get_star_masses(raw_snap)
		xs = get_xs_array(raw_snap)
		rc = as_dict['RC']
		t = as_dict['TTOT']*as_dict['TSCALE']
		print t
		try:
			for i in range(len(xs)):
    				r_2=0
    				for j in (0,1,2):
        				r_2 = r_2 + xs[i,j]**2
    				r=r_2**0.5
    				if r<2*rc:
        				f.write("%f %.12f %ld\n" %(r*as_dict['RSCALE'],m[i]*as_dict['ZMBAR'],i))
					inside_mass_all.append(m[i]*as_dict['ZMBAR'])
					
    				else:
        				f1.write("%f %.12f %ld\n" %(r*as_dict['RSCALE'],m[i]*as_dict['ZMBAR'],i))
					outside_mass_all.append(m[i]*as_dict['ZMBAR'])

		except IndexError:
			pass

	sorted_mass_in = sort(inside_mass_all)
	sorted_mass_in_log = []
	for i in range(len(sorted_mass_in)):
		sorted_mass_in_log.append(log10(sorted_mass_in[i]))
	no_bins_in = len(sorted_mass_in)/float(bin_mem_in)
	#print len(sorted_mass_in), no_bin_mem, width

	min, min_log = [], []
	pin, pin_log = [], []
	dpin, dpin_log = [], []
	win, win_log = 0., 0.
	#print range(0,len(sorted_mass_in), int(width))

	for i in range(no_bins_in):
		
		min.append((sorted_mass_in[int(bin_mem_in*i)]+sorted_mass_in[int(bin_mem_in*(i+1))])/2.)
		min_log.append( 0.5 * (sorted_mass_in_log[int(bin_mem_in*i)] + sorted_mass_in_log[int(bin_mem_in*(i+1))]) )

		pin.append(float(bin_mem_in)/(sorted_mass_in[int(bin_mem_in*(i+1))]-sorted_mass_in[int(bin_mem_in*i)]))
		pin_log.append(float(bin_mem_in)/(sorted_mass_in_log[int(bin_mem_in*(i+1))]-sorted_mass_in_log[int(bin_mem_in*i)]))

		dpin.append((float(bin_mem_in))**0.5/(sorted_mass_in[int(bin_mem_in*(i+1))]-sorted_mass_in[int(bin_mem_in*i)]))
		dpin_log.append((float(bin_mem_in))**0.5/(sorted_mass_in_log[int(bin_mem_in*(i+1))]-sorted_mass_in_log[int(bin_mem_in*i)]))

		win = win + float(bin_mem_in)/(sorted_mass_in[int(bin_mem_in*(i+1))]-sorted_mass_in[int(bin_mem_in*i)])
		win_log = win_log + float(bin_mem_in)/(sorted_mass_in_log[int(bin_mem_in*(i+1))]-sorted_mass_in_log[int(bin_mem_in*i)])

		print "m_l= %f m_u= %f m_mid= %f dn/dm= %f norm= %f i= %d" %(sorted_mass_in[int(bin_mem_in*i)], sorted_mass_in[int(bin_mem_in*(i+1))], min[i], pin[i], win, i)
		print "log_m_l= %f log_m_u= %f log_m_mid= %f dn/dlogm= %f log_norm= %f i= %d" %(sorted_mass_in[int(bin_mem_in*i)], sorted_mass_in[int(bin_mem_in*(i+1))], min[i], pin[i], win, i)
	
	insidefile=open(insidefile,'w')
	insidefile.write("#MF at t= %f +- %f, rc= %f +- %f, #snaps for averaging= %d, #stars/bin= %d\n#1.m_bin (MSun) 2.p_bin 3.dp_bin 4.log_m_bin(m/MSun) 5.dn/dlogm 6.err(dn/dlogm)\n" %(t_myr, err_t_myr, rc_pc, err_rc_pc, no_of_snaps, bin_mem_in))
	
	for i in range(len(pin)):
		pin[i] = pin[i]/float(len(sorted_mass_in))
		pin_log[i] = pin_log[i]/float(len(sorted_mass_in))

		dpin[i] = dpin[i]/float(len(sorted_mass_in))
		dpin_log[i] = dpin_log[i]/float(len(sorted_mass_in))

		insidefile.write("%.12f %.12f %.12f %.12f %.12f %.12f\n" %(min[i], pin[i], dpin[i], min_log[i], pin_log[i], dpin_log[i]))
	insidefile.close()
	

	sorted_mass_out = sort(outside_mass_all)
	no_bins_out = len(sorted_mass_out)/float(bin_mem_out)
	#print len(sorted_mass_in), bin_mem_out, width

	mout = []
	pout = []
	dpout = []
	wout = 0.
	#print range(0,len(sorted_mass_in), int(width))

	for i in range(no_bins_out):
		
		mout.append((sorted_mass_out[int(bin_mem_out*i)]+sorted_mass_out[int(bin_mem_out*(i+1))])/2.)
		pout.append(float(bin_mem_out)/(sorted_mass_out[int(bin_mem_out*(i+1))]-sorted_mass_out[int(bin_mem_out*i)]))
		dpout.append((float(bin_mem_out))**0.5/(sorted_mass_out[int(bin_mem_out*(i+1))]-sorted_mass_out[int(bin_mem_out*i)]))
		wout = wout + float(bin_mem_out)/(sorted_mass_out[int(bin_mem_out*(i+1))]-sorted_mass_out[int(bin_mem_out*i)])

		#print sorted_mass_out[int(bin_mem_out*i)], sorted_mass_out[int(bin_mem_out*(i+1))], mout[i], pout[i], wout, i
	
	outsidefile=open(outsidefile,'w')
	outsidefile.write("#MF at t= %f +- %f, rc= %f +- %f, #snaps for averaging= %d, #stars/bin= %d\n#1.m_bin (MSun) 2.p_bin 3.dp_bin\n" %(t_myr, err_t_myr, rc_pc, err_rc_pc, no_of_snaps, bin_mem_out))
	
	for i in range(len(pout)):
		pout[i] = pout[i]/float(len(sorted_mass_out))
		dpout[i] = dpout[i]/float(len(sorted_mass_out))
		outsidefile.write("%.12f %.12f %.12f\n" %(mout[i], pout[i], dpout[i]))
	outsidefile.close()

	
	import gracePlot
	gpl=gracePlot.gracePlot()
	gpl.hold()
	gpl.plot(min, pin, styles=1)
	gpl.plot(min_log, pin_log, styles=1)
		
	f.close()
	f1.close()

	return



def in_out_logequal_bins(no_of_snaps,start_snap,filename, insidefile):
	"""input:
		no of snaps to take for averaging
		starting snap
		filename containing the snaps
		no of members wanted in each bin: bins are divided with equal number of stars in each bin
	output:
		plots the mass function within 2rc of the cluster
		also prints out the masses and the time etc. in two different files
			r<2rc and r>2rc stellar masses
		insidefile: MF data inside r=2rc
		outsidefile: MF data outside r=2rc
	 """

	t_snaps = []
	rc_snaps = []
	for i in range(start_snap,start_snap+no_of_snaps,1):
		raw_snap = read_nbody6_snapshot_raw(filename,i)
		as_dict = get_as_dict(raw_snap)
		t_snaps.append(as_dict['TTOT']*as_dict['TSCALE'])
		rc_snaps.append(as_dict['RC']*as_dict['RSCALE'])
	t_myr = mean(t_snaps)
	err_t_myr = std(t_snaps)
	rc_pc = mean(rc_snaps)
	err_rc_pc = std(rc_snaps)
	
	inside_masses='inside_'+str(t_myr)+'.dat'
	f=open(inside_masses, 'w')
	outside_masses='inside_'+str(t_myr)+'.dat'
	f1=open(outside_masses, 'w')

	f.write("#t = %f +- %f Myr rc = %f +- %f pc\n" %(t_myr, err_t_myr, rc_pc, err_rc_pc))
	f.write("#1. r(pc) 2.m(Msun) 3.id\n")
	f1.write("#t = %f +- %f Myr rc = %f +- %f pc\n" %(t_myr, err_t_myr, rc_pc, err_rc_pc))
	f1.write("#1. r(pc) 2.m(Msun) 3.id\n")

	inside_mass_all=[]
	outside_mass_all=[]
	for k in range(start_snap,start_snap+no_of_snaps,1):
		raw_snap = read_nbody6_snapshot_raw(filename,k)
		as_dict = get_as_dict(raw_snap)
		m = get_star_masses(raw_snap)
		xs = get_xs_array(raw_snap)
		rc = as_dict['RC']
		t = as_dict['TTOT']*as_dict['TSCALE']
		print t
		try:
			for i in range(len(xs)):
    				r_2=0
    				for j in (0,1):
        				r_2 = r_2 + xs[i,j]**2
    				r=r_2**0.5
    				if r<2*0.3:
        				f.write("%f %.12f %ld\n" %(r*as_dict['RSCALE'],m[i]*as_dict['ZMBAR'],i))
					inside_mass_all.append(m[i]*as_dict['ZMBAR'])
					
    				else:
        				f1.write("%f %.12f %ld\n" %(r*as_dict['RSCALE'],m[i]*as_dict['ZMBAR'],i))
					outside_mass_all.append(m[i]*as_dict['ZMBAR'])

		except IndexError:
			pass

	sorted_mass_in = sort(inside_mass_all)
	sorted_mass_in_log = []
	for i in range(len(sorted_mass_in)):
		sorted_mass_in_log.append(log10(sorted_mass_in[i]))
	no_bins_in = 100
	binwidth = abs((sorted_mass_in_log[-1] - sorted_mass_in_log[0])/100.)
	#print len(sorted_mass_in), no_bin_mem, width

	min_log = []
	nin_log, dnin_log = [], []
	for i in range(no_bins_in):
		
		#min.append((sorted_mass_in[int(bin_mem_in*i]+sorted_mass_in[int(bin_mem_in*(i+1))])/2.)
		min_log.append( 0.5*binwidth + sorted_mass_in_log[i] )

		#pin.append(float(bin_mem_in)/(sorted_mass_in[int(bin_mem_in*(i+1))]-sorted_mass_in[int(bin_mem_in*i)]))
		n = 0
		for k in range(len(sorted_mass_in_log)):
			if sorted_mass_in_log[k] < sorted_mass_in_log[i]+binwidth and sorted_mass_in_log[k] >= sorted_mass_in_log[i]:
				n += 1
		nin_log.append(n)
		dnin_log.append(n**0.5)		

	insidefile=open(insidefile,'w')
	insidefile.write("#MF at t= %f +- %f, rc= %f +- %f, #snaps for averaging= %d, #stars/bin= %d\n#1.log_m_bin(m/MSun) 2.dn/dlogm 3.err(dn/dlogm)\n" %(t_myr, err_t_myr, rc_pc, err_rc_pc, no_of_snaps, no_bins_in))
	
	for i in range(len(min_log)):
	#	pin[i] = pin[i]/float(len(sorted_mass_in))
	#	pin_log[i] = pin_log[i]/float(len(sorted_mass_in))
#
#		dpin[i] = dpin[i]/float(len(sorted_mass_in))
#		dpin_log[i] = dpin_log[i]/float(len(sorted_mass_in))

		insidefile.write("%.12f %.12f %.12f\n" %(min_log[i], nin_log[i], dnin_log[i]))
	insidefile.close()
	

	return


#########################################################################################################

def in_out_logequal_bins_onesnap(snapno,filename, annulus, data_mrange, norm_mrange, nbin):
	"""input:
		no of snaps to take for averaging
		starting snap
		filename containing the snaps
		no of members wanted in each bin: bins are divided with equal number of stars in each bin
	output:
		plots the mass function within 2rc of the cluster
		also prints out the masses and the time etc. in two different files
			r<2rc and r>2rc stellar masses
		insidefile: MF data inside r=2rc
		outsidefile: MF data outside r=2rc
	 """

	raw_snap = read_nbody6_snapshot_raw(filename,snapno)
	as_dict = get_as_dict(raw_snap)
	m = get_star_masses(raw_snap)
	xs = get_xs_array(raw_snap)
	rc = as_dict['RC']
	t = as_dict['TTOT']*as_dict['TSCALE']
	print 't', t, 'rc', rc
	outhist, inhist = {}, {}

	#for i in range(len(m)):
	#	m[i] = log10(m[i]*as_dict['ZMBAR'])
	for i in range(len(data_mrange)):
		data_mrange[i] = log10(data_mrange[i])
	for i in range(len(norm_mrange)):
		norm_mrange[i] = log10(norm_mrange[i])
	

	#divide mass in inside and outside
	filestring = 'MF_t'+str(t)+'_rc'+str(rc)

	#3D
	infilestring = filestring+'in_3D'
	outfilestring = filestring+'out_3D'
	massin_3D, massout_3D = mass_in_out(as_dict, m, xs, (0,1,2), annulus, infilestring, outfilestring)
	
	#norm = find_norm(m, binwidth, range)
	outhist['3D'] = histogram_equal_bins(massout_3D, nbin, data_mrange, norm_mrange)
	inhist['3D'] = histogram_equal_bins(massin_3D, nbin, data_mrange, norm_mrange)

	#print outhist['3D'], inhist['3D']

	#2D
	#massin_2D, massout_2D = [], []
	histout, histin = {}, {}
	#infilestring = inmassfilestring+'_2D'
	#outfilestring = outmassfilestring+'_2D'

	for k in (0,1,2):
		proj = (k%3, (k+1)%3)
		infilestring = filestring+'in_2D_'+str(proj[0])+str(proj[1])
		outfilestring = filestring+'out_2D'+str(proj[0])+str(proj[1])
		massin_2D, massout_2D = mass_in_out(as_dict, m, xs, proj, annulus, infilestring, outfilestring)
	
		histout[k] = histogram_equal_bins(massout_2D, nbin, data_mrange, norm_mrange)
		histin[k] = histogram_equal_bins(massin_2D, nbin, data_mrange, norm_mrange)
		print len(histout[k]), k, proj
	
	outhist['2D'] = zeros((len(histout[0]), 3))
	inhist['2D'] = zeros((len(histout[0]), 3))

	for i in range(len(histout[0])):
		mpoint = histout[0][i,0]
		n_ave = (histout[0][i,1] + histout[1][i,1] + histout[2][i,1])/3.
		n_poiss = (histout[0][i,2] + histout[1][i,2] + histout[2][i,2])/3.
		n_proj_err = (max(histout[0][i,1], histout[1][i,1], histout[2][i,1]) - min(histout[0][i,1], histout[1][i,1], histout[2][i,1]) )/2.
		tot_err_n = n_poiss+n_proj_err

		outhist['2D'][i,0], outhist['2D'][i,1], outhist['2D'][i,2] = mpoint, n_ave, tot_err_n

		mpoint = histin[0][i,0]
		n_ave = (histin[0][i,1] + histin[1][i,1] + histin[2][i,1])/3.
		n_poiss = (histin[0][i,2] + histin[1][i,2] + histin[2][i,2])/3.
		n_proj_err = (max(histin[0][i,1], histin[1][i,1], histin[2][i,1]) - min(histin[0][i,1], histin[1][i,1], histin[2][i,1]) )/2.
		tot_err_n = n_poiss+n_proj_err

		inhist['2D'][i,0], inhist['2D'][i,1], inhist['2D'][i,2] = mpoint, n_ave, tot_err_n
	
	inhistfile = filestring+'hist_in.dat'
	outhistfile = filestring+'hist_out.dat'
	inhistfile = open(inhistfile, 'w')
	outhistfile = open(outhistfile, 'w')
	for i in range(nbin):
		for j in range(3):
			inhistfile.write("%f " %( inhist['3D'][i,j], ))
			outhistfile.write("%f " %( outhist['3D'][i,j], ))
		for j in range(1,3):
			inhistfile.write("%f " %(inhist['2D'][i,j]))
			outhistfile.write("%f " %(outhist['2D'][i,j]))
		inhistfile.write("\n")
		outhistfile.write("\n")

	inhistfile.close()
	outhistfile.close()
	return inhist, outhist


def mass_in_out(as_dict, m, xs, proj, annulus, infilestring, outfilestring):
	""""""
	print m
	infile = infilestring+'.dat'
	outfile = outfilestring+'.dat'
	infile = open(infile,'w')
	outfile = open(outfile, 'w')
	inmass, outmass = [], []
	for i in range(len(xs)):
    		r_2=0.
    		for j in proj:
        		r_2 = r_2 + xs[i,j]**2
    		r_pc=(r_2**0.5)*as_dict['RSCALE']
		#print r_pc
		if len(annulus)==2:
			if r_pc>annulus[0] and r_pc<=annulus[1]:
				inmass.append(log10(m[i]*as_dict['ZMBAR']))
				infile.write("%f\n" %( log10(m[i]*as_dict['ZMBAR']) ))
			elif r_pc>annulus[1]:
				outmass.append(log10(m[i]*as_dict['ZMBAR']))
				outfile.write("%f\n" %( log10(m[i]*as_dict['ZMBAR']) ))
		elif len(annulus)==1:
    			if r_pc<=annulus[0]:
				inmass.append(log10(m[i]*as_dict['ZMBAR']))
			elif r_pc>annulus[0]:
				outmass.append(log10(m[i]*as_dict['ZMBAR']))

	print 'done in out'
	for i in range(len(inmass)):
		if inmass[i]>1.:
			print inmass[i]
	infile.close()
	outfile.close()
	return inmass, outmass


def histogram_equal_bins(m, binwidth, mrange, normrange):
	""""""
	
	nbin = (mrange[1]-mrange[0])/float(binwidth)
	print binwidth, mrange, len(m)

	histdata = zeros((nbin, 3))
	for i in range(int(nbin)):
		ni = 0
		for j in range(len(m)):
			if m[j]>= mrange[0]+i*binwidth and m[j]< mrange[0]+(i+1)*binwidth:
				ni += 1
				#print j
		#print ni
		masspoint = mrange[0]+(i+0.5)*binwidth
		histdata[i,0], histdata[i,1], histdata[i,2] = float(masspoint), float(ni), float(ni)**0.5
	
	norm = find_norm(m, binwidth, normrange[0], normrange[1])
	for i in range(len(histdata)):
		histdata[i,1] = histdata[i,1]/norm
		histdata[i,2] = histdata[i,2]/norm

	#print histdata
	return histdata

def find_norm(m, binwidth, m1, m2):
	"""Finds the normalization for equal bin histograms: Normalizes so that the area within 
	mass range (m1, m2) is 1
	"""
	n = 0
	for dummy in range(len(m)):
		if m[dummy]>=m1 and m[dummy]<m2:
			n += 1
	norm = n*binwidth
	return norm


def compare_mf_obs_data(obsfile, out3file, snapno, binwidth, filestring):
	snap = read_nbody6_snapshot_raw(out3file, snapno)
	as_dict = get_as_dict(snap)
	t_snap = "%.2f" %(get_time(snap))

	obsnormedfile = filestring+'_MF_obs.snapno_'+str(snapno)+'.t_'+str(t_snap)+'.dat'
	simnormedfile_in = filestring+'_MF_sim_in.snapno_'+str(snapno)+'.t_'+str(t_snap)+'.dat'
	simnormedfile_out = filestring+'_MF_sim_out.snapno_'+str(snapno)+'.t_'+str(t_snap)+'.dat'
	inmassfile = filestring+'_masses_sim_in.snapno_'+str(snapno)+'.t_'+str(t_snap)+'.dat'
	outmassfile = filestring+'_masses_sim_out.snapno_'+str(snapno)+'.t_'+str(t_snap)+'.dat'

	m = get_mass_array(snap)
	xs = get_xs_array(snap)
	proj = (0,1)
	annulus = (0.19, 0.35)
	#annulus = (0.195,)
	mrange = (-1., 2.)
	#normrange = (0., 1.7)
	normrange = (0.3, 1.8)
	inmass, outmass = mass_in_out(as_dict, m, xs, proj, annulus, inmassfile, outmassfile)
	inhistdata = histogram_equal_bins(inmass, binwidth, mrange, normrange)
	outhistdata = histogram_equal_bins(outmass, binwidth, mrange, normrange)

	obsdata = loadtxt(obsfile)
	norm = sum(obsdata[3:,1])*0.1
	print 'norm', norm
	obshistdata = zeros((len(obsdata), 6))
	for i in range(len(obsdata)):
		obsm = obsdata[i,0]
		obsn = obsdata[i,1]/norm
		obsdnp = obsdata[i,2]/norm
		obsdnm = obsdata[i,3]/norm
		pobsdn = obsdnp + obsn
		mobsdn = obsn - obsdnm
		obshistdata[i,0], obshistdata[i,1], obshistdata[i,2], obshistdata[i,3],  obshistdata[i,4], obshistdata[i,5] = obsm, obsn, obsdnp, obsdnm, mobsdn, pobsdn

	f1=open(obsnormedfile, 'w')
	f1.write("#1.log10(m/msun) 2.p 3.dp+ 4.dp- 5.p-dp- 6.p+dp+\n")
	for i in range(len(obshistdata)):
		for j in range(len(obshistdata[i])):
			f1.write("%f " %(obshistdata[i,j]))
		f1.write("\n")
	f1.close()
	f2=open(simnormedfile_in, 'w')
	f2.write("#1.log10(m/msun) 2.p 3.dp\n")
	f3=open(simnormedfile_out, 'w')
	f3.write("#1.log10(m/msun) 2.p 3.dp\n")
	for i in range(len(inhistdata)):
		for j in range(len(inhistdata[i])):
			f2.write("%f " %(inhistdata[i,j]))
			f3.write("%f " %(outhistdata[i,j]))
		f2.write("\n")
		f3.write("\n")

	f1.close()
	f2.close()
	f3.close()
	
	return obshistdata, inhistdata, outhistdata


def choose_mf_part(obsfile, simfile):
	sim_array, obs_array = [], []
	obs = loadtxt(obsfile)
	sim = loadtxt(simfile)
	for i in range(1,len(obs)):
		try:
			for j in range(len(sim)):
				if obs[i,0] == sim[j,0]:
					sim_array.append(sim[j,1])
					obs_array.append(obs[i,1])
					print obs[i,0], sim[j,0], obs[i,1], sim[j,1]
					raise StopIteration()
		except StopIteration:
			pass
	
	return obs_array, sim_array


#########################################################################################################

def get_seek(infile):
	"""takes file fort.83 and finds the positions for the snapshots from where the 
	corresponding time and the number of stars are to be read.  """

	positions=[]
	indata = open(infile, 'r')
	indata.seek(0)
	try:
		while indata:
			line = indata.readline().split()
			if len(line)==2:
				#print 'success'
				positions.append([indata.tell(), int(line[0]), float(line[1])])
				for i in range(int(line[0])):
					indata.readline()
			elif len(line)==0:
				raise StopIteration()
			else:
				print 'bad bad algorithm'
	except StopIteration:
		pass
	return positions



def get_xs_L_T(snap, positions, infile):
	"""snapfile: OUT3
	positions: f.tell() pointers for the start of a snap in fort.83
	infile: fort.83
	outfile: printing stuff
	"""
	xs = get_xs_array(snap)
	names = get_star_names(snap)
	masses = get_mass_array(snap)
	time = get_time(snap)
	as_dict = get_as_dict(snap)
	rc_snap = as_dict['RC']  #in code units
	snapdict = {}
	for i in range(len(names)):
		snapdict[names[i]] = {'X': xs[i,0],  #in code units
				'Y': xs[i,1],
				'Z': xs[i,2],
				'M': masses[i]
				}

	indata = open(infile, 'r')
	pos_l_t = {}
	#dtype = [('name', int), ('kw', int), ('r/rc', float), ('mass1', float), ('mass2', float), ('logL', float), ('lograd', float), ('logT', float), ('x/rc', float), ('y/rc', float), ('z/rc', float), ('2Dr/rc', float), ('3Dr/rc', float)]

	for i in range(len(positions)):
		#pos_l_t = {}
		#dummy = 0
		#print positions[i]
		if (float(positions[i][2])<float(time) and float(positions[i+1][2])>float(time)) or float(positions[i][2]==float(time)):
			#print i, float(positions[i][2]), float(time), float(positions[i+1][2])
			indata.seek(positions[i][0])
			for j in range(positions[i][1]):
			#for j in range(100):
				line = indata.readline()
				#print line
				line_array = line.split()
				id = int(line_array[0])
				#print id
				if snapdict.has_key(id):
					pos_l_t[id] = {'kw': int(line.split()[1]),
						#'r': float((line.split()[2]))*as_dict['RSCALE'],
						#'r/rc': float((line.split()[2]))/float(rc_snap),
						'mass1': float((line_array[3])),
						'mass2': float(snapdict[id]['M'])*as_dict['ZMBAR'],
						'logL': float((line_array[4])),
						'lograd': float((line_array[5])),
						'logT': float((line_array[6])),
						'logRstar': float((line_array[5])),
						#'x/rc': float(snapdict[id]['X'])/float(rc_snap),
						'x': float(snapdict[id]['X'])*as_dict['RSCALE'], 
						#'y/rc': float(snapdict[id]['Y'])/float(rc_snap),
						'y': float(snapdict[id]['Y'])*as_dict['RSCALE'],
						#'z/rc': float(snapdict[id]['Z'])/float(rc_snap),
						'z': float(snapdict[id]['Z'])*as_dict['RSCALE'],
						'2Dr_xy': ((float(snapdict[id]['X'])**2. + float(snapdict[id]['Y'])**2.)**0.5) * as_dict['RSCALE'],
						'2Dr_xz': ((float(snapdict[id]['X'])**2. + float(snapdict[id]['Z'])**2.)**0.5) * as_dict['RSCALE'],
						'2Dr_zy': ((float(snapdict[id]['Z'])**2. + float(snapdict[id]['Y'])**2.)**0.5) * as_dict['RSCALE'],
						'3Dr': ((float(snapdict[id]['X'])**2. + float(snapdict[id]['Y'])**2. + float(snapdict[id]['Z'])**2.)**0.5) * as_dict['RSCALE']
						}
					


					#print pos_l_t[id]
					#if pos_l_t[id]['3Dr/rc']<1.:
					#	dummy += 1
					#	print dummy
				else:
					print 'problem id', id
	return pos_l_t


def sorted_l_r(data_dict, s):
	"""data containing L, T, 2Dr
	s1: tuple containing the projection coordinates
	e.g., to suppress z s=(x,y)"""
	
	dtype = [('id', int), ('kw', int), ('2Dr', float), ('3Dr', float), ('L', float)]
	a = []
	for i in data_dict.keys():
		r_proj_square = 0
		for j in range(len(s)):
			r_proj_square = r_proj_square + (data_dict[i][s[j]])**2.
		r_proj = r_proj_square**0.5
		a.append( (int(i), data_dict[i]['kw'], r_proj, data_dict[i]['3Dr'], 10.**data_dict[i]['logL']) )
	
	a=array(a, dtype=dtype)
	a=sort(a, order='2Dr')

	return a

def get_r_hl(a):  
	"""takes the structured array obtained from sorted_l_r and calculates the 
	projected radius at half light  """
	L_r_dist = []
	L_tot = sum(a[:]['L'])
	try:
		sum_L = 0
		for i in range(len(a)):
			sum_L += a[i]['L']
			if sum_L >= 0.5*L_tot:
				rs = (a[i-1]['2Dr'], a[i]['2Dr'])
				Ls = (sum(a[:i+1]['L']), sum(a[:i]['L'])) 
				raise StopIteration()
	except StopIteration:
		rc_obs = linear_interpolate(Ls[0], Ls[1], rs[0], rs[1], L_tot/2.)
		#print rs, Ls, L_tot, L_tot/2., i, sum_L
		pass
	
	return rc_obs

def linear_interpolate(x1, x2, y1, y2, x):
	"""interpolates linearly for y value
	xs: tuple with (x1, x2)
	ys: tuple with (y1, y2)"""

	y = y1 + (x - x1)*(y2 - y1)/(x2 - x1)

	return y

def avg_r_hl(no_of_snaps, start_snap, snapname, filename):
	"""positions: create array using get_seek(fort.83)
	no. of snaps over which the averaging will take place
	start snap: the starting timestep from where n_snaps snaps will be taken
	snapname: file containing the snapshots
	"""
	positions = get_seek(filename)
	rc = []
	for i in range(start_snap,start_snap+no_of_snaps,1):
		snap=read_nbody6_snapshot_raw(snapname, i)
		time = get_time(snap)
		print time
		pos_l_t = get_xs_L_T(snap, positions, filename)
		rc_tmp = []
		a = sorted_l_r(pos_l_t, ('x', 'y'))
		rc_tmp.append(get_r_hl(a))
		a = sorted_l_r(pos_l_t, ('x', 'z'))
		rc_tmp.append(get_r_hl(a))
		a = sorted_l_r(pos_l_t, ('z', 'y'))
		rc_tmp.append(get_r_hl(a))
		#print rc_tmp
		rc.append( [mean(rc_tmp), std(rc_tmp)] )
	rc = reshape(rc, (no_of_snaps, 2))
	rc_obs = [mean(rc[:,0]), std(rc[:,0])]

	return rc_obs, rc


def hrdiag(infile, outfile, time):
	"""takes fort.83 of Nbody6 run and creates a file with luminosity and T_eff at a given time t"""

	infile = open(infile, 'r')
	infile.seek(0)
	line = infile.readline()
	snapinfo = {'nstar': float(line.split()[0]),
			'time': float(line.split()[1])
			}

	for i in range(snapinfo['nstar']):
		readline




	print snapinfo



def mass_loss(infile, star_id):
	indata = open(infile, 'r')
	indata.seek(0)
	finished = 0
	time = []
	m = []
	L = []
	while finished==0:
		line = indata.readline()
		if len(line) == 0:
			finished = 1
		elif len(line.split()) == 2:
			time.append( float(line.split()[1]) )
		elif line.split()[0] == star_id:
			m.append( float(line.split()[3]) )
			L.append( 10.**(float(line.split()[4])) )
	print len(time), len(m), len(L)

	return time, m, L



def create_3D_data(i, positions, dtype):
	snap=read_nbody6_snapshot_raw('OUT3', i)
    	pos_dict=get_xs_L_T(snap, positions,'fort.83')
	print 'came here'
    	a=[]
    	for j in pos_dict.keys():
        	a.append( (int(j), pos_dict[j]['x'], pos_dict[j]['y'], pos_dict[j]['z'], pos_dict[j]['mass2'], 10**pos_dict[j]['logRstar'], pos_dict[j]['logL'], 10**pos_dict[j]['logT'], pos_dict[j]['kw']) )
    	a=array(a, dtype=dtype)
    	sorted_a = sort(a, order='id')

	return sorted_a


def create_in_out_logequal_bins_mf(filename, snapno, binrange, binwidth, normrange, annulus):
	#extract simulation data
	raw_snap = read_nbody6_snapshot_raw(filename, snapno)
	as_dict = get_as_dict(raw_snap)
	m = get_star_masses(raw_snap)
	xs = get_xs_array(raw_snap)
	rc = as_dict['RC']
	t = as_dict['TTOT']*as_dict['TSCALE']
	print 't', t, 'rc', rc
	outhist, inhist = {}, {}

	#now separate inside and outside mass for 3D in 2rc
	filestring = 'MF_t'+str(t)+'_rc'+str(rc)
	outfilestring3D = filestring+'out_3D'
	infilestring3D = filestring+'in_3D'

	massin3D, massout3D = mass_in_out(as_dict, m, xs, (0,1,2), annulus, infilestring3D, outfilestring3D)
	inhist3Dfile=infilestring3D+'_hist.dat'
	outhist3Dfile=outfilestring3D+'_hist.dat'
	simhistdata3D_in = create_simhistdata(massin3D, (-1., 2.), 0.1, (0., 1.8), inhist3Dfile)
	simhistdata3D_out = create_simhistdata(massout3D, (-1., 2.), 0.1, (0., 1.8), outhist3Dfile)

	#now separate inside and outside mass for 2D annulus between rc and 2rc => 0.195 and 0.35
	outfilestring2D = filestring+'out_2D'
	infilestring2D = filestring+'in_2D'
	massin2D, massout2D = mass_in_out(as_dict, m, xs, (0,1), annulus, infilestring2D, outfilestring2D)
	inhist2Dfile=infilestring2D+'_hist.dat'
	outhist2Dfile=outfilestring2D+'_hist.dat'
	simhistdata2D_in = create_simhistdata(massin2D, (-1., 2.), 0.1, (0., 1.8), inhist2Dfile)
	simhistdata2D_out = create_simhistdata(massout2D, (-1., 2.), 0.1, (0., 1.8), outhist2Dfile)





def create_simhistdata(m_array, binrange, binwidth, normrange, outputfile):
	f1 = open(outputfile, 'w')
	f1.write("#1.log10(m/msun) 2.normalized dn/dlogm 3.normalized d(dn/dlogm)\n")
	nbins = int((binrange[1]-binrange[0])/binwidth)
	#ms = log10(m_array)
	norm = find_norm(m_array, binwidth, normrange[0], normrange[1])
	dndlogm, mbins = histogram(m_array, bins=nbins, range=binrange, normed=False)
	simhistdata = zeros((len(dndlogm), 3))
	mmids = zeros(len(mbins)-1)
	for i in range(1,len(mbins)):
		mmids[i-1] = (mbins[i-1]+mbins[i])/2.
	if len(mmids) == len(dndlogm):
		for i in range(len(dndlogm)):
			temp_m, temp_dndlogm, temp_ddndlogm = mmids[i], dndlogm[i]/norm, (dndlogm[i])**0.5/norm
			simhistdata[i,0], simhistdata[i,1], simhistdata[i,2] = temp_m, temp_dndlogm, temp_ddndlogm
			f1.write("%f %f %f\n" %(temp_m, temp_dndlogm, temp_ddndlogm))
	else:
		print "error in getting histogram, check dimensions"
	
	return simhistdata	



def in_out_logequal_bins_onesnap(snapno,filename, annulus, data_mrange, norm_mrange, nbin):
	"""input:
		no of snaps to take for averaging
		starting snap
		filename containing the snaps
		no of members wanted in each bin: bins are divided with equal number of stars in each bin
	output:
		plots the mass function within 2rc of the cluster
		also prints out the masses and the time etc. in two different files
			r<2rc and r>2rc stellar masses
		insidefile: MF data inside r=2rc
		outsidefile: MF data outside r=2rc
	 """

	#for i in range(len(m)):
	#	m[i] = log10(m[i]*as_dict['ZMBAR'])
	for i in range(len(data_mrange)):
		data_mrange[i] = log10(data_mrange[i])
	for i in range(len(norm_mrange)):
		norm_mrange[i] = log10(norm_mrange[i])
	

	#divide mass in inside and outside
	filestring = 'MF_t'+str(t)+'_rc'+str(rc)

	#3D
	infilestring = filestring+'in_3D'
	outfilestring = filestring+'out_3D'
	massin_3D, massout_3D = mass_in_out(as_dict, m, xs, (0,1,2), annulus, infilestring, outfilestring)
	
	#norm = find_norm(m, binwidth, range)
	outhist['3D'] = histogram_equal_bins(massout_3D, nbin, data_mrange, norm_mrange)
	inhist['3D'] = histogram_equal_bins(massin_3D, nbin, data_mrange, norm_mrange)

	#print outhist['3D'], inhist['3D']

	#2D
	#massin_2D, massout_2D = [], []
	histout, histin = {}, {}
	#infilestring = inmassfilestring+'_2D'
	#outfilestring = outmassfilestring+'_2D'

	for k in (0,1,2):
		proj = (k%3, (k+1)%3)
		infilestring = filestring+'in_2D_'+str(proj[0])+str(proj[1])
		outfilestring = filestring+'out_2D'+str(proj[0])+str(proj[1])
		massin_2D, massout_2D = mass_in_out(as_dict, m, xs, proj, annulus, infilestring, outfilestring)
	
		histout[k] = histogram_equal_bins(massout_2D, nbin, data_mrange, norm_mrange)
		histin[k] = histogram_equal_bins(massin_2D, nbin, data_mrange, norm_mrange)
		print len(histout[k]), k, proj
	
	outhist['2D'] = zeros((len(histout[0]), 3))
	inhist['2D'] = zeros((len(histout[0]), 3))

	for i in range(len(histout[0])):
		mpoint = histout[0][i,0]
		n_ave = (histout[0][i,1] + histout[1][i,1] + histout[2][i,1])/3.
		n_poiss = (histout[0][i,2] + histout[1][i,2] + histout[2][i,2])/3.
		n_proj_err = (max(histout[0][i,1], histout[1][i,1], histout[2][i,1]) - min(histout[0][i,1], histout[1][i,1], histout[2][i,1]) )/2.
		tot_err_n = n_poiss+n_proj_err

		outhist['2D'][i,0], outhist['2D'][i,1], outhist['2D'][i,2] = mpoint, n_ave, tot_err_n

		mpoint = histin[0][i,0]
		n_ave = (histin[0][i,1] + histin[1][i,1] + histin[2][i,1])/3.
		n_poiss = (histin[0][i,2] + histin[1][i,2] + histin[2][i,2])/3.
		n_proj_err = (max(histin[0][i,1], histin[1][i,1], histin[2][i,1]) - min(histin[0][i,1], histin[1][i,1], histin[2][i,1]) )/2.
		tot_err_n = n_poiss+n_proj_err

		inhist['2D'][i,0], inhist['2D'][i,1], inhist['2D'][i,2] = mpoint, n_ave, tot_err_n
	
	inhistfile = filestring+'hist_in.dat'
	outhistfile = filestring+'hist_out.dat'
	inhistfile = open(inhistfile, 'w')
	outhistfile = open(outhistfile, 'w')
	for i in range(nbin):
		for j in range(3):
			inhistfile.write("%f " %( inhist['3D'][i,j], ))
			outhistfile.write("%f " %( outhist['3D'][i,j], ))
		for j in range(1,3):
			inhistfile.write("%f " %(inhist['2D'][i,j]))
			outhistfile.write("%f " %(outhist['2D'][i,j]))
		inhistfile.write("\n")
		outhistfile.write("\n")

	inhistfile.close()
	outhistfile.close()
	return inhist, outhist


