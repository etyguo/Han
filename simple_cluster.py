#!/usr/bin/python

import numpy as np
import imf

from matplotlib import pyplot as plt

from amuse.units import nbody_system
from amuse.units import units
from amuse.community.hermite0.interface import Hermite
import logging

from spherically_symmetric_powerlaw import new_powerlaw_model, clump_radius, gasveldisp
#logging.basicConfig(level=logging.DEBUG)

smoothing_length = 0.0 | nbody_system.length ** 2

G = 6.63784e-11
cm_to_m = 1. / 100.
msun_to_g = 1.988e33
pc_to_m = 3.08567758128e16
gcm_to_msunpc = 5.03e-34 * (3.086e18)**2

def print_log(time, gravity, particles, total_energy_at_t0):
    kinetic_energy = gravity.kinetic_energy
    potential_energy = gravity.potential_energy
    total_energy_at_this_time = kinetic_energy + potential_energy
    print "time                    : " , time
    print "energy error            : " , (total_energy_at_this_time - total_energy_at_t0) / total_energy_at_t0


def simulate_small_cluster(
        number_of_stars = 1000.,
        end_time = 40 | nbody_system.time,
        number_of_workers = 1
    ):
    #convert_nbody = nbody_system.nbody_to_si(number_of_stars | units.MSun, 1. | units.parsec)
    particles = new_powerlaw_model(number_of_stars)
    particles.scale_to_standard()

    gravity = Hermite(number_of_workers = number_of_workers)
    gravity.parameters.epsilon_squared = 0.15 | nbody_system.length ** 2
    
    gravity.particles.add_particles(particles)

    from_model_to_gravity = particles.new_channel_to(gravity.particles)
    from_gravity_to_model = gravity.particles.new_channel_to(particles)

    time = 0.0 | end_time.unit
    total_energy_at_t0 = gravity.kinetic_energy + gravity.potential_energy

    positions_at_different_times = []
    positions_at_different_times.append(particles.position)
    times = []
    times.append(time)
    velocities_at_different_times = []
    velocities_at_different_times.append(particles.velocity)

    print "evolving the model until t = " + str(end_time)
    while time < end_time:
        time +=  end_time / 3.0

        gravity.evolve_model(time)
        from_gravity_to_model.copy()

        positions_at_different_times.append(particles.position)
        times.append(time)
        print_log(time, gravity, particles, total_energy_at_t0)

    gravity.stop()
    unit_var = units.m, units.s, units.kg
    print 'units_parsec', units.parsec
    return times, positions_at_different_times, velocities_at_different_times, unit_var, particles


def adjust_spines(ax,spines, ticks):
    for loc, spine in ax.spines.iteritems():
        if loc in spines:
            spine.set_position(('outward',10)) # outward by 10 points
            spine.set_smart_bounds(True)
        else:
            spine.set_color('none') # don't draw spine

    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
        ax.yaxis.set_ticks(ticks)
    else:
        ax.yaxis.set_ticks([])

    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
        ax.xaxis.set_ticks(ticks)

    else:
        ax.xaxis.set_ticks([])


def plot_positions(times, positions_at_different_times):
    plt.ion()
    figure = plt.figure()
    plot_matrix_size = np.ceil(np.sqrt(len(positions_at_different_times))).astype(int)
    number_of_rows = len(positions_at_different_times) / plot_matrix_size
    figure.subplots_adjust(wspace = 0.15, hspace = 0.15)

    for index, (time, positions) in enumerate(zip(times, positions_at_different_times)):
        subplot = figure.add_subplot(plot_matrix_size, plot_matrix_size, index + 1)

        subplot.scatter(
            positions[...,0].value_in(nbody_system.length),
            positions[...,1].value_in(nbody_system.length),
            s = 1,
            edgecolors = 'red',
            facecolors = 'red'
        )

        subplot.set_xlim(-4.0,4.0)
        subplot.set_ylim(-4.0,4.0)
        subplot.set_aspect(1.)

        title = 'time = {0:.2f}'.format(time.value_in(nbody_system.time))

        subplot.set_title(title)#, fontsize=12)
        spines = []
        if index % plot_matrix_size == 0:
            spines.append('left')

        if index >= ((number_of_rows - 1)*plot_matrix_size):
            spines.append('bottom')


        adjust_spines(subplot, spines,np.arange(-4.0,4.1, 1.0))

        if index % plot_matrix_size == 0:
            subplot.set_ylabel('y')

        if index >= ((number_of_rows - 1)*plot_matrix_size):
            subplot.set_xlabel('x')
            
    plt.savefig('simple_cluster_powerlaw.png', clobber=True)


def write_init_file(m, pos, vel, fname):
    """
    Writes the initial masses, positions, and velocities of a cluster to a file.
    """
    f = open(fname, 'w')

    for i in np.arange(np.size(m)):
        f.write(str(m[i]) + '\t' + str(pos[i,0]) + '\t' + str(pos[i,1]) + '\t' + str(pos[i,2]) + '\t' + str(vel[i,0]) + '\t' + '\t' + str(vel[i,1]) + '\t' + str(vel[i,2]) + '\n')

    f.close()

    return

def gasvel_potential(sigma_cl, mstars, eps, k_rho, k_P=1., f_g=1., phi_Pc=4./3., phi_Pbar=1.32, R=None):        
    #R = 0.05 * np.power(A / (k_P**2 * eps**2 * phi_Pc * phi_Pbar), 1./4.) * np.power(np.sum(mstars) / 30., 1./2.) * sigma_cl**(-1./2.) * pc_to_m
    if R is None:
        R = clump_radius(mstars, sigma_cl, eps, k_rho, k_P, f_g, phi_Pc, phi_Pbar) * pc_to_m
    
    #return G / (2. * (5. - 2. * k_rho)) / R * (np.sum(mstars) * msun_to_g / 1000.)**2
    return - 3. / 5. * (1. - k_rho / 3.) / (1. - 2. * k_rho / 5.) * G * (np.sum(mstars) * msun_to_g / 1000.)**2. / R
    

def gasvel_ke(mstars, sigma_cl, eps, phi_b, k_rho, k_P=1., f_g=1., phi_Pc=4./3., phi_Pbar=1.32, R=None):
    """
    Calculates the kinetic energy expected for a cluster with stars moving using a mass-averaged velocity dispersion.

    mstars: array_like
        Masses of the stars in solar masses
    sigma: float
        Mass-averaged velocity dispersion in km/s
    """
    sigma = gasveldisp(mstars, sigma_cl, eps, phi_b, k_rho, k_P, f_g, phi_Pc, phi_Pbar)# * 7./6.

    print(sigma)
    
    #R = 0.05 * np.power(A / (k_P**2 * eps**2 * phi_Pc * phi_Pbar), 1./4.) * np.power(np.sum(mstars) / 30., 1./2.) * sigma_cl**(-1./2.)
    if R is None:
        R = clump_radius(mstars, sigma_cl, eps, k_rho, k_P, f_g, phi_Pc, phi_Pbar)
    
    rho_s = (3. - k_rho) / (4. * np.pi) / R**3.

    #print(sigma)
    #print(rho_s)
    #print(R_s)
    
    return 3./2. * (np.sum(mstars) * msun_to_g / 1000.) * (sigma * 1000.)**2.
    #return 2. * np.pi * (rho_s * np.sum(mstars) * msun_to_g / 1000) * (np.sqrt(3.) * sigma * 1000)**2. * R**3. / (5. - 2. * k_rho)

def radii(pos):
    return np.sqrt((pos[ : , 0].value_in(units.parsec))**2 + (pos[ : , 1].value_in(units.parsec))**2 +(pos[ : , 2].value_in(units.parsec))**2)

def convert_to_nbody(mass, cluster):
    """
    Converts a cluster to nbody units using the prescription of Heggie & Mathieu 1986
    """
    Gun = G | units.m**3 * units.kg**-1. * units.s**-2.
    
    mtot = np.sum(mass)

    ke = cluster.kinetic_energy()

    pe = cluster.potential_energy()

    etot = ke + pe

    if etot > 0 | units.m**2 * units.kg * units.s**-2:
        etot *= -1.

    l = -Gun * mtot**2. / (4 * etot)
    t = Gun * mtot**(5./2.) / (-4. * etot)**(3./2.) 

    return mass / mtot, cluster.position / l, cluster.velocity / l * t

def initialize_mcluster(m_cluster, eps, phib, sigma_cl=0.1, seed=42, k_rho=1.5, massfunc='kroupa', single_m=1.,  plotdir='./', velreduc=1.0):
    """
    Same as initialize_ncluster, but keeps total mass of cluster constant instead of N. See below.
    """
    np.random.seed(seed)

    epsrange = np.arange(0.2, 1.01, .01)

    rmax = np.zeros(np.size(eps))
    
    if massfunc == 'kroupa':
        masses = imf.make_cluster(m_cluster)
    elif massfunc is None:
        masses = np.zeros(m_cluster) + single_m

    plt.ion()
    
    for j in np.arange(np.size(phib)):
        plt.figure(1)
        plt.clf()
        plt.figure(10)
        plt.clf()
        plt.figure(11)
        plt.clf()
        plt.figure(12)
        plt.clf()

        anlpot = gasvel_potential(sigma_cl, masses, epsrange, k_rho)
        anlkin = gasvel_ke(masses, sigma_cl, epsrange, phib[j], k_rho)
        anlkinplot =  gasvel_ke(masses, sigma_cl, eps, phib[j], k_rho)

        for i in np.arange(np.size(eps)):
            cluster = new_powerlaw_model(masses, k_rho, eps[i], phib[j], sigma_cl, velreduc)

            write_init_file(masses, cluster.position.value_in(units.parsec), cluster.velocity.value_in(units.km / units.s), plotdir + 'seed' + str(seed) + '.eps' + str(eps[i]) + '.phib' + str(phib[j]) + '_fort.10')

            pot = cluster.potential_energy()
            kin = cluster.kinetic_energy()

            rmax[i] = np.max(radii(cluster.position))

            plt.figure(10)
            plt.scatter(eps[i], np.abs(pot.value_in(units.kg * units.m**2 * units.s**-2)), color='k')
            plt.figure(11)
            plt.scatter(eps[i], np.abs(kin.value_in(units.kg * units.m**2 * units.s**-2)), color='k')
            plt.figure(12)
            plt.scatter(eps[i], anlkinplot[i] / kin.value_in(units.kg * units.m**2 * units.s**-2))

            #print('Potential energy = ' + str(pot))
            #print('Kinetic energy = ' + str(kin))
            #print('Cluster Rs = ' + str(cluster.R_s[0]))

            plt.figure(1)
            plt.scatter(eps[i], np.abs(kin / pot), color='k')


        print('Clump radius - Rmax = ' + str(clump_radius(masses, sigma_cl, eps, k_rho, k_P=1., f_g=1., phi_Pc=4./3., phi_Pbar=1.32) - rmax))
        #print('Potential from analytic = ' + str(anlpot))
        #print('Kinetic from analytic = ' + str(anlkin))

        plt.figure(1)
        plt.plot(epsrange, np.abs(anlkin / anlpot), 'k-')
        plt.ylabel('|KE/GPE|')
        plt.xlabel(r'$\epsilon$')
        plt.savefig(plotdir + 'seed' + str(seed) + '.phib' + str(phib[j]) + '_virratio.png', dpi=200, clobber=True)

        plt.figure(10)
        plt.plot(epsrange, np.abs(anlpot), 'k-')
        plt.ylabel('GPE')
        plt.xlabel(r'$\epsilon$')
        plt.savefig(plotdir + 'seed' + str(seed) + '.phib' + str(phib[j]) + '_pot.png', dpi=200, clobber=True)
        
        plt.figure(11)
        plt.plot(epsrange, np.abs(anlkin), 'k-')
        plt.ylabel('KE')
        plt.xlabel(r'$\epsilon$')
        plt.savefig(plotdir + 'seed' + str(seed) + '.phib' + str(phib[j]) + '_kin.png', dpi=200, clobber=True)

        plt.figure(12)
        plt.ylabel('Analytic KE / KE')
        plt.xlabel(r'$\epsilon$')
        plt.savefig(plotdir + 'seed' + str(seed) + '.phib' + str(phib[j]) + '_anlkinratio.png', dpi=200, clobber=True)
        

    return cluster

def initialize_ncluster(n_cluster, eps, phib, sigma_cl=0.1, seed=42, k_rho=1.5, massfunc='kroupa', single_m=1.,  plotdir='./', velreduc=1.0, nbody=False):
    """
    Initializes a cluster given a number of stars, star formation efficiency, and megnetic potantial term following Tan & McKee's model of a gas clump. 
    """
    #Seed the RGN
    np.random.seed(seed)

    #Make an array of different star formation efficiencies for the same realization
    epsrange = np.arange(0.2, 1.01, .01)

    #Make an array for the farthest star from t he center of the cluster, R
    rmax = np.zeros(np.size(eps))

    #CHeck mass function
    if massfunc == 'kroupa':
        #If Kroupa, use Adam Ginsburg's code to sample the IMF
        #Note I modified his code to make N the inputer parameter as opposed to the total stellar mass.
        #You can still input the total stellar mass by using initialize_mcluster instead of this function
        masses = imf.make_ncluster(n_cluster)
    elif massfunc is None:
        masses = np.zeros(n_cluster) + single_m

    plt.ion()

    #Loop through different values for magnetic field term
    for j in np.arange(np.size(phib)):
        plt.figure(1)
        plt.clf()
        plt.figure(10)
        plt.clf()
        plt.figure(11)
        plt.clf()
        plt.figure(12)
        plt.clf()

        #Calculate the analytical energies given the cluster properties
        anlpot = gasvel_potential(sigma_cl, masses, epsrange, k_rho)
        anlkin = gasvel_ke(masses, sigma_cl, epsrange, phib[j], k_rho)
        anlkinplot =  gasvel_ke(masses, sigma_cl, eps, phib[j], k_rho)

        #Write the energies to a file
        qfile = open(plotdir + 'seed' + str(seed) + '.phib' + str(phib[j]) + '_virratios.dat', 'w')
        qfile.write('Epsilon\tQ\n')

        #Loop through different values for star formation efficiency
        for i in np.arange(np.size(eps)):
            #Get positions and velocities for masses
            cluster = new_powerlaw_model(masses, k_rho, eps[i], phib[j], sigma_cl, velreduc)

            #Check if output in Nbody units or not
            #Then write the fort.10 file
            if nbody is not True:
                write_init_file(masses, cluster.position.value_in(units.parsec), cluster.velocity.value_in(units.km / units.s), plotdir + 'seed' + str(seed) + '.eps' + str(eps[i]) + '.phib' + str(phib[j]) + '_fort.10')
            else:
                newm, newx, newv = convert_to_nbody(masses | units.MSun, cluster)
                write_init_file(newm, newx, newv, plotdir + 'seed' + str(seed) + '.eps' + str(eps[i]) + '.phib' + str(phib[j]) + '_fort.10')

            #Numerically calculate the energies
            pot = cluster.potential_energy()
            kin = cluster.kinetic_energy()

            #Write these to the energy file
            qfile.write(str(eps[i]) + '\t' + str(np.abs(kin/pot)) + '\n')

            #Find the radius of the cluster
            rmax[i] = np.max(radii(cluster.position))

            #Plot things
            #THis is mostly to check that the cluster is behaving as the analytic model predicts it should
            plt.figure(10)
            plt.scatter(eps[i], np.abs(pot.value_in(units.kg * units.m**2 * units.s**-2)), color='k')
            plt.figure(11)
            plt.scatter(eps[i], np.abs(kin.value_in(units.kg * units.m**2 * units.s**-2)), color='k')
            plt.figure(12)
            plt.scatter(eps[i], anlkinplot[i] / kin.value_in(units.kg * units.m**2 * units.s**-2))

            #print('Potential energy = ' + str(pot))
            #print('Kinetic energy = ' + str(kin))
            #print('Cluster Rs = ' + str(cluster.R_s[0]))

            plt.figure(1)
            plt.scatter(eps[i], np.abs(kin / pot), color='k')


            #print('Clump radius - Rmax = ' + str(clump_radius(masses, sigma_cl, eps, k_rho, k_P=1., f_g=1., phi_Pc=4./3., phi_Pbar=1.32) - rmax))
        #print('Potential from analytic = ' + str(anlpot))
        #print('Kinetic from analytic = ' + str(anlkin))

        #Plot more things to check
        plt.figure(1)
        plt.plot(epsrange, np.abs(anlkin / anlpot), 'k-')
        plt.ylabel('|$\mathcal{T} / \mathcal{W}$|')
        plt.xlabel(r'$\epsilon$')
        plt.savefig(plotdir + 'seed' + str(seed) + '.phib' + str(phib[j]) + '_virratio.png', dpi=200, clobber=True)

        plt.figure(10)
        plt.plot(epsrange, np.abs(anlpot), 'k-')
        plt.ylabel('GPE')
        plt.xlabel(r'$\epsilon$')
        plt.savefig(plotdir + 'seed' + str(seed) + '.phib' + str(phib[j]) + '_pot.png', dpi=200, clobber=True)
        
        plt.figure(11)
        plt.plot(epsrange, np.abs(anlkin), 'k-')
        plt.ylabel('KE')
        plt.xlabel(r'$\epsilon$')
        plt.savefig(plotdir + 'seed' + str(seed) + '.phib' + str(phib[j]) + '_kin.png', dpi=200, clobber=True)

        plt.figure(12)
        plt.ylabel('Analytic KE / KE')
        plt.xlabel(r'$\epsilon$')
        plt.savefig(plotdir + 'seed' + str(seed) + '.phib' + str(phib[j]) + '_anlkinratio.png', dpi=200, clobber=True)
        
    qfile.close()

    return cluster

#The single mass case
#For just one cluster
#cluster = initialize_ncluster(1000., np.arange(0.2, 1.1, .1), [2.8], sigma_cl=0.1, seed=100, massfunc=None, plotdir='gasvels/single_single/sigma0.1/astrophysical_input/', velreduc=1.0, nbody=False)
#For different random seeds
for i in np.arange(1, 12):
    cluster = initialize_ncluster(1000, np.arange(0.2, 1.1, .1), [2.8], sigma_cl=0.1, seed=i, massfunc=None, plotdir='gasvels/single_single/sigma0.1/astrophysical_input/', velreduc=1.0, nbody=False)

#Same for the Kroupa case
#Note that the N=2400 is so that the number of stars will be equal, and the total mases will be ~1000 Msuns
for i in np.arange(1, 12):
    cluster = initialize_ncluster(2400, np.arange(0.2, 1.1, .1), [2.8], sigma_cl=0.1, seed=i, massfunc='kroupa', plotdir='gasvels/kroupa_single/sigma0.1/const_n/astrophysical/', velreduc=1.0, nbody=False)
#Low and high epse; just one realization
#cluster = initialize_ncluster(2375, [0.25, 0.5, 0.75], [2.8], sigma_cl=0.1, seed=1, massfunc='kroupa', plotdir='gasvels/kroupa_single/sigma0.1/const_n/astrophysical/', velreduc=1.0, nbody=False)

#The Kroupa binary case (binary fraction = 0.5)
for i in np.arange(1, 12):
    cluster = initialize_ncluster(1200, np.arange(0.2, 1.1, .1), [2.8], sigma_cl=0.1, seed=i, massfunc='kroupa', plotdir='gasvels/kroupa_binary/binfrac0.5/sigma0.1/const_n/astrophysical/', velreduc=1.0, nbody=False)

