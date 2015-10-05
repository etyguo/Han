"""
Generates a spherically symmetirc, power-law density distributed cluster of stars

This model contains a function used to create a spherically symmetric cluster of stars. The density of mass in this model follows the form:
rho = rho_s * (R_s / r)**(-k_rho)

Returns the cluster as an AMUSE cluster object. Note that while AMUSE is used in this version of the code, it really isn't necessary. The only it's used for is dealing with the units. If you want to change things to not use AMUSE, you just have to modify it to deal with the units manually when it returns the cluster.
"""

import numpy as np

from amuse.units import nbody_system
from amuse import datamodel
from amuse.units import units
__all__ = ["new_powerlaw_sphere", "new_powerlaw_model"]

gcm_to_msunpc = 5.03e-34 * (3.086e18)**2

#Note, equations governing the radius and velocity dispersion from Tan et al. 2013
def clump_radius(mstars, sigma_cl, eps, k_rho, k_P=1., f_g=1., phi_Pc=4./3., phi_Pbar=1.32):
    """
    Calculates the expected radius of the clump given the conditions in it.
    """
    A = (3. - k_rho) * (k_rho - 1.) * f_g
    
    #return 0.05 * np.power(A / (k_P**2 * eps**2 * phi_Pc * phi_Pbar), 1./4.) * np.power(np.sum(mstars) / 30., 1./2.) * sigma_cl**(-1./2.)
    return 0.071 * np.power(A / (k_P**2 * phi_Pc * phi_Pbar), 1./4.) * np.power(np.sum(mstars) / eps / 60., 1./2.) * sigma_cl**(-1./2.)
    #return np.sqrt(np.sum(mstars) / eps / (np.pi * sigma_cl *gcm_to_msunpc ))
    
def gasveldisp(masses, sigma_cl, eps, phi_b, k_rho, k_P=1., f_g=1., phi_Pc=4./3., phi_Pbar=1.32):
    """
    Calculates the mass-averaged velocity dispersion in the clump given the mass and properties.
    """
    A = (3. - k_rho) * (k_rho - 1.) * f_g
    
    return 1.91 * 2. * (3. - k_rho) / (8. - 3. * k_rho) * np.power((phi_Pc * phi_Pbar) / (A * k_P**2 * phi_b**4), 1./8.) * np.power(np.sum(masses) / eps / 60., 1./4.) * np.power(sigma_cl, 1./4.)


class MakePowerLawSphere(object):
    def __init__(self, masses, k_rho, eps, phi_B, sigma_cl=0.1, velreduc=1.0, convert_nbody=None,
            do_scale = False, random_state = None):
        #Parameters of model from Tan et al. 2013
        k_P = 1.
        f_g = 1.
        phi_Pc = 4./3.
        phi_Pbar = 1.32
        A = (3. - k_rho) * (k_rho - 1.) * f_g
        
        self.mass = masses
        self.number_of_particles = np.size(masses)
        self.velreduc = velreduc
        self.convert_nbody = convert_nbody
        self.R_s = clump_radius(masses, sigma_cl, eps, k_rho, k_P, f_g, phi_Pc, phi_Pbar)
        #Set density at R_s so that mass enclosed is 1
        self.rho_s = (3. - k_rho) / (4. * np.pi) / (self.R_s)**3.
        self.k_rho = k_rho
        self.do_scale = do_scale
        if not random_state == None:
            print "DO NOT USE RANDOM STATE"
        self.random_state = None

        self.sigma_s = gasveldisp(self.mass, sigma_cl, eps, phi_B, k_rho, k_P, f_g, phi_Pc, phi_Pbar) * 7./6.

        print(self.sigma_s)

    def calculate_radius(self):
        return np.power(((3. - self.k_rho) * np.random.uniform(0.0, 1.0, (self.number_of_particles,1))) / (4 * np.pi * self.rho_s * (self.R_s**self.k_rho)), (1. / (3. - self.k_rho)))

    def new_positions_spherical(self):
        r = self.calculate_radius()
        theta = np.arccos(np.random.uniform(-1.0,1.0, (self.number_of_particles,1)))
        phi = np.random.uniform(0.0, 2 * np.pi, (self.number_of_particles,1))
    
        return (r, theta, phi)

    def gas_velocities(self, radii):
        v = np.zeros((self.number_of_particles, 3))
        
        for i, r in enumerate(radii):
            for j in np.arange(3):
                v[i, j] = np.random.normal(0., self.sigma_s * np.power(r / self.R_s, (2. - self.k_rho) / 2.)) * self.velreduc
                #print(v[i, j])

        return (v[:,0], v[:,1], v[:,2])
    
    def new_velocities_spherical_coordinates(self, radii):
        velocity = np.sqrt((4 * np.pi * self.rho_s * self.R_s**self.k_rho) * radii**(2 - self.k_rho) / (3 - self.k_rho))
        theta = np.arccos(np.random.uniform(-1.0,1.0, (self.number_of_particles,1)))
        phi = np.random.uniform(0.0,2. * np.pi, (self.number_of_particles,1))
        return (velocity,theta,phi)

    def coordinates_from_spherical(self, radius, theta, phi):
        x = radius * np.sin( theta ) * np.cos( phi )
        y = radius * np.sin( theta ) * np.sin( phi )
        z = radius * np.cos( theta )
        return (x,y,z)

    def new_model(self):
        m = self.mass
        r, theta, phi = self.new_positions_spherical()
        position = np.hstack(self.coordinates_from_spherical(r, theta, phi))
        vx, vy, vz = self.gas_velocities(r)
        velocity = np.hstack((vx, vy, vz))

        return (m, position, velocity)

    @property
    def result(self):
        masses, positions, velocities = self.new_model()
        result = datamodel.Particles(self.number_of_particles)
        result.mass = units.MSun.new_quantity(masses)
        result.position = units.parsec.new_quantity(positions)
        result.velocity = (units.km / units.s).new_quantity(velocities)
        result.radius = 0. | units.parsec
        result.R_s = self.R_s
        result.rho_s = self.rho_s

        result.move_to_center()
        if self.do_scale:
            result.scale_to_standard()

        if not self.convert_nbody is None:
            result = datamodel.ParticlesWithUnitsConverted(result, self.convert_nbody.as_converter_from_si_to_generic())
            result = result.copy_to_memory()
        
        return result

def new_powerlaw_model(mcluster, *list_arguments, **keyword_arguments):
    uc = MakePowerLawSphere(mcluster, *list_arguments, **keyword_arguments)
    
    return uc.result

new_powerlaw_sphere = new_powerlaw_model
