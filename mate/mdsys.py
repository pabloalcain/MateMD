"""
MateMD class heriarchy
"""

import random as R
import numpy as np
from itertools import combinations_with_replacement as cwr
from itertools import combinations

class MDSys(object):
    """Base class for molecular dynamics system. Only one type of
    particles.

    """
    def __init__(self, box, particles, interaction, dt):
        """Constructor. Needs a box to link to and particles. It takes the
        dimension from the dimension of the box. If the dimensions of the
        particles is different to that of the box, raises a value error.

        Parameters
        ----------

        box : Box instance
            Box in which the simulation will be carried
        particles: Particles instance
            Particles in the simulation
        interaction: Dictionary
            Dictionary with the interaction between types.
        dt: float
            Timestep

        """
        if particles.dim != box.dim:
            raise ValueError(("The dimension of the particles must be "
                              "equal to that of the box"))

        self.dt = float(dt)
        self.dim = box.dim
        self.box = box
        self.particles = particles
        self.interaction = interaction
        self.potential = 0

        self.part_list = {}
        self.pair_list = {}
        self.precondition()


    def precondition(self):
        """Get the indices through which we should run, and check if anything
        is missing. Populate the dictionaries part_list (that has the
        list of particle indices for each type) and pair_list (the
        list of particle pair indices for each particle pair).

        """

        self.part_list = {}
        for j in self.particles.names:
            self.part_list[j] = \
            [i for i, _x in enumerate(self.particles.name) if _x == j]

        #Create all pair_list arrays
        self.pair_list = {}
        for pair_name in cwr(self.particles.names, 2):
            self.pair_list[tuple(sorted(pair_name))] = []

        #Populate pair_list
        for pair in combinations(xrange(self.particles.N), 2):
            name_1 = self.particles.name[pair[0]]
            name_2 = self.particles.name[pair[1]]
            pair_name = tuple(sorted((name_1, name_2)))
            self.pair_list[pair_name].append((pair[0], pair[1]))



    def run(self, n):
        """Run the simulation

        Parameters
        ----------

        n : integer
            Steps to be integrated

        """

        for _ in xrange(n):
            print _
            self.first_step()
            self.forces()
            self.final_step()
            self.boundary()

    def first_step(self):
        """First step of Verlet integration
        """

        #Loop over all types
        for i in self.part_list.keys():
            mass = self.particles.mass[i]
            #Loop over all particles of this type
            for index in self.part_list[i]:
                dv = 0.5 * self.dt / mass * self.particles.f[index]
                self.particles.v[index] += dv
                dx = self.dt * self.particles.v[index] 
                self.particles.x[index] += dx
                
    def final_step(self):
        """Final step of Verlet integration
        """

        #Loop over all types
        for i in self.part_list.keys():
            mass = self.particles.mass[i]
            #Loop over all particles of this type
            for index in self.part_list[i]:
                dv = 0.5 * self.dt / mass * self.particles.f[index]
                self.particles.v[index] += dv


    def forces(self):
        """Calculate all forces
        """

        #Loop over all types
        length = [0] * self.dim
        self.potential = 0
        self.particles.f.fill(0.0)
        for i in range(self.dim):
            length[i] = self.box.shape[i][1] - self.box.shape[i][0]
        for i in self.pair_list.keys():
            F = self.interaction[i].force
            V = self.interaction[i].energy
            #Loop over all particles of this type
            for index in self.pair_list[i]:
                _ii = index[0]
                _jj = index[1]
                r1 = self.particles.x[_ii]
                r2 = self.particles.x[_jj]
                r = r1 - r2
                for k in range(self.dim):
                    while r[k] > 0.5 * length[k]:
                        r[k] -= length[k]
                    while r[k] < - 0.5 * length[k]:
                        r[k] += length[k] 
                rabs = np.linalg.norm(r)
                if rabs > self.interaction[i].rcut: continue
                self.particles.f[_ii] += F(rabs)/rabs * r
                self.particles.f[_jj] -= F(rabs)/rabs * r
                self.potential += V(rabs)

    def boundary(self):
        length = [0] * self.dim
        for i in range(self.dim):
            length[i] = self.box.shape[i][1] - self.box.shape[i][0]
        for i in range(self.dim):
            for j in range(self.particles.N):
                while self.particles.x[j, i] > self.box.shape[i][1]:
                    self.particles.x[j, i] -= 0.5 * length[i]
                while self.particles.x[j, i] < self.box.shape[i][0]:
                    self.particles.x[j, i] += 0.5 * length[i]

    def kinetic_energy(self):
        ke = 0.0
        for i in self.part_list.keys():
            mass = self.particles.mass[i]
            v2 = 0
            for index in self.part_list[i]:
                v2 += np.dot(self.particles.v[index], self.particles.v[index])
            ke += 0.5 * mass * v2
        return ke

        

class Particles(object):

    """All the information about particles in the molecular dynamics
    system. This particle will hold a dictionary for the values of
    each name-depending variable. For example, the dictionary mass
    will always exist.

    """

    def __init__(self, box, dim=3, N=0):
        """Constructor. Creates particles and sets the default value for the
        dimension of the particles in it.

        Parameters
        ----------

        box : Box instance
            Box in which to create the particles

        dim : integer
            Dimension of the particles
            Default is ``3``

        N : integer
            Particles to create, with random position and zero
            velocity, named 'A'
            Default is ``0``
        """

        self.dim = dim
        self.box = box
        self.x = None
        self.v = None
        self.f = None
        self.name = []

        self.names = []
        self.mass = {}
        self.values = {'mass': self.mass}

        self.N = 0
        if N < 0:
            raise ValueError("Tried to create {0} particles".format(N))
        if N > 0:
            self.add_particles(N)

    def add_particles(self, N, x=None, v=None, name='A'):

        """Add particles to the system.

        Parameters
        ----------
        N : integer
            Number of particles to add to the system
        x : None or array_like
            Position of the particles. If ``None``, particles
            positions will be random inside the box.
            Default is ``None``
        v : None or array_like
            Velocities of the particles. If ``None``, particles
            velocities will be set to zero.
            Default is ``None``
        part : string
            Name of particles to be added.
            Default is ``'A'``

        Returns
        -------
        None

        """
        dim = self.dim

        if x == None:
            x = np.zeros((N, dim))
            for _i in xrange(N):
                for _j in xrange(dim):
                    lo = self.box.shape[_j][0]
                    hi = self.box.shape[_j][1]
                    x[_i, _j] = R.uniform(lo, hi)
        elif np.shape(x) != (N, dim):
            raise ValueError(("Shape of x is {0}, but needed {1} to fill all "
                              "the values").format(np.shape(x), (N, dim)))

        if v == None:
            v = np.zeros((N, dim))
        elif np.shape(v) != (N, dim):
            raise ValueError(("Shape of v is {0}, but needed {1} to fill all "
                              "the values").format(np.shape(v), (N, dim)))

        if self.x == None:
            self.x = x
        else:
            self.x = np.vstack((self.x, x))

        if self.v == None:
            self.v = v
        else:
            self.v = np.vstack((self.v, v))

        if self.f == None:
            self.f = np.zeros((N, dim))
        else:
            self.f = np.vstack((self.v, np.zeros((N, dim))))

        if name not in self.names:
            self.names.append(name)
            self.mass[name] = 1.0
            
        for _ in xrange(N):
            self.name.append(name)
            
        self.N += N

class Box(object):
    """Box class, defining the boundary of the simulation

    """
    def __init__(self, shape):
        """Constructor. Sets the dimensions of the box.

        Parameters
        ----------
        shape: list
            List of shape (dim, 2), that writes the shape of the box. First
            column is low values, second is high values. One row per dimension.

        """
        dim = np.shape(shape)[0]

        if dim > 3 or dim < 1:
            raise ValueError("Tried to set {0} dimensions".format(dim))

        self.dim = dim
        self.shape = shape

class Potential(object):
    """Potential class. Stores the tables of potential and force.

    """
    def __init__(self, pot, rcut, npoints=5000):
        """Constructor: Builds the potential and force LUT.

        Parameters
        ----------
        pot : function
            Potential energy between the particles

        rcut : float
            Cut off radius of the potential

        npoints : integer
            Points with which to interpolate
            Default is ``5000``

        """
        self.rcut = rcut
        self.r = np.linspace(0, self.rcut, npoints + 1)[1:]
        self.V = pot(self.r)
        self.F = - np.diff(self.V)/np.diff(self.r)
        self.F.resize(npoints)

    def energy(self, r):
        """ Energy at a distance r, from the LUT.

        Parameters
        ----------
        r : float
            Distance between the particles

        Returns
        -------
        V : float
            Potential energy

        """
        return np.interp(r, self.r, self.V)

    def force(self, r):
        """ Force at a distance r, from the LUT.

        Parameters
        ----------
        r : float
            Distance between the particles

        Returns
        -------
        F : float
            Force

        """
        return np.interp(r, self.r, self.F)
