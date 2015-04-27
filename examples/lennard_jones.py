"""
Example of Lennard-Jones method
"""
import mdsys
import numpy as np
from itertools import product

#Set the box
L = [0, 5] 
box = mdsys.Box([L, L, L])

#Start particles instance
part = mdsys.Particles(box, dim = 3)

#Create 125 particles in a lattice
nL = 5
npart = nL**3
x = np.zeros((npart, 3))
for (i, j, k) in product(xrange(nL), xrange(nL), xrange(nL)):
    x[i * nL**2 + j * nL + k, :] = [i, j, k]
v = 0.5*np.random.rand(npart, 3)
part.add_particles(npart, x, v)
    
#Set interaction
lennard = lambda r: 4 * (r**(-12.0) - r**(-6.0))
lj = mdsys.Potential(lennard, 2.5, 5000)
interactions = {('A', 'A'): lj}

#Create system
system = mdsys.MDSys(box, part, interactions, 1e-4)

#Run
fp = open("thermo.dat", "w")
print>>fp, "#Potential Energy, Kinetic Energy, Total Energy"
for i in range(1000):
    print i
    system.run(1)
    system.particles.kinetic_energy()
    pe = system.particles.potential
    ke = system.particles.kinetic
    system.dump()
    print>>fp, pe, ke, pe + ke
