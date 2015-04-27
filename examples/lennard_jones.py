"""
Example of Lennard-Jones method
"""

import mdsys

#Set the box
box = mdsys.Box([[-10, 10], [-10, 10], [-10, 10]])

#Create 40 particles
part = mdsys.Particles(box, N = 40)

#Set interaction
lennard = lambda r: 4 * (r**(-12.0) - r**(-6.0))
lj = mdsys.Potential(lennard, 2.5)
interactions = {}
interactions[('A', 'A')] = lj

#Create system
system = mdsys.MDSys(box, part, interactions, 0.01)

system.run(1000)


