# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 11:49:22 2025

@author: ethan
"""

#imports
import numpy as np
from ReadFile import Read
import astropy.units as u


def ParticleInfo(filename, particle_type, particle_number):
    """This function will compute and return position, velocity, and mass data for a 
       user selected particle.
       Inputs: filename is the name of the data file the information
                   is being retrieved from
               particle_type is the type of particle being considered. It
                   accepts the following values: 1 = dark matter, 2 = disk stars
                   3 = bulge stars
                particle_number is the index of the particle selected for data retreival

       Outputs: velocity (astropy units km/s) is the calculated magnitude of the particle's velocity
                   measured from coordinate system at center of Milky Way
                distance (astropy units kpc) is the calculated magnitude of the distance of the particle
                    from the MW's center of mass
                mass (astropy units solar masses) is the mass of the particle
    """
    
    #extract file data
    time,total_part,data = Read(filename)
    
    #create index to separate for user inputted particle type
    index = np.where(data['type'] == particle_type)
    
    #split data into lists and assign units to each
    #apply index so only user inputted particle type are included in lists
    #x,y,z are cartesian coordinate positions, vx,vy,vz are the velocity components
    #m is the mass data in 10**10 solar masses
    x = data['x'][index]*u.kpc
    y = data['y'][index]*u.kpc
    z = data['z'][index]*u.kpc
    vx = data['vx'][index]*(u.km/u.s)
    vy = data['vy'][index]*(u.km/u.s)
    vz = data['vz'][index]*(u.km/u.s)
    m = data['m'][index]*u.M_sun*10**(10)
    
    #calculate magnitude of distance 
    distance = np.sqrt(x[particle_number-1]**2 + y[particle_number-1]**2 + z[particle_number-1]**2)
    
    #calculate magnitude of velocity
    velocity = np.sqrt(vx[particle_number-1]**2 + vy[particle_number-1]**2 + vz[particle_number-1]**2)

    #round velocity and distance values to 3 decimal places
    distance = np.around(distance,3)
    velocity = np.around(velocity,3)
    
    #retrieve the mass
    mass = m[particle_number-1]
    
    #return distance, velocity, and mass data for particle
    return distance, velocity, mass

#distance, velocity, mass = ParticleInfo('MW_000.txt',2,100)

#print the final values
#print('Distance:',distance)
#print('Velocity:', velocity)
#print('Mass:', mass)

#print the distance value converted to lightyears
#print('Distance in lightyears:',distance.to(u.lyr))

