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
    #takes a filename, particle type and particle number as inputs
    #and returns the magnitude of the distance in kpc, magnitude of
    #the particle velocity in km/s, and the mass in solar masses.
    #particle number is the specific particle being considered, filename
    #is the name of the data file, and particle type is the type of particle
    #being considered (1,2, or 3)=(dark matter, disk stars, bulge stars).
    #distance and velocity are calculated as follows d = sqrt(x^2 + y^2 + z^2)
    #mass is extracted directly from the file
    
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

