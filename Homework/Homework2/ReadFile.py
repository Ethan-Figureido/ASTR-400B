# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 11:21:02 2025

@author: ethan
"""

#imports
import numpy as np
import astropy.units as u

def Read(filename):
  #This function takes a filename (text file) as an input, opens and reads the file.
  #The file data is then extracted and returned as three variables: time,
  #total particles, and a data variable that includes data for positions, velocities
  #masses, and particle types for each particle. 

  file = open(filename,'r') #open the file

  line1 = file.readline() #read first line
  label_1,value_1 = line1.split() #split line along spaces. First variable is the label (time) the second is the value
  time = float(value_1)*u.Myr #store in units of Myr

  line2 = file.readline() #read second line
  label_2,value_2 = line2.split() #repeat as above, this time for the total number of particles
  total_part = float(value_2) #store total number of particles (no astropy unit for this as far as I know)
  
  file.close() #close the file
  
  #retreive remainder of file data
  data = np.genfromtxt(filename,dtype=None,names=True,skip_header=3)
  
  #return time, total particles, and data
  return time, total_part, data
  