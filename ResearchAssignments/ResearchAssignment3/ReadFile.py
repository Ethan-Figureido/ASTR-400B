# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 11:21:02 2025

@author: ethan
"""

#imports
import numpy as np
import astropy.units as u

def Read(filename):
  """ This function will open and read out a given data file and returns the time,
      total number of particles as variables, and returns particle info as a data array
      
      Input: filename is the name of the given data file
      
      Outputs: time (astropy units Myr) is the time
               total_part is the total number of particles in the data file
               data is an array including the following data for each particle in
                 the data file: particle type, mass, x, y, z, vx, vy, and vz
  """

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
  


#time, part, data = Read('MW_000.txt')
