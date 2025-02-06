#imports
import numpy as np #to work with numpy arrays
from ReadFile import Read #to read out file data
import pandas as pd  #to create table

def ComponentMass(filename,p_type):
    """This function returns the total mass
    of a desired galaxy component.
    
    Inputs: filename is the name of the file where the data is 
                extracted from
            p_type is the particle type, which are as follows:
                (1,2,3) = (halo type, disk type, bulge type). It 
                defines the desired component.
        
    Output: total_mass (astropy units solMass) is the
                total mass of the desired galaxy component.
    """
    
    #read in file data
    time, total_particles, data = Read(filename)
    
    #create index to separate for user inputted particle type
    index = np.where(data['type'] == p_type)
    
    #retrieve mass data for all particles of requested particle type
    #mass data in units of M_sun
    m = data['m'][index]*10**(10) 
    
    total_mass = sum(m)  #sum all entries in mass list
    
    total_mass = np.round(total_mass,3) #round to 3 decimal places
    
    return total_mass
    
#all of the following code is for creating a table for the requested values
#and exporting it
#creating table of mass data for 3 input files

#calculating total mass of each galaxy
MW_total_mass = ComponentMass("MW_000.txt",1) + ComponentMass("MW_000.txt",2) + ComponentMass("MW_000.txt",3)
M31_total_mass = ComponentMass("M31_000.txt",1) + ComponentMass("M31_000.txt",2) + ComponentMass("M31_000.txt",3)
M33_total_mass = ComponentMass("M33_000.txt",1) + ComponentMass("M33_000.txt",2)

#total mass of the local group
LG_mass = MW_total_mass + M31_total_mass + M33_total_mass

#finding Baryon Fraction for each galaxy
MW_mass_fraction = (ComponentMass("MW_000.txt",2)+ComponentMass("MW_000.txt",3))/MW_total_mass
M31_mass_fraction = (ComponentMass("M31_000.txt",2)+ComponentMass("M31_000.txt",3))/M31_total_mass
M33_mass_fraction = (ComponentMass("M33_000.txt",2))/M33_total_mass

#creating dictionary for table
mass_table = {'Galaxy Name':['MW','M31','M33'],
              'Halo Mass (10^12 M_sun)':[round(ComponentMass("MW_000.txt",1)/(10**12),3),
                                         round(ComponentMass("M31_000.txt",1)/(10**12),3),
                                         round(ComponentMass("M33_000.txt",1)/(10**12),3)],
              'Disk Mass (10^12 M_sun)':[round(ComponentMass("MW_000.txt",2)/(10**12),3),
                                         round(ComponentMass("M31_000.txt",2)/(10**12),3),
                                         round(ComponentMass("M33_000.txt",2)/(10**12),3)],
              'Bulge Mass (10^12 M_sun)':[round(ComponentMass("MW_000.txt",3)/(10**12),3),
                                          round(ComponentMass("M31_000.txt",3)/(10**12),3),
                                          "none"],
              'Total Galaxy Mass (10^12 M_sun)':[round(MW_total_mass/(10**12),3), 
                                                 round(M31_total_mass/(10**12),3),
                                                 round(M33_total_mass/(10**12),3)],
              'Baryon Fraction':[round(MW_mass_fraction,3),
                                 round(M31_mass_fraction,3),
                                 round(M33_mass_fraction,3)],
              'Total Local Group Mass (10^12 M_sun)':[LG_mass/(10**12),'-','-']}

#add total mass and stellar mass fraction (disk+bulge/total) columns to this table.

mass_table = pd.DataFrame(mass_table)
#print(mass_table)


#mass_table.to_html('400B HW3 Mass Table.html')  #export table to html, which you can then save as a pdf
