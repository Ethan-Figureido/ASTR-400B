

# Homework 6 Template
# G. Besla & R. Li




# import modules
import numpy as np
import astropy.units as u
from astropy.constants import G

# import plotting modules
import matplotlib.pyplot as plt
import matplotlib

# my modules
from ReadFile import Read
# Step 1: modify CenterOfMass so that COM_P now takes a parameter specifying 
# by how much to decrease RMAX instead of a factor of 2
from CenterOfMassSingle import CenterOfMass




def OrbitCOM(galaxy, start, end, n):
    """Function that loops over all the desired snapshots to compute the COM pos and vel as a function of time.
    
    Inputs: galaxy is the name of the galaxy i.e. "MW"
    
            start is the number of the first snapshot to be read in
            
            end is the number of the last snapshot to be read in
            
            n is an integer indicating the intervals over which 
                this function will return the COM
          
    Outputs: Orbit_galaxyname.txt. A text file that stores the 
                pos, time, and vel values of the COM of the input galaxy.
    """
    
    # compose the filename for output
    fileout = f"Orbit_{galaxy}.txt"
    
    # set tolerance and VolDec for calculating COM_P in CenterOfMass
    # for M33 that is stripped more, use different values for VolDec
    volDec = 2
    volDec_M33 = 4
    delta = 0.1
    
    # generate the snapshot id sequence 
    # it is always a good idea to also check if the input is eligible (not required)
    snap_ids = np.arange(start, end, step=n)
    
    #check if array is empty
    if len(snap_ids) == 0:
        print("snap_ids empty. Please input valid start and end values.")
        return
    
    # initialize the array for orbital info: t, x, y, z, vx, vy, vz of COM
    orbit = np.zeros([len(snap_ids),7])
    
    # a for loop 
    for i, snap_id in enumerate(snap_ids): # loop over files
        # compose the data filename (be careful about the folder)
        # add a string of the filenumber to the value "000"
        ilbl = '000' + str(snap_id)
        #remove all but last 3 digits
        ilbl = ilbl[-3:]
        filename = "%s_"%(galaxy) + ilbl +'.txt'
        
        # Initialize an instance of CenterOfMass class, using disk particles
        gal = CenterOfMass(filename, 2)
        
        # Store the COM pos and vel. Remember that now COM_P required VolDec
        #check for M33 case:
        if galaxy == "M33":
            com_pos = gal.COM_P(volDec_M33, delta)
            com_vel = gal.COM_V(com_pos[0], com_pos[1], com_pos[2])
        else:
            com_pos = gal.COM_P(volDec, delta)
            com_vel = gal.COM_V(com_pos[0], com_pos[1], com_pos[2])
    
        #store t, pos, and vel values. divide t by 1000 to get units of Gyr
        orbit[i][0] = (gal.time/1000).value #store time without units
        orbit[i][1:4] = com_pos.value #store pos without units
        orbit[i][4:] = com_vel.value #store vel without units
        
        # print snap_id to see the progress
        print(snap_id)
        
    # write the data to a file
    # we do this because we don't want to have to repeat this process 
    # this code should only have to be called once per galaxy.
    np.savetxt(fileout, orbit, fmt = "%11.3f"*7, comments='#',
               header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                      .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))
    #return
    return


# Recover the orbits and generate the COM files for each galaxy
# read in 800 snapshots in intervals of n=5
# Note: This might take a little while - test your code with a smaller number of snapshots first! 

#OrbitCOM("MW", 0, 800, 5)
#OrbitCOM("M31", 0, 800, 5)
#OrbitCOM("M33", 0, 800, 5)

# Read in the data files for the orbits of each galaxy that you just created
# headers:  t, x, y, z, vx, vy, vz
# using np.genfromtxt
#MW_COM = np.genfromtxt('Orbit_MW.txt',dtype=None,names=True)
#M33_COM = np.genfromtxt('Orbit_M33.txt',dtype=None,names=True)
#M31_COM = np.genfromtxt('Orbit_M31.txt',dtype=None,names=True)

# function to compute the magnitude of the difference between two vectors 
# You can use this function to return both the relative position and relative velocity for two 
# galaxies over the entire orbit  

def VectorMag(galdata_1, galdata_2):
    """This function computes the difference between
    two vectors and returns the magnitude of that vector.
    Computes this for every snapshot in the data
    
    Inputs: galdata_1 is COM data for the first galaxy
            
            galdata_2 is COM data for the second galaxy
            
        
    Outputs: vel_vectormag is an array of magnitudes of velocity vector difference of
                galaxy 1 and 2 
                
             pos_vectormag is an array of magnitude of the position vector difference
                of galaxy 1 and 2
    """
    #define position and velocity data for galaxy 1
    x_1 = galdata_1['x']
    y_1 = galdata_1['y']
    z_1 = galdata_1['z']
    vx_1 = galdata_1['vx']
    vy_1 = galdata_1['vy']
    vz_1 = galdata_1['vz']
    
    #define position and velocity data for galaxy 2
    x_2 = galdata_2['x']
    y_2 = galdata_2['y']
    z_2 = galdata_2['z']
    vx_2 = galdata_2['vx']
    vy_2 = galdata_2['vy']
    vz_2 = galdata_2['vz']
    
    #compute vector differences for pos and vel
    x_f = x_1 - x_2
    y_f = y_1 - y_2
    z_f = z_1 - z_2
    vx_f = vx_1 - vx_2
    vy_f = vy_1 - vy_2
    vz_f = vz_1 - vz_2
    
    #compute magnitudes
    pos_vectormag = np.sqrt(x_f**2 + y_f**2 + z_f**2)
    vel_vectormag = np.sqrt(vx_f**2 + vy_f**2 + vz_f**2)
      
    return pos_vectormag, vel_vectormag


# ALL OF THE BELOW CODE IS RELATED TO THE PLOTTING AND QUESTION ANSWERING OF
# THE HOMEWORK


#extract time array. it is the same for each body so we use MW for all
#time = MW_COM['t'] 

# Determine the magnitude of the relative position and velocities 

# of MW and M31

#MW_M31_p, MW_M31_v = VectorMag(MW_COM, M31_COM)

#print(f"Differences in position for MW and M31: {MW_M31_p}")
#print(f"Differences in velocity for MW adn M31: {MW_M31_v}")

# of M33 and M31 

#M33_M31_p, M33_M31_v = VectorMag(M33_COM, M31_COM)

#print(f"Differences in position for M33 and M31: {M33_M31_p}")
#print(f"Differences in velocity for M33 and M31: {M33_M31_v}")

# Plot the Orbit of the galaxies 
#################################

#plt.plot(time, MW_M31_p, color='red',label='MW-M31')
#plt.xlabel('Time (Gyr)')
#plt.ylabel('Separation (kpc)')
#plt.legend()

#plt.plot(time, M33_M31_p, color='blue',label='M33-M31')
#plt.xlabel('Time (Gyr)')
#plt.ylabel('Separation (kpc)')
#plt.legend()


# Plot the orbital velocities of the galaxies 
#################################

#plt.plot(time, MW_M31_v, color='red',label='MW-M31')
#plt.xlabel('Time (Gyr)')
#plt.ylabel('Relative Velocity (km/s)')
#plt.legend()

#plt.plot(time, M33_M31_v, color='blue',label='M33-M31')
#plt.xlabel('Time (Gyr)')
#plt.ylabel('Relative Velocity (km/s)')
#plt.legend()


# Questions
#################################

#1. 
#The troughs of the wavelike separataion graph indicate close encounters
#for the two bodies. The MW-M31 graph has 3 troughs, indicating three 
#close encounters for the MW and M31 in the future.

#2.
#The points of smallest separation are peaks of relative velocity in the 
#relative velocity graph. This makes sense, as bodies in orbit experience a 
#decrease in velocity as separation increases. We see that demonstrated on the 
#separation and relative velocity graphs, with peaks in the relative velocity 
#graphs matching with troughs in the separation graphs.

#3.
#The galaxy merger has been completed when the separation between M31 and the MW
#has decreased to zero. We see that on the graph at a time of ~6.2 Gyr. Once M31
#and the MW have merged, M33's orbit continues to spiral into the new merged
#galaxy, indicated by lowering amplitude of separation in the M33-M31 graph.

#4. BONUS

#apocenter 1 ~ 115 kpc
#apocenter 2 ~ 90
#period ~ 2 Gyr

#Decay Rate: alpha = (apocenter1 - apocenter2)/period

#alpha = 12.5 kpc/Gyr

#Assuming this rate is constant, and that M33 is initially at a distance of 
#75 kpc, we can find the time it will take for M33 to merge with the combined 
#MW+M31 remnant.

#time = 75/alpha = 75/12.5 = 6 Gyr

#So, it will take ~6 Gyr for M33 to merge with the combined MW+M31 remnant.



















