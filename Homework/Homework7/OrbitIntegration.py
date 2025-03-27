
# # Homework 7 Template
# 
# Rixin Li & G . Besla
# 
# Make edits where instructed - look for "****", which indicates where you need to 
# add code. 




# import necessary modules
# numpy provides powerful multi-dimensional arrays to hold and manipulate data
import numpy as np
# matplotlib provides powerful functions for plotting figures
import matplotlib.pyplot as plt
# astropy provides unit system and constants for astronomical calculations
import astropy.units as u
import astropy.constants as const
# import Latex module so we can display the results with symbols
from IPython.display import Latex

# **** import CenterOfMass to determine the COM pos/vel of M33
from CenterOfMass2 import CenterOfMass

# **** import the GalaxyMass to determine the mass of M31 for each component
import GalaxyMass

# M33AnalyticOrbit




class M33AnalyticOrbit:
    """ Calculate the analytical orbit of M33 around M31 """
    
    def __init__(self,filename): # **** add input
        """Class initialization method
        
        filename is the name of the file that will
            store the integrated orbit.
        """
        ### get the gravitational constant (the value is 4.498502151575286e-06)
        self.G = const.G.to(u.kpc**3/u.Msun/u.Gyr**2).value
        
        ### **** store the output file name
        self.filename = filename
        
        ### get the current pos/vel of M33 
        # **** create an instance of the  CenterOfMass class for M33 
        M33 = CenterOfMass('M33_000.txt', 2)
        # **** store the position VECTOR of the M33 COM (.value to get rid of units)
        M33_pos = M33.COM_P(volDec=4, delta=0.1).value #volDec value from Homework 6
        # **** store the velocity VECTOR of the M33 COM (.value to get rid of units)
        M33_vel = M33.COM_V(M33_pos[0], M33_pos[1], M33_pos[2]).value  #HERE
        
        ### get the current pos/vel of M31 
        # **** create an instance of the  CenterOfMass class for M31 
        M31 = CenterOfMass('M31_000.txt', 2)
        # **** store the position VECTOR of the M31 COM (.value to get rid of units)
        M31_pos = M31.COM_P(volDec=2, delta=0.1).value #volDec value from Homework 6
        # **** store the velocity VECTOR of the M31 COM (.value to get rid of units)
        M31_vel = M31.COM_V(M31_pos[0], M31_pos[1], M31_pos[2]).value
        
        ### store the DIFFERENCE between the vectors posM33 - posM31
        # **** create two VECTORs self.r0 and self.v0 and have them be the
        # relative position and velocity VECTORS of M33
        self.r = M33_pos - M31_pos
        self.v = M33_vel - M31_vel
        
        
        ### get the mass of each component in M31 
        
        ### disk
        # **** self.rdisk = scale length (no units)
        self.rdisk = 5 #kpc

        # **** self.Mdisk set with ComponentMass function. Remember to *1e12 to get the right units. Use the right ptype
        self.Mdisk = 0.12*1e12 #solar masses from assignment 3
        
        ### bulge
        # **** self.rbulge = set scale length (no units)
        self.rbulge = 1 #kpc

        # **** self.Mbulge  set with ComponentMass function. Remember to *1e12 to get the right units Use the right ptype
        self.Mbulge = 0.019*1e12 #solar masses
        
        # Halo
        # **** self.rhalo = set scale length from HW5 (no units)
        self.rhalo = 62 #kpc #from assignment 5
    
        # **** self.Mhalo set with ComponentMass function. Remember to *1e12 to get the right units. Use the right ptype
        self.Mhalo = 1.921*1e12 #solar masses
    
    
    def HernquistAccel(self, M, r_a, r): # it is easiest if you take as an input the position VECTOR 
        """This method computes the acceleration vector from a Hernquist
        potential using the following equation:
            
        a = -G*M/(rmag*(r_a + rmag)**2) *r
        
        (rmag is not an input, and is calculated via the input postion vector)
        
        Inputs:
            M is the total halo or bulge mass of the considered galaxy
            (units of solar masses)
            
            r_a is the corresponding scale length
            (units of kpc)
            
            r is the position vector of the galaxy, inputted as an array
            (units of kpc)
            
        Output:
            a is the acceleration vector from a Hernquist potential.
                It is an array with a[0] being the x component, a[1] the y,
                and a[2] the z. (units of kpc/Gyr**2)
        """
        
        ### **** Store the magnitude of the position vector
        rmag = np.sqrt(r[0]**2 + r[1]**2 + r[2]**2)
        
        ### *** Store the Acceleration
        Hern = -self.G*M/(rmag*(r_a + rmag)**2)*r #follow the formula in the HW instructions
        # NOTE: we want an acceleration VECTOR so you need to make sure that in the Hernquist equation you 
        # use  -G*M/(rmag *(ra + rmag)**2) * r --> where the last r is a VECTOR 
        
        return Hern
    
    
    
    def MiyamotoNagaiAccel(self, M, r):# it is easiest if you take as an input a position VECTOR  r 
        """This method computes the acceleration vector from a Miyamoto-
        Nagai 1975 profile using the following equation:
            
            acc = -G*M/(R**2 + B**2)**(1.5)
            
        Where:
            R = (r[0]**2 + r[1]**2)**(1/2)
            
            B = r_d + ()**(1/2)
        
        PARAMETERS
        ----------
            M: 'float'
                The mass of the galaxy in units of M_sun
            r: 'array'
                Array of the x, y, and z components of the position vector.
                Elements are floats with units of kpc
            
        OUTPUTS
        -------
            MNacc: 'array'
                Array of x, y, and z components of the Miyamoto-Nagai 
                acceleration.
        """
        r_d = self.rdisk
        z_d = self.rdisk/5.0
        
        ### Acceleration **** follow the formula in the HW instructions
        # AGAIN note that we want a VECTOR to be returned  (see Hernquist instructions)
        # this can be tricky given that the z component is different than in the x or y directions. 
        # we can deal with this by multiplying the whole thing by an extra array that accounts for the 
        # differences in the z direction:
        # multiply the whle thing by :   np.array([1,1,ZSTUFF]) 
        # where ZSTUFF are the terms associated with the z direction
        
        R = np.sqrt(r[0]**2 + r[1]**2)
        B = r_d + np.sqrt(r[2]**2 + z_d**2)
        
        #the Right Hand Side of eqn.4 from homework sheet
        RHS = np.array([1,1,B/np.sqrt(r[2]**2 + z_d**2)])
        
        #calculate acceleration
        MNacc = -self.G*M/(R**2 + B**2)**(1.5) *r*RHS
    
        return MNacc
        # the np.array allows for a different value for the z component of the acceleration
     
    
    def M31Accel(self, r): # input should include the position vector, r
        """This method sums all acceleration vectors from each galaxy component
        and returns it as a 3D vector of the total acceleration
        
        PARAMETERS
        ----------
            r: 'array'
                3D position vector of M31
        
        OUTPUTS
        -------
            acc: 'array'
                3D, summed acceleration vector for M31
        """

        ### Call the previous functions for the halo, bulge and disk
        # **** these functions will take as inputs variable we defined in the initialization of the class like 
        # self.rdisk etc. 
        
        #Determine acceleration for each component of M31
        disk_acc = self.MiyamotoNagaiAccel(self.Mdisk, r)
        bulge_acc = self.HernquistAccel(self.Mbulge, self.rbulge, r) 
        halo_acc = self.HernquistAccel(self.Mhalo, self.rhalo, r)
        
        #sum acceleration of each galaxy component
        acc = disk_acc + bulge_acc + halo_acc 
            
            # return the SUM of the output of the acceleration functions - this will return a VECTOR 
        return acc
    
    
    
    def LeapFrog(self, dt, r, v): # take as input r and v, which are VECTORS. Assume it is ONE vector at a time
        """This method integrates a function using a variant of the
        'Leap Frog' integration scheme
        
        PARAMETERS
        ----------
            dt: 'float'
                The time interval for integration. In our case, we want
                dt to be positive.
            r: 'array'
                The starting 3D position vector for the M33 COM position
                relative to M31.
            v: 'array'
                The starting velocity vector for M33 relative to M31.
        
        OUTPUTS
        -------
            rnew: 'array'
                3D position vector advanced one full time step from r.
                
            vnew: 'array'
                3D velocity vector advanced one full time step from v.
        """
        
        # predict the position at the next half timestep
        rhalf = r + v*dt/2
        
        #find acceleration at half timestep
        ahalf = self.M31Accel(rhalf)
        
        # predict the final velocity at the next timestep using the acceleration field at the rhalf position 
        vnew = v + ahalf*dt
        
        # predict the final position using the average of the current velocity and the final velocity
        # this accounts for the fact that we don't know how the speed changes from the current timestep to the 
        # next, so we approximate it using the average expected speed over the time interval dt. 
        rnew = rhalf + vnew*dt/2
        
        return rnew, vnew
    
    
    
    def OrbitIntegration(self, t0, dt, tmax):
         """This method computes the future orbit of M33 for 10 Gyr 
         in the future.
         
         PARAMETERS
         ----------
             t0: 'float'
                 The starting time for integration
             dt: 'float'
                 The time interval of integration. In this case we want it
                 to be positive.
             tmax: 'float'
                 The final time for integration
         
         OUTPUTS
         -------
         """

         # initialize the time to the input starting time
         t = t0
        
         # initialize an empty array of size :  rows int(tmax/dt)+2  , columns 7
         orbit = np.zeros((int(tmax/dt)+2, 7))
        
         # initialize the first row of the orbit
         orbit[0] = t0, *tuple(self.r), *tuple(self.v)
         # this above is equivalent to 
         # orbit[0] = t0, self.r0[0], self.r0[1], self.r0[2], self.v0[0], self.v0[1], self.v0[2]
        
        
         # initialize a counter for the orbit.  
         i = 1 # since we already set the 0th values, we start the counter at 1
        
         # start the integration (advancing in time steps and computing LeapFrog at each step)
         while (t <= tmax):  # as long as t has not exceeded the maximal time 
            
             # **** advance the time by one timestep, dt
             t += dt 
             # **** store the new time in the first column of the ith row
             orbit[i,0] = t
            
             # ***** advance the position and velocity using the LeapFrog scheme
             # remember that LeapFrog returns a position vector and a velocity vector  
             # as an example, if a function returns three vectors you would call the function and store 
             # the variable like:     a,b,c = function(input)
             r_step, v_step = self.LeapFrog(dt, orbit[i-1, 1:4], orbit[i-1, 4:])
         
    
             # ****  store the new position vector into the columns with indexes 1,2,3 of the ith row of orbit
             # TIP:  if you want columns 5-7 of the Nth row of an array called A, you would write : 
             # A[n, 5:8] 
             # where the syntax is row n, start at column 5 and end BEFORE column 8
             orbit[i, 1:4] = r_step
             orbit[i, 4:] = v_step
            
             # ****  store the new position vector into the columns with indexes 1,2,3 of the ith row of orbit
            
            
             # **** update counter i , where i is keeping track of the number of rows (i.e. the number of time steps)
             i += 1
        
        
         # write the data to a file
         np.savetxt(self.filename, orbit, fmt = "%11.3f"*7, comments='#', 
                   header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                   .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))
        
         # there is no return function



#%%

dt = 0.1 #Gyr  
tmax = 10 #Gyr
t0 = 0 #Gyr

M33_AO = M33AnalyticOrbit('M33_Analytic_Orbit.txt')

#M33_AO.OrbitIntegration(t0, dt, tmax) #integrate orbit COMPLETE

# ANALYSIS
filename = 'M33_Analytic_Orbit.txt'
#retrieve data from file
data = np.genfromtxt(filename,dtype=None,names=True)

pos_data = np.array([data['x'], data['y'], data['z']])    #select position data
vel_data = np.array([data['vx'], data['vy'], data['vz']]) #select velocity data
integration_time = data['t']                              #select time data

#Calculate velocity and position magnitudes
pos_mag = np.sqrt(pos_data[0]**2 + pos_data[1]**2 + pos_data[2]**2)
vel_mag = np.sqrt(vel_data[0]**2 + vel_data[1]**2 + vel_data[2]**2)

#get orbit COM data from homework 6
M33_COM = np.genfromtxt('Orbit_M33.txt',dtype=None,names=True)
M31_COM = np.genfromtxt('Orbit_M31.txt',dtype=None,names=True)

#import to compute COM vector magnitude differences
from OrbitCOM import VectorMag

M33_M31_p, M33_M31_v = VectorMag(M33_COM, M31_COM) #compute vector magnitude differences
time = M33_COM['t'] #extract time array for homework 6 orbit array

#plotting
#position
#plt.plot(time, M33_M31_p, color='red',label='M33-M31(OrbitCOM)')
#plt.plot(integration_time, pos_mag, color='blue', label='(Orbit Integration)')
#plt.xlabel('Time (Gyr)')
#plt.ylabel('Separation (kpc)')
#plt.legend()

#velocity
#plt.plot(integration_time, vel_mag, color='blue', label='(Orbit Integration)')
#plt.plot(time, M33_M31_v, color='red',label='M33-M31 (OrbitCOM)')
#plt.xlabel('Time (Gyr)')
#plt.ylabel('Relative Velocity (km/s)')
#plt.legend()


#QUESTION 2

'''Both the velocity and position plots are in approximate agreement until
t~2 Gyr. From here, the HW7 position plot takes a larger orbit, with its
pericenter at ~470kpc. This is vastly different from the HW6 position plot, 
which show many smaller orbits of decreasing orbital distance.
The velocity plots show a similar trend. after 2 Gyr the HW7 plot shows a single,
large orbit over the 10Gyr simulation time, whereas the HW6 graph gives a
sinusoidal plot, indicating many orbits over the 10Gyr simulation time.
'''

#QUESTION 3
'''One thing missing from the HW7 system is the MW. The neglected effects of
the MW could explain the differences in the data.
'''

#QUESTION 4
'''If you wanted to include the MW in this code, you can just add its 
acceleration to the sum of accelerations in the M31Accel method. Find the 
component accelerations for the MW using the same Hernquist and Miyamoto-Nagai 
acceleration methods as used for M31, and sum them along with the M31 
acceleration values.
'''

