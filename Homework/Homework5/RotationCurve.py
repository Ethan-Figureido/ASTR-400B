#imports
from ReadFile import Read
import numpy as np
from CenterOfMass import CenterOfMass
import astropy.units as u
import matplotlib.pyplot as plt



class MassProfile:
    """Class to determine the mass distribution of a given galaxy
    at a given simulation snapshot
    """
    
    def __init__(self, galaxy, snap):
        """A class to determine the mass distribution of a galaxy
        at a given snap number.
        
        Inputs: galaxy is the galaxy name input as a string
                    e.g. "MW", "M31"
                    
                snap is the snapshot number. Int 0,1,2,3,...
        """
        # add a string of the filenumber to the value "000"
        ilbl = '000' + str(snap)
        #remove all but last 3 digits
        ilbl = ilbl[-3:]
        self.filename = "%s_"%(galaxy) + ilbl +'.txt'
        
        #retrieve snapshot data
        self.time, self.total_p, self.data = Read(self.filename)
        
        #read in positoin and mass data.
        #apply units only to position at this stage.
        self.x = self.data['x']*u.kpc
        self.y = self.data['y']*u.kpc
        self.z = self.data['z']*u.kpc
        self.m = self.data['m']*(1e10)
        
        #store galaxy name as global property
        self.gname = galaxy
        
        
    def MassEnclosed(self, ptype, radii):
        """Method that computes the mass enclosed by a 
        given radius of the COM position for a specified galaxy
        and specified component of the galaxy
        
        Inputs: ptype is particle type and specifies 
                    galaxy component (int halo = 1, disk = 2, or bulge = 3)
                
                radii is an array of radii (magnitudes)
            
        Outputs: mass_enclosed is an array of masses enclosed by
                    the respective radii inputs (units of
                    1e12*solMass)
        """
        #initialize mass array
        mass_enclosed = np.zeros(len(radii))
        
        #applying units to radii
        radii = radii*u.kpc
        
        #determine COM position of galaxy using disk particles
        COM = CenterOfMass(self.filename, 2) 
        galaxy_com = COM.COM_P()  #COM position as array of [x,y,z]
   
        #select for specified particle type, and shift to COM reference frame
        index = np.where(self.data['type'] == ptype) 
        x_pos = self.x[index] - galaxy_com[0]
        y_pos = self.y[index] - galaxy_com[1]
        z_pos = self.z[index] - galaxy_com[2]
        particle_mass = self.m[index]
        
        #iterate over radius array
        for i in range(len(radii)):
            enclosing_radius = radii[i] #radius being considered
            
            #find distance of each particle from COM
            part_radius = np.sqrt(x_pos**2 + y_pos**2 + z_pos**2)
            
            #select all particles within enclosed radius
            index2 = np.where(enclosing_radius > part_radius)
            enclosed_particles = particle_mass[index2]
            
            #sum mass of all enclosed particles, and append to mass_enclosed array
            mass_enclosed[i] += sum(enclosed_particles)
        
        #apply mass units (solMass) to mass raray
        mass_enclosed = mass_enclosed*u.solMass
        
        #return array of enclosed mass values
        return mass_enclosed
    
    
    def MassEnclosedTotal(self, r):
        """This method computes the total enclosed mass
        (bulge + halo + disk) within a given set of radii from the 
        center of mass.
        
        Inputs: r is an array of radii in kpc
            
        Outputs: tot_mass_enclosed is an array of the total 
                    enclosed masses for each radii in r
                    (astropy units 1e12*solMass)
        """
        #determine halo, disk, and bulge masses for given radii
        #check for M33 input. M33 does not have a bulge, so it requires
        #a different computation
        if self.gname == "M33": 
            halo_mass = self.MassEnclosed(1,r)
            disk_mass = self.MassEnclosed(2,r)
            
            tot_mass_enclosed = halo_mass + disk_mass
            
            return tot_mass_enclosed
        
        #for MW and M31 
        else:
            halo_mass = self.MassEnclosed(1,r)
            disk_mass = self.MassEnclosed(2,r)
            bulge_mass = self.MassEnclosed(3,r)
        
            #compute total mass array
            tot_mass_enclosed = halo_mass + disk_mass + bulge_mass
        
            return tot_mass_enclosed
        
        
    def HernquistMass(self, radius, a, Mhalo):
        """This method computes the mass enclosed within a given 
        radius using the following theoretical profile:
            
            rho = (Mhalo*a)/(2*pi*radius*(radius + a)**3)
            
            M = Mhalo*radius**2/(a + radius)**2
            
        Inputs: radius is the radius enclosing the considered mass
                
                a is the scale factor
                
                Mhalo is the halo mass
            
        Outputs: M (astropy units solMass) is the halo mass
                    computed via the above profile
                    
                rho is the density profile computed via
                    the above profile 
        """
        
        #compute Mass profile
        M = Mhalo*(radius**2/(a+radius)**2)
        
        #compute density profile
        #rho = (Mhalo*a)/(2*np.pi*radius*(radius + a)**3)
        
        return M
    
    
    def CircularVelocity(self, ptype, r):
        """This method computes the circular speed for a specified 
        particly type using the mass enclosed at each radius, and
        returns an array of circular speeds in km/s
        
        v = (G*M/r)**(1/2)
        
        Inputs: ptype is the particle type
                    (int 1=halo, 2=disk, 3=bulge)
                
                r is an array of radii in units of kpc
            
        Outputs: v_array is an array of circular speeds for 
                    the input radii in km/s
        """
        #import gravitational constant and adjust units for our purposes
        from astropy import constants as const
        G = const.G.to(u.kpc*u.km**2/u.s**2/u.Msun)
        
        #extract array of enclosed masses
        M = self.MassEnclosed(ptype,r)
        
        #apply units to r array
        r = r*u.kpc
        
        v_array = np.sqrt(G*M/r)
            
        #return final circular speed array
        return v_array
    
    
    def CircularVelocityTotal(self,r):
        """This method computes the total circular velocity of all
        galaxy components.
        
        Input: r is an array of radii magnitudes with units kpc
        
        Output: vel is an array of circular speed values with units
                of km/s. It is the total circular velocity created by
                all galaxy components
        """
        #import gravitational constant and adjust units for our purposes
        from astropy import constants as const
        G = const.G.to(u.kpc*u.km**2/u.s**2/u.Msun)
        
        #Compute total mass enclosed by radii 
        #no need to account for M33 lacking a bulge component, already
        #accounted for in MassEnclosedTotal method.
        M = self.MassEnclosedTotal(r)
        
        #apply units to radius array
        r = r*u.kpc
        
        vel = np.sqrt(G*M/r)
        
        #return final circular velocity array
        return vel
    
    def HernquistVCirc(self,radius,a,Mhalo):
        """This method computes the circular speed using the Herquist
        mass profile.
        
        Inputs: radius is the distance from the COM of the galaxy
                    where mass will be enclosed
                
                a is the scale factor
                
                Mhalo is the halo mass of the galaxy in solMass units
            
        Output: circ_vel is the circular speed in units of km/s
        """
        #import gravitational constant and adjust units for our purposes
        from astropy import constants as const
        G = const.G.to(u.kpc*u.km**2/u.s**2/u.Msun)
        
        #call the HernquistMass method to compute the Mass 
        M = self.HernquistMass(radius,a,Mhalo)
        
        #compute the circular velocity
        circ_vel = np.sqrt(G*M/radius)
        
        return circ_vel
    

           
#%%    
MW = MassProfile("MW",0) #initialize mass profile for MW
r = np.arange(0.25, 30.5, 1.5) #create array of radii

halo_mass = MW.MassEnclosed(1,r) #get enclosed halo masses at each radius
disk_mass = MW.MassEnclosed(2,r) #get enclosed disk masses at each radius
bulge_mass = MW.MassEnclosed(3,r) #get enclosed bulge masses at each radius
total_mass = MW.MassEnclosedTotal(r)

Hern_mass = MW.HernquistMass(r, 63, 1.975*1e12) #scale factor 63 matches halo mass for MW

#plotting mass profiles of components of the MW
fig = plt.figure(figsize=(7,7))
ax = plt.subplot(111)

plt.plot(r, halo_mass, color='black', linewidth=3, label='Halo')
plt.plot(r, disk_mass, color='red', linewidth=3, label='Disk')
plt.plot(r, bulge_mass, color='green', linewidth=3, label='Bulge')
plt.plot(r, total_mass, color='blue', linewidth=3, label='Total')
plt.plot(r, Hern_mass, color='yellow', linewidth=3, label='Hernquist a=63')

plt.xlabel('Radius (kpc)')
plt.ylabel('Mass (solMass)')
plt.title('MW Mass Profile')

plt.semilogy() #y axis is in log 

#ax.set_ylim(0,1e14) #extend y axis range

legend = ax.legend(loc='lower right', fontsize='x-large')

#1.975*1e12Msun MW
#1.921*1e12Msun M31
#0.187*1e12Msun M33

plt.savefig('MW Mass Profile.png')
#%%
M31 = MassProfile("M31",0) #initialize mass profile for M31
r = np.arange(0.25, 30.5, 1.5) #create array of radii

halo_mass = M31.MassEnclosed(1,r) #get enclosed halo masses at each radius
disk_mass = M31.MassEnclosed(2,r) #get enclosed disk masses at each radius
bulge_mass = M31.MassEnclosed(3,r) #get enclosed bulge masses at each radius
total_mass = M31.MassEnclosedTotal(r)

Hern_mass = M31.HernquistMass(r, 60, 1.921*1e12) #a = 60

#plotting mass profiles of components of M31
fig = plt.figure(figsize=(7,7))
ax = plt.subplot(111)

plt.plot(r, halo_mass, color='black', linewidth=3, label='Halo')
plt.plot(r, disk_mass, color='red', linewidth=3, label='Disk')
plt.plot(r, bulge_mass, color='green', linewidth=3, label='Bulge')
plt.plot(r, total_mass, color='blue', linewidth=3, label='Total')
plt.plot(r, Hern_mass, color='yellow', linewidth=3, label='Hernquist a=60')

plt.xlabel('Radius (kpc)')
plt.ylabel('Mass (solMass)')
plt.title('M31 Mass Profile')

plt.semilogy() #y axis is in log 

#ax.set_ylim(0,1e14) #extend y axis range

legend = ax.legend(loc='lower right', fontsize='x-large')

plt.savefig('M31 Mass Profile.png')

#%%
M33 = MassProfile("M33",0) #initialize mass profile for MW
r = np.arange(0.25, 30.5, 1.5) #create array of radii

halo_mass = M33.MassEnclosed(1,r) #get enclosed halo masses at each radius
disk_mass = M33.MassEnclosed(2,r) #get enclosed disk masses at each radius
total_mass = M33.MassEnclosedTotal(r)

Hern_mass = M33.HernquistMass(r, 25, 0.187*1e12) #a = 25

#plotting mass profiles of components of M33
fig = plt.figure(figsize=(7,7))
ax = plt.subplot(111)

plt.plot(r, halo_mass, color='black', linewidth=3, label='Halo')
plt.plot(r, disk_mass, color='red', linewidth=3, label='Disk')
plt.plot(r, total_mass, color='blue', linewidth=3, label='Total')
plt.plot(r, Hern_mass, color='yellow', linewidth=3, label='Hernquist a=25')

plt.xlabel('Radius (kpc)')
plt.ylabel('Mass (solMass)')
plt.title('M33 Mass Profile')

plt.semilogy() #y axis is in log 

#ax.set_ylim(0,1e14) #extend y axis range

legend = ax.legend(loc='lower right', fontsize='x-large')

plt.savefig('M33 Mass Profile.png')

#%%
MW = MassProfile("MW",0) #initialize mass profile for MW
r = np.arange(0.25, 30.5, 1.5) #create array of radii

halo_vel = MW.CircularVelocity(1,r) 
disk_vel = MW.CircularVelocity(2,r) 
bulge_vel = MW.CircularVelocity(3,r)
total_vel = MW.CircularVelocityTotal(r)
halo_mass = MW.MassEnclosed(1,r)

Hern_vel = MW.HernquistVCirc(r, 63, 1.975*1e12) 

#plotting mass profiles of components of the MW
fig = plt.figure(figsize=(7,7))
ax = plt.subplot(111)

plt.plot(r, halo_vel, color='black', linewidth=3, label='Halo')
plt.plot(r, disk_vel, color='red', linewidth=3, label='Disk')
plt.plot(r, bulge_vel, color='green', linewidth=3, label='Bulge')
plt.plot(r, total_vel, color='blue', linewidth=3, label='Total')
plt.plot(r, Hern_vel, color='yellow', linewidth=3, label='Hernquist a=63')

plt.xlabel('Radius (kpc)')
plt.ylabel('Circular Velocity (km/s)')
plt.title('MW Rotation Curve')

plt.semilogy() #y axis is in log 

#ax.set_ylim(0,1e14) #extend y axis range

legend = ax.legend(loc='lower right', fontsize='x-large')

plt.savefig('MW Rotation Curve.png')

#%%
M31 = MassProfile("M31",0) #initialize mass profile for MW
r = np.arange(0.25, 30.5, 1.5) #create array of radii

halo_vel = M31.CircularVelocity(1,r) 
disk_vel = M31.CircularVelocity(2,r) 
bulge_vel = M31.CircularVelocity(3,r) 
total_vel = M31.CircularVelocityTotal(r)
halo_mass = M31.MassEnclosed(1,r)

Hern_vel = M31.HernquistVCirc(r, 60, 1.921*1e12) 

#plotting Rotation Curve for M31
fig = plt.figure(figsize=(7,7))
ax = plt.subplot(111)

plt.plot(r, halo_vel, color='black', linewidth=3, label='Halo')
plt.plot(r, disk_vel, color='red', linewidth=3, label='Disk')
plt.plot(r, bulge_vel, color='green', linewidth=3, label='Bulge')
plt.plot(r, total_vel, color='blue', linewidth=3, label='Total')
plt.plot(r, Hern_vel, color='yellow', linewidth=3, label='Hernquist a=60')

plt.xlabel('Radius (kpc)')
plt.ylabel('Circular Velocity (km/s)')
plt.title('M31 Rotation Curve')

plt.semilogy() #y axis is in log 

#ax.set_ylim(0,1e14) #extend y axis range

legend = ax.legend(loc='lower right', fontsize='x-large')

plt.savefig('M31 Rotation Curve.png')

#%%
M33 = MassProfile("M33",0) #initialize mass profile for MW
r = np.arange(0.25, 30.5, 1.5) #create array of radii

halo_vel = M33.CircularVelocity(1,r) 
disk_vel = M33.CircularVelocity(2,r) 
total_vel = M33.CircularVelocityTotal(r)
halo_mass = M33.MassEnclosed(1,r)

Hern_vel = MW.HernquistVCirc(r, 25, 0.187*1e12)

#plotting rotation curve of M33
fig = plt.figure(figsize=(7,7))
ax = plt.subplot(111)

plt.plot(r, halo_vel, color='black', linewidth=3, label='Halo')
plt.plot(r, disk_vel, color='red', linewidth=3, label='Disk')
plt.plot(r, total_vel, color='blue', linewidth=3, label='Total')
plt.plot(r, Hern_vel, color='yellow', linewidth=3, label='Hernquist a=25')

plt.xlabel('Radius (kpc)')
plt.ylabel('Circular Velocity (km/s)')
plt.title('M33 Rotation Curve')

plt.semilogy() #y axis is in log 

#ax.set_ylim(0,1e14) #extend y axis range

legend = ax.legend(loc='lower right', fontsize='x-large')

plt.savefig('M33 Rotation Curve.png')