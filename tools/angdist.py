import numpy as np
import math
#Constants
c = 2.9979e8 #speed of light
G = 6.6726e-11 #Gravitational constant in si
m_sol = 1.989e30 #Solar mass in kilograms
rad_arcs = 206264.806 #Conversion factor for going from radians to arcseconds
Mpc_to_m = 3.08567758e22


class constants:
    
	'''
	class of cosmology constants
	'''
	def __init__(self, omega_m, omega_l, H_0):
        
		self.omega_m = omega_m 					#matter density
		self.omega_l = omega_l					#dark energy density / lambda
		self.omega_k = 1.0 - omega_l - omega_m 	#curvature density, for flat normally zero
		self.H = H_0						#Hubble parameter today
		self.h = H_0/100.					#dimensionless Hubble parameter
		self.d_h = 299792.458/H_0			#Hubble distance, c/H_0 with c in units of km * s^-1


        
class distances:
    
	'''
	class of different distances in cosmology
	'''
    
	def codist(self, z1, z2, cosmo):
        
		'''
		This function calculates the comoving distance d_c in units of Mpc.
		z1, z2 :: redshifts with z1<z2
		cosmo  :: class of cosmology constants (omega_m, omega_l, h, etc.)
		d_c    :: comoving distance
		'''
		dz = 0.0001
		z = z1 + 0.5 * dz
		d_c = 0.0
		while(z < z2):
			f = 1.0 + z
			demon = cosmo.omega_m * f**3 + cosmo.omega_k * f**2 + cosmo.omega_l
			d_c += dz / np.sqrt(demon)
			z += dz
		return d_c*cosmo.d_h
    
	def tcodist(self, z1, z2, d_c, cosmo):
        
		'''
		This function calculates the transverse comoving distance d_M in units of Mpc.
		z1, z2	:: redshifts with z1<z2
		cosmo	:: class of cosmology constants (omega_m, omega_l, h, etc.)
		d_c		:: comoving distance
		d_M		:: transverse comoving distance
		'''
		d_M = 0.0
		f = np.sqrt(np.abs(cosmo.omega_k))
		if(cosmo.omega_k > 0.0):
			d_M = cosmo.d_h/f*np.sinh(f*d_c/cosmo.d_h)
			#print 'Omega_k > 0'
		elif(cosmo.omega_k < 0.0):
			d_m = cosmo.d_h/f*np.sin(f*d_c/cosmo.d_h)
			#print 'Omega_k < 0'
		else:
			d_M = d_c
		return d_M
    
	def angdist(self, d_M, z, d_A):
        
		'''
		This function calculates the angular diameter distance d_M in units of Mpc.
		z1		:: redshift z
		d_M		:: transverse comoving distance
		d_A		:: angular diameter distance
		'''
		d_A = 0.0
		d_A = d_M/(1.0+z)
		return d_A



def getAngularDiameterDistances(cosmo,z_o,z_l,z):
        N=len(z)

        # Define cosmology.
#         cosmo = constants(0.3, 0.7, 70)

#         #Check Curvature is near zero
#         print ("Omega_k: ", cosmo.omega_k)

#         # Define redshifts
#         z_o = 0.0  # observer redshift
#         z_l = 0.1832  # lens redshift for Abell 2744

        # Initialize angular diameter arrays
        d_ol = np.zeros(N)
        d_ls = np.zeros(N)
        d_os = np.zeros(N)
        scrit = np.zeros(N)

        d_c, d_M, d_A = 0.0, 0.0, 0.0  # Initialize distances
        dist = distances()  # create cosmology.distance instance

        # For loop over images to calculate the angular diameter distances and critical surface densities.
        for i in range(N):
            z_s = z[i]

            # Calculate the observer to lens ang dist.
            d_c = dist.codist(z_o, z_l, cosmo)
            dco=d_c
            d_M = dist.tcodist(z_o, z_l, d_c, cosmo)
            d_A = dist.angdist(d_M, z_l, d_A)

            d_ol[i] = d_A * Mpc_to_m  # Convert from Mpc to m

            # Calculate the lens to source ang dist.
            d_c = dist.codist(z_l, z_s, cosmo)
            d_M = dist.tcodist(z_l, z_s, d_c, cosmo)
            d_A = dist.angdist(d_M, z_s, d_A)

            d_ls[i] = d_A * Mpc_to_m  # Convert from Mpc to m

            # Calculate the observer to source ang dist.
            d_c = dist.codist(z_o, z_s, cosmo)
            d_M = dist.tcodist(z_o, z_s, d_c, cosmo)
            d_A = dist.angdist(d_M, z_s, d_A)

            d_os[i] = d_A * Mpc_to_m  # Convert from Mpc to m

            #Calculate critical surface density
            scrit[i] = c**2/(4*np.pi*G)*d_os[i]/d_ol[i]/d_ls[i]

        # print(scrit,z,d_ol)
        return d_ol, d_os, d_ls, scrit

