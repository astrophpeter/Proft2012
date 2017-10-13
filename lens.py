from astropy.time import Time
import numpy as np
import matplotlib.pyplot as plt

class Lens:

        #rough start and end time for the gaia mission - scaled in
	#Barycentric coordinate time (TCB)

	GAIA_START = '2014-01-01'
	GAIA_END = '2019-01-01'
 
        #min angular with for lens to pass in front of source [Arcseconds]
	min_angular_width = 0.7  

	#mas to degree conversion
	mas_to_deg = (1.0 / 3600)*10**(-6)

	def __init__ (self,id,ra_0,dec_0,pmra,pmdec,epoch_0,scale_in='tcb',format_in='jyear'):
		"""
		lens object constructor

		Args:
			id (long): Lens id as specified by the
		                         'source_id' column in Gaia
					 TGAS catalogue

			ra_0 (double): Barycentric right ascension 
                                       of the source in ICRS at the
                                       reference epoch Angle [degrees]

			dec_0 (double): Barycentric declination of 
                                        the source in ICRS at the 
                                        reference epoch Angle [degress]  
		

                        pmra (double): Proper motion in right ascension 
                                       direction Angular Velocity[mas/year]

                        pmdec (double): Proper motion in declination 
                                        direction Angular Velocity[mas/year]

			epoch_0 (double) : Reference epoch in Julian year 
                                           with Barycentric coordinate time (TCB)
		"""
		self.id = id
		self.ra_0 = ra_0
		self.dec_0 = dec_0
		self.pmra = pmra
		self.pmdec = pmdec
		self.epoch_0 = epoch_0

                #Time(epoch_0,scale=scale_in,format=format_in)

	def getId(self):
		return self.id

	def get_eq_coords(self,ra,dec):
		return [dec * np.cos(np.deg2rad(ra)),ra]
        

	def get_eq_coords_at_epoch(self,epoch):
		ra = self.ra_0 + (epoch - self.epoch_0)* self.pmra * self.mas_to_deg
		dec = self.dec_0 + (epoch - self.epoch_0) * self.pmdec * self.mas_to_deg
		return self.get_eq_coords(ra,dec)

	def datetime_to_jyTCB(self,date):
		time = Time(date,scale='tcb')
		time.format = 'jyear'
		return time.value
                
	def get_GAIA_start_end_coords(self):
		start= self.get_eq_coords_at_epoch(self.datetime_to_jyTCB(self.GAIA_START))
		end = self.get_eq_coords_at_epoch(self.datetime_to_jyTCB(self.GAIA_END))
		return [[start[0],end[0]],[start[1],end[1]]]	
		
		
lens1 = Lens(12,1,2,1500000,1500000,2012.0)


ans = lens1.getId()
x = lens1.get_GAIA_start_end_coords()
init = lens1.get_eq_coords_at_epoch(2012.0)
plt.scatter(x[0],x[1])
plt.scatter(init[0],init[1])
plt.text(1.996, 1.01, r'Lens with $\delta_{0}=1,\alpha_{0}=2,\mu_{\delta}=\mu_{\alpha}=150$[mas/yr]')
plt.ylabel(r'$\delta$')
plt.xlabel(r'$\alpha\cos\delta$')
plt.text(1.999,1,r't$_{0}$')
plt.text(2.001,1,r' t$_{GAIA-START}$')
plt.text(2.003,1.0030,r't$_{GAIA-END}$')
plt.show()
#print(y)
print(x)

print(ans)


