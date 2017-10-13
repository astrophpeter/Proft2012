from astropy.time import Time



class Lens:

        #rough start and end time for the gaia mission - scaled in
	#Barycentric coordinate time (TCB)

	Gaia_start_time = Time('2014-01-01',scale='tcb')
	Gaia_end_time = Time('2019-01-01' ,scale='tcb')

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
		self.epoch_0 = Time(epoch_0,scale=scale_in,format=format_in)

	def getId(self):
		return self.id

#        def getRefEpoch(format):
		
lens1 = Lens(12,1,2,3,4,5)

ans = lens1.getId()

print(ans)


