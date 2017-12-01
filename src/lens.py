from astropy.time import Time
import numpy as np
import matplotlib.pyplot as plt
#import doctest
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
import math

class Lens:

        #rough start and end time for the gaia mission - scaled in
	#Barycentric coordinate time (TCB)

	GAIA_START = '2013-01-01'
	GAIA_END = '2025-01-01'
 
        #min angular with for lens to pass in front of source [Arcseconds]
	min_angular_width = 0.7  

	#mas to degree conversion
	mas_to_deg = 0.0000002777777

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
		#set all properties for lens
		self._id = id
		self._ra_0 = ra_0
		self._dec_0 = dec_0
		self._pmra = pmra
		self._pmdec = pmdec
		self._epoch_0 = epoch_0
		
		#calulate the search window for this lens
		#done on initialisation to speed things up.
		self._lens_box = self.get_lens_box()

        
	def getId(self):
		"""
		Returns the ID of a lens as specified in the
		source_id column of GAIA TGAS Table.
               
		Returns:
			id (long) : source id from source_id
			column in GAIA TGAS Table.
		
		Units Tests:
		>>> unittestinglens.getId()
		123456789
		"""
		
		return self._id

	def get_eq_coords(self,ra,dec):
		"""
		Returns a speficed ra,dec postion in degress
		in equatorial coordinates [ra*cos(dec),dec] 
	
		Args: 
			ra (double) : Barycentric right ascension
				      [Degrees]
			dec (double): Barycentric declination
				      [Degrees]
		Reuturns:
			Eq_coord (array(double)) : Equatorial 
					coordinates [ra*cos(dec),dec]
		
		Unit Tests:
		>>> math.isclose(unittestinglens.get_eq_coords(30.0,60.0)[0],15.0)
		True
	
		>>> math.isclose(unittestinglens.get_eq_coords(30.0,60.0)[1],60.0)
                True
		
		"""	
		return [ra* np.cos(np.deg2rad(dec)),dec]
        

	def get_eq_coords_at_epoch(self,epoch):
		"""
		Returns the equtorial coordinates of the lens
		and the time epoch.
		Args:
			epoch (double) : Time at which to get coordinates
					 Julian year with Barycentric 
                                         coordinate time (TCB)
		Reutrns:
			Eq_coord (array(double)) : Coordinates of the lens
						   at time = epoch in the form
						   [ra*cos(dec),dec][degrees]
		Unit Tests:
		>>> math.isclose(unittestinglens.get_eq_coords_at_epoch(2025.0)[0],15.000444444)
		True
		>>> math.isclose(unittestinglens.get_eq_coords_at_epoch(2025.0)[1],60.000444444)
                True
		>>> math.isclose(unittestinglens.get_eq_coords_at_epoch(1956.0)[0],14.99852778)
		True
		>>> math.isclose(unittestinglens.get_eq_coords_at_epoch(1956.0)[1],59.99852779)
		True
		>>> math.isclose(unittestinglens.get_eq_coords_at_epoch(2017.0)[0],15.00022222)
		True
		>>> math.isclose(unittestinglens.get_eq_coords_at_epoch(2017.0)[1],60.00022222)
                True
		"""
		
		
                #compute ra*cos(dec) at time epoch 
		#This is where I suspect the issues lie - had a dicussion with Leigh about
		# this and we agreed it was fine however.
		#note self._pmra already includes the cos(dec) factor. 
		raCosDec = (self._ra_0*np.cos(np.deg2rad(self._dec_0))) + (epoch - self._epoch_0)* self._pmra * self.mas_to_deg
		
                #compute the dec at time epoch
		dec = self._dec_0 + (epoch - self._epoch_0) * self._pmdec * self.mas_to_deg
		
		return [raCosDec,dec]

	

	def datetime_to_jyTCB(self,date):
		"""
		Returns the value of a time in Julian years from a date string
		of the format 'YYYY-MM-DD' scaled in (TCB)
		Args: 
			date (string) : string format of the date 'YYYY-MM-DD'
		Returns:
			julianYear (double) : the date in julian years
		
		Units Tests:
		#>>> math.isclose(unittestinglens.datetime_to_jyTCB('2009-01-01 00:00:00'),2009.0)
		True
		#>>> math.isclose(unittestinglens.datetime_to_jyTCB('2016-07-02'),2016.5)
		True
		"""
		
		time = Time(date,scale='tcb')
		time.format = 'jyear'
		
		return time.value
	                
	def get_GAIA_start_end_coords(self):
		"""
		Returns the position of the lens at the start and end times of the 
		GAIA missin defined by the class constants GAIA_START GAIA_END
		Returns:
			coords (array(double)) : Cordinates of the lens at the
						 start and end of the gaia mission
						[[start],[end]] -
						[[ra * cos(dec),dec] ,[ra*cos(dec),dec]] 
						[Degrees]
		"""

		start= self.get_eq_coords_at_epoch(self.datetime_to_jyTCB(self.GAIA_START))
		end = self.get_eq_coords_at_epoch(self.datetime_to_jyTCB(self.GAIA_END))
		
		return [start,end]	

	def get_lens_box(self):
		"""
		Returns a search box of with defined by the class constant 
		min_angular_width and with length of the GAIA Mission.
		This used as a search window for other catalogues to
		find sources the lens passes close to.
		Returns:
			box (array(double)) : Four corners of the search
					      window box [[x1,y1],[x2,y2],
					      [x3,y3],[x4,y4]] with format	
					      [ra*cos(dec),dec][degrees]   
		"""

		start_end = self.get_GAIA_start_end_coords()
                
		start = start_end[0]
		end = start_end[1]
		
		#some coordinate geometry in the (ra*cos(dec),dec) plane to get the
		#the coners of the box
		hypt = np.sqrt((start[0] - end[0])**(2) + (start[1] - end[1])**(2))

		dX = abs((1.0/3600)*self.min_angular_width * (start[1] - end[1]) / hypt)
		dY = abs((1.0/3600)*self.min_angular_width * (start[0] - end[0]) / hypt)		
               
		#put the corners of the box in the right place.
		if ((start[0] < end[0]) and (start[1] < end[1])):
			return [[start[0] - dX,start[1] + dY],[start[0] + dX,start[1] - dY],[end[0] + dX,end[1] - dY],[end[0] -dX,end[1] + dY]]
		elif ((start[0] < end[0]) and (start[1] >  end[1])):
			return [[start[0] - dX,start[1] - dY],[start[0] + dX,start[1] + dY],[end[0] + dX,end[1] + dY],[end[0] -dX,end[1] - dY]]
		elif ((start[0] > end[0]) and (start[1] >  end[1])):
			return [[start[0] - dX,start[1] + dY],[start[0] + dX,start[1] - dY],[end[0] + dX,end[1] - dY],[end[0] -dX,end[1] + dY]]
		else:
			return [[start[0] + dX,start[1] + dY],[start[0] - dX,start[1] - dY],[end[0] - dX,end[1] - dY],[end[0] + dX,end[1] + dY]]
	
	
	def is_coord_in_box(self,ra,dec):
		"""
		Checks if a source with coordinates
		ra, dec is in the search window for the lens.
		Args:
			ra (double) : Barycentric right ascension
                                      of a source [Degrees]
			dec (double) :Barycentric declination
                                      of a source [Degrees]
		Returns:
			inBox (bool) : Return true if source is in
				       len's search box, false otherwise.
		
		Unit Tests:
		>>> unittestinglens.is_coord_in_box(30.00069760,60.00024)
		True
		>>> unittestinglens.is_coord_in_box(30.00069760,60.00124)
                False
		"""

		point = Point(ra*np.cos(np.deg2rad(dec)),dec)
		box1 = self._lens_box
		
		polygon = Polygon([(box1[0][0],box1[0][1]),(box1[1][0],box1[1][1]),(box1[2][0],box1[2][1]),(box1[3][0],box1[3][1])])
		
		return polygon.contains(point)

	def get_angular_separation_at_epoch(self,epoch,source_ra,source_dec):
		"""
		Returns the agular separation between the lens and source
		at time epoch assuming the source does not move.
		Implemented using the distance formula in (2.1) in
		Candidate Lensing Events for Gaia - doc.
		Args: 
			epoch (double):  Time at which to find separation epoch 
					 in Julian year with Barycentric 
					 coordinate time (TCB)
			
			source_ra (double) : Barycentric right ascension of source
                                             [Degrees] 
			source_dec (double) : Barycentric declination of source
                                             [Degrees]
		Returns:
			angular_sep (double) : Angular separation of source and 	
				     		lens at time = epoch [mas]
		Units Tests:
		testing this the usual way is tricky as mas is so precise, I think 1 d.p of mas
		precision is fine here 
		>>> math.isclose(np.round(unittestinglens.get_angular_separation_at_epoch(2017.0,30.00069760,60.00024),1),90.4)
		True
		
		
		
		"""
		cosSourceDec = np.cos(np.deg2rad(source_dec))	
		coords = self.get_eq_coords_at_epoch(epoch)

		return (np.sqrt(((source_dec - coords[1])**2) + (((source_ra * cosSourceDec) - (coords[0]))**2))) / self.mas_to_deg 
		
	def get_time_of_closest_app(self,source_ra,source_dec):
		"""
		Returns the time of closest approach between a 
		lens and source assuming tha the source does not
		move. Implemented using formula (1) in 
		Candidate Lensing Events for Gaia doc.
		
		Args:
			source_ra (double) : Barycentric right ascension of source
                                             [Degrees]
                        source_dec (double) : Barycentric declination of source
                                             [Degrees]
		Returns:
	
			time (double) : time of closest approach [Julian Years]
		Unit Tests:
		>>> math.isclose(np.round(unittestinglens.get_time_of_closest_app(30.00069760,60.00024),2),2017.64)
		True
	
		"""
		
		#convert \mu_{\apha *} from [mas/yr] to [deg/yr]  		
		#note this term already includes the cos(dec)
		pmRaDeg = self._pmra * self.mas_to_deg
		
		#convert \mu_{\dec} from [mas/yr] to [deg/yr]
		pmDecDeg = self._pmdec * self.mas_to_deg

                #pre compute some frequently used quantities 
		cosSourceDec = np.cos(np.deg2rad(source_dec))
		cosLensDec = np.cos(np.deg2rad(self._dec_0))
	
		#calulate time of closest approach from the lens ref epoch
		#_epoch_0
		top = - (((pmDecDeg) * (self._dec_0 - source_dec)) + (pmRaDeg) * ((self._ra_0 * cosLensDec) - (source_ra * cosSourceDec)))
		bottom = (pmDecDeg)**2 + (pmRaDeg)**2
	
		#return the time of closest approach 
		return self._epoch_0 + top / bottom 

	def get_dist_of_closest_app(self,source_ra,source_dec):
		"""
		Returns the distance of closest apprach in [mas]
		Assuming the source is stationary. 
		
		Args: 
			source_ra (double) : Barycentric right ascension of source
                                             [Degrees]
                                             
                        source_dec (double) : Barycentric declination of source
                                             [Degrees] 
		
		Returns: distance (double) : Min Angular separation of source and
                                             lens [mas]
		"""

		
		time_closest = self.get_time_of_closest_app(source_ra,source_dec)
		return self.get_angular_separation_at_epoch(time_closest,source_ra,source_dec)

	def get_time_closest_app_sourcepm(self,source_ra,source_dec,source_pmra,source_pmdec):
		"""
		Calculates the time of closest approach of
		a source and takes into account the proper
		motion of the source.
		Args: 
			lens1 (lens object) : the lens 
		
			source_ra (double) : Right acension of the
			             source [Degress] at 2015.0
			source_dec (double) : Declination of the 
				     source [Degrees] at 2015.0
			source_epoch (double) : Reference epoch of 
					the source [Julian Yrs]
			source_pmra (double) ; Proper motion in ra
					of source [mas/yr]
			source_pmdec (double) : Proper motion in dec
					of source [mas/yr]
		Returns:
			time (double) : time of closest approach 
					[Julian yrs]
		"""
		#convert \mu_{\apha *} from [mas/yr] to [deg/yr]
        	#note this term already includes the cos(dec)
		LpmRaDeg = self._pmra * self.mas_to_deg
		SpmRaDeg = source_pmra * self.mas_to_deg

		#convert \mu_{\dec} from [mas/yr] to [deg/yr]
		LpmDecDeg = self._pmdec * self.mas_to_deg
		SpmDecDeg = source_pmdec * self.mas_to_deg

		#pre compute some frequently used quantities
		cosSourceDec = np.cos(np.deg2rad(source_dec))
		cosLensDec = np.cos(np.deg2rad(self._dec_0))

		top = - (((LpmDecDeg-SpmDecDeg) * (self._dec_0 - source_dec)) + (LpmRaDeg-SpmRaDeg) * ((self._ra_0 * cosLensDec) - (source_ra * cosSourceDec)))
		bottom = (LpmDecDeg-SpmDecDeg)**2 + (LpmRaDeg-SpmRaDeg)**2

		return 2015.0 + (top / bottom)

	def get_dist_closest_app_sourcepm(self,source_ra,source_dec,source_pmra,source_pmdec):
		"""
		Args:
                	lens1 (lens object) : the lens
                	source_ra (double) : Right acension of the
                                     source [Degress] at 2015.0
                	source_dec (double) : Declination of the
                                     source [Degrees] at 2015.0
                	source_pmra (double) ; Proper motion in ra
                                        of source [mas/yr]
                	source_pmdec (double) : Proper motion in dec
                                        of source [mas/yr]
        	Returns:
                	distance (double) : distance of closest approach
                        	        [Julian yrs]
		"""	

		time = self.get_time_closest_app_sourcepm(source_ra,source_dec,source_pmra,source_pmdec)
		Lcoords = self.get_eq_coords_at_epoch(time)
	
		#compute the source position at time of closest approach
		SraCosDec = (source_ra*np.cos(np.deg2rad(source_dec))) + (time - 2015.0)* source_pmra * self.mas_to_deg  
		Sdec = source_dec + (time - 2015.0) * source_pmdec * self.mas_to_deg
	
		return np.sqrt((Sdec-Lcoords[1])**2+(SraCosDec - Lcoords[0])**2) / self.mas_to_deg


	#def get_time_closest_app_lens(self, otherlens):
	#	"""
	#	Returns the time of cloesest approach between a lens
        #		and another source where the source (otherlens)
#		which takes into account proper motion.
#		
#		Args:
#			otherlens (lens object) : Source built using another
#						  lens object so pm can be acconted
#						  for 
#		Returns:
#			time (double) : time of closest approach [Julian Years]
#		"""
#		
#				
#		#convert \mu_{\apha *} from [mas/yr] to [deg/yr]
#                #note this term already includes the cos(dec)
#                pmRaDegSelf = self._pmra * self.mas_to_deg
#
#                #convert \mu_{\dec} from [mas/yr] to [deg/yr]
#                pmDecDegSelf = self._pmdec * self.mas_to_deg
#
#                
#		#convert \mu_{\apha *} from [mas/yr] to [deg/yr]
#                #note this term already includes the cos(dec)
#                pmRaDegOther = otherlens._pmra * self.mas_to_deg
#
#                #convert \mu_{\dec} from [mas/yr] to [deg/yr]
#                pmDecDegOther = otherlens._pmdec * self.mas_to_deg
#
#
#		#pre compute some frequently used quantities
#                cosOtherDec = np.cos(np.deg2rad(otherlens._dec_0))
#                cosSelfDec = np.cos(np.deg2rad(self._dec_0))
#
#		#find the postion of the other lens at the refernce epoch of
#		# this lens.
#
#		coord = otherlens.get_eq_coord_at_epoch(self._epoch_0)
#		
#
#		top = - (((pmDecDeg) * (self._dec_0 - source_dec)) + (pmRaDeg) * ((self._ra_0 * cosLensDec) - (source_ra * cosSourceDec)))
	
	
                
		
#doctest.testmod(extraglobs={'unittestinglens':Lens(123456789,30.0,60.0,100.0,100.0,2009.0)})

