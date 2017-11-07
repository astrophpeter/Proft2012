from astropy.time import Time
import numpy as np
import matplotlib.pyplot as plt
import doctest
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from scipy.optimize import minimize_scalar

class Lens:

        #rough start and end time for the gaia mission - scaled in
	#Barycentric coordinate time (TCB)

	GAIA_START = '2013-01-01'
	GAIA_END = '2025-01-01'
 
        #min angular with for lens to pass in front of source [Arcseconds]
	min_angular_width = 0.7  

	#mas to degree conversion
	mas_to_deg = (1.0 / 3600)*10**(-3)

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
		"Candidate Lensing Events for Gaia" doc.

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
		"""
		cosSourceDec = np.cos(np.deg2rad(source_dec))	
		coords = self.get_eq_coords_at_epoch(epoch)

		return (np.sqrt(((source_dec - coords[1])**2) + (((source_ra * cosSourceDec) - (coords[0]))**2))) / self.mas_to_deg 
		
	def get_time_of_closest_app(self,source_ra,source_dec):
		"""
		Returns the time of closest approach between a 
		lens and source assuming tha the source does not
		move. Implemented using formula (1) in 
		"Candidate Lensing Events for Gaia" doc.
		
		Args:
			source_ra (double) : Barycentric right ascension of source
                                             [Degrees]

                        source_dec (double) : Barycentric declination of source
                                             [Degrees]

		Returns:
	
			time (double) : time of closest approach [Julian Years]

		"""
		
		#convert \mu_{\apha *} from [mas/yr] to [deg/yr]  		
		#note this term already includes the cos(dec)
		pmRaDeg = self._pmra * self.mas_to_deg
		
		#convert \mu_{\dec *} from [mas/yr] to [deg/yr]
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



#doctest.testmod(extraglobs={'testlens':Lens(0,0,0,0,0,0)})



