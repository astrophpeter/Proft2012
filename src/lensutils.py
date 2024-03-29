import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import numpy as np
import aplpy
from astropy.wcs import WCS
from astropy.io import fits, ascii
from astroquery.skyview import SkyView
import sqlutil as sqlutil 
import doctest
import math


def plot_mwd(RA,Dec,org=0,title='Mollweide projection', projection='mollweide',filename='Mollweide-Projection'):
	''' RA, Dec are arrays of the same length.
	RA takes values in [0,360), Dec in [-90,90],
	which represent angles in degrees.
	org is the origin of the plot, 0 or a multiple of 30 degrees in [0,360).
	title is the title of the figure.
	projection is the kind of projection: 'mollweide', 'aitoff', 'hammer', 'lambert'
	'''
	coord = SkyCoord(ra=RA,dec=Dec,unit='deg',frame='icrs')
	l = coord.galactic.l.deg
	b = coord.galactic.b.deg
	RA = l	
	Dec = b
	x = np.remainder(RA+360-org,360) # shift RA values
	ind = x>180
	x[ind] -=360    # scale conversion to [-180, 180]
	x=-x    # reverse the scale: East to the left
	tick_labels = np.array([150, 120, 90, 60, 30, 0, 330, 300, 270, 240, 210])
	tick_labels = np.remainder(tick_labels+360+org,360)
	fig = plt.figure(figsize=(10, 5))
	ax = fig.add_subplot(111, projection=projection, axisbg ='white')
	ax.scatter(np.radians(x),np.radians(Dec))  # convert degrees to radians
	ax.set_xticklabels(tick_labels)     # we add the scale on the x axis
	ax.set_title(title)
	ax.title.set_fontsize(15)
	ax.set_xlabel("l")
	ax.xaxis.label.set_fontsize(12)
	ax.set_ylabel("b")
	ax.yaxis.label.set_fontsize(12)
	ax.grid(True)
	plt.savefig('out.png')


def plt_len_traj(lens,sourceRa,sourceDec):
	"""
	Creates and saves a plot of the trajectory of
	of a lens as it goes past closest approah of
	a source assumed stationary

	Args:
		lens (lens object) : the lens star
	
		sourceRa (double) : Source Right acession
				    [Degrees]
		
		sourceDec (double) : Source Declination 
				    [Degrees]

	Returns:
		file (.png) : save in current directoy with
			      name filename.

	"""
	
	#filename 
	lensidstring = 'trject_test' + str(lens.getId()) + '.png'

	#get two point on lens trjectory so line can be drawn 
	pos1 = lens.get_eq_coords_at_epoch(2009.0)
	pos2 = lens.get_eq_coords_at_epoch(2025.0)
        
	#get the search window for the lens 
	box = np.transpose(lens.get_lens_box())	
	
	#add the first coordinates of search window corners
	#to the end of the lost for plotting purposes to be 
	#able to shade search window with matplotlib 
	boxjX = np.append(box[0],box[0][0])
	boxjY = np.append(box[1],box[1][0])
        
	#find the time and sitance of closest approah of
	# lens to source
	timeCl = lens.get_time_of_closest_app(sourceRa,sourceDec)
	posClose = lens.get_eq_coords_at_epoch(timeCl)
	
	# shade background search window
	plt.fill_between(boxjX,boxjY,alpha=0.3,label='Background source search window')
	
	#set plot limints close to edge of search window 
	plt.xlim(min(box[0]) - 0.0001,max(box[0]) + 0.0001)
	plt.ylim(min(box[1]) - 0.0001,max(box[1]) + 0.0001)
	
	plt.ylabel(r'$\delta$')
	plt.xlabel(r'$\alpha\cos\delta$')
        
	plt.title(r'Lens with $\mu_{\alpha *}$=%.1f and $\mu_{\delta}$=%.1f [mas yr$^{-1}$]'%(lens._pmra,lens._pmdec))
	plt.scatter(sourceRa*np.cos(np.deg2rad(sourceDec)),sourceDec,label='source')
	plt.plot([pos1[0],pos2[0]],[pos1[1],pos2[1]],'r--',label='lens trajectory')
	plt.scatter(posClose[0],posClose[1],label='lens at closest app')
	
	plt.legend()	
	
	plt.savefig(lensidstring)
	
	plt.clf()
	
def plt_lens_env(lens, sourceRa, sourceDec,lensMag,sourceMag,sourceId):
	"""
	Creates and saves a plot of the source and lens 
	environment. Usea DSS sky cut out images. In the 
	plot the lens blue with trajectory in green. The 
	source is identified in red. It is assumed source
	is stationary.

	Ags:
		lens (lens object): the lens star

		sourceRa (double): Source Right acession
                                   [Degrees]

		sourceDec (double) : Source Declination
                                    [Degrees]

	Returns:
		file (.png) : Output file saved with the id
			      of the lens.
	"""

	#define cutout image size [arsec]
	imsize = 2

	##### 2MASS ##### - not currently being used - DSS more complete. 
	#filter_2m = 'j' # 2mass filter to use

	#metadata_url = 'http://irsa.ipac.caltech.edu/ibe/search/twomass/allsky/allsky?POS={0},{1}'.format(sourceRa,sourceDec)
	#metadata = ascii.read(metadata_url, Reader=ascii.ipac.Ipac)
    	
	# selecting the first obs in the requested filter in case there is more than one
	#target_obs = metadata[metadata['filter']==filter_2m][0]
    
	#params = { 'ordate': target_obs['ordate'],
	#	'hemisphere': target_obs['hemisphere'],
	#	'scanno': target_obs['scanno'],
	#	'fname': target_obs['fname'],
	#	'ra': sourceRa,
	#	'de': sourceDec,
	#	'imsize': imsize }

	params = {'ra': sourceRa,'de': sourceDec,'imsize': imsize }
        
	#DSS Search
	url = "http://stdatu.stsci.edu/cgi-bin/dss_search?v=poss2ukstu_red&r={ra}&d={de}&e=J2000&h={imsize}&w={imsize}&f=fits".format(**params)
	hdulist = fits.open(url)

	print(url) 
	
	#Find the time the image was taken YYYY-MM-DD
	timeString  = hdulist[0].header['DATE-OBS'][:10]

	# Get the position of the lens at the time the image was taken
	time = lens.datetime_to_jyTCB(timeString)
	coord = lens.get_eq_coords_at_epoch(time)
    
	#Get ra and dec of lens at the image time (raCos(dec) is converted into ra)
	raLensimag = coord[0] / np.cos(np.deg2rad(lens._dec_0))
	decLensimag = coord[1]
    
	# Find the postition of the lens in the future so its trajectory can be plotted
	timeend = lens.datetime_to_jyTCB('2040-01-01')
	coordend = lens.get_eq_coords_at_epoch(timeend)
	raLensimagend = coordend[0] / np.cos(np.deg2rad(lens._dec_0))
	decLensimagend = coordend[1]
	
	# Find the postition of the lens at gaia epoch 2015.0
	coordGaiaEpoch = lens.get_eq_coords_at_epoch(2015.0)
	raLensimagGaiaEpoch = coordGaiaEpoch[0] / np.cos(np.deg2rad(lens._dec_0))
	decLensimagGaiaEpoch = coordGaiaEpoch[1]
    
	# Find the line of the lens trajectory between the time when the image was taken to
	# to some time in the future. 
	line = np.array([[raLensimag,raLensimagend],[decLensimag,decLensimagend]])
    
	# Plot the image, with reverse gray color map.
	fig = aplpy.FITSFigure(hdulist[0])
	fig.show_colorscale(cmap='gray_r')
	
	#Plot the positions of the source, lens and lens tragectory on the image
	fig.show_markers(sourceRa,sourceDec,marker='o',edgecolor='r',label='source')
	fig.show_markers(raLensimag,decLensimag,marker='o',edgecolor='b',label='lens at image time (1989-11-22)')
	fig.show_lines([line],color='g',linestyle='--',label='lens-trajectory')

	#PLot gaia source positions also in the image
	pos = get_gaia_source_pos(sourceRa,sourceDec,2)
	fig.show_markers(pos[0],pos[1],marker='v',edgecolor='magenta',label='Gaia source')

	#Plot lens at gaia epoch on the image
	fig.show_markers(raLensimagGaiaEpoch,decLensimagGaiaEpoch,marker='>',edgecolor='b')

	#Find Time and distance of Cloeset approach
	timeCl = lens.get_time_of_closest_app(sourceRa,sourceDec)
	distCl = lens.get_angular_separation_at_epoch(timeCl,sourceRa,sourceDec)

        # get some info about the lens system to display on the plot
	Diststr = 'Distance of closest Approach: ' + str("%.2f" % distCl) + ' [mas]'
	Timestr = 'Time of closest Approach: ' + str("%.2f" %timeCl) + ' [Jyr]'
	Lensstr = 'TGAS Lens Id (Blue): ' + str(lens.getId()) 
	Sourcestr = 'PPMXL Source (Red) Id:' + str(sourceId)
	SourceMag = 'PPMXL source mag [j,h,k,b1,b2,r1,r2]: '
	LensMag = 'TGAS Lens Gmag: ' + str("%.2f" %lensMag)
    
	# Add some extra info about the source and lens to the plot.
	fig.add_label(0.75,0.98,Diststr,relative=True)
	fig.add_label(0.75,0.94,Timestr,relative=True)
	fig.add_label(0.75,0.90,Lensstr,relative=True)
	fig.add_label(0.75,0.86,Sourcestr,relative=True)
	fig.add_label(0.80,0.82,SourceMag,relative=True)
	fig.add_label(0.70,0.78,str(sourceMag),relative=True)
	fig.add_label(0.85,0.74,LensMag,relative=True)

	
	filename = 'GaiasourceMatchfuture/TGAS_' + str(lens.getId()) + '.png'
	fig.save(filename,dpi=200)


def get_gaia_source_pos(centerRa,centerDec,size):
	"""
	Get all the gaia source positions are epoch 2015.0
	in a search radius of size.


	Args: 
		centerRA (double) : Right acession of the center of 
				    the search radius [Degrees] 
		
		centerDec (souble) : Declination of the center of the 
			             the search radius [Degrees]

		size (double) : size of the search radius [arminutes]

	Returns: 
		pos (np.array(double)) : Position of gaia sources [
					[Ra,Ra],[dec,dec]..]


	"""


	query = 'select ra,dec from gaia_dr1.gaia_source where ra BETWEEN ' + str(centerRa) + ' - 0.016 and ' + str(centerRa)  + '+ 0.016 AND dec BETWEEN ' + str(centerDec)  + '- 0.016 and ' + str(centerDec) + '+ 0.016'

	ra, dec = sqlutil.get(query,db='wsdb',host='cappc127.ast.cam.ac.uk', user='peter_mcgill', password='Ln3g.wsk')

	 
	return [ra,dec]

def calc_centroid_shift(lensMass,lensDist,impactAngle):
	"""
	calculates the centriod shit of the lens
	from the lens mass and distance

	Args:
		lensMass (double) : mass of the lens
			            [Msol]

		lensDist (double) : Distance to the lens
					[parces]
		
		impactAngle (double) : Angular distance of closest	
					approch between lens and 
					source [mas]

	Return: 
		Centriod shift (double) : Shift of the centriod

					   [mas]
	Unit Tests:
	>>> calc_centroid_shift(0.5,120.0,241)
	0.15
	
	"""
	#speed of light in units [cm/s]
	c = 2.99 * (10 ** 10)
	
	#G in units of [cm^3/gs^2]
	G = 6.672 * (10 **(-8))
 
	#Calculate Enstien Radius
	enstienR = 90.0 * np.sqrt((lensMass) / (lensDist))
		
	print(enstienR)
	
	#calculate mu
	mu = impactAngle / enstienR
	
	print(mu)

	return (mu * enstienR) / ((mu**2) + 2)


def get_closes_gaia_source_match(centerRa,centerDec,size):
        """
        Get all the gaia source information for all gaia sources
	within a small search radius. And have low proper motions
	 


        Args:
                centerRA (double) : Right acession of the center of
                                    the search radius [Degrees]

                centerDec (souble) : Declination of the center of the
                                     the search radius [Degrees]

                size (double) : size of the search radius [arminutes]

        Returns:
                pos (np.array(double)) : Position of gaia sources [
                                        [[id, id,..],[Ra,Ra],[dec,dec],[gmag,gmag...]..]


        """


        query = 'select ra,dec, source_id, phot_g_mag_mean from gaia_dr1.gaia_source where ra BETWEEN ' + str(centerRa) + ' - 0.026 and ' + str(centerRa)  + '+ 0.026 AND dec BETWEEN ' + str(centerDec)  + '- 0.026 and ' + str(centerDec) + '+ 0.026 AND hypot(pmra,pmdec) < 50.0'

        ra, dec, id, gmag = sqlutil.get(query,db='wsdb',host='cappc127.ast.cam.ac.uk', user='peter_mcgill', password='Ln3g.wsk')


        return [ra,dec,id,gmag]	

def lens_on_color_mag(TGASid,filename):
	"""
	Plots the lens with TGASid on a color magnitude diagram
	with the rest of the TGAS sources - uses 2mass [j-k] colour
	and absolute G magnitude, saves the plot with file name.

	Args:
		TGASid (long) : indentified for lens id 
				in TGAS source_id column

		filename (string) : output file name	

	Returns:

		file (.png) : output file of the plot with
			      filename.

	"""

	#Get color mag info for lens
	querystringLens = 'select j_m,k_m,parallax,phot_g_mean_mag from gaia_dr1_aux.gaia_source_2mass_xm where source_id=' + str(TGASid)
	JmagLens,kMagLens,parallaxLens,gMagLens = sqlutil.get(querystringLens,db='wsdb',host='cappc127.ast.cam.ac.uk', user='peter_mcgill', password='Ln3g.wsk')

	#Calculate colour and absolute mag for lens 
	absMagLens = gMagLens + 5 * (np.log10(parallaxLens / 1000.0) + 1)
	colorLens = JmagLens - kMagLens


	#Get colour mag for 10000 source form tgas
	querystring2mass= 'select j_m,k_m,parallax,phot_g_mean_mag from gaia_dr1_aux.gaia_source_2mass_xm where parallax > 0.0001 and j_m > 0.0001 and k_m > 0.0001 limit 100000'
	Jmag,kMag,parallax,gMag = sqlutil.get(querystring2mass,db='wsdb',host='cappc127.ast.cam.ac.uk', user='peter_mcgill', password='Ln3g.wsk')


	#Calculate colour and absolute mags
	color = Jmag-kMag
	absmag = gMag + 5 * (np.log10(parallax / 1000.0) + 1)

	#plot 2dhist with scatter of lens overlaid
	plt.hist2d(color,absmag,bins=300,range=[[-0.1,1.25],[-2,13.5]],cmap='plasma')
	plt.scatter([colorLens],[absMagLens],color='white',label='TGAS Lens')
	plt.gca().invert_yaxis()
	plt.xlabel(r'J-K [2MASS]')
	plt.ylabel('Absolute G [mag]')
	plt.legend()
	plt.savefig(filename,dpi=200)
	plt.clf()
	
	

doctest.testmod()
