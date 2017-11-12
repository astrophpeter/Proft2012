import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import numpy as np

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


def plot_len_traj(lens,sourceRa,sourceDec):
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
def plt_lens_env(lens,sourceRa, sourceDec):
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

	##### 2MASS ##### - not currently being used - DSS more complete. 
	filter_2m = 'j' # 2mass filter to use

	metadata_url = 'http://irsa.ipac.caltech.edu/ibe/search/twomass/allsky/allsky?POS={0},{1}'.format(sourceRa[0],sourceDec[0])
	metadata = ascii.read(metadata_url, Reader=ascii.ipac.Ipac)
    	
	# selecting the first obs in the requested filter in case there is more than one
	target_obs = metadata[metadata['filter']==filter_2m][0]
    
	params = { 'ordate': target_obs['ordate'],
		'hemisphere': target_obs['hemisphere'],
		'scanno': target_obs['scanno'],
		'fname': target_obs['fname'],
		'ra': sourceRa,
		'de': sourceDec,
		'imsize': imsize }
        
	#DSS Search
	url = "http://stdatu.stsci.edu/cgi-bin/dss_search?v=poss2ukstu_red&r={ra}&d={de}&e=J2000&h={imsize}&w={imsize}&f=fits".format(**params)
	hdulist = fits.open(url)
	
	#Find the time the image was taken YYYY-MM-DD
	timeString  = hdulist[0].header['DATE-OBS'][:10]

    	# Get the position of the lens at the time the image was taken
	time = lens.datetime_to_jyTCB(timeString)
	coord = lens.get_eq_coords_at_epoch(time)
    
	#Get ra and dec of lens at the image time (raCos(dec) is converted into ra)
	raLensimag = coord[0] / np.cos(np.deg2rad(lensDec[0]))
	decLensimag = coord[1]
    
	# Find the postition of the lens in the future so its trajectory can be plotted
	timeend = testlens1.datetime_to_jyTCB('2040-01-01')
	coordend = testlens1.get_eq_coords_at_epoch(timeend)
    	raLensimagend = coordend[0] / np.cos(np.deg2rad(lensDec[0]))
	decLensimagend = coordend[1]
    
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
    
	# Add some extra info about the source and lens to the plot.
	fig.add_label(0.7,0.98,Diststr,relative=True)
	fig.add_label(0.7,0.94,Timestr,relative=True)
	fig.add_label(0.7,0.90,Lensstr,relative=True)
	fig.add_label(0.7,0.86,Sourcestr,relative=True)
    
	filename = 'source_env/TGAS_' + str(lens.getId()) + '.png'
	fig.save(filename,dpi=200)

def get_cut_out(lens, source):



	
