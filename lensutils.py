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



def get_cut_out(lens, source):



	
