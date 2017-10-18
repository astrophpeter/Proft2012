import sqlutil as sqlutil 
import matplotlib.pyplot as plt
import numpy as np
from lens import Lens
from matplotlib import rcParams
rcParams['axes.titlepad'] = 20 



#Getting High proper motion lenses.
lensRa, lensDec, id,pmra,pmdec,ref_epoch = sqlutil.get('select ra, dec, source_id,pmra,pmdec,ref_epoch from gaia_dr1.tgas_source where POWER(pmra,2) + POWER(pmdec,2) > 1000000',
                       db='wsdb',host='cappc127.ast.cam.ac.uk', user='peter_mcgill', password='Ln3g.wsk')

#backgroud catalogue
sourceRa, sourceDec = sqlutil.get('select ra, dec from ppmxl.main limit 1000000',
                       #db='wsdb',host='cappc127.ast.cam.ac.uk', user='peter_mcgill', password='Ln3g.wsk')


for i range(0,len(lensRa)):
	lens = Lens(id[i],lensRa[i],lensDec[i],pmra[i],pmdec[i],ref_epoch[i])
	box = lens.get_lens_box()
	
	
	
