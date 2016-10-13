from __future__ import division

import os
import sys
import lal
import numpy as np
from lal import gpstime

from lal.lal import MSUN_SI as LAL_MSUN_SI   # kg -- mass of the Sun
from lal.lal import C_SI as LAL_C_SI         # m s^-1 -- speed of light

from glue.ligolw import table as ligolw_table
from glue.ligolw import utils as ligolw_utils
from glue.ligolw import lsctables

from lalinference.bayestar import ligolw as ligolw_bayestar

import numpy as np
import healpy as hp
from lalinference import fits

file1 = sys.argv[1]
(path1, filename1) = os.path.split(file1)

file2 = sys.argv[2]
(path2, filename2) = os.path.split(file2)

xmldoc1 = ligolw_utils.load_filename(file1,
                                     contenthandler=ligolw_bayestar.LSCTablesContentHandler)

xmldoc2 = ligolw_utils.load_filename(file2,
                                     contenthandler=ligolw_bayestar.LSCTablesContentHandler)

sim_inspiral_table = ligolw_table.get_table(xmldoc1,
                                            lsctables.SimInspiralTable.tableName)

coinc_inspiral_table = ligolw_table.get_table(xmldoc2, 
                                              lsctables.CoincInspiralTable.tableName)

# Store infos in a .dat file.                                                                 
with open('summary.txt', 'w') as file :

    print >> file, "# GPS date mass1 mass2 dist SNR RA dec inclination skymap Egw"
    #print "# GPS mass1 mass2 dist SNR RA dec inclination"                                                                                 

    for coinc in coinc_inspiral_table:

        for sim_inspiral in sim_inspiral_table:
            if abs(sim_inspiral.geocent_end_time-coinc.end_time) < 2:
                break
            else:
                sim_inspiral=None

        if not sim_inspiral:
            continue

        end_time = lal.LIGOTimeGPS(sim_inspiral.geocent_end_time, sim_inspiral.geocent_end_time_ns)
        (RA,dec) = sim_inspiral.ra_dec
        event_id = str(coinc.coinc_event_id)
        
        # Compute Egw - Take BNS Egw estimation from "How loud are neutron star mergers ?" - S. Bernuzzi et al. (2016)
        Egw = 1.5e-2 * (sim_inspiral.mass1 + sim_inspiral.mass2) * LAL_MSUN_SI * (LAL_C_SI**2) * 1e7 # in erg

        # Compute search area
        skymap_filename = "{}/{}.toa_phoa_snr.fits.gz".format(path1,event_id[-1])

        skymap, metadata = fits.read_sky_map(skymap_filename, nest=None)
        nside = hp.npix2nside(len(skymap))

        # Convert sky map from probability to probability per square degree.
        probperdeg2 = skymap / hp.nside2pixarea(nside, degrees=True)
        deg2=hp.nside2pixarea(nside, degrees=True)

        indices = np.argsort(-skymap)
        region = np.empty(skymap.shape)
        region[indices] = 100 * np.cumsum(skymap[indices])
        mylist=region[indices]

        # Compute search areas
        area90 = area50 = 0.0
        for pix_prob in mylist:
            if pix_prob < 90.0:
                area90 += deg2
            if pix_prob < 50.0:
                area50 += deg2
    
        # Store data in output file.
        current_infos =  '{gps};"{date}";{mass1};{mass2};{dist};{snr};{RA};{dec};{inclination};{skymap};{egw};{area90};{area50}'.format(gps=end_time,
                                                                                                                      date=lal.gpstime.gps_to_utc(end_time).isoformat(' '),
                                                                                                                      mass1=sim_inspiral.mass1,
                                                                                                                      mass2=sim_inspiral.mass2,
                                                                                                                      dist=sim_inspiral.distance,                                      
                                                                                                                      snr=coinc.snr,
                                                                                                                      RA=RA,
                                                                                                                      dec=dec,
                                                                                                                      inclination=sim_inspiral.inclination,
                                                                                                                      skymap=skymap_filename,
                                                                                                                      egw=Egw,
                                                                                                                      area90=area90,
                                                                                                                      area50=area50)
        print >> file, current_infos
