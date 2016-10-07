import os 
import glob      # List all files in a dir.
import numpy
import shutil    # Copy a file in another dir.
import scipy.integrate
import matplotlib.pyplot as plt


# XXX TODO XXX

# XXX HOWTOUSE XX
# Remove old infos_*EM*.dat
# Open IPython console then : import analysis_tools_sim20160929; d = analysis_tools_sim20160929.read_GWinfos(); ab = analysis_tools_sim20160929.compute_A_BAND(d); f = analysis_tools_sim20160929.compute_F_75_2000(d, ab); analysis_tools_sim20160929.write_EMinfos(d, ab, f)


def read_GWinfos():
    """ 
    Reads the 'infos_recoveredGWsources.dat' file which contains the infos from simulated GW events.
    """
    filename = 'infos_recoveredGWsources_sim20160929.dat'
    lines = []

    print 'Reading {} content...'.format(filename)

    with open(filename, 'r') as datafile:

        for line in datafile:

            if line.startswith('#'): 
                continue

            GPStime, date, m1, m2, D, SNR, RA, dec, iota, skymap, Egw = line.split(';')

            m1 = float(m1)
            m2 = float(m2)
            D = float(D)
            SNR = float(SNR)
            RA = float(RA) 
            dec = float(dec)
            Egw = float(Egw)
            iota = float(iota)

            data = (GPStime, date, m1, m2, D, SNR, RA, dec, iota, skymap, Egw)
            lines.append(data)

    return numpy.array(lines,
                       dtype=[('GPStime', '|S20'), ('date', '|S30'), ('mass1', float), ('mass2', float),
                              ('distances', float), ('SNR', float), ('RA', float),
                              ('dec', float), ('inclination', float), ('file_skymap', '|S40'), ('Egw', float)])
 
def write_EMinfos(data, A_band, F_75_2000):
    """
    Create a EM data file of the recovered events named 'infos_recoveredEMsources.dat'
    """
    filename = 'infos_recoveredEMsources_sim20160929.dat'
    lines = []
    
    GPStime = data['GPStime']
    date = data['date']
    m1 = data['mass1']
    m2 = data['mass2']
    dist = data['distances']
    snr = data['SNR']
    RA = data['RA']
    dec = data['dec']
    iota = data['inclination']
    skymaps = data['file_skymap']
    egw = data['Egw']

    print 'Creating EM infos data file {}...'.format(filename)
    
    
    with open(filename, 'a') as datafile:

        print >> datafile, '# GPStime date mass1 mass2 distance SNR RA dec inclination skymap Egw A_band F_75_2000'
        print >> datafile, '# s    Msun Msun Mpc   rad rad rad   erg ph/s/cm2/keV erg/s/cm2'

        for gps_, date_, m1_, m2_, d_, snr_, RA_, dec_, iota_, sky_, egw_, aband_, f_  in zip(GPStime, date, m1, m2, dist, snr, RA, dec, iota, skymaps, egw, A_band, F_75_2000):

            print >> datafile, gps_+' '+date_+' '+str(m1_)+' '+str(m2_)+' '+str(d_)+' '+str(snr_)+' '+str(RA_)+' '+str(dec_)+' '+str(iota_)+' '+sky_+' '+str(egw_)+' '+str(aband_)+' '+str(f_)
       

def plot_skydistrib(data):
    """ 
    Plot sky distribution of the simulated GW events.

    Input : data, is a numpy recorded array produced by the read function.
    """
    RA = map(float, data['RA'] - numpy.pi)
    DEC = map(float, data['dec'])

    assert len(RA) == len(DEC) ,\
        "Non matching RA and DEC."

    print 'Plotting sky distribution...'

    plt.subplot(111, projection='mollweide')
    plt.plot(RA, DEC, '+k')
    plt.xlabel('Right ascension')
    plt.ylabel('Declination')
    plt.title('Skymap of the GW events simulated in sim_20160929')
    plt.grid()

    plt.show()
    #plt.savefig('skymap_plot.png')

def plot_snrdistrib(data):
    """                                                                                     
    Plot sky distribution of the simulated GW events.                                                                                                                                                                                                                 
    Input : data, is a numpy recorded array produced by the read function.                                
    """
    snr = data['SNR']

    print 'Printing the snr distribution...'

    plt.plot((4., 4.), (0., 180.), '--k', label='SNR threshold = 4') # Show the threshold in SNR in ECM simulation.
    plt.hist(snr, bins=30)
    plt.xlabel('SNR')
    plt.ylabel('counts')
    plt.title('Histogram showing the SNR distribution for sim_20160929')
    plt.legend()
    plt.grid()

    #plt.show()
    plt.savefig('snr_distrib.png')

def plot_inclinationdistrib(data):
    """                                                                                                                 
    Plot sky distribution of the simulated GW events.                                                                   

    Input : data, is a numpy recorded array produced by the read function.                                                   
    """
    iota = data['inclination']

    print 'Printing the inclination distribution...'

    plt.plot((0., 0.), (0., 60.), '--k', label='Lowe bound : 0 rad')                   # Draw vertical line at iota = 0 rad
    plt.plot((numpy.pi, 2*numpy.pi), (0., 60.), '--k', label='Upper bound : 2 pi  rad')  # Draw vertical line at iota = 2 pi rad
    plt.hist(iota, bins=30)
    plt.xlabel('iota (rad)')
    plt.ylabel('counts')
    plt.title('Histogram showing the inclination distribution for sim_20160929')
    plt.legend()
    plt.grid()

    #plt.show()
    plt.savefig('inclination_distrib.png')

def extract_skymaps_recoveredGWevents():
    """                                                                                                                   
    Gather all skymaps associated to recovered events into a repository named 'skymaps_sim20160929'.                         
    """ 
    print 'Extracting skymaps...'

    # Create final directory which will contain all the skymaps.
    dirname = 'skymaps_sim20160929'
    if not os.path.exists(dirname):
        os.mkdir(dirname)

    dst = '/afs/in2p3.fr/home/p/pbacon/private/INTEGRAL_work/sim_20160929/skymaps_sim20160929/'

    # Loop over /mdcX.
    for directory in glob.glob('./mdc*'):
        src = glob.glob('/afs/in2p3.fr/home/p/pbacon/private/INTEGRAL_work/sim_20160929'+directory[1:]+'/*.fits.gz')

        # Loop over all skymaps in a given mdc.
        for s in src:
       
            # Copy skymap to the desired directory.
            src_dir = '/afs/in2p3.fr/home/p/pbacon/private/INTEGRAL_work/sim_20160929'+directory[1:]+'/'
            dst_dir = dst
            src_file = s
            shutil.copy(src_file, dst_dir)

            # Rename the copied skymap to identify its belonging mdc.
            dst_file = dst+s.split('/')[-1]
            new_dst_filename = dst_dir+directory[2:]+'_'+s.split('/')[-1]
            os.rename(dst_file, new_dst_filename)
        

    # Count numbers of skymaps.
    N = len(glob.glob(dst+'*'))
    print 'Found {} skymaps.'.format(N)

def BAND_spectrum_WN(E, alpha, beta, Epeak):
    """
    BAND spectrum specified in Eq (2) of "The FERMI GBM GRB spectral catalog : Four years of data" - D. Gruber et al. (2014)
    withour A_band factor. 

    Inputs : 
    -----------------------------------
    * E : energy (in kev)
    * alpha : alpha coefficient in BAND spectrum
    * beta : beta coefficient in BAND spectrum
    * Epeak : peak energy in BANd spectrum (in keV)
    """

    Ethreshold = (alpha - beta) * Epeak / (alpha + 2.)
    
    if (E >= 0) and (E < Ethreshold):
        return numpy.power(E/100., alpha) * numpy.exp(-((alpha + 2.)*E/Epeak))
    elif (E >= 0) and (E >= Ethreshold):
        return numpy.power(E/100., beta) * numpy.exp(beta - alpha) * numpy.power(((alpha - beta)*Epeak)/(100.*(alpha + 2.)), alpha - beta)
    else:
        print 'Non positive energies are unphysical quantities ! : found E = {} keV.'.format(E)


def func_NE(E, alpha, beta, Epeak):
    """
    Function BAND spectrum
    """
    return BAND_spectrum_WN(E, alpha, beta, Epeak)

def func_E_NE(E, alpha, beta, Epeak):
    """
    Function E * N(E), where N(E) is the BAND spectrum.
    """
    return E * BAND_spectrum_WN(E, alpha, beta, Epeak)

def compute_A_BAND(data):
    """
    Compute the A_band factor that should be placed with the BAND spectrum.
    """
    print 'Computing A band factor...'

    Egw = data['Egw']                # in erg
    D = data['distances'] * 3.086e24 # in cm

    # Compute the total EM energy emitted during the coalescence, provided
    # the conversion factor is ETA = 10 %.
    ETA = .10
    Eem = ETA * Egw

    # EM luminosity of the events provided the emission lasts DURATION seconds.
    DURATION = 1  # in s
    Lem = Eem / DURATION # in erg/s

    # Computing EM flux 
    Fem = Lem / (4*numpy.pi*D*D) # in erg/s/cm^{-2}

    # Computing f_total[keV.erg]
    LOWER_E = 1        # in keV
    UPPER_E = 1e4      # in keV
    ALPHA_BAND = - 0.5 
    BETA_BAND = - 2.25
    E_PEAK = 800       # in keV

    f_total, f_total_err = scipy.integrate.quad(func_E_NE, LOWER_E, UPPER_E, (ALPHA_BAND, BETA_BAND, E_PEAK)) # in keV.erg

    # Return A_band (in ph/s/cm2/keV)
    return  Fem / f_total
    
def compute_F_75_2000(data, A_band_arr):
    """
    Compute the f_75_2000 factor, then deduce F_75_2000 from A_band. 
    """
    print 'Computing F_75_2000...'

    Egw = data['Egw']                # in erg
    A_band = A_band_arr               # in ph/s/cm2/keV

    # Computing f_75_2000. The 75-2000 keV is the energy range of INTEGRAL.
    LOWER_E = 75        # in keV                                                                                                                         
    UPPER_E = 2e3      # in keV                                                                                                                          
    ALPHA_BAND = - 0.5
    BETA_BAND = - 2.25
    E_PEAK = 800       # in keV
    
    f_75_2000, f_75_2000_err = scipy.integrate.quad(func_E_NE, LOWER_E, UPPER_E, (ALPHA_BAND, BETA_BAND, E_PEAK)) # in keV.erg
    
    # Return F_75_2000
    return A_band * f_75_2000

def mapping_iota(iota0):
    """
    Maps ECM and Barbara s conventions on inclination angle.
    """
    if ((0. <= iota0) and (iota0 <= numpy.pi/2)):
        return abs(iota0)
    elif ((numpy.pi/2 < iota0) and (iota0 <= 3*numpy.pi/2)):
        return abs(numpy.pi - iota0)
    elif ((3*numpy.pi/2 <= iota0) and (iota0 <= 2*numpy.pi)):
        return abs(2*numpy.pi - iota0)
    else:
        return 0.
