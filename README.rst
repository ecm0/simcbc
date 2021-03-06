
Instructions to run the complete simulation
===========================================

1. Generate a simulated galaxy catalog

nohup octave -q generate_galaxies.m

Produces galaxies.txt and galaxies.mat

2. Generate a catalog of mergers at defined R_merger rate [23.5 Myr^-3 per MW galaxy]
   [note: uses the Glade catalog for mergers with distance < 200 Mpc]

nohup octave -q associate_merger_to_gal.m

Produces mergers.txt and mergers.mat

3. Generate XML tables for further processing with Bayestar/little_hope

ipython generate_xml_tables.py mergers.txt

Produces mdc.xml and psd.xml

Generates merger time, inclination, phase at coalescence, polarization angle

Select events with at least one single-detector SNR > threshold [= 4]

Power spectral densities are taken from the "Prospects" paper
http://arxiv.org/abs/1304.0670

aLIGO (mid-stage, low-sensitivity curve -- 80 Mpc)
adV (early-stage, low-sensitivity curve -- 20 Mpc)

4. Run little_hope on laptop
   
bayestar_sim_to_tmpltbank mdc.xml -o templates.xml
bayestar_littlehope --detector H1 --detector L1 --min-triggers 2 --snr-threshold 4.0 --reference-psd psd.xml --template-bank templates.xml --waveform "TaylorF2threePointFivePN" mdc.xml -o coinc.xml
bayestar_localize_coincs --waveform "TaylorF2threePointFivePN" --f-low 30 coinc.xml
bayestar_plot_allsky 0.toa_phoa_snr.fits.gz --contour 90 --radec long lat -o skymap.png

5. Run little_hope at CC Lyon

   Generates templates for each MDC:
   for d in `seq 1 YY`; do ./my_bayestar_sim_to_tmpltbank sim_XXXXX/mdc$d.xml -o templates$d.xml; done
   
   Set up the environment [see .bashrc below]
   Install healpy (v 1.9.0) [see below]

   Store all MDC and template files produced at step 3. into a single folder sim_XXXX/orig_files
   Generate batch jobs with:   
        createjobs.py sim_XXXX/orig_files sim_XXXX "-l s_rss=10G" 3
   This creates folders sim_XXXX/mdcXX and a script sim_XXXX/jobsubmission.sh
   Submit jobs to the queue with ./jobsubmission.sh
   To check that jobs are in the queue: qstat
   Note: skymap computation does not work at CC Lyon (skymap files are
   produced but are empty files).
   To create skymaps, run ./createskymaps.sh on your computer
   To create a sim summary, run ./create_datafile_GWsources_infos.sh


export PATH="$HOME/local/bin:${THRONG_DIR}/bin:$PATH"
export LD_LIBRARY_PATH="$HOME/local/lib:${THRONG_DIR}/lib:$LD_LIBRARY_PATH"
export PKG_CONFIG_PATH="${THRONG_DIR}/lib/pkgconfig"

export PATH="/usr/local/python/python-2.7/bin/:/afs/in2p3.fr/throng/virgo/virgoApp/lalsuite/v6r30p3/amd64_sl6/bin:$PATH"
export LD_LIBRARY_PATH="/afs/in2p3.fr/home/throng/virgo/virgoApp/lalsuite/v6r30p3/amd64_sl6/lib:$LD_LIBRARY_PATH"
export PYTHONPATH="/afs/in2p3.fr/home/e/ecm/local/lib/python2.7/site-packages:/usr/local/python/python-2.7/lib/site-packages:/afs/in2p3.fr/home/throng/virgo/virgoApp/lalsuite/v6r30p3/amd64_sl6/lib/python2.7/site-packages"

curl -O https://pypi.python.org/packages/source/h/healpy/healpy-1.9.0.tar.gz
tar xzf healpy-1.9.0.tar.gz
cd healpy-1.9.0
python setup.py install --user
