from __future__ import division

import sys
import math
import os
from shutil import copy

from glue.ligolw import ligolw
from glue.ligolw import ilwd
from glue.ligolw import lsctables
from glue.ligolw import utils as ligolw_utils

LITTLEHOPE_HOME = os.getcwd()
LITTLEHOPE_OPTS2 = '--detector H1 --detector L1 --detector V1 --min-triggers 2 --snr-threshold 4.0 --trigger-window 1 --reference-psd {home}/psds_2016-17.xml --waveform "TaylorT4threePN"'.format(home=LITTLEHOPE_HOME)
LITTLEHOPE_OPTS3 = '--detector H1 --detector L1 --min-triggers 2 --snr-threshold 4.0 --trigger-window .1 --reference-psd {home}/psds_2016-17.xml --waveform "TaylorT4threePN"'.format(home=LITTLEHOPE_HOME)
CMD_LITTLEHOPE = '{home}/my_bayestar_littlehope {opts} --template-bank {simdir}/{template_file} {simdir}/{mdc_file} -o coinc.xml\n'

LOCALCOINCS_OPTS = '--waveform "TaylorT4threePN" --f-low 30'
CMD_LOCALCOINCS = 'bayestar_localize_coincs {opts} coinc.xml\n'
CMD_PLOT = 'for f in *.fits.gz; do bayestar_plot_allsky $f --contour 90 --radec 0.0 0.0 -o ${f%.*}.png; done\n'
CMD_JOBSUB = '/opt/sge/bin/lx24-amd64/qsub -N {} -o {} -e {} {} {}'

if __name__ == "__main__":

    # parse input args
    if len(sys.argv) <= 1:
        sys.exit()

    indir = sys.argv[1]
    
    # list contents of input dir
    if os.path.isdir(indir):
        files = os.listdir(indir)
    else:
        sys.exit()

    if len(sys.argv) > 2:
        out_prefix = sys.argv[2]
    else:
        out_prefix = ""

    if len(sys.argv) > 3:
        jobsub_opts = sys.argv[3]
    else:
        jobsub_opts = ""

    if len(sys.argv) > 4:
        num_ifo = sys.argv[4]
    else:
        num_ifo = 2

    # loop on coincX.xml
    submission_cmds = []
    select_mdc_xml = lambda f: f.lower().startswith('mdc') and f.lower().endswith('.xml')
    for my_file in filter(select_mdc_xml, files):
        
        # create job dirs
        simdir,_ = os.path.splitext(my_file)
        simdir = "{}/{}".format(prefix,simdir)
        if not os.path.exists(simdir):
            os.makedirs(simdir)

        copy(indir + '/' + my_file, simdir + '/')

        # retrieve corresponding template filename
        tmpl_file = my_file.replace("mdc","templates")
        copy(indir + '/' + tmpl_file, simdir + '/')
        
        fullpath_simdir = os.path.abspath(simdir)

        if num_ifo == 2:
            options = LITTLEHOPE_OPTS2
        else:
            options = LITTLEHOPE_OPTS3

        # create job scripts
        batch_filename = "{}/batch.sh".format(simdir)
        with open(batch_filename, "w") as batch_script:
            batch_script.write("#!/usr/bin/env bash\n")
            batch_script.write("cd {simdir}\n".format(simdir=fullpath_simdir))
            batch_script.write(CMD_LITTLEHOPE.format(mdc_file=my_file,
                                                     template_file=tmpl_file,
                                                     simdir=fullpath_simdir,
                                                     home=LITTLEHOPE_HOME,
                                                     opts=options))
            batch_script.write(CMD_LOCALCOINCS.format(simdir=fullpath_simdir,
                                                      opts=LOCALCOINCS_OPTS))
            # batch_script.write(CMD_PLOT)
            batch_script.write('echo "hello word!"')

        os.chmod(batch_filename, 0744)
        
        # create script with submission list
        submission_cmds.append(CMD_JOBSUB.format(simdir, 
                                                 "{}/{}.out".format(fullpath_simdir,simdir), 
                                                 "{}/{}.err".format(fullpath_simdir,simdir), 
                                                 jobsub_opts,
                                                 batch_filename))

    outfile = open('jobsubmission.sh', 'w')
    print>>outfile, "#!/usr/bin/env bash"
    for line in submission_cmds:
        print>>outfile, line

    os.chmod('jobsubmission.sh', 0744)
