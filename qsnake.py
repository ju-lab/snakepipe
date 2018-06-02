#!/usr/bin/env python3
# From https://github.com/broadinstitute/viral-ngs/blob/master/pipes/Broad_LSF/run-pipe.sh
# usage: snakemake -p [final_output] --cluster 'python qsnake.py' -j 24
import os
import sys
import re
from snakemake.utils import read_job_properties

LOGDIR = 'logs'
jobscript = sys.argv[1]
props = read_job_properties(jobscript)

# # set up job name, project name
jobname = "{rule}_job{jobid}".format(rule=props["rule"], jobid=props["jobid"])

# cmdline = f'qsub -N {jobname} '

# # log file output
# if "-N" not in props["params"].get("LSF", ""):

cmdline = 'qsub '
cmdline += f"-o {LOGDIR}/{jobname}.out -e {LOGDIR}/{jobname}.err "

# pass memory and cpu resource request to LSF
ncpus = props['threads']
mem = props.get('resources', {}).get('mem_mb')
if mem:
    cmdline += f'-l nodes=1:ppn={ncpus} -l mem={mem}M '
    
else:
    cmdline += f'-l nodes=1:ppn={ncpus} '

# # figure out job dependencies
# dependencies = set(sys.argv[1:-2])
# if dependencies:
#     cmdline += "-W depend=afterok'{}' ".format(" && ".join(dependencies))

# the actual job
cmdline += jobscript

# call the command
os.system(cmdline)