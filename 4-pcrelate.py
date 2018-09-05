#! /usr/bin/env python2.7

"""PC-Relate"""

import TopmedPipeline
import sys
import os
from argparse import ArgumentParser
from copy import deepcopy

description = """
PC-Relate
"""
parser = ArgumentParser(description=description)
parser.add_argument("config_file", help="configuration file")
#parser.add_argument("--cluster_type", default="UW_Cluster",
#                    help="type of compute cluster environment [default %(default)s]")
#parser.add_argument("--cluster_file", default=None,
#                    help="json file containing options to pass to the cluster")
parser.add_argument("--verbose", action="store_true", default=False,
                    help="enable verbose output to help debug")
#parser.add_argument("-e", "--email", default=None,
#                    help="email address for job reporting")
#parser.add_argument("--print_only", action="store_true", default=False,
#                    help="print cluster commands without submitting")
#parser.add_argument("--version", action="version",
#                    version="TopmedPipeline "+TopmedPipeline.__version__,
#                    help="show the version number and exit")
args = parser.parse_args()

configfile = args.config_file
#cluster_file = args.cluster_file
#cluster_type = args.cluster_type
#email = args.email
#print_only = args.print_only
verbose = args.verbose

#version = "--version " + TopmedPipeline.__version__

#cluster = TopmedPipeline.ClusterFactory.createCluster(cluster_type, cluster_file, verbose)

pipeline = os.path.dirname(os.path.abspath(sys.argv[0]))
driver = os.path.join(pipeline, "runRscript.sh")

configdict = TopmedPipeline.readConfig(configfile)
for subdir in ["config","data", "log", "plots"]:
    if not os.path.exists(configdict['output_file'] + '/' + subdir):
        os.mkdir(configdict['output_file'] + '/' + subdir)
#configdict = TopmedPipeline.directorySetup(configdict, subdirs=["config", "data", "log", "plots"])

job = "pcrelate"

rscript = os.path.join(pipeline, "R", job + ".R")

config = deepcopy(configdict)
config["out_prefix"] = configdict['output_file'] + '/data/' + configdict["data_prefix"]
configfile =  configdict['output_file'] + '/config/' + configdict["config_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(config, configfile)



cmd =  " ".join(['bsub -q short', #"-R 'rusage[mem=45000]'",
#                 '-L /bin/bash',
#                 'R -q --vanilla --args',
#                 config_dir,'<'+rscript+'> /dev/null']
                    'Rscript', rscript, configfile])
print(cmd)
os.system(cmd)


#jobid = cluster.submitJob(job_name=job, cmd=driver, args=[rscript, configfile, version], email=email, print_only=print_only)

job = "kinship_plots"

rscript = os.path.join(pipeline, "R", job + ".R")

config = deepcopy(configdict)
config["kinship_file"] =  configdict['output_file'] + '/data/' +configdict["data_prefix"] + "_pcrelate.gds"
config["kinship_method"] = "pcrelate"
config["out_file_all"] =  configdict['output_file'] + '/plots/' +configdict["plots_prefix"] + "_kinship_all.pdf"
config["out_file_cross"] =  configdict['output_file'] + '/plots/' +configdict["plots_prefix"] + "_kinship_cross.pdf"
config["out_file_study"] =  configdict['output_file'] + '/plots/' +configdict["plots_prefix"] + "_kinship_study.pdf"
configfile = configdict['output_file'] + '/config/' +configdict["config_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(config, configfile)

if os.path.isfile(config['kinship_file']):
    cmd =  " ".join(['bsub -q short', #"-R 'rusage[mem=45000]'",
#                 '-L /bin/bash',
#                 'R -q --vanilla --args',
#                 config_dir,'<'+rscript+'> /dev/null']
                    'Rscript', rscript, configfile])
    print(cmd)
    os.system(cmd)
else:
    print('kinship file has not completed yet; plots not generated')
#jobid = cluster.submitJob(job_name=job, cmd=driver, args=[rscript, configfile, version], holdid=[jobid], email=email, print_only=print_only)


#cluster.submitJob(job_name="cleanup", cmd=os.path.join(pipeline, "cleanup.sh"), holdid=[jobid], print_only=print_only)
