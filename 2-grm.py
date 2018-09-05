#! /usr/bin/env python2.7

"""GRM"""

import TopmedPipeline
import sys
import os
from argparse import ArgumentParser
from copy import deepcopy

description = """
Genetic Relationship Matrix (GRM)
"""
parser = ArgumentParser(description=description)
parser.add_argument("config_file", help="configuration file")
parser.add_argument("-c", "--chromosomes", default="1-22",
                    help="range of chromosomes [default %(default)s]")
#parser.add_argument("--cluster_type", default="UW_Cluster",
#                    help="type of compute cluster environment [default %(default)s]")
#parser.add_argument("--cluster_file", default=None,
#                    help="json file containing options to pass to the cluster")
parser.add_argument("--verbose", action="store_true", default=False,
                    help="enable verbose output to help debug")
#parser.add_argument("-n", "--ncores", default="1-8",
#                    help="number of cores to use; either a number (e.g, 1) or a range of numbers (e.g., 1-4) [default %(default)s]")
#parser.add_argument("-e", "--email", default=None,
#                    help="email address for job reporting")
#parser.add_argument("--print_only", action="store_true", default=False,
#                    help="print cluster commands without submitting")
#parser.add_argument("--version", action="version",
#                    version="TopmedPipeline "+TopmedPipeline.__version__,
#                    help="show the version number and exit")
args = parser.parse_args()

configfile = args.config_file
chromosomes = args.chromosomes
#cluster_file = args.cluster_file
#cluster_type = args.cluster_type
#ncores = args.ncores
#email = args.email
#print_only = args.print_only
verbose = args.verbose

#version = "--version " + TopmedPipeline.__version__

#cluster = TopmedPipeline.ClusterFactory.createCluster(cluster_type, cluster_file, verbose)

pipeline = os.path.dirname(os.path.abspath(sys.argv[0]))
#driver = os.path.join(pipeline, "runRscript.sh")

configdict = TopmedPipeline.readConfig(configfile)
for subdir in ["config","data", "log", "plots"]:
    if not os.path.exists(configdict['output_file'] + '/' + subdir):
        os.mkdir(configdict['output_file'] + '/' + subdir)

#configdict = TopmedPipeline.directorySetup(configdict, subdirs=["config", "data", "log"])

chrom_string = TopmedPipeline.parseChromosomes(chromosomes)
chrom_list = chrom_string.split(' ')
job = "grm"

rscript = os.path.join(pipeline, "R", job + ".R")

config = deepcopy(configdict)
config["out_file"] = configdict['output_file'] + '/data/' + configdict["data_prefix"] + "_grm_chr .gds"
configfile = configdict['output_file'] + '/config/' + configdict["config_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(config, configfile)

#jobid = cluster.submitJob(job_name=job, cmd=driver, args=["-c", rscript, configfile, version], request_cores=ncores, array_range=chromosomes, email=email, print_only=print_only)
for chrom in chrom_list:
    if os.path.isfile(config["gds_file"].replace('chr ', 'chr' + chrom)) == True and os.path.isfile(config["out_file"].replace('chr ', 'chr'+ chrom)) == False:
        cmd =  " ".join(['bsub -q big -R "rusage[mem=20000]" -n 4', #"-R 'rusage[mem=45000]'",
#                 '-L /bin/bash',
#                 'R -q --vanilla --args',
#                 config_dir,'<'+rscript+'> /dev/null']
                    'Rscript', rscript,  configfile, '--chromosome ' + chrom])
        print(cmd)
        os.system(cmd)
    elif os.path.isfile(config["gds_file"].replace('chr ', 'chr' + chrom)) == False:
        print('gds file for chromosome' + chrom + 'is missing')


job = "grm_combine"

rscript = os.path.join(pipeline, "R", job + ".R")

config = dict()
config["chromosomes"] = TopmedPipeline.parseChromosomes(chromosomes)
config["in_file"] = configdict['output_file'] + '/data/' + configdict["data_prefix"] + "_grm_chr .gds"
config["out_file"] = configdict['output_file'] + '/data/' + configdict["data_prefix"] + "_grm.gds"
configfile = configdict["config_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(config, configfile)

exit_flag = False
if os.path.isfile(config["out_file"]) == True:
    exit_flag = True
if exit_flag == True:
    print('Merged GRM file has already been created')
    exit()
for chrom in chrom_list:
    if os.path.isfile(config["in_file"].replace('chr ', 'chr'+ chrom)) == False:
        exit_flag = True
if exit_flag == True:
    print('Some chromosome files are not ready yet; run script again later to merge')
    exit()
else:
    cmd =  " ".join(['bsub -q big', 
                    'Rscript', rscript,  configfile])
    os.system(cmd)
#jobid = cluster.submitJob(job_name=job, cmd=driver, args=[rscript, configfile, version], holdid=[jobid], email=email, print_only=print_only)


#cluster.submitJob(job_name="cleanup", cmd=os.path.join(pipeline, "cleanup.sh"), holdid=[jobid], print_only=print_only)
