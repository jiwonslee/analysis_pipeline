#! /usr/bin/env python2.7

"""Convert VCF to GDS"""
import TopmedPipeline
import sys
import os
from argparse import ArgumentParser
from copy import deepcopy

description = """
Convert VCF to GDS with the following steps:
1) Convert per-chromosome VCF files to GDS
2) Merge genotypes from per-chromosomes GDS files into a combined file
3) Assign unique variant id from merged file to per-chromosome GDS files
"""

parser = ArgumentParser(description=description)
parser.add_argument("config_file", help="configuration file")
#parser.add_argument("--cluster_type", default="UW_Cluster",
#                    help="type of compute cluster environment [default %(default)s]")
#parser.add_argument("--cluster_file", default=None,
#                    help="json file containing options to pass to the cluster")
parser.add_argument("--verbose", action="store_true", default=False,
                    help="enable verbose output to help debug")
parser.add_argument("-c", "--chromosomes", default="1-23",
                    help="range of chromosomes [default %(default)s]")
#parser.add_argument("-n", "--ncores", default="1-8",
#                    help="number of cores to use; either a number (e.g, 1) or a range of numbers (e.g., 1-4) [default %(default)s]")
#parser.add_argument("-e", "--email", default=None,
#                    help="email address for job reporting")
#parser.add_argument("--print_only", action="store_true", default=False,
#                    help="print qsub commands without submitting")
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
#config_dir = os.path.join(pipeline, configfile)
config_dir = configfile

configdict = TopmedPipeline.readConfig(configfile)

for subdir in ['config', 'log']:
    if not os.path.exists(configdict['output_file'] + '/' + subdir):
        os.mkdir(configdict['output_file'] + '/' + subdir)

#configdict = TopmedPipeline.directorySetup(configdict, subdirs=["config", "log"])

job = "vcf2gds"

rscript = os.path.join(pipeline, "R", job + ".R")

# parsing bcf files relies on streaming bcftools output, so can't run in parallel
if os.path.splitext(configdict["vcf_file"])[1] == ".bcf":
    ncores = None


chrom_string = TopmedPipeline.parseChromosomes(chromosomes)
chrom_list = chrom_string.split(' ')

for chrom in chrom_list:
    if os.path.isfile(configdict['gds_file'].replace('chr ', 'chr'+chrom)) == False:
        cmd =  " ".join(['bsub -q big -n 4', "-R 'rusage[mem=45000]'",
                 'Rscript', rscript,
                 config_dir, '--chromosome ' + chrom])
        print(cmd)
        os.system(cmd)
#jobid = cluster.submitJob(job_name=job, cmd=driver, args=["-c", rscript, configfile], array_range=chromosomes, request_cores=ncores)

job = "merge_gds"

rscript = os.path.join(pipeline, "R", job + ".R")

configdict["chromosomes"] = TopmedPipeline.parseChromosomes(chromosomes)
configfile = configdict['output_file'] + '/config/'+ configdict["config_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(configdict, configfile)
#config_dir = os.path.join(pipeline, configfile)
chrom_list = configdict["chromosomes"].split(" ")
flag = False
for chrom in chrom_list:
    if os.path.isfile(configdict['gds_file'].replace('chr ', 'chr'+chrom)) == False:
        flag = True
    if os.path.isfile(configdict["merged_gds_file"]) == True:
        flag = True
if flag == False:
    cmd =  " ".join(['bsub -q big', #"-R 'rusage[mem=45000]'",
#                 '-L /bin/bash',
#                 'R -q --vanilla --args',
#                 config_dir,'<'+rscript+'> /dev/null']
                    'Rscript', rscript, configfile])
    print(cmd)
    os.system(cmd)
else:
    print('Run this script again to merge gds files')


job = "unique_variant_ids"

rscript = os.path.join(pipeline, "R", job + ".R")

if os.path.isfile(configdict["merged_gds_file"]) == True:
    cmd =  " ".join(['bsub -q big', #"-R 'rusage[mem=45000]'",
#                 '-L /bin/bash',
#                 'R -q --vanilla --args',
#                 config_dir,'<'+rscript+'> /dev/null']
                    'Rscript', rscript, configfile])
    print(cmd)
    os.system(cmd)
#cluster.submitJob(job_name="cleanup", cmd=os.path.join(pipeline, "cleanup.sh"), holdid=[jobid], print_only=print_only)

