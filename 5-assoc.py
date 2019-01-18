#! /usr/bin/env python2.7

"""Association tests"""

import TopmedPipeline
import sys
import os
import subprocess
from time import localtime, strftime
from argparse import ArgumentParser
from copy import deepcopy

description = """
Association tests
"""

default_segment_length = "10000"

parser = ArgumentParser(description=description)
parser.add_argument("assoc_type", choices=["single", "window", "aggregate"],
                    help="type of association test")
parser.add_argument("config_file", help="configuration file")
parser.add_argument("-c", "--chromosomes", default="1-23",
                    help="range of chromosomes [default %(default)s]")
parser.add_argument("--segment_length", default=default_segment_length,
                    help="segment length in kb [default %(default)s]")
parser.add_argument("--n_segments", default=None,
                    help="number of segments for the entire genome (overrides segment_length)")
#parser.add_argument("-e", "--email", default=None,
#                    help="email address for job reporting")
#parser.add_argument("--print_only", action="store_true", default=False,
#                    help="print qsub commands without submitting")
parser.add_argument("--verbose", action="store_true", default=False,
                    help="enable verbose output to help debug")
#parser.add_argument("--version", action="version",
#                    version="TopmedPipeline "+TopmedPipeline.__version__,
#                    help="show the version number and exit")
args = parser.parse_args()

assoc_type = args.assoc_type
configfile = args.config_file
chromosomes = args.chromosomes
segment_length = args.segment_length
n_segments = args.n_segments
#email = args.email
#print_only = args.print_only
verbose = args.verbose
#version = "--version " + TopmedPipeline.__version__

pipeline = os.path.dirname(os.path.abspath(sys.argv[0]))

configdict = TopmedPipeline.readConfig(configfile)

for subdir in ["config","data", "results", "log", "plots", "report"]:
    if not os.path.exists(configdict['output_file'] + '/' + subdir):
        os.mkdir(configdict['output_file'] + '/' + subdir)

# check type of association test - single-variant unrelated is handled differently
no_pcrel = "pcrelate_file" not in configdict or configdict["pcrelate_file"] == "NA"
no_grm = "grm_file" not in configdict or configdict["grm_file"] == "NA"
single_unrel = assoc_type == "single" and no_pcrel and no_grm

#######job: null model#############
if not single_unrel:
    job = "null_model"
    assocScript = "assoc_" + assoc_type
    if "test_type" in configdict:
        job = job + '_saige' if configdict["test_type"].lower() == "saige" else job
        assocScript = assocScript + '_saige' if configdict["test_type"].lower() == "saige" else assocScript
    logname = configdict['output_file'] + '/log/' + configdict['config_prefix'] + '_' + job + '.log' 
    rscript = os.path.join(pipeline, "R", job + ".R")
    config = deepcopy(configdict)
    config["out_file"] = configdict['output_file'] + '/data/' + configdict["data_prefix"] + "_null_model.RData"
    config["out_phenotype_file"] = configdict['output_file'] + '/data/' + configdict["data_prefix"] + "_phenotypes.RData"
    configfile = configdict['output_file'] + '/config/' +configdict["config_prefix"] + "_" + job + ".config"
    TopmedPipeline.writeConfig(config, configfile)
    if os.path.isfile(config["out_file"]) == False:
        cmd =  " ".join(["bsub -q big -R 'rusage[mem=50000]'","-o", logname, 'Rscript', rscript,  configfile])
        print(cmd)
        os.system(cmd)
else:
    assocScript = "assoc_single_unrel"

chrom_list = TopmedPipeline.parseChromosomes(chromosomes).split(" ")
########## for aggregate tests, generate variant list################
if assoc_type == "aggregate":
    job = "aggregate_list"
    rscript = os.path.join(pipeline, "R", job + ".R")
    config = deepcopy(configdict)
    config["out_file"] = configdict['output_file'] + '/data/' + configdict["data_prefix"] + "_" + job + "_chr .RData"
    configfile = configdict['output_file'] + '/config/' +configdict["config_prefix"] + "_" + job + ".config"
    TopmedPipeline.writeConfig(config, configfile)
    log_prefix = configdict['output_file'] + '/log/' + configdict['config_prefix'] + '_' + job
    for chrom in chrom_list:
        if os.path.isfile(config["out_file"].replace('chr ', 'chr'+chrom)) == False:
            cmd =  " ".join(['bsub -q short', "-o", log_prefix + ".log", 'Rscript', rscript,  configfile, "--chromosome " + chrom])
            os.system(cmd)

##################### define segments#########################
if segment_length == default_segment_length and n_segments is None:
    build = configdict.setdefault("genome_build", "hg38")
    segment_file = os.path.join(pipeline, "segments_" + build + ".txt")
    print("Using default segment file for build " + build + " with segment_length " + default_segment_length + " kb")
else:
    job = "define_segments"

    rscript = os.path.join(pipeline, "R", job + ".R")

    config = deepcopy(configdict)
    config["out_file"] = configdict['output_file'] + '/data/' + configdict["config_prefix"] + "_segments.txt"
    configfile = configdict['output_file'] + '/config/' +configdict["config_prefix"] + "_" + job + ".config"
    TopmedPipeline.writeConfig(config, configfile)

    segment_file = configdict['output_file'] + '/data/' +config["out_file"]
    logname = configdict['output_file'] + '/log/' + configdict['config_prefix'] + '_' + job + '.log'

    # run and wait for results
    print("Defining segments...")
#    log_file = job + "_" + strftime("%Y-%m-%d-%H-%M-%S", localtime()) + ".log"
    if n_segments is not None:
        cmd =  " ".join(['bsub -q short',  '-o', logname, 'Rscript', rscript,  configfile,  "--n_segments " + n_segments])
        os.system(cmd)
#        cmd = [driver, rscript, configfile, "--n_segments " + n_segments]
    else:
        cmd =  " ".join(['bsub -q short',  '-o', logname, 'Rscript', rscript,  configfile,  "--segment_length " + segment_length])
        os.system(cmd)
#        cmd = [driver, rscript, configfile, "--segment_length " + segment_length]
#    cluster.runCmd(job_name=job, cmd=cmd, logfile=log_file)

############# set up config for association test ##############
config = deepcopy(configdict)
config["assoc_type"] = assoc_type
config["null_model_file"] = configdict['output_file'] + '/data/' +configdict["data_prefix"] + "_null_model.RData"
if not single_unrel:
    config["phenotype_file"] = configdict['output_file'] + '/data/' +configdict["data_prefix"] + "_phenotypes.RData"
if assoc_type == "aggregate":
    config["aggregate_variant_file"] = configdict['output_file'] + '/data/' +configdict["data_prefix"] + "_aggregate_list_chr .RData"
config["out_prefix"] = configdict['output_file'] + '/results/' +configdict["data_prefix"] + "_" + assocScript
config["segment_file"] = segment_file
configfile = configdict['output_file'] + '/config/' +configdict["config_prefix"] + "_" + assocScript + ".config"
TopmedPipeline.writeConfig(config, configfile)

############### get segments for each chromosome##################
chrom_list = TopmedPipeline.parseChromosomes(chromosomes).split(" ")
segment_list = TopmedPipeline.getChromSegments(segment_file, chrom_list)
print(segment_list)
segment_str = ["-".join([str(i) for i in s]) for s in segment_list]
print(segment_str)
segments = dict(zip(chrom_list, segment_str))
print(segments)
################## run association tests ####################
#hold_combine = []
if not single_unrel:
    if os.path.isfile(configdict['output_file'] + '/data/' +configdict["data_prefix"] + "_null_model.RData") == False:
        print('null model output unfinished')
        quit()
check_file = configdict['output_file'] + '/log/' + configdict["data_prefix"] + "_" + assocScript + '.log'
for chromosome in chrom_list:
    run = True
    job_assoc = assocScript + "_chr" + chromosome
    rscript = os.path.join(pipeline, "R", assocScript + ".R")
    chrom_segments = TopmedPipeline.chromosomeRangeToList(segments[chromosome])
    chrom_segments = [str(x) for x in chrom_segments]
    if os.path.isfile(check_file) == True:
        run = False
    elif os.path.isfile(check_file) == False:
        if assoc_type == 'aggregate':
            if os.path.isfile(configdict['output_file'] + '/data/' +configdict["data_prefix"] + "_aggregate_list_chr22.RData")== False:
                print('Aggregate list files have not been generated yet')
                exit()
    if run == True:
        for segment in chrom_segments:
            logname = configdict['output_file'] + '/log/' + configdict['config_prefix'] + '_' + job_assoc + '_' + segment + '.log'
            cmd=" ".join(['bsub -q big -R ', "'rusage[mem=24000]'", '-o', logname, 'Rscript', rscript,  configfile, "--chromosome " + chromosome, '--segment ' + segment])
            print(cmd)
            os.system(cmd)

if run == True:
    output_log = open(check_file, 'w')
    output_log.write(assocScript+'.R running on all chromosome segments')
    output_log.close()
    print('First run - running segments of chromosomes. Run script again to combine segments for chromosome files')
    exit()
else:
    print('Second run - running combine script to join chromosomes')


for chromosome in chrom_list:
    combScript = "assoc_combine"
    job_comb = combScript + "_chr" + chromosome
    rscript = os.path.join(pipeline, "R", combScript + ".R")
    log = os.path.join(configdict['output_file'], 'log', configdict["data_prefix"]+ '_' + assocScript + "_chr" + chromosome + '.log')
#    args = [rscript, configfile, "--chromosome " + chromosome, version]
    check_file = configdict['output_file'] + '/results/' + configdict["data_prefix"] + "_" + assocScript + "_chr .RData"
    if not os.path.isfile(check_file.replace('chr ', 'chr'+chromosome)):
        cmd=" ".join(['bsub -q short', '-o', log, 'Rscript', rscript,  configfile, "--chromosome " + chromosome])
        print(cmd)
        os.system(cmd)

############job: assoc_plots################
job = "assoc_plots"

rscript = os.path.join(pipeline, "R", job + ".R")

config = deepcopy(configdict)
config["assoc_file"] = configdict['output_file'] + '/results/' + configdict["data_prefix"] + "_" + assocScript + "_chr .RData"
config["assoc_type"] = assoc_type
config["chromosomes"] = TopmedPipeline.parseChromosomes(chromosomes)
config["out_file_manh"] = configdict['output_file'] + '/plots/' +configdict["plots_prefix"] + "_manh.png"
config["out_file_qq"] =  configdict['output_file'] + '/plots/' +configdict["plots_prefix"] + "_qq.png"
configfile =  configdict['output_file'] + '/config/' +configdict["config_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(config, configfile)
logname = configdict['output_file'] + '/log/' + configdict['config_prefix'] + '_' + job + '.log' 

exit_flag = False
chrom_list = config['chromosomes'].strip().split()
for chrom in chrom_list:
    if os.path.isfile(config["assoc_file"].replace('chr ', 'chr' + chrom)) == False:
        exit_flag = True

if exit_flag == True:
    print('association files not generated yet')
    exit()

if os.path.isfile(config["out_file_qq"]) == False:
    cmd =  " ".join(['bsub -q big ', '-R', "'rusage[mem=11000]'", '-o', logname, 'Rscript', rscript,  configfile])
    print(cmd)
    os.system(cmd)
    print('Run again to generate report after plots are complete')
    exit()


###########job: assoc_report (analysis report)############
job = "assoc_report"

rscript = os.path.join(pipeline, "R", job + ".R")

config = deepcopy(configdict)
config["out_file"] = configdict['output_file'] + '/report/' +configdict["out_prefix"] + "_analysis_report"
configfile = configdict['output_file'] + '/config/' +configdict["config_prefix"] + "_" + job + ".config"
TopmedPipeline.writeConfig(config, configfile)
logname = configdict['output_file'] + '/log/' + configdict['config_prefix'] + '_' + job + '.log' 
cmd =  " ".join(['bsub -q short', "-o", logname, 'Rscript', rscript,  configfile])
print(cmd)
os.system(cmd)


