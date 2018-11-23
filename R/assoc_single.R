library(argparser)
library(TopmedPipeline)
library(SeqVarTools)
library(Biobase)
library(GENESIS)
sessionInfo()

argp <- arg_parser("Association test - single variant")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--chromosome", help="chromosome (1-24 or X,Y)", type="character")
argp <- add_argument(argp, "--segment", help="segment number", type="integer")
argp <- add_argument(argp, "--version", help="pipeline version number")
argv <- parse_args(argp)
cat(">>> TopmedPipeline version ", argv$version, "\n")
config <- readConfig(argv$config)
chr <- intToChr(argv$chromosome)
segment <- argv$segment

required <- c("gds_file",
              "null_model_file",
              "phenotype_file")
optional <- c("mac_threshold"=5, # takes precedence
              "maf_threshold"=0.001,
              "out_prefix"="assoc_single",
              "pass_only"=TRUE,
              "segment_file"=NA,
              "test_type"="score",
              "variant_include_file"=NA,
              "output_file" = NA)
config <- setConfigDefaults(config, required, optional)
print(config)
###insert report output destination here for Rmd to work###
if (!is.na(config["output_file"])){
  params_path = paste0(file.path(config["output_file"], "report", basename(argv$config)), ".assoc_single.params")
} else {
  params_path = paste0(basename(argv$config), ".assoc_single.params")
}

writeConfig(config, params_path)

## gds file can have two parts split by chromosome identifier
gdsfile <- config["gds_file"]
varfile <- config["variant_include_file"]
if (!is.na(chr)) {
    bychrfile <- grepl(" ", gdsfile) # do we have one file per chromosome?
    gdsfile <- insertChromString(gdsfile, chr)
    varfile <- insertChromString(varfile, chr)
}

gds <- seqOpen(gdsfile)


# get phenotypes
annot <- getobj(config["phenotype_file"])

# createSeqVarData object
seqData <- SeqVarData(gds, sampleData=annot)
# get null model
nullModel <- getobj(config["null_model_file"])

# get samples included in null model
sample.id <- nullModel$sample.id

if (!is.na(segment)) {
  filterBySegment(seqData, segment, config["segment_file"])
}

if (!is.na(varfile)) {
  filterByFile(seqData, varfile)
}

## if we have a chromosome indicator but only one gds file, select chromosome
if (!is.na(chr) && !bychrfile) {
  filterByChrom(seqData, chr)
}

if (as.logical(config["pass_only"])) {
  filterByPass(seqData)
}

## MAC/MAF filtering
mac.min <- as.numeric(config["mac_threshold"])
maf.min <- as.numeric(config["maf_threshold"])
if (!is.na(mac.min)) {
  filterByMAC(seqData, sample.id, mac.min=mac.min)
} else {
  filterByMAF(seqData, sample.id, maf.min=maf.min)
}

checkSelectedVariants(seqData)
###change from 'variant.id' to 'annotation/id' results in subsequent error:
#Error in getSnpIndex(data = genoData, snp.include, chromosome) :
#  None of the SNPs in snp.include are in the provided data
#Calls: assocTestMM -> getSnpIndex
#Execution halted
#variant.id <- seqGetData(gds, "variant.id")
##variant.id <- seqGetData(gds, "annotation/id")
#seqResetFilter(gds, verbose=FALSE)



test <- switch(tolower(config["test_type"]),
               score="Score",
               wald="Wald")

assoc <- assocTestMM(seqData, nullModel, test=test)#, snp.include=variant.id)

## make output consistent with aggregate tests
message('formatting results...')
assoc <- formatAssocSingle(seqData, assoc)
message('adding MAC...')
#save(assoc, file=constructFilename(config["out_prefix"], chr, segment))
assoc <- addMAC(assoc, "single")
message('Saving results...')

save(assoc, file=constructFilename(config["out_prefix"], chr, segment))

seqClose(seqData)

# mem stats
ms <- gc()
cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")
