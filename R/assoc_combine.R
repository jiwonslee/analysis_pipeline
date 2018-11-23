library(argparser)
library(TopmedPipeline)
sessionInfo()

argp <- arg_parser("Association test - combine files")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--chromosome", help="chromosome (1-24 or X,Y)", type="character")
argp <- add_argument(argp, "--version", help="pipeline version number")
argv <- parse_args(argp)
cat(">>> TopmedPipeline version ", argv$version, "\n")
config <- readConfig(argv$config)
chr <- intToChr(argv$chromosome)

required <- c("assoc_type",
              "out_prefix")
optional <- c("output_file" = NA)
config <- setConfigDefaults(config, required, optional)
print(config)

if (!is.na(config["output_file"])){
  params_path = paste0(file.path(config["output_file"], "report", basename(argv$config)), ".assoc_combine.params")
} else {
  params_path = paste0(basename(argv$config), ".assoc_combine.params")
}
writeConfig(config, params_path)

file.pattern <- paste0(basename(config["out_prefix"]), "_chr", chr, "_seg[[:digit:]]+.RData")
files <- list.files(path=dirname(config["out_prefix"]), pattern=file.pattern, full.names=TRUE)

assoc <- combineAssoc(files, config["assoc_type"])

save(assoc, file=constructFilename(config["out_prefix"], chr))

# delete segment files
unlink(files)

# mem stats
ms <- gc()
cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")
