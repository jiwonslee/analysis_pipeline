library(argparser)
library(TopmedPipeline)
library(Biobase)
library(GENESIS)
library(gdsfmt)
sessionInfo()

argp <- arg_parser("Null model for association tests")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--version", help="pipeline version number")
argv <- parse_args(argp)
cat(">>> TopmedPipeline version ", argv$version, "\n")
config <- readConfig(argv$config)

required <- c("outcome",
              "phenotype_file")
optional <- c("gds_file"=NA, # required for conditional variants
              "pca_file"=NA,
              "pcrelate_file"=NA,
              "grm_file"=NA,
              "binary"=FALSE,
              "conditional_variant_file"=NA,
              "covars"=NA,
              "group_var"=NA,
              "inverse_normal"=TRUE,
              "n_pcs"=3,
              "out_file"="null_model.RData",
              "out_phenotype_file"="phenotypes.RData",
              "rescale_variance"="marginal",
              "resid_covars"=TRUE,
              "sample_include_file"=NA,
              "output_file" = NA,
              "meta" = FALSE)
config <- setConfigDefaults(config, required, optional)
print(config)
if (!is.na(config["output_file"])){
  params_path = paste0(file.path(config["output_file"], "report", basename(argv$config)), ".null_model.params")
} else {
  params_path = paste0(basename(argv$config), ".null_model.params")
}

writeConfig(config, params_path)

# get phenotypes
phen <- getPhenotypes(config)
annot <- phen[["annot"]]
outcome <- phen[["outcome"]]
covars <- phen[["covars"]]
group.var <- phen[["group.var"]]
sample.id <- phen[["sample.id"]]

save(annot, file=config["out_phenotype_file"])

if (as.logical(config["binary"])) {
  stopifnot(all(annot[[outcome]] %in% c(0,1,NA)))
  family <- binomial
} else {
  family <- gaussian
}

message('getting grm or kinship matrix...')
# kinship matrix or GRM
grm <- getGRM(config, sample.id)

# print model
random <- if (!is.na(config["pcrelate_file"])) "kinship" else if (!is.na(config["grm_file"])) "GRM" else NULL
model.string <- modelString(outcome, covars, random, group.var)
message("Model: ", model.string)
message(length(sample.id), " samples")

###rewrite null model function so that rescaling is not done separately##
if (!is.null(grm)) {
  ## fit null model allowing heterogeneous variances among studies
  nullmod <- fitNullMM(annot, outcome=outcome, covars=covars,
                       covMatList=grm, scan.include=sample.id,
                       family=family, group.var=group.var)
  ## if we need an inverse normal transform, take residuals and refit null model
  if (as.logical(config["inverse_normal"])) {
 #   if (is.null(group.var)) {
    annot <- addInvNorm(annot, nullmod, outcome, covars)
  #  } else {
   #   groups <- unique(annot[[group.var]])
    #  ## inverse-normal transform residuals from each study separately (mean=0, var=1)
  #    resid.group <- do.call(rbind, lapply(groups, function(g) {
  #      samp.g <- intersect(nullmod$scanID, annot$sample.id[annot[[group.var]] == g])
   #     resid.g <- nullmod$resid.marginal[nullmod$scanID %in% samp.g]
    #    resid.norm <- rankNorm(resid.g)
        ## rescale the inverse-normal residuals to have study-specific variances =
        ## kinship variance component + study-specific residual
    #     if (config["rescale_variance"] == "varcomp") {
    #       resid.scale <- nullmod$varComp["V_A"] + nullmod$varComp[paste0("V_", g)]
    #       resid.norm <- resid.norm * sqrt(resid.scale)
    #     } else if (config["rescale_variance"] == "marginal") {
    #       resid.scale <- var(resid.g)
    #       resid.norm <- resid.norm * sqrt(resid.scale)
    #     }
    #     data.frame(sample.id=samp.g, resid.norm, stringsAsFactors=FALSE)
    #   }))
    #   annot$resid.norm <- resid.group$resid.norm[match(annot$sample.id, resid.group$sample.id)]
    # }
    # 
    ## fit null model again with these residuals as outcome and allowing heterogeneous variances
    if (as.logical(config['meta'])){
      message('meta=TRUE; rescaling rank normalized residuals by sd of residuals...')
      scale <- sqrt(var(nullmod$resid.marginal))
      annot$resid.norm <- annot$resid.norm*scale
    }
    resid.covars <- if (config["resid_covars"]) covars else NULL
    nullmod <- fitNullMM(annot, outcome="resid.norm", covars=resid.covars,
                         covMatList=grm, scan.include=sample.id,
                         family=family, group.var=group.var)
  }
} else {
  message("No kinship file specified, assuming samples are unrelated.")
  
  nullmod <- fitNullReg(annot, outcome=outcome, covars=covars,
                        scan.include=sample.id, family=family)
  
  if (as.logical(config["inverse_normal"])) {
    annot <- addInvNorm(annot, nullmod, outcome, covars)
    resid.covars <- if (config["resid_covars"]) covars else NULL
    nullmod <- fitNullReg(annot, outcome="resid.norm", covars=resid.covars,
                          scan.include=sample.id, family=family)
  }
  
}

save(nullmod, file=config["out_file"])

# mem stats
ms <- gc()
cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")





## fit null model allowing heterogeneous variances among studies - much too slow. reverting back to fitNullMM
#nullmod <- fitNullModel(annot, outcome=outcome, covars=covars,
#                        cov.mat=grm, sample.id=sample.id,
#                        family=family, group.var=group.var)
#nullmod <- fitNullMM(annot, outcome=outcome, covars=covars,
#                        covMatList=grm, scan.include=sample.id,
#                        family=family, group.var=group.var)

#####note: new GENESIS functions are faulty. reverting back to old version####
## if we need an inverse normal transform, take residuals and refit null model
#if (as.logical(config["inverse_normal"]) & !as.logical(config["binary"])) {
#  norm.option <- "all"
#  if (is.null(group.var)) {
#    rescale <- "none"
#  } else {
#    if (config["rescale_variance"] == "varcomp") {
#      rescale <- "model"
#    } else if (config["rescale_variance"] == "marginal") {
#      rescale <- "residSD"
#    } else{
#      print('No rescale_variance parameter provided; setting default to config parameter rescale_variance = varcomp / nullModelInvNorm function parameter rescale=model')
#      rescale <- 'model'
#    }
#  }
  
#  nullmod <- nullModelInvNorm(nullmod, cov.mat=grm,
#                              norm.option=norm.option,
#                              rescale=rescale, group.var=group.var)
#}

#save(nullmod, file=config["out_file"])

# mem stats
#ms <- gc()
#cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")
