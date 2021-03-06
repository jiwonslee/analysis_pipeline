
#' @importFrom dplyr "%>%" arrange_ select_
.arrange_chr_pos <- function(df, chr="chr", pos="pos") {
    df$chr.factor <- factor(df[[chr]], levels=c(1:22, "X", "Y"))
    df %>%
        arrange_("chr.factor", pos) %>%
        select_("-chr.factor")
}

.index_chr_pos <- function(x, chr="chr", pos="pos") {
    df <- do.call(rbind, lapply(x, function(xx) xx[1,c(chr,pos)]))
    df$chr.factor <- factor(df[[chr]], levels=c(1:22, "X", "Y"))
    order(df$chr.factor, df[[pos]])
}


#' Combine association test results
#' 
#' Combine association test results from multiple files into a single object.
#' 
#' Useful for combining per-segment results into a single file per chromosome.
#' 
#' @param files Vector of file names with association test results
#' @param assoc_type Type of association test ("single", "aggregate", "window")
#' @return Association test object
#'
#' @importFrom dplyr "%>%" distinct_ filter_ group_by_
#' @export
combineAssoc <- function(files, assoc_type) {
    stopifnot(assoc_type %in% c("single", "aggregate", "window"))
    x <- lapply(unname(files), getobj)
    if (assoc_type == "single") {
        assoc <- do.call(rbind, x) %>%
            .arrange_chr_pos()
    } else if (assoc_type  == "aggregate") {
        assoc <- x[[1]][c("param", "nsample")]
        varInfo <- do.call(c, lapply(x, function(y) y$variantInfo))
        # get index to put units in order by chr, pos
        index <- .index_chr_pos(varInfo)
        assoc$results <- do.call(rbind, lapply(x, function(y) y$results))[index,]
        assoc$variantInfo <- varInfo[index]
    } else if (assoc_type == "window") {
        assoc <- x[[1]][c("param", "window", "nsample")]
        assoc$results <- do.call(rbind, lapply(x, function(y) y$results)) %>%
            .arrange_chr_pos(pos="window.start") %>%
            distinct_() %>%
            group_by_("chr", "window.start", "window.stop") %>%
            filter_(~(n.site == max(n.site)), ~(!duplicated(n.site))) %>%
            as.data.frame()
        assoc$variantInfo <- do.call(rbind, lapply(x, function(y) y$variantInfo)) %>%
            .arrange_chr_pos() %>%
            distinct_()
    }
    assoc
}

#' Get association test results
#'
#' Return association test results in a standard format
#'
#' Read association test results in multiple files and combine all into a single
#' data frame with standard column names.
#'
#' If a single aggregate unit contains variants from multiple chromosomes, each chromosome will have its own row in the output.
#' 
#' @inheritParams combineAssoc
#' @return data.frame including standard columns ("chr", "pos", "start", "end", "stat", "pval")
#'
#' @importFrom dplyr "%>%" filter_ group_by_ left_join mutate_ n rename_ select_ summarise_
#' @export
getAssoc <- function(files, assoc_type) {
    stopifnot(assoc_type %in% c("single", "aggregate", "window"))
    assoc <- do.call(rbind, lapply(unname(files), function(f) {
        x <- getobj(f)
        if (assoc_type  == "aggregate") {
            tmp <- x$results %>%
                mutate_(group_id=~(1:n())) %>%
                filter_(~(n.site > 0))
            group.info <- do.call(rbind, lapply(tmp$group_id, function(g) {
                x$variantInfo[[g]] %>%
                    group_by_("chr") %>%
                    summarise_(start=~min(pos),
                               end=~max(pos),
                               pos=~(floor((min(pos) + max(pos))/2))) %>%
                    mutate_(group_id=g) %>%
                    as.data.frame()
            }))
            x <- left_join(tmp, group.info, by="group_id")
        } else if (assoc_type == "window") {
            x <- filter_(x$results, ~(n.site > 0), ~(dup == 0)) %>%
                mutate_(pos=~(floor((window.start + window.stop)/2))) %>%
                rename_(start="window.start", end="window.stop")
        } else {
            x <- mutate_(x, start="pos", end="pos")
        }
        x
    }))
    
    if ("pval_0" %in% names(assoc)) {
        ## SKAT
        pval.col <- if ("pval_SKATO" %in% names(assoc)) "pval_SKATO" else "pval_0"
        assoc <- select_(assoc, "chr", "pos", "start", "end", pval.col) %>%
            rename_(pval=pval.col)
    } else {
        ## burden or single
        assoc <- select_(assoc, "chr", "pos", "start", "end", ~ends_with("stat"), ~ends_with("pval"))
        names(assoc)[5:6] <- c("stat", "pval")
    }
    assoc <- filter_(assoc, ~(!is.na(pval))) %>%
        mutate_(chr=~factor(chr, levels=c(1:22, "X")))
    assoc
}


#' Format single-variant assocation test results
#'
#' Return association test results for single-variant tests in a standard format
#'
#' Ensures that single-variant association test results from different sources
#' (\code{\link[GENESIS]{assocTestMM}}, \code{\link[SeqVarTools]{regression}})
#' will all have a standard format.
#'
#' @param seqData A \code{\link[SeqArray]{SeqVarGDSClass}} object
#'   (needed to get chromosome and position if not present)
#' @param assoc data.frame with assocation test results
#' @return data.frame including standard columns ("variantID", "chr", "pos", "n", "MAF", "minor.allele")
#' 
#' @import SeqArray
#' @importFrom dplyr "%>%" mutate_ mutate
#' @export
formatAssocSingle <- function(seqData, assoc) {

    names(assoc)[names(assoc) %in% c("snpID", "variant.id")] <- "variantID"
    names(assoc) <- sub(".Stat", ".stat", names(assoc), fixed=TRUE)
    names(assoc) <- sub(".Pval", ".pval", names(assoc), fixed=TRUE)
    
    seqSetFilter(seqData, variant.id=assoc$variantID, action="push+set", verbose=FALSE)
    assoc$pos <- seqGetData(seqData, "position")
    assoc$snpid <- seqGetData(seqData, "annotation/id")
    assoc$alleles <- seqGetData(seqData, "allele") ##split by comma and add in ref alt later
    assoc <- assoc %>% mutate(effect.allele=sapply(strsplit(alleles, ","), "[", 2))
    if (!("chr" %in% names(assoc))) {
        assoc$chr <- seqGetData(seqData, "chromosome")
    }
    if ("n0" %in% names(assoc)) {
        assoc$n <- rowSums(assoc[,c("freq0", "freq1")], na.rm=TRUE)
        assoc$freq <- rowMeans(assoc[,c("freq0", "freq1")], na.rm=TRUE)
    }
    if (!("MAF" %in% names(assoc))) {
        assoc$MAF <- pmin(assoc$freq, 1 - assoc$freq)
        assoc$minor.allele <- ifelse(assoc$freq > 0.5, "ref", "alt")
    }
    seqSetFilter(seqData, action="pop", verbose=FALSE)
    
    init.cols <- c("variantID", "chr", "pos", "snpid", "alleles", "n", "MAF", "minor.allele", "effect.allele")
    cols <- setdiff(names(assoc), c(init.cols, "freq"))
    assoc <- assoc[,c(init.cols, cols)]

    assoc
}


#' Omit known hits from an association test data frame
#'
#' @param assoc data.frame with assocation test results (including columns chr, pos)
#' @param hits data.frame with known hits (including columns chr, pos)
#' @param flank Number of kb on either side of each known hit to exclude
#' @return data.frame with assocation test results not in regions around known hits
#'
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors queryHits
#' @export
omitKnownHits <- function(assoc, hits, flank=500) {
    stopifnot(all(c("chr", "pos") %in% names(hits)))
    assoc.gr <- GRanges(seqnames=assoc$chr, ranges=IRanges(start=assoc$pos, end=assoc$pos))
    hits.gr <- GRanges(seqnames=hits$chr, ranges=IRanges(start=hits$pos-(flank*1000),
                                              end=hits$pos+(flank*1000)))
    ol <- findOverlaps(assoc.gr, hits.gr)
    assoc[-queryHits(ol),]
}



#' Add MAC column to association test output
#'
#' @param assoc results from \code{\link{assocTestSingle}} or \code{\link{assocTestAggregate}}
#' @param assoc_type Type of association test ("single", "aggregate", "window")
#' @return \code{assoc} with a "MAC" column added to the results data.frame
#'
#' @export
addMAC <- function(assoc, assoc_type) {
  mac <- function(x) {
    round(2 * x$n * pmin(x$MAF, 1-x$MAF))
  }
  if (assoc_type == "single") {
    assoc$MAC <- mac(assoc)
  } else if (assoc_type %in% c("aggregate", "window")) {
    assoc$results$MAC <- sapply(assoc$variantInfo, function(x) sum(mac(x)))
  }
  assoc
}



