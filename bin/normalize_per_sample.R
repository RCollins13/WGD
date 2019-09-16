#!/usr/bin/env Rscript

# Copyright (c) 2019 Talkowski Laboratory
# Contact: Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# Normalize a bincov matrix by per-sample medians


### Set master parameters
options(stringsAsFactors=F, scipen=1000)


####################
### HELPER FUNCTIONS
####################
# Read & format binCov file
read.bincov <- function(bincov.in){
  bincov <- read.table(bincov.in, header=T, sep="\t", comment.char="")
  colnames(bincov)[1] <- "chr"
  bincov[, -c(1:3)] <- apply(bincov[, -c(1:3)], 2, as.numeric)
  return(bincov)
}

# Read & format medians
read.medians <- function(medians.in){
  medians <- read.table(medians.in, header=T, sep="\t", comment.char="")
  samples <- colnames(medians)
  medians <- as.numeric(medians[1, ])
  names(medians) <- samples
  return(medians)
}

# Normalize bincov by sample medians
normalize.bincov <- function(bincov, medians){
  bincov <- bincov[, c(1:3, which(colnames(bincov) %in% names(medians)))]
  norm.vals <- as.data.frame(sapply(4:ncol(bincov), function(i){
    vals <- bincov[, i]
    sample <- colnames(bincov)[i]
    med <- unique(as.numeric(medians[which(names(medians) == sample)]))
    if(length(med)==1){
      return(round(2*vals/med, 3))
    }
  }))
  colnames(norm.vals) <- colnames(bincov)[-c(1:3)]
  norm.bincov <- as.data.frame(cbind(bincov[, 1:3],
                                     norm.vals))
  return(norm.bincov)
}


################
### RSCRIPT BLOCK
################
require(optparse, quietly=T)
### Get command-line arguments & options
option_list <- list()
args <- parse_args(OptionParser(usage="%prog BINCOV MEDIANS OUTFILE", 
                                option_list=option_list), 
                   positional_arguments=TRUE)
opts <- args$options

### Checks for appropriate positional arguments
if(length(args$args) != 2){
  stop("Incorrect number of required positional arguments\n")
}

### Writes args & opts to vars
bincov.in <- args$args[1]
medians.in <- args$args[2]
OUTFILE <- args$args[3]

# #Dev parameters (local)
# bincov.in <- "~/scratch/gnomAD_v2_SV_PCRPLUS_Q1_batch_1.binCov.bed.gz"
# medians.in <- "~/scratch/gnomAD_v2_SV_PCRPLUS_Q1_batch_1.binCov_medians.txt"
# OUTFILE <- "~/scratch/gnomAD_v2_SV_PCRPLUS_Q1_batch_1.binCov.norm.bed"

### Read input data
bincov <- read.bincov(bincov.in)
medians <- read.medians(medians.in)

### Normalize bincov
norm.bincov <- normalize.bincov(bincov, medians)

### Write out
colnames(norm.bincov)[1] <- c("#chr")
write.table(norm.bincov, OUTFILE, col.names=T, row.names=F, sep="\t", quote=F)
