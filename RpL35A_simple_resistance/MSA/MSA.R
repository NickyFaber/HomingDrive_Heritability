rm(list = ls())

########################################
########## Environment set-up ##########
########################################

setwd("/Users/s2018147/Library/CloudStorage/OneDrive-WageningenUniversity&Research/PhD/03_Genetic_diversity/r1_r2_PCRs/Plasmidsaurus sequencing")

library(msa, quietly = TRUE)
library(tidyverse, quietly = TRUE)
library(patchwork, quietly = TRUE)

mySequences <- readDNAStringSet("consensus_haplotypes.fasta")
mySequences

myAlignment <- msa(mySequences, order="input", 
                        method = "Muscle")
myAlignment
print(myAlignment, show="complete")

smallSubset <- myAlignment
largeSubset <- myAlignment
everything <- myAlignment

msaPrettyPrint(smallSubset, y = c(323,417),
               output="pdf", showNames="left", showConsensus="none", showNumbering="none",
               shadingMode="identical", shadingColors="greens", 
               showLogo="none", showLegend=FALSE, 
               askForOverwrite=FALSE, verbose=FALSE)
msaPrettyPrint(largeSubset, y = c(323,512),
               output="pdf", showNames="left", showConsensus="none", showNumbering="none",
               shadingMode="identical", shadingColors="greens", 
               showLogo="none", showLegend=FALSE, 
               askForOverwrite=FALSE, verbose=FALSE)
msaPrettyPrint(everything, 
               output="pdf", showNames="left", showConsensus="none", showNumbering="none",
               shadingMode="identical", shadingColors="greens", 
               showLogo="none", showLegend=FALSE, 
               askForOverwrite=FALSE, verbose=FALSE)












