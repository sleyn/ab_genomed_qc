#!/usr/bin/env Rscript
print("Load libraries")
library(seqinr)
library(tibble)
library(dplyr, warn.conflicts = FALSE)
library(tidyr)
library(readr)
library(stringr)

# Read alignment file and calculate distances
args <- commandArgs(trailingOnly = TRUE)
print(paste("Read file", args[1]))
aln <- read.alignment(args[1], "fasta")
print("Calculate distance matrix")
dist_matrix <- dist.alignment(aln, "identity")

# Look at distances from the reference sequence
print("Extract distances to the reference sequence")
target_dist <- as.matrix(dist_matrix)["reference", ]
distances <- tibble(Distance = target_dist, Protein = names(target_dist)) %>%
  mutate(Genome = str_match(Protein, "fig\\|(\\d+\\.\\d+)\\.")[, 2]) %>%
  arrange(Distance) %>%
  group_by(Genome) %>%
  filter(row_number() == 1, Protein != "reference") %>%
  write_tsv(args[2])
