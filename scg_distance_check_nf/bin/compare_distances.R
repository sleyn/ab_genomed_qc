#!/usr/bin/env Rscript
print("Load libraries")
library(readr)
library(dplyr, warn.conflicts = FALSE)
library(stringr)
library(ggplot2)

# Read alignment file and calculate distances
# first argument - directory with distance tables
# second argument - output folder
args <- commandArgs(trailingOnly = TRUE)
dist_files <- Sys.glob(file.path(args[1], "*"))

read_distance <- function(dist_file){
  pgfam = str_remove(basename(dist_file), '.dist')
  return(
        read_tsv(
                  dist_file, 
                  show_col_types = F, 
                  col_types = cols(.default = col_character(), Distance = col_double())
                 ) %>% 
           select(-Protein) %>% 
           rename({{pgfam}} := Distance)
         )
}

print('Read distance tebles')
dist_combined <- dist_files %>% purrr::map(read_distance) %>% purrr::reduce(full_join, by="Genome")
#dist_combined <- dist_combined %>% select(-PGF_00067467)

print('Calculate sum, mean, sd and number of NAs')
dist_combined <- dist_combined %>% 
  rowwise() %>% 
  mutate(
    sum = sum(c_across(starts_with("PGF_")), na.rm = T),
    mean = mean(c_across(starts_with("PGF_")), na.rm = T),
    sd = sd(c_across(starts_with("PGF_")), na.rm = T),
    nas = sum(is.na(c_across(starts_with("PGF_"))))
    )

print('Filter all genomes that have number of NAs > 75% quantile')
dist_combined <- dist_combined %>% filter(nas <= quantile(dist_combined$nas, 0.75))

print('Plot mean vs sd counts')
p_mean_sd <- dist_combined %>% ggplot(aes(mean, sd)) + 
  geom_hex(aes(fill = log(..count..))) +
  xlab('Mean distance') +
  ylab('Mean standard deviation') +
  theme(legend.position = "none")

ggsave(file.path(args[2], 'dist_mean_vs_sd.pdf'), p_mean_sd)

print('Calculate distance p-values assuming normal distribution of distances')
dist_combined_pvalue <- dist_combined %>% select(-c('sum', 'mean', 'sd', 'nas')) %>% 
  tidyr::pivot_longer(-Genome, names_to = 'pgfam', values_to = 'distance') %>% 
  filter(!is.na(distance)) %>% 
  group_by(pgfam) %>% 
  mutate(pvalue = pnorm(distance, mean = mean(distance), sd = sd(distance))) %>% 
  select(-distance) %>% 
  tidyr::pivot_wider(names_from = pgfam, values_from = pvalue)

dist_combined_pvalue <- dist_combined_pvalue %>% rowwise() %>% 
  mutate(
    sum = sum(c_across(starts_with("PGF_")), na.rm = T),
    mean = mean(c_across(starts_with("PGF_")), na.rm = T),
    sd = sd(c_across(starts_with("PGF_")), na.rm = T)
  )

p_pvalue <- dist_combined_pvalue %>% 
  ggplot(aes(mean)) + 
  geom_histogram() +
  xlab('Mean p-value') +
  ggtitle('Mean p-value assuming normal distribution of\ndistances in each PGfam')

ggsave(file.path(args[2], 'dist_pvalue.pdf'), p_pvalue)

print('Save distance and p-value tables.')
dist_combined %>% write_tsv(file.path(args[2], 'distance_combined.tsv'))
dist_combined_pvalue %>% filter(mean <= 0.6) %>% write_tsv(file.path(args[2], 'distance_pvalue_combined.tsv'))
