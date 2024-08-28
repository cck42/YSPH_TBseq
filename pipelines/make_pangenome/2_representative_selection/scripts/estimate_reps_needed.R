#!/usr/bin/env Rscript

library(argparse)
library(tidyverse)
library(ggplot2)
library(dplyr)

# argparse initialization
parser <- ArgumentParser(description= 'Create reps needed to account for max diversity plot')

parser$add_argument('--input', '-i', help= 'Input reps accounting for diversity file')
parser$add_argument('--output', '-o', help= 'Output plot')

xargs<- parser$parse_args()

# read in csv
diversity_dat <- read.csv(xargs$input)

# plot with ggplot
divfig <- diversity_dat %>%
  ggplot() +
  geom_point(aes(x = Representatives, y = Diversity_covered)) +
  theme_light() +
  labs(x = "Number of representative genomes", y = "% of total tree diversity covered")

# save plot
jpeg(filename = xargs$output)
print(divfig)
dev.off()
