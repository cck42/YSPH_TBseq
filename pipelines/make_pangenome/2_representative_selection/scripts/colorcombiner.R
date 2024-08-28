#!/usr/bin/env Rscript

library(dplyr)
library("argparse")

parser <- ArgumentParser(description = 'Reformat output from Parnas')

parser$add_argument('--representatives-file','-r',help='File containing a list of the chosen representatives') #need to put better explanations here
parser$add_argument('--colors-file','-c',help='Cluster groups for each strain')
parser$add_argument('--output','-o',help='output')

xargs <- parser$parse_args()

representatives <- read.delim(xargs$representatives_file,sep="\t",header=FALSE,col.names=c('strain'))
colors <- read.delim(xargs$colors_file,sep="\t",header=FALSE,col.names=c('strain','group'))

#Label each representative with its group
reps_labeled <- representatives %>%
  left_join(colors,by='strain')

colnames <- list(
  "strain"="strain.x",
  group="group",
  "representative"="strain.y"
)
all_relabeled <- colors %>% 
  left_join(reps_labeled, by='group') %>%
  rename(!!!colnames) %>%
  select(strain,representative)

n_per_rep <- all_relabeled %>%
  group_by(representative) %>%
  summarise(number=n())

write.table(all_relabeled,xargs$output,quote=FALSE,row.names=FALSE,sep="\t")





