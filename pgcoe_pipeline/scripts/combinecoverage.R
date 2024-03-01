

library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(argparse)

statsdirs <- c('/Users/chaneykalinich/Documents/TB-sequencing/github/tb-seq/tb_pipeline/results/Yale-TB001_coveragestats_sub1.0.csv',
               '/Users/chaneykalinich/Documents/TB-sequencing/github/tb-seq/tb_pipeline/results/Yale-TB002_coveragestats_sub1.0.csv',
               '/Users/chaneykalinich/Documents/TB-sequencing/github/tb-seq/tb_pipeline/results/Yale-TB003_coveragestats_sub1.0.csv',
               '/Users/chaneykalinich/Documents/TB-sequencing/github/tb-seq/tb_pipeline/results/Yale-TB004_coveragestats_sub1.0.csv',
               '/Users/chaneykalinich/Documents/TB-sequencing/github/tb-seq/tb_pipeline/results/Yale-TB005_coveragestats_sub1.0.csv',
               '/Users/chaneykalinich/Documents/TB-sequencing/github/tb-seq/tb_pipeline/results/Yale-TB006_coveragestats_sub1.0.csv',
               '/Users/chaneykalinich/Documents/TB-sequencing/github/tb-seq/tb_pipeline/results/Yale-TB007_coveragestats_sub1.0.csv',
               '/Users/chaneykalinich/Documents/TB-sequencing/github/tb-seq/tb_pipeline/results/Yale-TB008_coveragestats_sub1.0.csv',
               '/Users/chaneykalinich/Documents/TB-sequencing/github/tb-seq/tb_pipeline/results/Yale-TB009_coveragestats_sub1.0.csv')


statsheaders <- c('sample','subsample','rname','startpos','endpos','numreads','covbases','coverage','meandepth','meanbaseq','meanmapq')

for sampkle