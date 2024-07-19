library(tidyverse)
library(readr)

tb003_mykrobe <- read.csv('/Users/chaneykalinich/Documents/PGCoE/github/pgcoe_pipeline/pgcoe_pipeline/results/combined.csv')
refdata_mykrobe <- read.csv('/Users/chaneykalinich/Documents/PGCoE/github/tb-comparison-data/results/combined.csv')
srr_sample <- read.csv('/Users/chaneykalinich/Documents/PGCoE/github/tb-comparison-data/data/matched_ids.csv')
pgcoe_datasheet <- read.csv('/Users/chaneykalinich/Documents/PGCoE/github/tb-comparison-data/data/PGCOE_Datasheet-M.tuberculosis_Seq.csv') %>%
  filter(NGS_Run_ID=='TB003') %>%
  mutate(Original_ID=as.factor(Original_ID)) 


tb003_refsamples <-levels(pgcoe_datasheet$Original_ID)
tb003_refdata_mykrobe <- refdata_mykrobe %>%
  filter(sample %in% tb003_refsamples)


tb003_ref_mykrobe_combined <- tb003_mykrobe %>%
  left_join(pgcoe_datasheet,by=c(sample='Seq_ID')) %>%
  filter(NGS_primers=='TBv1.1') %>%
  left_join(refdata_mykrobe,by=c(Original_ID='sample',drug='drug')) %>%
  select(sample,Original_ID,drug,susceptibility.x,susceptibility.y,"variants..dna_variant.AA_variant.ref_kmer_count.alt_kmer_count.conf...use...format.json.for.more.info..x",
         "variants..dna_variant.AA_variant.ref_kmer_count.alt_kmer_count.conf...use...format.json.for.more.info..y",phylo_group.x,phylo_group.y,species.x,species.y,lineage.x,lineage.y,
         phylo_group_per_covg.x,phylo_group_per_covg.y,species_per_covg.x,species_per_covg.y,phylo_group_depth.x,phylo_group_depth.y,species_depth.x,species_depth.y,
         CT,Template_dilution) %>%
  mutate(agree=if_else(susceptibility.x==susceptibility.y,'YES','NO'))

ctvcov_withres <- tb003_ref_mykrobe_combined %>%
  mutate(CT=replace_na(CT,40)) %>%
  ggplot(aes(x=CT,y=phylo_group_per_covg.x)) +
    geom_point(aes(color=agree)) +
  facet_wrap('drug')
plot(ctvcov_withres)
