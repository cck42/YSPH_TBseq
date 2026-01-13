rm(list = ls())

library("tidyverse")
library(RColorBrewer)
library(patchwork)
library(ggtext)


# data import  ------------------------------------------------------------

qpcr_results <- read.table("data/qpcr_data.tsv",sep="\t",header=T) %>%
  select(c(Sample_ID,starting_quantity))

metadataTB <- read.csv("data/M.tuberculosis.csv") %>% 
  rename("Sample_Dilution"="Template_dilution") %>% 
  select(c("Seq_ID", "Original_ID", "Sample_Type", 
           "Sample_Dilution",  "Primer_Conc",  
           "NGS_Prep_Method", "CT","starting_quant")) %>% 
  rename(start_quant=starting_quant) %>%
  rename_with(tolower) %>%
  mutate("species"="M.tuberculosis") 

covtabTB1 <- read.table("data/TB.coverage.tsv",sep="\t",header=T) %>% 
  rename_with(tolower)


covtabTB2 <- read.table("data/TB005-6_combinedcoverage.txt",sep="\t",header=T) %>% 
  rename_with(tolower) %>%
  rename(seq_id = sample)


# TB amp vs unamp ---------------------------------------------------------
covtab <- rbind(covtabTB1,covtabTB2) %>% 
  mutate(seq_id=gsub("Yale-","",seq_id),
         subsample = as.numeric(subsample)) %>%
  filter(subsample==1)

metatab <- metadataTB %>% 
  mutate(seq_id=gsub("Yale-","",seq_id)) %>%
  filter(primer_conc %in% c("100uM","200uM","0")) %>% 
  filter(!is.na(start_quant)) %>% 
  filter(!is.na(sample_dilution)) %>%
  filter(!original_id %in% c("111-2565-18")) %>%
  filter(sample_type %in% c("Culture Isolate","DNA")) %>%
  filter(!(is.na(ct) & ngs_prep_method=='COVIDseq')) #remove these CDC dilutions which were re-sequenced

dilutiontab <- merge(metatab,covtab,by="seq_id") %>% mutate(group=paste(original_id,ngs_prep_method))

dilutiontab$ngs_prep_method <- fct_recode(dilutiontab$ngs_prep_method,
                                          "Unamplified"="Hybrid CovidSeq without amplification",
                                          "Amplicon"="Hybrid COVIDseq",
                                          "Amplicon"="Mpox/Hybrid CovidSeq",
                                          "Amplicon"="COVIDseq")

hascf <- dilutiontab %>% filter(ngs_prep_method=="Unamplified") %>% pull(original_id) %>% unique()
dilutiontab <- dilutiontab %>% filter(original_id %in% hascf)

cov_sqtab <- dilutiontab %>% select(original_id,sample_dilution,ngs_prep_method,coverage) %>%
  arrange(sample_dilution) %>% 
  filter(sample_dilution > 1e-09) %>% 
  pivot_wider(names_from = c(sample_dilution,ngs_prep_method), values_from = coverage)

write.table(cov_sqtab,"amp_unamp_coverage_table.tsv",sep="\t",row.names=F,quote=F)

dep_sqtab <- dilutiontab %>% select(original_id,sample_dilution,ngs_prep_method,meandepth) %>%
  arrange(sample_dilution) %>% 
  filter(sample_dilution > 1e-09) %>% 
  pivot_wider(names_from = c(sample_dilution,ngs_prep_method), values_from = meandepth)

write.table(cov_sqtab,"amp_unamp_depth_table.tsv",sep="\t",row.names=F,quote=F)

